"""
PyMOL & Mol* visualisation utilities for MD minimisation outputs.

Provides reusable functions for loading, aligning, colouring, and
generating interactive Web viewer files (MVSJ) from structures
produced by the Amber MD pipeline (§9).  Designed to be copied into
each project via ``cp -r templates/MD/*`` and imported by
project-specific entry-point scripts.

Functions
---------
process_final_min(project, experiment, job)
    Load stripped PDB files from the final-minimisation stage, align
    them, apply a publication-quality colour scheme, and save two
    PyMOL session files (.pse): one base-coloured and one with
    constraint-coloured residues.

load_constraints(project, experiment)
    Parse contact constraints from the experiment's YAML file and
    return deduplicated residue-pair tuples.

color_constrained_residues(constraints)
    Colour constrained residues in the current PyMOL session by
    connected-component palette.  Automatically falls back to
    residue-number-only selection when chain IDs are absent.

find_minimized_pdbs(job_dir)
    Discover stripped PDB files in a job's task directories.

generate_minimized_viewer(project, experiment, job)
    Generate a multi-snapshot Mol* viewer (.mvsj) for minimised
    structures, with DNA bases coloured by type. Each snapshot is
    titled "model N" (the short identifier derived from the task
    directory name).

Dependencies
------------
- PyMOL (with Python API, e.g. ``pymol-open-source``)
- PyYAML (``pyyaml``)
- Gemmi (``gemmi``)
- MolViewSpec (``molviewspec``)

Usage
-----
Import from a project script::

    from pymol_utils import process_final_min
    process_final_min("CSS", "Seq1_C0R_constrained", "J1234567")

Or use individual functions::

    from pymol_utils import load_constraints, color_constrained_residues

    constraints = load_constraints("CSS", "Seq1_C0R_constrained")
    if constraints:
        color_constrained_residues(constraints)

Generate a Mol* viewer::

    from pymol_utils import generate_minimized_viewer
    generate_minimized_viewer("CSS", "CSS1_free_constrained", "J1129505")
"""

import base64
import json
import os
import tempfile

import gemmi
import molviewspec as mvs
from molviewspec.nodes import States, GlobalMetadata
import yaml
# PyMOL is imported lazily inside functions that use it (process_final_min,
# color_constrained_residues) so that generate_minimized_viewer and other
# non-PyMOL utilities can be imported without a PyMOL installation.


# =====================================================================
#  Constraint palette
# =====================================================================

# Muted / pastel colour palette for constraint components (PyMOL names).
# Prioritises softer tones for a less aggressive, more aesthetic look.
CONSTRAINT_COLORS = [
    "palecyan", "palegreen", "paleyellow", "lightpink", "lightblue",
    "lightorange", "lightmagenta", "slate", "teal", "violet",
    "salmon", "lime", "skyblue", "wheat", "olive",
    "deepteal", "aquamarine", "raspberry", "darksalmon", "pink",
    "tan", "silver", "forest", "ruby", "hotpink",
]


# =====================================================================
#  Constraint loading and PyMOL colouring
#  (adapted from boltz2_utils/process_boltz_results.py)
# =====================================================================


def load_constraints(project, experiment):
    """
    Load unique contact constraints from the YAML file for a constrained
    experiment.

    Only loads constraints if *experiment* contains the ``"_constrained"``
    suffix.

    Parameters
    ----------
    project : str
        Project directory name (e.g. ``"CSS"``).
    experiment : str
        Experiment name (e.g. ``"Seq1_C0R_constrained"``).

    Returns
    -------
    list of tuple or None
        Each element is ``((chain_i, res_i), (chain_j, res_j))``, with
        canonical ordering so that ``(A,1)-(A,46)`` and ``(A,46)-(A,1)``
        are treated as the same constraint (deduplicated).
        Returns ``None`` if the experiment is not constrained or the YAML
        file cannot be found/loaded.
    """
    # Only attempt loading if the experiment name indicates constraints
    if "_constrained" not in experiment:
        print("Experiment name does not indicate constraints — skipping YAML loading.")
        return None

    # 1. Locate YAML file for this experiment
    yaml_path = f"{project}/yaml/{experiment}.yaml"
    if not os.path.exists(yaml_path):
        print(f"YAML file not found: {yaml_path}")
        return None

    print(f"Loading constraints from: {yaml_path}")
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    # 2. Extract raw constraint list
    constraints_raw = data.get("constraints", [])
    if not constraints_raw:
        print("No constraints section found in YAML file.")
        return None

    # 3. Deduplicate using frozenset so order doesn't matter,
    #    then convert to sorted tuples for predictable ordering
    unique_pairs = set()
    for entry in constraints_raw:
        contact = entry.get("contact", {})
        token1 = contact.get("token1")
        token2 = contact.get("token2")
        if token1 is None or token2 is None:
            continue
        if len(token1) != 2 or len(token2) != 2:
            continue
        pair = frozenset({(token1[0], token1[1]), (token2[0], token2[1])})
        unique_pairs.add(pair)

    constraints_list = [tuple(sorted(p)) for p in unique_pairs]
    print(f"Extracted {len(constraints_list)} unique constraint pairs.")
    return constraints_list


def color_constrained_residues(constraints):
    """
    Colour constrained residues in the current PyMOL session.

    Builds a graph from the constraint pairs, finds connected components
    (resolving ambiguity transitively), and assigns each component a
    unique colour from :data:`CONSTRAINT_COLORS`.

    Automatically detects whether the loaded PDB structures have chain IDs.
    If the chains referenced in the constraints do not exist in the session
    (e.g. PDB files with blank chain IDs), the chain filter is dropped and
    residues are selected by residue number only.

    Parameters
    ----------
    constraints : list of tuple
        List of ``((chain_i, res_i), (chain_j, res_j))`` pairs as returned
        by :func:`load_constraints`.
    """
    from pymol import cmd

    # 1. Build adjacency graph: each residue node -> neighbours
    graph = {}
    for (node1, node2) in constraints:
        graph.setdefault(node1, set()).add(node2)
        graph.setdefault(node2, set()).add(node1)

    if not graph:
        print("No residue nodes to colour.")
        return

    # 2. Detect whether the referenced chains actually exist in the session.
    #    PDB files from MD simulations often have blank chain IDs, while the
    #    YAML constraints reference chain "A".  In that case we drop the
    #    chain filter and select by residue number alone.
    referenced_chains = {ch for (ch, _) in graph.keys()}
    available_chains = set(cmd.get_chains())
    use_chain_filter = referenced_chains.issubset(available_chains)
    if not use_chain_filter:
        print(
            f"Referenced chain(s) {referenced_chains} not found in PDB "
            f"(available: {available_chains}). "
            "Falling back to residue-number-only selection."
        )

    # 3. Find connected components via BFS
    visited = set()
    components = []
    for node in graph:
        if node in visited:
            continue
        queue = [node]
        component = set()
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            component.add(current)
            for neighbour in graph.get(current, []):
                if neighbour not in visited:
                    queue.append(neighbour)
        components.append(component)

    print(f"Found {len(components)} connected component(s) among constrained residues.")

    # 4. Assign a colour from the palette to each component
    for idx, component in enumerate(components):
        color = CONSTRAINT_COLORS[idx % len(CONSTRAINT_COLORS)]

        # Group residues by chain for a compact selection string
        chain_groups = {}
        for (chain, resnum) in component:
            chain_groups.setdefault(chain, []).append(resnum)

        selections = []
        for chain, resids in chain_groups.items():
            resi_str = "+".join(str(r) for r in sorted(resids))
            if use_chain_filter:
                selections.append(f"(chain {chain} and resi {resi_str})")
            else:
                selections.append(f"(resi {resi_str})")

        # Apply colour
        sel_expr = " or ".join(selections)
        cmd.color(color, sel_expr)

        res_list = sorted(component)
        res_desc = "; ".join(f"{ch}:{r}" for (ch, r) in res_list)
        print(f"  Component {idx + 1} -> {color}: {res_desc}")


# =====================================================================
#  Structure loading and colouring
# =====================================================================


def process_final_min(project, experiment, job):
    """
    Load, align, and colour final-minimised structures in PyMOL.

    Scans the output directory for the given *project/experiment/job*
    triple, loads every ``*_stripped.pdb`` file found under each task's
    ``final_min/pdb/`` sub-directory, aligns all models, and applies
    a layered colour scheme.  Two session files are saved: one with
    base-coloured nucleic acids and one with white nucleic acids and
    constraint-coloured residues.

    Parameters
    ----------
    project : str
        Project directory name (e.g. ``"CSS"``).
    experiment : str
        Experiment name, e.g. ``"Seq1_C0R_constrained"``.
    job : str
        Job identifier (e.g. ``"J1234567"``).
    """
    from pymol import cmd, util

    final_min_path = f"{project}/MD/pmemd/out/{experiment}/{job}"

    # --- Discover task directories ------------------------------------
    tasks = [
        d for d in os.listdir(final_min_path)
        if os.path.isdir(os.path.join(final_min_path, d))
    ]
    print(f"Found {len(tasks)} tasks in {final_min_path}:\n{tasks}\n{'=' * 60}\n")

    # --- Initialise PyMOL and load structures ------------------------
    cmd.reinitialize()

    for task in tasks:
        task_path = os.path.join(final_min_path, task, "final_min/pdb")
        print(f"Processing task: {task} in path: {task_path}")

        if not os.path.exists(task_path):
            print(f"  Path does not exist — skipping")
            continue

        pdb_files = [
            f for f in os.listdir(task_path) if f.endswith("_stripped.pdb")
        ]
        print(f"  Found {len(pdb_files)} PDB file(s): {pdb_files}")

        for pdb_file in pdb_files:
            pdb_path = os.path.join(task_path, pdb_file)
            obj_name = f"{task}_{pdb_file[:-4]}"
            task_number = task.split("_")[-1]

            cmd.load(pdb_path, obj_name)
            # Rename to a short model identifier (e.g. "model_1")
            cmd.set_name(obj_name, f"model_{task_number}")
            print(f"  Loaded: {pdb_path} -> model_{task_number}")

    # --- Global PyMOL settings ---------------------------------------
    cmd.space("pymol")
    cmd.set("cartoon_ring_mode", 3)

    # --- Hide solvent and ions ---------------------------------------
    cmd.hide("everything", "resn WAT")
    cmd.hide("everything", "resn Na\\+")
    cmd.hide("everything", "resn Cl\\-")

    cmd.set("sphere_transparency", 0.3)

    # --- Align all loaded objects to the first one -------------------
    cmd.select("none")
    cmd.alignto()

    # --- Colour nucleic acids by base identity -----------------------
    # Covers standard DNA, RNA, and modified/deleted termini variants.
    cmd.color("wheat",     "resn DG+DG5+DG3+G+G3+G5+DI")   # Guanine
    cmd.color("palecyan",  "resn DA+DA5+DA3+A+A5+A3")        # Adenine
    cmd.color("palegreen", "resn DT+DT5+DT3+T+T5+T3")        # Thymine
    cmd.color("lime",      "resn DU+U+DU3+U3+DU5+U5")        # Uracil
    cmd.color("lightpink", "resn DC+DC5+DC3+C+C5+C3")        # Cytosine

    cmd.set("cartoon_nucleic_acid_color", "grey90")
    cmd.set("cartoon_discrete_colors", 1)

    # --- Colour protein in a single muted tone -----------------------
    # ``polymer.protein`` is a PyMOL built-in that selects all amino
    # acid residues regardless of chain.
    if cmd.count_atoms("polymer.protein") > 0:
        cmd.color("warmpink", "polymer.protein")

    # --- Colour small molecules by element ---------------------------
    util.cbay("br. organic")

    cmd.select("none")
    cmd.hide("everything", "hydro")
    cmd.set("sphere_scale", 0.5)

    # --- Save session 1: base-coloured nucleic acids -----------------
    session_path = os.path.join(final_min_path, f"{experiment}_final_min.pse")
    cmd.save(session_path)
    print(f"Saved PyMOL session to {session_path}")

    # --- Session 2 prep: white nucleic acids + constraints -----------
    # Reset nucleic acid colour to white before adding constraint
    # colouring so that only the constrained residues stand out.
    cmd.color("white", "polymer.nucleic")

    constraints = load_constraints(project, experiment)
    if constraints:
        print(f"Loaded {len(constraints)} unique constraint pairs.")
        color_constrained_residues(constraints)
    else:
        print("No constraints found — skipping constraint coloring.")

    # --- Save session 2: constraint-coloured -------------------------
    constraints_path = os.path.join(
        final_min_path, f"{experiment}_final_min_constraints.pse"
    )
    cmd.save(constraints_path)
    print(f"Saved constraint session to {constraints_path}")


# =====================================================================
#  DNA base colours (Mol* web viewer)
# =====================================================================

# Matches the PyMOL colour scheme used in process_final_min.
_BASE_COLORS = {
    "A": "#ADD8E6",  # lightblue
    "C": "#DDA0DD",  # plum
    "G": "#D2B48C",  # tan
    "T": "#90EE90",  # lightgreen
}

_NON_NUCLEIC_COLOR = "#E0E0E0"  # light grey for non-DNA residues

# Amber (and PDB) nucleic residue name → single-letter base.
_RESIDUE_BASE = {
    "DA": "A", "DA5": "A", "DA3": "A",
    "DC": "C", "DC5": "C", "DC3": "C",
    "DG": "G", "DG5": "G", "DG3": "G",
    "DT": "T", "DT5": "T", "DT3": "T",
    "DI": "G",
    "A": "A", "A5": "A", "A3": "A",
    "C": "C", "C5": "C", "C3": "C",
    "G": "G", "G5": "G", "G3": "G",
    "U": "T", "U5": "T", "U3": "T",
}


def _resolve_base(resname):
    """Return single-letter base code (A, C, G, T) or None."""
    return _RESIDUE_BASE.get(resname.strip())


# =====================================================================
#  MVSJ viewer for minimised structures
# =====================================================================


def find_minimized_pdbs(job_dir):
    """Return sorted list of stripped PDB paths under *job_dir*.

    Scans ``task_*/final_min/pdb/`` sub-directories for files ending
    with ``_stripped.pdb``.
    """
    paths = []
    for task in sorted(os.listdir(job_dir)):
        task_pdb = os.path.join(job_dir, task, "final_min", "pdb")
        if not os.path.isdir(task_pdb):
            continue
        for f in sorted(os.listdir(task_pdb)):
            if f.endswith("_stripped.pdb"):
                paths.append(os.path.join(task_pdb, f))
    return paths


# Terminal 5'/3' residue name mapping (Amber naming → standard gemmi names)
# so that gemmi classifies them as DNA polymer instead of unknown.
_TERMINAL_RENAME = {
    "DG5": "DG", "DG3": "DG",
    "DA5": "DA", "DA3": "DA",
    "DT5": "DT", "DT3": "DT",
    "DC5": "DC", "DC3": "DC",
}

# Calm uniform color for DNA cartoon
_DNA_COLOR = "#e1e9ef"

# Flashy color for ligand carbon atoms
_LIGAND_CARBON_COLOR = "#f5e507"

# Element colors (CPK-like) for ligand atoms — carbon excluded (uses _LIGAND_CARBON_COLOR)
_ELEMENT_COLORS = {
    "H":  "#FFFFFF",
    "O":  "#FF0000",
    "N":  "#3050F8",
    "P":  "#FFA500",
    "S":  "#FFD700",
    "F":  "#00FF00",
    "Cl": "#00FF00",
    "Br": "#A52A2A",
    "I":  "#800080",
    "Mg": "#20B2AA",
    "Na": "#B0B0FF",
}


def _write_minimized_mvsj(pdb_paths, output_path, verbose=True, title=None):
    """Core: build and write a multi-snapshot MVSJ from the given PDBs.

    Parameters
    ----------
    pdb_paths : list of str
        Paths to stripped PDB files, one per task.
    output_path : str
        Destination path for the .mvsj file.
    verbose : bool
    title : str | None
        Display title for the MVSJ metadata. Defaults to the output
        basename without extension.
    """
    snapshots = []

    for pdb_path in pdb_paths:
        structure = gemmi.read_structure(str(pdb_path))

        # Fix blank chain IDs so the annotation matches the embedded PDB
        for model in structure:
            for chain in model:
                if not chain.name.strip():
                    chain.name = "A"

        # Rename terminal residues so gemmi classifies them as DNA
        for model in structure:
            for chain in model:
                for residue in chain:
                    new_name = _TERMINAL_RENAME.get(residue.name)
                    if new_name:
                        residue.name = new_name

        # Remove ions (Na+, Cl-, etc.) — they are distracting in the viewer
        ion_names = {"NA", "CL", "Na+", "Cl-", "NA+", "CL-", "NA2+"}
        for model in structure:
            for chain in list(model):
                for i in range(len(chain) - 1, -1, -1):
                    if chain[i].name.strip() in ion_names:
                        del chain[i]
                if len(chain) == 0:
                    model.remove_chain(chain)

        # Build atom-level colour annotation for non-polymer residues
        # (ligands, MG, etc.): carbon atoms → flashy color, others → element color
        lig_annotation = []
        for model in structure:
            for chain in model:
                chain_id = chain.name.strip()
                for residue in chain:
                    info = gemmi.find_tabulated_residue(residue.name)
                    is_nucleic = (
                        info is not None
                        and info.kind in (gemmi.ResidueKind.DNA, gemmi.ResidueKind.RNA)
                    )
                    if is_nucleic:
                        continue
                    for atom in residue:
                        elem = atom.element.name
                        color = (
                            _LIGAND_CARBON_COLOR
                            if elem == "C"
                            else _ELEMENT_COLORS.get(elem, "#CCCCCC")
                        )
                        lig_annotation.append({
                            "label_asym_id": chain_id,
                            "label_seq_id": residue.seqid.num,
                            "label_atom_id": atom.name.strip(),
                            "color": color,
                        })

        if lig_annotation:
            ann_bytes = json.dumps(lig_annotation).encode()
            ann_b64 = base64.b64encode(ann_bytes).decode()
            ann_uri = f"data:application/json;base64,{ann_b64}"

        # Write fixed structure to temp PDB and embed as data URI
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            structure.write_pdb(tmp_path)
            with open(tmp_path, "rb") as fh:
                pdb_bytes = fh.read()
        finally:
            os.unlink(tmp_path)

        pdb_b64 = base64.b64encode(pdb_bytes).decode()
        pdb_uri = f"data:chemical/x-pdb;base64,{pdb_b64}"

        # --- Build MVSJ snapshot ---
        b = mvs.create_builder()
        b.canvas(custom={
            "molstar_postprocessing": {
                "enable_outline": True,
                "enable_ssao": False,
            },
        })
        ds = b.download(url=pdb_uri)
        ps = ds.parse(format="pdb")
        ms = ps.model_structure()

        # Polymer (DNA) → cartoon with uniform calm color
        poly = ms.component(selector="polymer")
        poly_rep = poly.representation(
            type="cartoon",
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        poly_rep.color(color=_DNA_COLOR)

        # Non-polymer (ligands, MG) → ball-and-stick with atom-level coloring
        lig = ms.component(selector="ligand")
        lig_rep = lig.representation(
            type="ball_and_stick",
            size_factor=0.5,
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        if lig_annotation:
            lig_rep.color_from_uri(
                uri=ann_uri,
                format="json",
                schema="atom",
                field_name="color",
            )

        # Snapshot title = short model identifier (e.g. "model 0")
        task_dir = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(pdb_path))))
        snap = b.get_snapshot(
            title=f"model {task_dir.split('_')[-1]}",
            linger_duration_ms=3000,
            transition_duration_ms=500,
        )
        snapshots.append(snap)

    mvsj_basename = os.path.basename(output_path)
    states = States(
        metadata=GlobalMetadata(title=title or os.path.splitext(mvsj_basename)[0]),
        snapshots=snapshots,
    )
    states_dict = states.model_dump(exclude_none=True)
    with open(output_path, "w") as fh:
        json.dump(states_dict, fh, indent=2)
    if verbose:
        print(f"  Wrote MD minimised MVSJ: {output_path}")


def generate_minimized_viewer(project, experiment, job, verbose=True):
    """Generate a multi-snapshot .mvsj file from Amber MD minimised PDBs.

    Colours DNA bases by type matching the PyMOL session colours:
    - Adenine → lightblue
    - Cytosine → plum
    - Guanine → tan
    - Thymine → lightgreen

    Non-nucleic residues (ions) are shown in light grey.

    The MVSJ file is written to the job directory with the name::

        {experiment}_{job}_final_min.mvsj

    Parameters
    ----------
    project : str
        Project directory name (e.g. ``"CSS"``, ``"PQ4"``).
    experiment : str
        Experiment name (e.g. ``"CSS1_free_constrained"``).
    job : str
        Job identifier (e.g. ``"J1129505"``).
    verbose : bool, default True
        If False, suppresses status prints.
    """
    job_dir = os.path.join(project, "MD", "pmemd", "out", experiment, job)
    pdb_paths = find_minimized_pdbs(job_dir)

    if not pdb_paths:
        if verbose:
            print(f"  No stripped PDBs found in {job_dir}")
        return

    output_path = os.path.join(job_dir, f"{experiment}_{job}_final_min.mvsj")

    # Neat display title (e.g. "CSS1 · HCY bound") instead of the raw
    # {experiment}_{job}_final_min directory-style name
    if "free" in experiment:
        label = "Free"
    elif "HCY" in experiment:
        label = "HCY bound"
    elif "PQ" in experiment:
        label = "PQ bound"
    else:
        label = experiment
    mvs_title = f"{experiment.split('_')[0]} \u00B7 {label}"

    _write_minimized_mvsj(pdb_paths, output_path, verbose, title=mvs_title)


# =====================================================================
#  CLI entry point
# =====================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate Mol* viewer (.mvsj) for Amber MD minimised structures.",
    )
    parser.add_argument("project", help="Project directory name (e.g. CSS, PQ4)")
    parser.add_argument("experiment", help="Experiment name (e.g. CSS1_free_constrained)")
    parser.add_argument("job", help="Job identifier (e.g. J1129505)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress status output")

    args = parser.parse_args()
    generate_minimized_viewer(args.project, args.experiment, args.job, verbose=not args.quiet)
