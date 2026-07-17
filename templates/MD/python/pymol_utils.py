"""
PyMOL visualisation utilities for MD minimisation outputs.

Provides reusable functions for loading, aligning, and colouring
structures produced by the Amber MD pipeline (§9).  Designed to be
copied into each project via ``cp -r templates/MD/*`` and imported
by project-specific entry-point scripts.

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

Dependencies
------------
- PyMOL (with Python API, e.g. ``pymol-open-source``)
- PyYAML (``pyyaml``)

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
"""

import os

import yaml
from pymol import cmd, util


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
