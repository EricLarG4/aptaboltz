"""
Boltz-2 result post-processing pipeline.

Loads Boltz-2 prediction outputs (CIF, JSON, NPZ) for one or more
experimental runs and produces:

- Aligned, coloured PyMOL sessions (plain / constraint-coloured / pLDDT-coloured)
- Extracted ligand SDF files and chirality-annotated grid images (PNG)
- Aggregated confidence metrics, PAE, PDE, and pLDDT tables (CSV)

Entry point for batch processing: process_single_experiment()
Entry point for a single job directory: process_experiment()
"""

import yaml
from natsort import natsort
from pymol import cmd
import gemmi
import os
import numpy as np
import pandas as pd
import json
from rdkit import Chem as chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D as chemdraw
from rdkit.Chem.Draw import MolsToGridImage as chemdrawgrid

try:
    rdDepictor.SetPreferCoordGen(True)
    HAVE_COORDGEN = True
except ImportError:
    HAVE_COORDGEN = False


# ── Package-relative path to base.py (PyMOL styling script) ────────────
# Resolved at import time so it works regardless of how boltz2_utils is
# installed (editable or regular pip install).
_BASE_SCRIPT = os.path.join(os.path.dirname(__file__), "base.py")


# =====================================================================
#  Single-experiment prediction processing
# =====================================================================


def process_experiment(
    project,
    model,
    experiment,
    LGD_DICT=None,
    NAME_MAP=None,
    ligand_processing=True,
    generate_mvsj=True,
):
    """
    Load Boltz-2 CIF predictions into PyMOL, align structures, apply
    colouring (DNA, RNA, protein, ligand, constraints, pLDDT), and
    extract ligand conformers for chirality analysis.

    Parameters
    ----------
    project : str
        Project directory name.
    model : str
        Experiment name, e.g. "Seq1_C0R_constrained".
    experiment : str
        Job directory name (e.g. "J1121077").
    LGD_DICT : dict, optional
        Maps aptamer names -> CCD codes.  Must be provided when
        *ligand_processing* is True.
    NAME_MAP : dict, optional
        Maps CCD codes -> display names.  Must be provided when
        *ligand_processing* is True for experiments with a ligand.
    ligand_processing : bool
        If True, extract ligand SDF, fetch reference, and draw chirality grid.
    """
    print(f"Processing experiment: {project}/{model}/{experiment}\n")

    # --- File paths ---
    cif_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{cif_path}"

    # 1. List and sort CIF files (best model first via natural sort)
    cif_files = [f for f in os.listdir(full_path) if f.endswith(".cif")]
    cif_files = natsort.natsorted(cif_files)
    print(f"Found CIF files: {cif_files}\n\n")

    # 2. Initialise PyMOL and load all CIF models
    cmd.reinitialize()
    for cif_file in cif_files:
        print(f"Loading CIF file: {full_path + cif_file}")
        cmd.load(full_path + cif_file)

    # 3. Apply base styling (cartoon colours, solvent/ion hiding),
    #    align loaded models, and save as multi-object CIF + PSE
    cmd.load(_BASE_SCRIPT)
    cmd.save(f"{experiment_path}{model}_{experiment}.cif", state=0)
    print(f"Saved multi-object file: {experiment_path}{model}_{experiment}.cif\n")
    cmd.deselect()
    cmd.save(f"{experiment_path}{model}_{experiment}.pse")
    print(f"Saved PyMOL session: {experiment_path}{model}_{experiment}.pse\n")

    # 4. Colour DNA chain(s) white for legibility
    cmd.color("white", "chain A")
    if cmd.count_atoms("chain B") > 0:
        cmd.color("white", "chain B")
        cmd.util.cbaw("chain B")
    else:
        # Fallback: colour any non-A chain white
        all_chains = cmd.get_chains()
        for ch in all_chains:
            if ch != "A":
                cmd.color("white", f"chain {ch}")
                cmd.util.cbaw(f"chain {ch}")

    # 5. If constraints YAML exists, load and colour constrained residues
    #    with a distinct palette (pastel colours)
    constraints = load_constraints(project, model, experiment)
    if constraints:
        print(f"Loaded {len(constraints)} unique constraint pairs.")
        color_constrained_residues(constraints)
        cmd.save(f"{experiment_path}{model}_{experiment}_constraints.pse")
    else:
        print("No constraints found — skipping constraint coloring.")

    # 6. Colour by pLDDT (confidence) using a red→blue spectrum,
    #    then save as a separate session file
    cmd.spectrum(
        "b", "tv_red yelloworange palecyan density", "all", minimum=50, maximum=90
    )
    cmd.save(f"{experiment_path}{model}_{experiment}_plddt.pse")
    print(
        f"Saved pLDDT colored PyMOL session: {experiment_path}{model}_{experiment}_plddt.pse\n"
    )

    # 7. Mol* web viewer generation (MVSJ files for quarto-molstar)
    if generate_mvsj:
        from boltz2_utils.generate_web_viewer import (
            generate_plddt_viewer,
            generate_constraint_viewer,
        )

        mvsj_base = f"{experiment_path}{model}_{experiment}"

        # pLDDT viewer (per-model MVSJs with individual CIFs)
        generate_plddt_viewer(
            f"{experiment_path}{model}_{experiment}.cif",
            experiment_path,
        )

        # Constraint viewer (only for constrained experiments)
        if constraints:
            yaml_path = f"{project}/yaml/{model}.yaml"
            generate_constraint_viewer(
                f"{experiment_path}{model}_{experiment}.cif",
                yaml_path,
                f"{mvsj_base}_constraints.mvsj",
            )

        print()

    # 8. Ligand extraction and chirality visualisation
    if ligand_processing:
        # 7a. Determine aptamer and ligand from model name
        apt = model.split("_")[0]
        lgd = model.split("_")[1]

        # Skip if no ligand present
        if lgd.lower() == "free":
            print(f"No ligand to extract for aptamer {apt} (ligand: {lgd})")
            return

        name = NAME_MAP[lgd]

        # 7b. Select ligand atoms by residue name (CCD code); fallback to LIG1
        cmd.select(f"resn {lgd}")
        print(f"Number of selected atoms: {cmd.count_atoms('sele')}")
        if cmd.count_atoms("sele") == 0:
            print(
                f"No atoms found for ligand {lgd} in the structure. Switching to LIG1 selection."
            )
            lgd_selector = "resn LIG1"
        else:
            lgd_selector = f"resn {lgd}"
        print(f"Extracting {lgd_selector} ({name}) for aptamer {apt}")
        lgd_path = f"{experiment_path}{name}.sdf"
        print(f"Extracting ligand: {name}, Path: {lgd_path}")
        cmd.save(lgd_path, selection=lgd_selector, state=1)
        cmd.reinitialize()

        # 7c. Fetch reference ligand from PDB (by CCD code) as comparison
        cmd.fetch(lgd)
        if cmd.get_object_list() == []:
            print(
                f"{lgd} is not a valid PDB ligand code. Skipping reference ligand extraction."
            )
        else:
            cmd.save(f"{experiment_path}{lgd}.sdf")
        cmd.reinitialize()

        # 7d. Load all SDFs into RDKit mols; sanitise, remove Hs,
        #     assign stereochemistry, and straighten depiction
        mols = {}
        sdf_files = [f for f in os.listdir(experiment_path) if f.endswith(".sdf")]
        print(f"sdf files found: {sdf_files}")
        for file in sdf_files:
            mol_name = file.split(".")[0]
            mol = chem.SDMolSupplier(
                f"{experiment_path}{file}", removeHs=False, sanitize=True
            )[0]
            if mol is None:
                # Retry without sanitisation if the first attempt failed
                mol = chem.SDMolSupplier(
                    f"{experiment_path}{file}", removeHs=False, sanitize=False
                )[0]
            else:
                mol = chem.RemoveHs(mol)
            chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            rdDepictor.StraightenDepiction(mol)
            mol.RemoveAllConformers()
            mols[mol_name] = mol

        mols = dict(sorted(mols.items()))

        # 7e. Identify chiral centres in each conformer and reference;
        #     discard centres marked as unassigned ("?")
        chiral_centers = []
        for mol in mols.values():
            chiral_center = chem.FindMolChiralCenters(mol, includeUnassigned=True)
            chiral_center = [
                (idx, stereo) for idx, stereo in chiral_center if stereo != "?"
            ]
            chiral_index = [idx for idx, _ in chiral_center]
            chiral_centers.append(chiral_index)

        # 7f. Draw a grid image with molecule names as legends and
        #     chiral centres highlighted in gold
        template_drawer = chemdraw.MolDraw2DCairo(500, 500)
        grid_opts = template_drawer.drawOptions()
        grid_opts.addStereoAnnotation = True
        grid_opts.highlightColour = (1, 0.843, 0, 0.75)  # gold, ~transparent
        grid_opts.highlightBondWidthMultiplier = True
        grid_opts.annotationFontScale = 1
        grid_opts.atomNoteColour = (1, 0.678, 0, 1)  # pure orange, opaque
        grid_opts.bondLineWidth = 3.5
        grid_opts.legendFontSize = 30
        grid_opts.singleColourWedgeBonds = True
        grid_opts.unspecifiedStereoIsUnknown = False
        grid_opts.additionalAtomLabelPadding = 0.05
        grid_image = chemdrawgrid(
            mols=list(mols.values()),
            molsPerRow=len(mols),
            subImgSize=(500, 500),
            legends=list(mols.keys()),
            highlightAtomLists=chiral_centers,
            drawOptions=grid_opts,
        )
        grid_image.save(f"{experiment_path}{name}.png")
        print(f"Wrote {name}.png in {experiment_path}")


# =====================================================================
#  Constraint loading and PyMOL colouring
# =====================================================================


# Muted / pastel colour palette for constraint components (PyMOL colour names)
# Prioritises softer tones for a less aggressive, more aesthetic look.
CONSTRAINT_COLORS = [
    "palecyan", "palegreen", "paleyellow", "lightpink", "lightblue",
    "lightorange", "lightmagenta", "slate", "teal", "violet",
    "salmon", "lime", "skyblue", "wheat", "olive",
    "deepteal", "aquamarine", "raspberry", "darksalmon", "pink",
    "tan", "silver", "forest", "ruby", "hotpink",
]


def load_constraints(project, model, experiment):
    """
    Load unique contact constraints from the YAML file for a constrained
    experiment.

    Only loads constraints if *model* contains the ``"_constrained"`` suffix
    (the naming convention used by :func:`~.generate_boltz2_yaml.generate_yaml_file`).

    Returns
    -------
    list of tuple
        Each element is ``((chain_i, res_i), (chain_j, res_j))``, with
        canonical ordering so that ``(A,1)-(A,46)`` and ``(A,46)-(A,1)``
        are treated as the same constraint (deduplicated).
        Returns ``None`` if the experiment is not constrained or the YAML
        file cannot be found/loaded.
    """
    # Only attempt loading if the model name indicates constraints
    if "_constrained" not in model:
        print("Model name does not indicate constraints — skipping YAML loading.")
        return None

    # 1. Locate YAML file for this model
    yaml_path = f"{project}/yaml/{model}.yaml"
    if not os.path.exists(yaml_path):
        print(f"YAML file not found: {yaml_path}")
        return None

    print(f"Loading constraints from: {yaml_path}")
    constraints = _parse_constraint_pairs(yaml_path)
    if not constraints:
        print("No constraints section found in YAML file.")
        return None

    # 3. Deduplicate using frozenset so order doesn't matter,
    #    then convert to sorted tuples for predictable ordering
    constraints_list = list(constraints.keys())
    print(f"Extracted {len(constraints_list)} unique constraint pairs.")
    return constraints_list


def _parse_constraint_pairs(yaml_path):
    """
    Parse unique contact constraints from a Boltz-2 YAML file.

    Returns
    -------
    dict
        Mapping from canonical ``((chain_i, res_i), (chain_j, res_j))``
        tuples (order-insensitive, deduplicated) to the restraint's
        ``max_distance``.
    """
    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    unique_pairs = {}
    for entry in data.get("constraints", []):
        contact = entry.get("contact", {})
        token1 = contact.get("token1")
        token2 = contact.get("token2")
        if token1 is None or token2 is None:
            continue
        if len(token1) != 2 or len(token2) != 2:
            continue
        pair = tuple(sorted({(token1[0], token1[1]), (token2[0], token2[1])}))
        unique_pairs[pair] = float(contact.get("max_distance", 4.0))

    return unique_pairs


def color_constrained_residues(constraints):
    """
    Colour constrained residues in the current PyMOL session.

    Builds a graph from the constraint pairs, finds connected components
    (resolving ambiguity transitively), and assigns each component a
    unique colour from :data:`CONSTRAINT_COLORS`.

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

    # 2. Find connected components via BFS
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

    # 3. Assign a colour from the palette to each component
    for idx, component in enumerate(components):
        color = CONSTRAINT_COLORS[idx % len(CONSTRAINT_COLORS)]

        # Group residues by chain for a compact selection string
        chain_groups = {}
        for (chain, resnum) in component:
            chain_groups.setdefault(chain, []).append(resnum)

        selections = []
        for chain, resids in chain_groups.items():
            resi_str = "+".join(str(r) for r in sorted(resids))
            selections.append(f"(chain {chain} and resi {resi_str})")

        # Apply colour (restricted to nucleic/DNA atoms to avoid ligand atoms)
        sel_expr = " or ".join(selections)
        cmd.color(color, sel_expr)

        res_list = sorted(component)
        res_desc = "; ".join(f"{ch}:{r}" for (ch, r) in res_list)
        print(f"  Component {idx + 1} -> {color}: {res_desc}")


def _model_heavy_atoms(block):
    """
    Extract per-residue heavy-atom coordinates from one mmCIF model block.

    Returns
    -------
    dict
        Mapping from ``(chain, res_seq)`` to an (n, 3) float array of heavy
        atom coordinates. Chain is the auth chain id and res_seq the label
        sequence number (1-based), matching the residue tokens in the
        constraint YAML files.
    """
    atoms = {}
    table = block.find_mmcif_category("_atom_site.")
    for row in table:
        if row.str(2) in ("H", "D"):
            continue
        chain = row.str(16)
        try:
            seq = int(row.str(8))
        except (TypeError, ValueError):
            continue
        key = (chain, seq)
        atoms.setdefault(key, []).append(
            (float(row.str(10)), float(row.str(11)), float(row.str(12)))
        )
    return {k: np.asarray(v) for k, v in atoms.items()}


def _min_heavy_atom_distance(atom_arrays, pair):
    """
    Minimum Euclidean distance (Å) between two residues' heavy atoms.

    Parameters
    ----------
    atom_arrays : dict
        Output of :func:`_model_heavy_atoms`.
    pair : tuple
        ``((chain_i, res_i), (chain_j, res_j))``.
    """
    (chain_i, res_i), (chain_j, res_j) = pair
    a_i = atom_arrays.get((chain_i, res_i))
    a_j = atom_arrays.get((chain_j, res_j))
    if a_i is None or a_j is None or len(a_i) == 0 or len(a_j) == 0:
        return float("inf")
    diff = a_i[:, None, :] - a_j[None, :, :]
    return float(np.linalg.norm(diff, axis=-1).min())


def _pair_label(pair):
    """
    Compact human-readable label for a constraint pair, e.g. ``A14-A24``.
    """
    (chain_i, res_i), (chain_j, res_j) = pair
    return f"{chain_i}{res_i}-{chain_j}{res_j}"


def process_constraint_verification(project, model, experiment, verbose=True):
    """
    Compute per-model satisfaction of contact restraints for one Boltz-2 job.

    For each unique contact restraint ``(i, j)`` extracted from the
    experiment's YAML file and for every model in the job's CIF file, the
    minimum heavy-atom distance between residues i and j is compared against
    the restraint's ``max_distance``.

    Verification is ambiguity-aware: restraint sets may list a residue with
    more than one possible partner (e.g. both ``A16-A22`` and ``A16-A23``),
    but a base-paired residue can only contact one partner at a time, so not
    all listed restraints can be satisfied simultaneously. A model is
    therefore marked as fully verified when **every constrained residue takes
    part in at least one satisfied restraint** (i.e. no constrained residue
    is left without any satisfied contact). Unsatisfied restraints whose
    endpoints are both covered by other satisfied restraints (explained
    ambiguity, e.g. ``A16-A23`` when T16 pairs A22) are not treated as
    failures.

    Only experiments whose model name contains the ``"_constrained"`` suffix
    are processed (the naming convention used by
    :func:`~.generate_boltz2_yaml.generate_yaml_file`).

    Results are written to
    ``{project}/{model}/{experiment}/{model}_{experiment}_constraints.csv``
    with one row per model:

    - model: 0-based model index (matches the sample numbering used for MD
      forwarding)
    - n_constraints / n_satisfied: unique restraint counts
    - all_verified: 1 if every constrained residue is covered by at least
      one satisfied restraint, else 0
    - failed_constraints: semicolon-separated labels of restraints that are
      unsatisfied while at least one of their endpoints is not covered
      ("—" when there are none)
    - worst_min_distance: largest min-distance observed across restraints

    Parameters
    ----------
    project : str
    model : str
    experiment : str
    verbose : bool
        Print progress messages to stdout (default True).
    """
    if "_constrained" not in model:
        if verbose:
            print("Model name does not indicate constraints — skipping constraint verification.")
        return None

    yaml_path = f"{project}/yaml/{model}.yaml"
    if not os.path.exists(yaml_path):
        if verbose:
            print(f"YAML file not found: {yaml_path}")
        return None
    constraints = _parse_constraint_pairs(yaml_path)
    if not constraints:
        if verbose:
            print("No contact constraints found in YAML file.")
        return None

    experiment_path = f"{project}/{model}/{experiment}/"
    cif_files = [f for f in os.listdir(experiment_path) if f.endswith(".cif")]
    if not cif_files:
        if verbose:
            print(f"No CIF files found in {experiment_path}")
        return None
    cif_path = os.path.join(experiment_path, cif_files[0])
    if verbose:
        print(f"Reading models from: {cif_path}")

    # Residues -> incident constraints, for ambiguity-aware coverage
    residue_constraints = {}
    for pair in constraints:
        residue_constraints.setdefault(pair[0], []).append(pair)
        residue_constraints.setdefault(pair[1], []).append(pair)

    doc = gemmi.cif.read(cif_path)
    rows = []
    for model_idx, block in enumerate(doc):
        atom_arrays = _model_heavy_atoms(block)
        satisfied = set()
        distances = {}
        for pair, max_distance in constraints.items():
            dist = _min_heavy_atom_distance(atom_arrays, pair)
            distances[pair] = dist
            if dist <= max_distance:
                satisfied.add(pair)

        # A residue is covered if at least one of its restraints is satisfied
        covered = {
            res
            for res, pairs in residue_constraints.items()
            if any(p in satisfied for p in pairs)
        }
        all_verified = len(covered) == len(residue_constraints)

        # Report only unsatisfied restraints with an uncovered endpoint
        failed = [
            pair
            for pair in constraints
            if pair not in satisfied
            and (pair[0] not in covered or pair[1] not in covered)
        ]
        failed.sort(key=_pair_label)

        rows.append(
            {
                "model": model_idx,
                "n_constraints": len(constraints),
                "n_satisfied": len(satisfied),
                "all_verified": int(all_verified),
                "failed_constraints": "; ".join(_pair_label(p) for p in failed)
                if failed
                else "—",
                "worst_min_distance": round(max(distances.values()), 3),
            }
        )

    result_df = pd.DataFrame(rows)
    out_path = f"{experiment_path}{model}_{experiment}_constraints.csv"
    result_df.to_csv(out_path, index=False)

    if verbose:
        n_ok = int(result_df["all_verified"].sum())
        print(
            f"Constraint verification for {project}/{model}/{experiment}: "
            f"{n_ok}/{len(result_df)} models fully verified. "
            f"Saved: {out_path}"
        )
    return out_path


# =====================================================================
#  Confidence, PAE, PDE, and pLDDT CSV aggregation
# =====================================================================


def process_confidence(project, model, experiment):
    """
    Aggregate confidence JSON files from one Boltz-2 job into a single CSV.

    JSON dict values (e.g. ``pair_chains_iptm``) that contain sub-dicts
    like ``{'0': val, '1': val}`` are kept as-is — the R-side processor
    handles expansion.

    Parameters
    ----------
    project : str
    model : str
    experiment : str
    """
    print(f"Processing confidence for experiment: {project}/{model}/{experiment}")

    # 1. Locate and sort JSON files
    json_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{json_path}"
    json_files = [f for f in os.listdir(full_path) if f.endswith(".json")]
    json_files = natsort.natsorted(json_files)
    print(f"Found JSON files: {json_files}")

    # 2. Load each JSON as a single-row DataFrame (wrapping in list to
    #    avoid pd.read_json expanding nested dicts into multiple rows)
    dataframes = []
    for json_file in json_files:
        json_file_path = os.path.join(full_path, json_file)
        with open(json_file_path) as f:
            data = json.load(f)
        df = pd.DataFrame([data])

        # Flatten dict-valued columns like pair_chains_iptm
        for col in ["pair_chains_iptm", "pair_chains_pae"]:
            if col in df.columns:
                df[col] = df[col].apply(
                    lambda x: list(x.values())[0] if isinstance(x, dict) else x
                )
        dataframes.append(df)

    # 3. Concatenate and save
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(
            f"{experiment_path}{model}_{experiment}_confidence.csv", index=False
        )
        print(f"Saved combined CSV file: {experiment_path}{model}_{experiment}.csv")


def process_pae(project, model, experiment):
    """
    Aggregate PAE (predicted aligned error) NPZ matrices into one CSV.

    Parameters
    ----------
    project : str
    model : str
    experiment : str
    """
    print(f"Processing PAE for experiment: {project}/{model}/{experiment}")

    # 1. Locate PAE NPZ files
    pae_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{pae_path}"
    pae_files = [
        f
        for f in os.listdir(full_path)
        if f.endswith(".npz") and f.startswith("pae_input_model_")
    ]
    pae_files = natsort.natsorted(pae_files)
    print(f"Found PAE files: {pae_files}")

    # 2. Load .npz, extract 'pae' array, append model number as column
    dataframes = []
    for pae_file in pae_files:
        pae_file_path = os.path.join(full_path, pae_file)
        with np.load(pae_file_path) as data:
            pae_matrix = data["pae"]
            df = pd.DataFrame(pae_matrix)
            df["model"] = pae_file.split("_")[-1].split(".")[0]
            dataframes.append(df)

    # 3. Concatenate and save
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(
            f"{experiment_path}{model}_{experiment}_pae.csv", index=False
        )
        print(
            f"Saved combined PAE CSV file: {experiment_path}{model}_{experiment}_pae.csv"
        )


def process_pde(project, model, experiment):
    """
    Aggregate PDE (predicted distance error) NPZ matrices into one CSV.

    Parameters
    ----------
    project : str
    model : str
    experiment : str
    """
    print(f"Processing PDE for experiment: {project}/{model}/{experiment}")

    # 1. Locate PDE NPZ files
    pde_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{pde_path}"
    pde_files = [
        f
        for f in os.listdir(full_path)
        if f.endswith(".npz") and f.startswith("pde_input_model_")
    ]
    pde_files = natsort.natsorted(pde_files)
    print(f"Found PDE files: {pde_files}")

    # 2. Load .npz, extract 'pde' array, append model number
    dataframes = []
    for pde_file in pde_files:
        pde_file_path = os.path.join(full_path, pde_file)
        with np.load(pde_file_path) as data:
            pde_matrix = data["pde"]
            df = pd.DataFrame(pde_matrix)
            df["model"] = pde_file.split("_")[-1].split(".")[0]
            dataframes.append(df)

    # 3. Concatenate and save
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(
            f"{experiment_path}{model}_{experiment}_pde.csv", index=False
        )
        print(
            f"Saved combined PDE CSV file: {experiment_path}{model}_{experiment}_pde.csv"
        )


def process_plddt(project, model, experiment):
    """
    Aggregate pLDDT (predicted local distance difference test) NPZ arrays
    into one CSV with one row per residue × model.

    Parameters
    ----------
    project : str
    model : str
    experiment : str
    """
    print(f"Processing pLDDT for experiment: {project}/{model}/{experiment}\n\n")

    # 1. Locate pLDDT NPZ files
    plddt_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{plddt_path}"
    plddt_files = [
        f
        for f in os.listdir(full_path)
        if f.endswith(".npz") and f.startswith("plddt_input_model_")
    ]
    plddt_files = natsort.natsorted(plddt_files)
    print(f"Found pLDDT files: {plddt_files}\n\n")

    # 2. Load .npz, extract 'plddt' array, append model number
    dataframes = []
    for plddt_file in plddt_files:
        plddt_file_path = os.path.join(full_path, plddt_file)
        with np.load(plddt_file_path) as data:
            plddt_array = data["plddt"]
            df = pd.DataFrame(plddt_array, columns=["pLDDT"])
            df["model"] = plddt_file.split("_")[-1].split(".")[0]
            dataframes.append(df)

    # 3. Concatenate and save
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(
            f"{experiment_path}{model}_{experiment}_plddt.csv", index=False
        )
        print(
            f"Saved combined pLDDT CSV file: {experiment_path}{model}_{experiment}_plddt.csv"
        )


# =====================================================================
#  Batch experiment processing (top-level entry point)
# =====================================================================


def process_single_experiment(
    project,
    experiment_name,
    LGD_DICT=None,
    NAME_MAP=None,
    ligand_processing=True,
    job_dirs=None,
    generate_mvsj=True,
):
    """
    Find J-directories for an experiment and delegate to process_boltz_results.

    Parameters
    ----------
    project : str
        Project name.
    experiment_name : str
        Experiment name (e.g. "Seq1_C0R_constrained").
    LGD_DICT : dict, optional
        Dictionary mapping aptamer names to ligand codes.  Must be
        provided when *ligand_processing* is True.
    NAME_MAP : dict, optional
        Dictionary mapping ligand codes to display names.  Must be
        provided when *ligand_processing* is True for experiments
        with a ligand.
    ligand_processing : bool, default=True
        Whether to extract and process ligand data.
    job_dirs : None, "last", str, or list of str, optional
        Controls which J-directories to process:

        - ``None`` (default): process **all** J-directories found.
        - ``"last"``: sort J-directories naturally and process only the
          highest-numbered (last) one.
        - A single string, e.g. ``"J1121077"``: process only that directory.
        - A list of strings, e.g. ``["J1121077", "J1123593"]``: process only
          those directories.
    """
    experiment_path = f"{project}/{experiment_name}/"

    # 1. Check experiment exists
    if not os.path.exists(experiment_path):
        print(f"Experiment does not exist: {experiment_path} | skipping")
        return

    # 2. Find all job (J*) directories within the experiment
    directories = [d for d in os.listdir(experiment_path) if d.startswith("J")]
    if not directories:
        print(f"No data found for {experiment_name} | skipping")
        return

    directories = natsort.natsorted(directories)
    print(f"Found directories for {experiment_name}: {directories}")

    # 3. Select which directories to process based on job_dirs argument
    if job_dirs is None:
        selected_dirs = directories
    elif job_dirs == "last":
        selected_dirs = [directories[-1]]
    elif isinstance(job_dirs, str):
        selected_dirs = [job_dirs]
    elif isinstance(job_dirs, list):
        selected_dirs = job_dirs
    else:
        print(f"Invalid job_dirs argument: {job_dirs} | falling back to all directories")
        selected_dirs = directories

    # 4. Validate selections actually exist on disk
    valid_dirs = [d for d in selected_dirs if d in directories]
    if not valid_dirs:
        print(f"None of the requested job directories exist in {directories} | skipping")
        return
    if len(valid_dirs) < len(selected_dirs):
        missing = set(selected_dirs) - set(directories)
        print(f"Warning: requested job directories not found: {sorted(missing)} | skipping those")

    # 5. Process each selected job directory
    for job_dir in valid_dirs:
        print(f"\n{'='*60}")
        print(f"Processing job directory: {job_dir}")
        print(f"{'='*60}\n")

        process_experiment(
            project,
            experiment_name,
            job_dir,
            LGD_DICT=LGD_DICT,
            NAME_MAP=NAME_MAP,
            ligand_processing=ligand_processing,
            generate_mvsj=generate_mvsj,
        )

        # 6. If prediction JSON directory exists, aggregate confidence/PAE/PDE/pLDDT
        json_path = os.path.join(
            experiment_path,
            job_dir,
            "boltz_results_input/predictions/input/",
        )
        if os.path.exists(json_path):
            process_confidence(project, experiment_name, job_dir)
            process_pae(project, experiment_name, job_dir)
            process_pde(project, experiment_name, job_dir)
            process_plddt(project, experiment_name, job_dir)
        else:
            print(f"No predictions input directory found in {job_dir}; "
                  f"skipping confidence/PAE/PDE/pLDDT processing.")