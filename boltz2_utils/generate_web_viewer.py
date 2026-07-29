"""Mol* Web Viewer generation for Boltz-2 prediction outputs.

Produces MolViewSpec JSON (.mvsj) files for interactive 3D visualisation
in the Mol* web viewer (https://molstar.org/viewer/), providing an
alternative to PyMOL-based visualisation.
"""

import json
import os
import base64
import tempfile

import gemmi
import yaml
import molviewspec as mvs
from molviewspec.nodes import States, GlobalMetadata, CategoricalPalette


# Pastel colours for constraint components (hex approximation of PyMOL names)
CONSTRAINT_COLORS = [
    "#AFEEEE",  # palecyan
    "#98FB98",  # palegreen
    "#FFFF99",  # paleyellow
    "#FFB6C1",  # lightpink
    "#ADD8E6",  # lightblue
    "#FFD27F",  # lightorange
    "#FFB3FF",  # lightmagenta
    "#C0C0C0",  # slate
    "#008080",  # teal
    "#EE82EE",  # violet
    "#FA8072",  # salmon
    "#00FF00",  # lime
    "#87CEEB",  # skyblue
    "#F5DEB3",  # wheat
    "#808000",  # olive
    "#008B8B",  # deepteal
    "#7FFFD4",  # aquamarine
    "#E30B5C",  # raspberry
    "#E9967A",  # darksalmon
    "#FFC0CB",  # pink
    "#D2B48C",  # tan
    "#C0C0C0",  # silver
    "#228B22",  # forest
    "#FF2400",  # ruby
    "#FF69B4",  # hotpink
]


# pLDDT gradient matching PyMOL "tv_red yelloworange palecyan density"
_PLDDT_GRADIENT = [
    (50.0, (1.0, 0.0, 0.0)),        # tv_red
    (63.3, (1.0, 0.647, 0.0)),      # yelloworange
    (76.6, (0.686, 0.933, 0.933)),  # palecyan
    (90.0, (0.878, 1.0, 1.0)),      # density
]

_PLDDT_MIN = 50.0
_PLDDT_MAX = 90.0


# =====================================================================
#  pLDDT viewer
# =====================================================================


def _plddt_to_color(value):
    """Map a single pLDDT score (0-100) to a hex colour via the gradient.

    Values outside [50, 90] are clamped to the nearest stop.
    """
    x = max(_PLDDT_MIN, min(value, _PLDDT_MAX))
    for i in range(len(_PLDDT_GRADIENT) - 1):
        x0, (r0, g0, b0) = _PLDDT_GRADIENT[i]
        x1, (r1, g1, b1) = _PLDDT_GRADIENT[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
            r = int((r0 + t * (r1 - r0)) * 255)
            g = int((g0 + t * (g1 - g0)) * 255)
            b = int((b0 + t * (b1 - b0)) * 255)
            return f"#{r:02x}{g:02x}{b:02x}"
    return f"#{int(_PLDDT_GRADIENT[-1][1][0]*255):02x}" \
            f"{int(_PLDDT_GRADIENT[-1][1][1]*255):02x}" \
            f"{int(_PLDDT_GRADIENT[-1][1][2]*255):02x}"


def generate_plddt_viewer(cif_path, output_dir, verbose=True):
    """Generate a multi-snapshot .mvsj file with per-model pLDDT colouring.

    Produces a single MVSJ file with one snapshot per model.  The Mol*
    web viewer displays a state gallery allowing the user to switch
    between models while preserving pLDDT colouring.

    Files are written to *output_dir* with the naming pattern::

        {cif_stem}_plddt.mvsj

    Parameters
    ----------
    cif_path : str
        Path to the PyMOL-aligned multi-model CIF file.
    output_dir : str
        Directory where the MVSJ is written (also expected to contain the CIF).
    verbose : bool, default True
        If False, suppresses status prints to stdout.
    """
    cif_basename = os.path.basename(cif_path)
    cif_stem = os.path.splitext(cif_basename)[0]

    doc = gemmi.cif.read_file(str(cif_path))

    snapshots = []

    for model_idx, block in enumerate(doc):
        # --- Extract per-atom b-factors and build colour annotation ---
        chain_col = block.find_values("_atom_site.label_asym_id")
        seq_col = block.find_values("_atom_site.label_seq_id")
        atom_col = block.find_values("_atom_site.label_atom_id")
        b_col = block.find_values("_atom_site.B_iso_or_equiv")

        annotation = []
        for chain, seq, atom, b in zip(chain_col, seq_col, atom_col, b_col):
            annotation.append({
                "label_asym_id": str(chain),
                "label_seq_id": int(seq),
                "label_atom_id": str(atom),
                "color": _plddt_to_color(float(b)),
            })

        ann_bytes = json.dumps(annotation).encode()
        ann_b64 = base64.b64encode(ann_bytes).decode()
        data_uri_ann = f"data:application/json;base64,{ann_b64}"

        # --- Write single-model CIF to temp file and embed as data URI ---
        single_doc = gemmi.cif.Document()
        single_doc.add_copied_block(block)
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            single_doc.write_file(tmp_path)
            with open(tmp_path, "rb") as fh:
                cif_bytes = fh.read()
        finally:
            os.unlink(tmp_path)

        cif_b64 = base64.b64encode(cif_bytes).decode()
        cif_data_uri = f"data:chemical/x-cif;base64,{cif_b64}"

        # --- Build per-model snapshot ---
        b = mvs.create_builder()
        b.canvas(custom={"molstar_postprocessing": {"enable_outline": True, "enable_ssao": False}})
        ds = b.download(url=cif_data_uri)
        ps = ds.parse(format="mmcif")
        ms = ps.model_structure()

        # Polymer (DNA) → cartoon with pLDDT coloring + illustrative style
        poly_comp = ms.component(selector="polymer")
        poly_rep = poly_comp.representation(
            type="cartoon",
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        poly_rep.color_from_uri(
            uri=data_uri_ann,
            format="json",
            schema="atom",
            field_name="color",
        )

        # Ligand (small molecules) → ball-and-stick with pLDDT coloring
        lig_comp = ms.component(selector="ligand")
        lig_rep = lig_comp.representation(type="ball_and_stick", size_factor=0.5, custom={"molstar_representation_params": {"ignoreLight": True}})
        lig_rep.color_from_uri(
            uri=data_uri_ann,
            format="json",
            schema="atom",
            field_name="color",
        )

        snap = b.get_snapshot(
            title=f"Model {model_idx}",
            linger_duration_ms=3000,
            transition_duration_ms=500,
        )
        snapshots.append(snap)

    mvsj_name = f"{cif_stem}_plddt.mvsj"
    mvsj_path = os.path.join(output_dir, mvsj_name)
    states = States(
        metadata=GlobalMetadata(title=f"{cif_stem} pLDDT"),
        snapshots=snapshots,
    )
    states_dict = states.model_dump(exclude_none=True)
    with open(mvsj_path, "w") as fh:
        json.dump(states_dict, fh, indent=2)
    if verbose:
        print(f"  Wrote multi-model pLDDT MVSJ: {mvsj_name}")

    
# =====================================================================


def generate_constraint_viewer(cif_path, yaml_path, output_path, verbose=True):
    """Generate a multi-snapshot .mvsj file with constraint-component colouring.

    Loads contact constraints from the Boltz-2 YAML file, finds connected
    components (BFS), assigns each component a colour from the pastel
    palette, and creates one snapshot per model.

    Unconstrained residues are set to white so the coloured components
    are visually distinct.

    Parameters
    ----------
    cif_path : str
        Path to the PyMOL-aligned multi-model CIF file.
    yaml_path : str
        Path to the Boltz-2 input YAML file.
    output_path : str
        Destination path for the .mvsj file.
    verbose : bool, default True
        If False, suppresses status prints to stdout.
    """
    constraints = _load_constraints(yaml_path, verbose=verbose)
    if not constraints:
        if verbose:
            print("  No constraints found -- skipping constraint MVSJ generation.")
        return

    components = _find_connected_components(constraints)

    annotation = []
    for idx, component in enumerate(components):
        color = CONSTRAINT_COLORS[idx % len(CONSTRAINT_COLORS)]
        for chain, resnum in component:
            annotation.append({
                "label_asym_id": chain,
                "label_seq_id": resnum,
                "color": color,
            })

    ann_bytes = json.dumps(annotation).encode()
    ann_b64 = base64.b64encode(ann_bytes).decode()
    data_uri = f"data:application/json;base64,{ann_b64}"

    doc = gemmi.cif.read_file(str(cif_path))
    snapshots = []

    for model_idx, block in enumerate(doc):
        # Write single-model CIF to temp file and embed as data URI
        single_doc = gemmi.cif.Document()
        single_doc.add_copied_block(block)
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            single_doc.write_file(tmp_path)
            with open(tmp_path, "rb") as fh:
                cif_bytes = fh.read()
        finally:
            os.unlink(tmp_path)

        cif_b64 = base64.b64encode(cif_bytes).decode()
        cif_data_uri = f"data:chemical/x-cif;base64,{cif_b64}"

        b = mvs.create_builder()
        b.canvas(custom={"molstar_postprocessing": {"enable_outline": True, "enable_ssao": False}})
        ds = b.download(url=cif_data_uri)
        ps = ds.parse(format="mmcif")
        ms = ps.model_structure()

        # Polymer → cartoon → white + constraint coloring + illustrative style
        poly_comp = ms.component(selector="all")
        poly_rep = poly_comp.representation(
            type="cartoon",
            custom={"molstar_representation_params": {"ignoreLight": True}},
        )
        poly_rep.color(color="white")
        poly_rep.color_from_uri(
            uri=data_uri,
            format="json",
            schema="residue",
            field_name="color",
        )

        # Ligand → ball-and-stick → element (CPK) coloring
        lig_comp = ms.component(selector="ligand")
        lig_rep = lig_comp.representation(type="ball_and_stick", size_factor=0.5, custom={"molstar_representation_params": {"ignoreLight": True}})
        lig_rep.color_from_source(
            schema="atom",
            category_name="_atom_site",
            field_name="type_symbol",
            palette=CategoricalPalette(colors="ElementSymbol"),
        )

        snap = b.get_snapshot(
            title=f"Model {model_idx}",
            linger_duration_ms=3000,
            transition_duration_ms=500,
        )
        snapshots.append(snap)

    mvsj_name = os.path.basename(output_path)
    states = States(
        metadata=GlobalMetadata(title=os.path.splitext(mvsj_name)[0]),
        snapshots=snapshots,
    )
    states_dict = states.model_dump(exclude_none=True)
    with open(output_path, "w") as fh:
        json.dump(states_dict, fh, indent=2)
    if verbose:
        print(f"  Wrote multi-model constraint MVSJ: {output_path}")


# =====================================================================
#  Constraint helpers (adapted from process_boltz_results.py)
# =====================================================================


def _load_constraints(yaml_path, verbose=True):
    """Load unique contact constraints from a Boltz-2 YAML file.

    Returns ``None`` if the file doesn't exist or the constraints section
    is empty/missing.
    """
    if not os.path.exists(yaml_path):
        if verbose:
            print(f"  YAML file not found: {yaml_path}")
        return None

    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    constraints_raw = data.get("constraints", [])
    if not constraints_raw:
        if verbose:
            print("  No constraints section found in YAML file.")
        return None

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

    return [tuple(sorted(p)) for p in unique_pairs]


def _find_connected_components(constraints):
    """Build an adjacency graph and return connected components (BFS).

    Parameters
    ----------
    constraints : list of tuple
        List of ``((chain_i, res_i), (chain_j, res_j))`` pairs.

    Returns
    -------
    list of set
        Each set contains ``(chain, resnum)`` tuples for one component.
    """
    graph = {}
    for node1, node2 in constraints:
        graph.setdefault(node1, set()).add(node2)
        graph.setdefault(node2, set()).add(node1)

    if not graph:
        return []

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

    return components
