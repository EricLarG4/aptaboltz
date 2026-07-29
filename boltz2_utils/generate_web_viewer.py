"""Mol* Web Viewer generation for Boltz-2 prediction outputs.

Produces MolViewSpec JSON (.mvsj) files for interactive 3D visualisation
via the quarto-molstar extension.
"""

import json
import os
import base64

import gemmi
import yaml
import molviewspec as mvs
from molviewspec.nodes import MVSJ


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


def generate_plddt_viewer(cif_path, output_dir):
    """Generate per-model .mvsj files with b-factor (pLDDT) colouring.

    Each model from the multi-model CIF gets its own individual CIF and
    MVSJ file.  pLDDT scores are read from ``_atom_site.B_iso_or_equiv``,
    averaged per residue, and pre-computed into hex colours embedded via
    a ``color_from_uri`` data URI.

    Files are written to *output_dir* with the naming pattern::

        {cif_stem}_model_{N}.cif
        {cif_stem}_plddt_model_{N}.mvsj

    Parameters
    ----------
    cif_path : str
        Path to the PyMOL-aligned multi-model CIF file.
    output_dir : str
        Directory where per-model CIFs and MVSJs are written.
    """
    cif_basename = os.path.basename(cif_path)
    cif_stem = os.path.splitext(cif_basename)[0]

    doc = gemmi.cif.read_file(str(cif_path))

    for model_idx, block in enumerate(doc):
        # --- Extract per-residue average b-factors for this model ---
        chain_col = block.find_values("_atom_site.label_asym_id")
        seq_col = block.find_values("_atom_site.label_seq_id")
        b_col = block.find_values("_atom_site.B_iso_or_equiv")

        residue_bfactors = {}
        for chain, seq, b in zip(chain_col, seq_col, b_col):
            key = (str(chain), int(seq))
            residue_bfactors.setdefault(key, []).append(float(b))

        # --- Compute annotation (pre-computed hex colours) ---
        annotation = []
        for (chain, seq), b_vals in residue_bfactors.items():
            avg_b = sum(b_vals) / len(b_vals)
            annotation.append({
                "label_asym_id": chain,
                "label_seq_id": seq,
                "color": _plddt_to_color(avg_b),
            })

        # --- Write individual (single-model) CIF ---
        model_cif_name = f"{cif_stem}_model_{model_idx}.cif"
        model_cif_path = os.path.join(output_dir, model_cif_name)
        single_doc = gemmi.cif.Document()
        single_doc.add_copied_block(block)
        single_doc.write_file(model_cif_path)

        # Read CIF bytes back for embedding in MVSJ
        with open(model_cif_path, "rb") as fh:
            cif_bytes = fh.read()
        cif_b64 = base64.b64encode(cif_bytes).decode()

        # --- Build MVSJ with data-URI annotation ---
        ann_bytes = json.dumps(annotation).encode()
        ann_b64 = base64.b64encode(ann_bytes).decode()
        data_uri_ann = f"data:application/json;base64,{ann_b64}"

        b = mvs.create_builder()
        ds = b.download(url=f"data:chemical/x-cif;base64,{cif_b64}")
        ps = ds.parse(format="mmcif")
        ms = ps.model_structure()
        comp = ms.component(selector="all")
        rep = comp.representation(type="cartoon")

        rep.color_from_uri(
            uri=data_uri_ann,
            format="json",
            schema="residue",
            field_name="color",
        )

        mvsj_name = f"{cif_stem}_plddt_model_{model_idx}.mvsj"
        mvsj_path = os.path.join(output_dir, mvsj_name)
        MVSJ(data=b.get_state()).dump(mvsj_path, indent=2)
        print(f"  Wrote pLDDT MVSJ:     {mvsj_name}")


# =====================================================================
#  Constraint viewer
# =====================================================================


def generate_constraint_viewer(cif_path, yaml_path, output_path):
    """Generate a .mvsj file with constraint-component colouring.

    Loads contact constraints from the Boltz-2 YAML file, finds connected
    components (BFS), assigns each component a colour from the pastel
    palette, and embeds the annotation as a data URI.

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
    """
    constraints = _load_constraints(yaml_path)
    if not constraints:
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

    # Embed the CIF as a data URI for self-contained MVSJ
    with open(cif_path, "rb") as fh:
        cif_bytes = fh.read()
    cif_b64 = base64.b64encode(cif_bytes).decode()

    b = mvs.create_builder()
    ds = b.download(url=f"data:chemical/x-cif;base64,{cif_b64}")
    ps = ds.parse(format="mmcif")
    ms = ps.model_structure()
    comp = ms.component(selector="all")
    rep = comp.representation(type="cartoon")

    rep.color(color="white")
    rep.color_from_uri(
        uri=data_uri,
        format="json",
        schema="residue",
        field_name="color",
    )

    MVSJ(data=b.get_state()).dump(output_path, indent=2)
    print(f"  Wrote constraint MVSJ: {output_path}")


# =====================================================================
#  Constraint helpers (adapted from process_boltz_results.py)
# =====================================================================


def _load_constraints(yaml_path):
    """Load unique contact constraints from a Boltz-2 YAML file.

    Returns ``None`` if the file doesn't exist or the constraints section
    is empty/missing.
    """
    if not os.path.exists(yaml_path):
        print(f"  YAML file not found: {yaml_path}")
        return None

    with open(yaml_path) as fh:
        data = yaml.safe_load(fh)

    constraints_raw = data.get("constraints", [])
    if not constraints_raw:
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
