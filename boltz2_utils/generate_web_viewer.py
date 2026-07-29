"""Mol* Web Viewer generation for Boltz-2 prediction outputs.

Produces MolViewSpec JSON (.mvsj) files for interactive 3D visualisation
via the quarto-molstar extension.
"""

import json
import os
import base64

import yaml
import molviewspec as mvs
from molviewspec.nodes import MVSJ, ContinuousPalette


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


# =====================================================================
#  pLDDT viewer
# =====================================================================


def generate_plddt_viewer(cif_path, output_path):
    """Generate a .mvsj file with b-factor (pLDDT) colouring.

    Uses ``color_from_source`` to read the ``_atom_site.B_iso_or_equiv``
    field and map it through a continuous gradient matching the PyMOL
    *tv_red -> yelloworange -> palecyan -> density* spectrum, with the
    range clamped to [50, 90].

    Parameters
    ----------
    cif_path : str
        Path to the PyMOL-aligned multi-model CIF file.
    output_path : str
        Destination path for the .mvsj file.
    """
    cif_filename = os.path.basename(cif_path)

    b = mvs.create_builder()
    ds = b.download(url=cif_filename)
    ps = ds.parse(format="cif")
    ms = ps.model_structure()
    comp = ms.component(selector="all")
    rep = comp.representation(type="cartoon")

    rep.color_from_source(
        schema="atom",
        category_name="atom_site",
        field_name="B_iso_or_equiv",
        palette=ContinuousPalette(
            colors=[
                ("#FF0000", 50.0),    # tv_red
                ("#FFA500", 63.3),    # yelloworange
                ("#AFEEEE", 76.6),    # palecyan
                ("#E0FFFF", 90.0),    # density
            ],
            mode="absolute",
            value_domain=(50.0, 90.0),
        ),
    )

    MVSJ(data=b.get_state()).dump(output_path, indent=2)
    print(f"  Wrote pLDDT MVSJ: {output_path}")


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

    cif_filename = os.path.basename(cif_path)

    b = mvs.create_builder()
    ds = b.download(url=cif_filename)
    ps = ds.parse(format="cif")
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
