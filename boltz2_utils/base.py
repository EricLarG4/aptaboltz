from pymol import cmd, util

def setup_base():
    cmd.space("pymol")
    cmd.set("cartoon_ring_mode", 3)

    # Hide solvent and ions
    cmd.select("water", "resn WAT")
    cmd.hide("everything", "water")

    cmd.set("sphere_transparency", 0.3)

    cmd.select("sodium", "name Na+")
    cmd.hide("everything", "sodium")

    cmd.select("chloride", "name Cl-")
    cmd.hide("everything", "chloride")

    cmd.select("none")
    cmd.alignto()

    # Colour nucleic acids by base identity (works for both DNA and RNA)
    cmd.color("wheat",     "resn DG+DG5+DG3+G+G3+G5+DI")
    cmd.color("palecyan",  "resn DA+DA5+DA3+A+A5+A3")
    cmd.color("palegreen", "resn DT+DT5+DT3+T+T5+T3")
    cmd.color("lime",      "resn DU+U+DU3+U3+DU5+U5")
    cmd.color("lightpink", "resn DC+DC5+DC3+C+C5+C3")

    cmd.set("cartoon_nucleic_acid_color", "grey90")
    cmd.set("cartoon_discrete_colors", 1)

    # Colour protein chains in a single muted, elegant colour
    # `polymer.protein` is a PyMOL built-in keyword that selects all
    # protein residues (amino acids) regardless of chain.
    if cmd.count_atoms("polymer.protein") > 0:
        cmd.color("warmpink", "polymer.protein")

    # Colour small molecules by element
    cmd.select("smallmol", "br. organic")
    util.cbay("smallmol")

    cmd.select("none")
    cmd.hide("everything", "hydro")
    cmd.set("sphere_scale", 0.5)


setup_base()
