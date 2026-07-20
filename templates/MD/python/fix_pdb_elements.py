#!/usr/bin/env python
"""Add PDB element symbols (cols 77-78) to Amber AC-format PDB files.

AmberTools antechamber/prepgen generate "PDB" files that are actually Amber
AC format — they lack the standard PDB element field in columns 77-78.
This script detects the element from the atom name (columns 13-16) and
inserts the correct right-justified element symbol.

Usage:
    fix_pdb_elements.py <input.pdb> [--output <output.pdb>]

If --output is omitted, the file is edited in-place.
"""

import sys
import os


def get_element(atom_name: str) -> str:
    """Determine the PDB element symbol from a left-justified 4-char atom name.

    Returns a right-justified 2-char element string (e.g. " C", "Cl", " H", " N").
    """
    name = atom_name.strip()
    if not name:
        return "  "

    # Strip trailing digits ("Cl1" -> "Cl", "C25" -> "C", "H1" -> "H")
    candidate = name.rstrip("0123456789")

    if not candidate:
        candidate = name[0]

    if len(candidate) == 1:
        return f" {candidate}"

    # Two or more characters
    if candidate[0].isupper() and candidate[1:].islower():
        # Multi-character element: Cl, Br, Na, Fe, etc.
        return candidate[:2]
    else:
        # PDB position labels (CA, CB, CD, CG, etc.) -> use first char
        return f" {candidate[0]}"


def fix_pdb(input_path: str, output_path: str | None = None) -> None:
    """Read an AC-format PDB, add element symbols, and write the result."""
    if output_path is None:
        output_path = input_path

    lines_out = []
    with open(input_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                if len(line) < 78:
                    line = line.rstrip("\n").ljust(78) + "\n"
                atom_name = line[12:16]
                element = get_element(atom_name)
                line = line[:76] + element + line[78:]
            lines_out.append(line)

    with open(output_path, "w") as fh:
        fh.writelines(lines_out)

    print(f"Fixed: {input_path} -> {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = None
    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        output_pdb = sys.argv[idx + 1]

    if not os.path.exists(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    fix_pdb(input_pdb, output_pdb)
