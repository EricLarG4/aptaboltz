#!/usr/bin/env python
"""
Template SLURM generator for Amber MD Preparation + Final Minimisation
=======================================================================

Reads the template PrepAndMin.slurm, replaces placeholders with your
project parameters, and writes a configured .slurm file ready for sbatch.

Usage (from the MD/ directory)::

    python generate_md_slurm.py \\
        --experiment CSS1_HCY_constrained \\
        --array 1,2,3,5,7 \\
        --scratchdir ${SLURM_SUBMIT_DIR}/scratch/boltz/projects/CSS/MD

Output:  ``./slurm/<experiment>_PrepAndMin.slurm``

All parameters can be overridden via CLI flags; sensible defaults are
provided for the most common cluster setup.  See ``--help`` for details.
"""

import argparse
import os
import pathlib

PLACEHOLDER_MAP = {
    "__JOB_NAME__": None,
    "__EXPERIMENT__": None,
    "__ARRAY__": None,
    "__SCRATCHDIR__": None,
    "__GPU_PARTITION__": None,
    "__RESIDUES__": None,
    "__SOLVENT__": None,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a configured PrepAndMin.slurm for an Amber MD run."
    )
    parser.add_argument(
        "--experiment",
        required=True,
        help="Experiment name (e.g. CSS1_HCY_constrained). "
        "Used for SBATCH -J, EXPERIMENT=, and the output filename.",
    )
    parser.add_argument(
        "--array",
        default="0-4",
        help="SLURM array task IDs (e.g. 0-4, 1,2,3,5,7). Default: 0-4",
    )
    parser.add_argument(
        "--scratchdir",
        required=True,
        help="Path to the scratch directory containing pmemd/ and leap/.",
    )
    parser.add_argument(
        "--gpu",
        default="gpu-h100,gpu-l40",
        help="GPU partition(s) for #SBATCH -p. Default: gpu-h100,gpu-l40",
    )
    parser.add_argument(
        "--residues",
        default="1-46",
        help="Residue mask for cpptraj align (e.g. 1-46). Default: 1-46",
    )
    parser.add_argument(
        "--solvent",
        default="WAT",
        help="Solvent residue name for cpptraj strip (e.g. WAT). Default: WAT",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    script_dir = pathlib.Path(__file__).resolve().parent
    template_path = script_dir / "slurm" / "PrepAndMin.slurm"

    if not template_path.exists():
        print(f"Error: template not found at {template_path}")
        raise SystemExit(1)

    content = template_path.read_text(encoding="utf-8")

    replacements = {
        "__JOB_NAME__": args.experiment,
        "__EXPERIMENT__": args.experiment,
        "__ARRAY__": args.array,
        "__SCRATCHDIR__": args.scratchdir,
        "__GPU_PARTITION__": args.gpu,
        "__RESIDUES__": args.residues,
        "__SOLVENT__": args.solvent,
    }

    for placeholder, value in replacements.items():
        content = content.replace(placeholder, value)

    found = [p for p in PLACEHOLDER_MAP if p in content]
    if found:
        print(f"Warning: unresolved placeholders remain: {found}")

    output_dir = pathlib.Path("slurm").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{args.experiment}_PrepAndMin.slurm"
    output_path.write_text(content, encoding="utf-8", newline="\n")

    print(f"Generated: {output_path}")


if __name__ == "__main__":
    main()
