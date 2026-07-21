#!/usr/bin/env python
"""
Template model_prep.py — Prepare Boltz-2 prediction CIFs for MD simulation.
============================================================================

OUTPUT FILE — Copy this file into your project's MD/ directory and edit the
placeholders below.

DESCRIPTION
  Loads Boltz-2 prediction .cif models into PyMOL, performs system-specific
  modifications (e.g. 5' phosphate removal + terminal hydroxyl capping for
  DNA aptamers), aligns the reference ligand PDB, reassigns chains, and
   saves ready-for-MD PDB files into pdb_for_md/<SEQ>_<OUTPUT_LGD>_constrained/.

  Two usage modes:
    1. Interactive (from within PyMOL after cd to MD/):
         run model_prep.py
         prep_model(seq="PLACEHOLDER_SEQ", lgd="PLACEHOLDER_LGD", ...)

    2. Command line (from bash / PowerShell):
          python model_prep.py [--seq PLACEHOLDER_SEQ] [--lgd PLACEHOLDER_LGD]
                               [--lgd-file placeholder_lgd_resp] [--unconstrained]
                               [--job J0000000] [--models 0-24]

     Use --lgd free (or any ligand name set to "free" in your Boltz-2 YAMLs)
     to skip all ligand-related operations for apo (ligand-free) predictions.

DEPENDENCIES
  - PyMOL (with Python API, e.g. pymol-open-source)
  - Boltz-2 prediction output in  ../<SEQ>_<LGD><suffix>/<JOB>/

CAUTION
  PyMOL operations (phosphate capping, terminal hydroxyl modification,
  ligand alignment, chain reassignment) are specific to **DNA aptamers**
  with a 5' terminal phosphate.  If your system uses RNA, protein, or a
  different terminal chemistry, or if the ligand name differs from the
  Boltz-2 prediction residue name, you MUST edit the marked blocks below.

USAGE OUTLINE
  1. Run Boltz-2 predictions on the cluster and copy results locally.
  2. Run ligand_prep.sh to generate ff/<lgd_file>.pdb.
  3. Edit the placeholders and the marked blocks below for your system.
  4. Run this script (from within the MD/ directory).
  5. Run leap.sh to solvate and ionise the prepared PDB models.

See also: README.md §8.7
"""

import os


def prep_model(
    seq="PLACEHOLDER_SEQ",
    lgd="PLACEHOLDER_LGD",
    lgd_file=None,
    constrained=True,
    job="J0000000",
    models=range(0, 25),
    output_lgd=None,
    use_pair_fit=False,
):
    """Load, modify, and save Boltz prediction models for MD.

    Parameters
    ----------
    seq : str
        Sequence identifier (e.g. "CSS1").         Used in directory paths.
    lgd : str
        Ligand identifier for directory path construction (e.g. "Piperaquine").
        When Boltz-2 was run with a SMILES ligand, the output residue is
        auto-detected as "LIG1" and renamed to output_lgd.
        Use "free" for apo (ligand-free) predictions.
    lgd_file : str or None
        Name of the reference ligand PDB file (without .pdb suffix),
        located in ff/{lgd_file}.pdb (from ligand_prep.sh).
        Ignored when lgd="free".
    output_lgd : str or None
        Desired three-letter residue name in the output PDB (e.g. "PQ").
        If None, defaults to the value of lgd.
    constrained : bool
        If True, appends "_constrained" to the run directory name.
    job : str
        Job identifier for locating Boltz-2 prediction files.
    models : iterable
        Model indices to process (e.g. range(0, 25)).
    use_pair_fit : bool
        If True, strip hydrogens from the reference ligand and use
        ``cmd.pair_fit`` (atom-type-based pairing) instead of
        ``cmd.align`` (name-based alignment).  Useful for SMILES-based
        predictions where atom names diverge between the QM reference
        and the Boltz prediction.
    """

    from pymol import cmd

    suffix = "_constrained" if constrained else ""
    has_ligand = lgd.lower() != "free"
    if output_lgd is None:
        output_lgd = lgd

    input_dir = f"../{seq}_{lgd}{suffix}/{job}/boltz_results_input/predictions/input"
    output_dir = f"pdb_for_md/{seq}_{output_lgd}{suffix}"

    os.makedirs(output_dir, exist_ok=True)

    for model in models:
        cif_path = f"{input_dir}/input_model_{models[model]}.cif"
        if has_ligand:
            lgd_path = f"ff/{lgd_file}.pdb"

        print(f"Processing input_model_{models[model]}==========")

        cmd.load(cif_path)
        print(f"Loaded {cif_path}")

        # ╔════════════════════════════════════════════════════════════════╗
        # ║  EDIT THIS BLOCK — System-specific PyMOL modifications        ║
        # ║                                                               ║
        # ║  The example below is for **DNA aptamers**: removes the 5'    ║
        # ║  phosphate and caps with a terminal hydroxyl.  For RNA,       ║
        # ║  proteins, or other systems, replace this block entirely.     ║
        # ╚════════════════════════════════════════════════════════════════╝

        # --- DNA-specific: remove 5' phosphate, cap terminal hydroxyl ---
        n_res = cmd.count_atoms("chain A and name P")
        cmd.remove("resi 1 and (name P or name OP1 or name OP2)")
        cmd.h_add("resi 1 and name O5'")
        cmd.alter("resi 1 and name H01", "name='H5T'")
        print("Terminal phosphate was removed and hydroxyl was capped.")

        # ────────────────────────────────────────────────────────────────
        #  END OF EDITABLE BLOCK
        # ────────────────────────────────────────────────────────────────

        if has_ligand:
            # Determine the actual ligand residue name in the prediction
            pred_lgd = output_lgd
            if cmd.count_atoms(f"input_model_{models[model]} and resn {output_lgd}") == 0 and cmd.count_atoms(f"input_model_{models[model]} and resn LIG1") > 0:
                pred_lgd = "LIG1"
                print(f"Ligand '{output_lgd}' not found; using 'LIG1' from prediction.")

            # Align the reference ligand PDB onto the predicted ligand
            cmd.load(lgd_path)
            if use_pair_fit:
                cmd.remove("resn UNL and h.")
                cmd.pair_fit("resn UNL", f"resn {pred_lgd}")
            else:
                cmd.align("resn UNL", f"resn {pred_lgd}")
            cmd.alter("resn UNL", f"resn = '{output_lgd}'")
            cmd.alter(f"{lgd_file} and resn {output_lgd}", f"resi={n_res + 1}")
            print(f"Ligand {output_lgd} aligned and relabeled to resi {n_res + 1}.")

            # Remove the original predicted ligand, reassign chain B
            cmd.remove(f"chain B and resn {pred_lgd}")
            cmd.alter(f"resn {output_lgd}", "chain='B'")
            cmd.alter("chain A", "segi=''")
            print(f"Original ligand removed. Chain B reassigned to {output_lgd}.")

            # Merge into a single object and clean up
            cmd.create("seq", f"input_model_{models[model]} or {lgd_file}")
            cmd.delete(f"input_model_{models[model]} or {lgd_file}")
            cmd.sort()
            print("Merged into a single object 'seq' and cleaned up.")
        else:
            # No ligand — rename the loaded model directly
            cmd.create("seq", f"input_model_{models[model]}")
            cmd.delete(f"input_model_{models[model]}")
            cmd.alter("chain A", "segi=''")
            cmd.sort()
            print("No ligand — model prepared without ligand operations.")

        out_path = f"{output_dir}/input_model_{model}.pdb"
        cmd.save(out_path, "seq", -1, "pdb")
        cmd.delete("seq")
        print(f"Saved modified model to {out_path}\n")

    print(f"========All models processed.======== \nOutput files are in {output_dir}.")

    readme_path = os.path.join(output_dir, "README.txt")
    with open(readme_path, "w") as readme_file:
        readme_file.write(f"Preparation Summary for {seq}_{output_lgd}{suffix}\n")
        readme_file.write(f"Job Identifier: {job}\n")
        if has_ligand:
            readme_file.write(f"Ligand Residue Name: {output_lgd}\n")
            readme_file.write(f"Ligand Reference File: {lgd_file}.pdb\n")
        readme_file.write(f"Constrained Run: {'Yes' if constrained else 'No'}\n")
        readme_file.write(f"Models Processed: {', '.join(map(str, models))}\n")
        readme_file.write("All models have been modified and saved in this directory.\n")
    print(f"README file created at {readme_path}.")


if __name__ == "__main__":
    import argparse
    import glob
    import re

    parser = argparse.ArgumentParser(
        description="Prepare Boltz-2 prediction CIFs for MD simulation."
    )
    parser.add_argument("--seq", default="PLACEHOLDER_SEQ", help="Sequence identifier")
    parser.add_argument("--lgd", default="PLACEHOLDER_LGD", help="Ligand identifier (directory path)")
    parser.add_argument(
        "--output-lgd",
        default=None,
        help="Desired output residue name (defaults to --lgd value)",
    )
    parser.add_argument(
        "--lgd-file",
        default=None,
        help="Reference ligand PDB (without extension, from ff/). Ignored when --lgd is 'free'.",
    )
    parser.add_argument(
        "--unconstrained",
        action="store_true",
        help="Omit the _constrained suffix",
    )
    parser.add_argument(
        "--job",
        default=None,
        help="Job identifier (default: latest J... directory)",
    )
    parser.add_argument(
        "--models",
        default="0-24",
        help="Model range, e.g. '0-24' or '0,1,2'",
    )
    parser.add_argument(
        "--use-pair-fit",
        action="store_true",
        help="Use pair_fit + hydrogen-stripping instead of align (SMILES ligands)",
    )
    args = parser.parse_args()

    suffix = "" if args.unconstrained else "_constrained"
    parent = f"../{args.seq}_{args.lgd}{suffix}"

    if args.job is None:
        job_dirs = glob.glob(f"{parent}/J*")
        job_nums = [
            int(re.sub(r".*J", "", d)) for d in job_dirs if re.search(r"J\d+$", d)
        ]
        if not job_nums:
            print(f"Error: no J... directories found in {parent}")
            exit(1)
        args.job = f"J{max(job_nums)}"
        print(f"Auto-detected latest job: {args.job}")

    if "-" in args.models:
        start, end = args.models.split("-", 1)
        models = range(int(start), int(end) + 1)
    else:
        models = [int(x) for x in args.models.split(",")]

    prep_model(
        seq=args.seq,
        lgd=args.lgd,
        lgd_file=args.lgd_file,
        constrained=not args.unconstrained,
        job=args.job,
        models=models,
        output_lgd=args.output_lgd,
        use_pair_fit=args.use_pair_fit,
    )
