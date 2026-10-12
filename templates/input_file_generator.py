"""
Template input file generator for Boltz-2.
============================================
INPUT FILE — Copy this file into your project directory.

Three mutually exclusive modes of operation:

1. Config file (recommended for complex projects):
       python input_file_generator.py --config config.json

2. CLI flags (quick runs):
       python input_file_generator.py --project my_project ^
           --seq Seq1 "GGGACGACGCCCGCATGTTCCATGGATAGTCTTGACTAGTCGTCCC" ^
           --ligand C0R "C0R" --free

3. Console script (if boltz2_utils is installed):
       boltz-generate-input --project my_project --seq Seq1 "ACG..." --free

WORKING DIRECTORY MATTERS
-------------------------
`project` is resolved relative to the CURRENT WORKING DIRECTORY, not to
this script's location and not to the repository root.  YAML files are
written to `{project}/yaml/` and SLURM files to `{project}/`, so running
this from inside the project directory creates a nested `my_project/`
inside itself.  Always launch from the directory that CONTAINS the project
folder:

    cd <repo root>
    python my_project/input_file_generator.py --config my_project/input_file_config.json

The same applies to `output_file_processing.py` (README §5) and to the
console scripts.

For full documentation see README.md §3.
"""

from boltz2_utils.cli.generate_input import main

if __name__ == "__main__":
    main()
