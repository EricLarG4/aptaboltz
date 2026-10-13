"""
Template output file processor for Boltz-2.
==============================================
OUTPUT FILE — Copy this file into your project directory.

Processes Boltz-2 prediction outputs: PyMOL alignment, ligand extraction,
confidence/PAE/PDE/pLDDT CSV aggregation.

Three mutually exclusive modes of operation:

1. CLI flags (recommended for quick runs):
       python output_file_processing.py --project my_project ^
           --model Seq1 --model Seq2 --model Seq3 ^
           --ligand-dict "{\"Ligand01\": \"L01\", \"free\": \"free\"}" ^
           --name-map "{\"L01\": \"Ligand01\"}"

2. Config file:
       python output_file_processing.py --config config.json

3. Console script (if boltz2_utils is installed):
       boltz-process-output --project my_project --model Seq1 --ligand-dict "..." --name-map "..."

WORKING DIRECTORY MATTERS
-------------------------
`project` is resolved relative to the CURRENT WORKING DIRECTORY, not to
this script's location and not to the repository root.  Each experiment is
looked up at `{project}/{model}_{ligand}{suffix}/`, so running this from
inside the project directory makes it search for a nested path that does
not exist:

    cd my_project
    python output_file_processing.py --config output_file_config.json
        -> looks for my_project/my_project/<exp>/  ... skips everything

The failure is SILENT: each experiment prints
"Experiment does not exist: ... | skipping" and the script exits 0, so it
looks like it succeeded.  Always launch from the directory that CONTAINS
the project folder:

    cd <repo root>
    python my_project/output_file_processing.py --config my_project/output_file_config.json

The same applies to the console scripts, and to `input_file_generator.py`
(README §3.4) — both resolve `project` against the cwd.

For full documentation see README.md §5.
"""

from boltz2_utils.cli.process_output import main

if __name__ == "__main__":
    main()
