"""
Template output file processor for Boltz-2.
==============================================
OUTPUT FILE — Copy this file into your project directory.

Processes Boltz-2 prediction outputs: PyMOL alignment, ligand extraction,
confidence/PAE/PDE/pLDDT CSV aggregation.

Three mutually exclusive modes of operation:

1. CLI flags (recommended for quick runs):
       python output_file_processing.py --project CSS ^
           --model CSS1 --model CSS2 --model CSS3 ^
           --ligand-dict "{\"Corticosteron\": \"C0R\", \"free\": \"free\"}" ^
           --name-map "{\"C0R\": \"Corticosteron\"}"

2. Config file:
       python output_file_processing.py --config config.json

3. Console script (if boltz2_utils is installed):
       boltz-process-output --project CSS --model CSS1 --ligand-dict "..." --name-map "..."

For full documentation see README.md §5.
"""

from boltz2_utils.cli.process_output import main

if __name__ == "__main__":
    main()
