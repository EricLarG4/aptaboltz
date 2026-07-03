"""
Template output file processor for Boltz-2.
==============================================
OUTPUT FILE — Copy this file into your project directory and edit the
placeholders below to process Boltz-2 prediction outputs (PyMOL alignment,
ligand extraction, confidence/PAE/PDE/pLDDT CSV aggregation).  This script
sources the shared utilities from `boltz2_utils/`.

See README.md §5 for full documentation.
"""

from boltz2_utils.process_boltz_results import process_single_experiment
from itertools import product

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — Edit these values for your project                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Dictionary mapping aptamer display names to their PDB/CCD ligand codes.
# Used for ligand identification in the structure and for fetching reference
# conformations from the PDB.
LGD_DICT = {
    "PLACEHOLDER_APTAMER": "PLACEHOLDER_CCD",
    "free": "free",
}

# Reverse mapping: CCD code -> display name (for file naming and legends).
NAME_MAP = {
    "PLACEHOLDER_CCD": "PLACEHOLDER_DISPLAY_NAME",
}

project = "PLACEHOLDER_PROJECT"                # ← Change to your project name
models = ["PLACEHOLDER_MODEL"]                 # ← List of model/sequence names
suffixes = ["_constrained", ""]                 # ← Typically include both


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  PROCESSING — Iterates over all (model × ligand × suffix) combinations   ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Process all model/ligand/suffix combinations
for model, lgd, suffix in product(models, LGD_DICT.values(), suffixes):
    print(f"Processing {model} with ligand {lgd} and suffix {suffix}")
    experiment_name = f"{model}_{lgd}{suffix}"

    # To process all job directories:        job_dirs=None
    # To process specific job directories:   job_dirs=["J1121077", "J1123593"]
    # To process a single job directory:     job_dirs="J1121077"
    # To process only the last (highest-numbered) directory:  job_dirs="last"
    process_single_experiment(
        project,
        experiment_name,
        LGD_DICT=LGD_DICT,
        NAME_MAP=NAME_MAP,
        job_dirs="last",
    )