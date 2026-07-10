"""
Template output file processor for Boltz-2.
==============================================
OUTPUT FILE — Copy this file into your project directory and edit the
placeholders below to process Boltz-2 prediction outputs (PyMOL alignment,
ligand extraction, confidence/PAE/PDE/pLDDT CSV aggregation).  This script
sources the shared utilities from `boltz2_utils/`.

See README.md §5 for full documentation.

CONFIGURATION NOTES
-------------------
LGD_DICT maps display names → CCD codes.
  Keys   = aptamer display names (used in PNG titles, legends).
  Values = CCD ligand codes matching the Boltz-2 YAML ligand definitions.
  The key "free" must be present with value "free" for apo structures.

NAME_MAP reverses LGD_DICT: CCD code → display name.
  This is used for file naming and chirality-grid annotations.
  Only the bound ligands need entries; "free" is handled automatically.

Example for 3 sequences, 3 ligands (2 CCD + 1 SMILES):
  LGD_DICT = {"Display1": "C0R", "Display2": "HCY", "Display3": "11DC", "free": "free"}
  NAME_MAP = {"C0R": "Display1", "HCY": "Display2", "11DC": "Display3"}

models: list of sequence names matching the SEQUENCES keys used during
        YAML generation (e.g. ["CSS1", "CSS2", "CSS3"]).

suffixes: YAML/SLURM suffixes to process.  Typically ["_constrained", ""]
          to process both constrained and unconstrained predictions.
"""

from boltz2_utils.process_boltz_results import process_single_experiment
from itertools import product

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — Edit these values for your project                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Dictionary mapping aptamer display names to their PDB/CCD ligand codes.
# Used for ligand identification in the structure and for fetching reference
# conformations from the PDB.
#
# Example:
#   LGD_DICT = {"Corticosteron": "C0R", "Cortisol": "HCY",
#               "11-deoxycortisol": "11DC", "free": "free"}
#
LGD_DICT = {
    "PLACEHOLDER_APTAMER": "PLACEHOLDER_CCD",
    "free": "free",
}

# Reverse mapping: CCD code → display name (for file naming and legends).
# Only bound ligands need entries.
#
# Example:
#   NAME_MAP = {"C0R": "Corticosteron", "HCY": "Cortisol", "11DC": "11-deoxycortisol"}
#
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

    # job_dirs controls which Boltz-2 job directories to process:
    #   None               — process ALL J... directories
    #   "last"             — process only the highest-numbered directory
    #   "J1121077"         — process only that specific directory
    #   ["J1121077", "J1123593"] — process only the listed directories
    #
    # Default "last" is usually what you want (latest run).
    process_single_experiment(
        project,
        experiment_name,
        LGD_DICT=LGD_DICT,
        NAME_MAP=NAME_MAP,
        job_dirs="last",
    )