"""
Template input file generator for Boltz-2.
============================================
INPUT FILE — Copy this file into your project directory and edit the
placeholders below to define DNA sequences, ligands, and generation
parameters for Boltz-2 YAML and SLURM files.

HOW TO USE
----------
1. Set `project` to project directory name (controls all local
   directory paths and the cluster output path).
   Define DNA sequences in SEQUENCES dict and set stem_length /
   max_distance for the stem base-pair constraints.  [STEP 1]
2. (Optional) Define `additional_pairs` for extra contact constraints.
   Three formats are accepted (see STEP 2 for full details and examples):
   - **Flat list** `[(res_i, res_j), ...]` — same pairs for all entries.
   - **Per-sequence dict** `{seq: [(chain, res_i, chain, res_j), ...]}`
   - **Per-seq + per-ligand dict** with `"*"` wildcard fallback.
   When resolved to a non-empty list, filenames include a `_constrained`
   suffix (e.g. `Seq1_C0R_constrained.yaml`).  Set to `None` to omit.
   [STEP 2]
3. Define LIGANDS — specify ligands by CCD code or SMILES string.
   [STEP 3]
4. [Automatic] The script builds the flat list of (sequence × ligand)
   entries.  No user action required.  [STEP 4]
5. Toggle `custom_msa` (True / False) depending on whether to use
   pre-computed MSA (.a3m) files.
   a. If True, place MSA files in {project}/msa/ and name them
      {seq_name}.a3m.
   b. When both `custom_msa` and `additional_pairs` are used, filenames
      combine both suffixes: e.g. `Seq1_C0R_custom_msa_constrained.yaml`.
   [STEP 5]
6. Run the script.
   It will:
   a. Create YAML files in {project}/yaml/  [STEP 6]
   b. Create SLURM submission scripts in {project}/  [STEP 7]
7. Send the slurm files and msa/ and yaml/ directories to the cluster.
8. Run Boltz-2 on the cluster. Commands depend on which sets you want:
   a. Run ALL jobs (standard, custom_msa, constrained):
      for f in scratch/boltz/projects/{project}/*.slurm; do sbatch "$f"; sleep 2; done
   b. Run only custom_msa jobs:
      for f in scratch/boltz/projects/{project}/*_custom_msa.slurm; do sbatch "$f"; sleep 2; done
   c. Run only constrained jobs:
      for f in scratch/boltz/projects/{project}/*_constrained.slurm; do sbatch "$f"; sleep 2; done
   d. Run only jobs that are NOT custom_msa:
      for f in scratch/boltz/projects/{project}/*.slurm; do [[ "$f" != *_custom_msa.slurm ]] && sbatch "$f" && sleep 2; done
   e. Run only jobs that are NOT constrained:
      for f in scratch/boltz/projects/{project}/*.slurm; do [[ "$f" != *_constrained.slurm ]] && sbatch "$f" && sleep 2; done
"""

from boltz2_utils.generate_boltz2_yaml import generate_boltz2_yamls
from boltz2_utils.generate_boltz2_slurm import generate_boltz2_slurm_files
from pathlib import Path


# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║  STEP 1 — Set project name and define DNA sequences                         ║
# ╚═════════════════════════════════════════════════════════════════════════════╝
# `project` is used throughout this script to build local directory paths
# (e.g. {project}/yaml/ for YAML output, {project}/msa/ for MSA files) and
# the cluster output path (scratch/boltz/projects/{project}/).
# Must match the project directory name.

project = "PLACEHOLDER_PROJECT"     # ← Change to your project directory name

# Molecule type: "dna" (default), "rna", or "protein".
# Controls the sequence-type key in the YAML output.
# Stem-contact constraints are only generated for "dna" and "rna".
molecule_type = "dna"

custom_msa = False                  # ← Toggle custom MSA usage (see STEP 5 below)

# Each key is a short name used in filenames (e.g. "Seq1").
# Each value is the nucleotide sequence (DNA/RNA).
#
#   TO MODIFY: Add/remove entries here as needed.
#   KEY FORMAT: A string that will appear in YAML/SLURM filenames.
#   VALUE:      The nucleic-acid sequence (5' → 3').
#
SEQUENCES = {
    "PLACEHOLDER_SEQ_NAME": "PLACEHOLDER_SEQUENCE",
}

stem_length = 8
max_distance = 4.0

# ╔═════════════════════════════════════════════════════════════════════════════╗
# ║  STEP 2 — (Optional) Define additional residue-pair contact constraints     ║
# ╚═════════════════════════════════════════════════════════════════════════════╝
# Define extra contact constraints between specific residue pairs.  The
# following formats are accepted (from simplest to most powerful):
#
# ── Format A: Flat list (backward-compatible) ──────────────────────────────
#   A simple list of (res_i, res_j) 1-based residue indices.  Both residues
#   are assigned to the DNA chain (default "A").  The same set of pairs is
#   applied to every (sequence × ligand) entry.
#
#     additional_pairs = [(10, 30)]                       # single pair
#     additional_pairs = [(1, 44), (2, 43), (3, 42)]     # manual stem pairs
#     additional_pairs = None                             # no extra pairs
#
# ── Format B: Per-sequence dictionary ──────────────────────────────────────
#   A dict mapping each sequence name to its own list of pairs.  Use the
#   same keys as SEQUENCES.  Useful when different sequences need different
#   constraints, regardless of the ligand.
#
#   Each pair can be a 2-tuple (res_i, res_j) — both on the DNA chain —
#   OR a 4-tuple (chain_i, res_i, chain_j, res_j) for explicit chains per
#   residue (e.g. cross-chain contacts between DNA chain "A" and ligand
#   chain "B").
#
#     additional_pairs = {
#         "Seq1": [("A", 5, "A", 36)],          # DNA–DNA contact, explicit chains
#         "Seq2": [("A", 3, "B", 10)],           # cross-chain: DNA → ligand
#     }
#
# ── Format C: Per-sequence + per-ligand dictionary (with wildcard) ─────────
#   A two-level dict: outer key = sequence name, inner key = ligand name.
#   This lets you specify different constraints for different conditions
#   (e.g. "free" vs. a specific bound ligand).
#
#   Wildcard: Use "*" as the ligand key to provide a fallback for any
#   ligands not explicitly listed.  This avoids repeating the same set
#   of pairs for multiple ligands.
#
#     additional_pairs = {
#         "Seq1": {
#             "free": [("A", 1, "A", 40)],       # only when no ligand
#             "C0R":  [("A", 5, "B", 10)],       # only with C0R bound
#             "*":    [("A", 2, "A", 39)],        # fallback for all other ligands
#         },
#     }
#
# What happens when the *_constrained* suffix is added to filenames?
#   The generated YAML and SLURM filenames automatically include the
#   "_constrained" suffix whenever additional_pairs resolves to a non-empty
#   list for a given entry — regardless of which format you used.
#
# TO MODIFY: Choose one of the formats above and assign it below.
additional_pairs = None                         # ← Set to your constraint pairs, or keep None


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 3 — Define ligands                                                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# Each ligand entry specifies how the ligand is identified in the YAML file.
#
#   Option A — By CCD code  (recommended for standard ligands):
#       {"ccd": "C0R", "smiles": None}
#     → Uses the CCD chemical-component dictionary entry (e.g. from the PDB).
#
#   Option B — By SMILES string  (for novel / non-standard ligands):
#       {"ccd": None, "smiles": "O=C4..."}
#     → The SMILES string is placed directly in the YAML file.
#
#   TO MODIFY: Add/remove ligand definitions here.
#              Keep exactly the two keys: "ccd" and "smiles".
#              One must be a non-None value; If both are provided,
#              ccd will be used if provided, and smiles will be ignored.
#
LIGANDS = {
    "free": {"ccd": None, "smiles": None},                                 # no ligand
    "PLACEHOLDER_LIGAND_NAME": {"ccd": "PLACEHOLDER_CCD", "smiles": None},  # CCD-based
    # "PLACEHOLDER_LIGAND_NAME": {                                            # SMILES-based
    #     "ccd": None,
    #     "smiles": "PLACEHOLDER_SMILES",
    # },
}


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 4 — Build the flat list of (sequence × ligand) entries               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# This nested loop creates one entry for every combination of sequence and
# ligand.  For example, 3 sequences × 3 ligands = 9 entries.
#
# Each entry is a dict with keys:
#   seq_name     : str  — short name of the sequence
#   sequence     : str  — nucleic-acid sequence itself
#   ligand_name  : str  — short name of the ligand (for filenames)
#   ccd          : str | None — CCD code (or None if SMILES is used)
#   smiles       : str | None — SMILES string (or None if CCD is used)
#
#   NO MODIFICATION NEEDED unless there is a need for per-entry overrides
#   (see optional: dna_chain, ligand_chain, cyclic, stem_length,
#    max_distance, msa).
#
entries = [
    {
        "seq_name": seq_name,
        "sequence": sequence,
        "ligand_name": lig_name,
        "ccd": lig["ccd"],
        "smiles": lig["smiles"],
    }
    for seq_name, sequence in SEQUENCES.items()
    for lig_name, lig in LIGANDS.items()
]

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 5 — Toggle custom MSA usage                                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# NOTE:
#   If custom_msa is True, each entry looks for an MSA file named
#   "{project}/msa/{seq_name}.a3m".  If the file exists, the path is added
#   to that entry's dict as the "msa" value.  If it does NOT exist, the
#   "msa" value is set to None, and the YAML file will not include an MSA
#   field.
#
#   To add custom MSA files:
#     Place them in  {project}/msa/  and name them  {seq_name}.a3m
#     (e.g. Seq1.a3m, Seq2.a3m).
#
#   The `project` variable is used in the path construction below.
#
if custom_msa:
    for entry in entries:
        msa_path_full = Path(f"{project}/msa/{entry['seq_name']}.a3m")
        msa_path_yaml = f"msa/{entry['seq_name']}.a3m"
        entry["msa"] = str(msa_path_yaml) if msa_path_full.exists() else None

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 6 — Generate YAML files (Boltz-2 input)                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# The output directory ({project}/yaml/) is created automatically.
#
# KEYWORD ARGUMENTS (adjust as needed):
#   output_dir    : str   — where to write the .yaml files
#   stem_length   : int   — number of Watson-Crick bp constraints at each end
#                           of the stem.  Set to 0 to omit constraints block.
#   max_distance  : float — upper-bound distance (Å) for each stem contact.
#
paths = generate_boltz2_yamls(
    entries,
    output_dir=f"{project}/yaml",
    stem_length=stem_length,
    max_distance=max_distance,
    additional_pairs=additional_pairs,
    molecule_type=molecule_type,
)

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 7 — Generate SLURM submission scripts                               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# One .slurm file is created per YAML file.
#
# The custom_msa toggle is automatically reflected in the SLURM scripts:
#   - If a YAML file has an "msa" field → the SLURM script OMITS the
#     --max_msa_seqs and --use_msa_server flags (because Boltz-2 will read
#     the MSA from the YAML file).
#   - If a YAML file does NOT have an "msa" field → the SLURM script
#     INCLUDES those flags, so Boltz-2 will generate an MSA on the fly
#     via the MSA server.
#
# KEYWORD ARGUMENTS (adjust as needed):
#   cluster_path : str — path relative to $SLURM_SUBMIT_DIR where output
#                        directories will be created on the cluster.
#                        This is built from the `project` variable:
#                        f"scratch/boltz/projects/{project}"
#   output_dir   : str — local directory for the generated .slurm files.
#                        Uses the `project` variable directly.
#
slurm_paths = generate_boltz2_slurm_files(
    paths,
    cluster_path=f"scratch/boltz/projects/{project}",
    output_dir=project,
)