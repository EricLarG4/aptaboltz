# Template process.R for Boltz-2 visualisation
# ==============================================
# OUTPUT FILE — Copy this file into your project directory and edit the
# placeholders below to generate confidence tables, PAE/PDE heatmaps, and
# pLDDT per-residue plots from the CSV files produced by
# output_file_processing.py.
#
# Prerequisites:
#   1. Run install_packages.R (once per machine) to set up the 'boltz'
#      conda environment and install required R packages.
#   2. Run output_file_processing.py first:
#         reticulate::source_python("PLACEHOLDER_PROJECT/output_file_processing.py")
#      This generates the confidence, PAE, PDE, and pLDDT CSV files
#      that this script reads.
#   3. The 'boltz' conda environment must be active (reticulate will
#      use it automatically after install_packages.R has been sourced).
#
# Dependencies:
#   - R packages:  defined in install_packages.R
#   - Python env:  defined in environment.yml (conda env 'boltz')

source('boltz_R_utils/install_packages.R')
source('boltz_R_utils/processor.R')

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — Edit these values for your project                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

project <- 'PLACEHOLDER_PROJECT'  # ← Change to your project name

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  1. Interactive confidence table (saved as HTML)                          ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

table_confidence(process_confidence(project))

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  2. PAE / PDE heatmaps                                                   ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# `ligand_number` is the residue index at which the ligand region begins
# (i.e. the number of DNA residues).  The matrix is split into DNA and
# Ligand regions at this boundary.  Set to NULL to plot the full matrix.
#
# Example: for a 46-mer DNA aptamer, the ligand starts at residue 47:
#   ligand_number = 47
# For a 30-mer RNA, ligand starts at residue 31:
#   ligand_number = 31

pae_dt <- process_pxe(project, type = "pae")
pde_dt <- process_pxe(project, type = "pde")

lapply(
  names(pae_dt),
  function(n) {
    plot_pxe(
      pae_dt[n],
      type = "pae",
      ligand_number = PLACEHOLDER_LIGAND_NUMBER  # ← Set to number of DNA residues
    )
  }
)

lapply(
  names(pde_dt),
  function(n) {
    plot_pxe(
      pde_dt[n],
      type = "pde",
      ligand_number = PLACEHOLDER_LIGAND_NUMBER  # ← Set to number of DNA residues
    )
  }
)

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  3. pLDDT per-residue line plots                                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

plot_plddt(project)