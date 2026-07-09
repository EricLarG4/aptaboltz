#!/bin/bash
#===============================================================================
# leap.sh — Amber tLEAP solvation and ionisation for MD simulations
#
# DESCRIPTION
#   Creates a single tLEAP input file that loads force fields, the ligand
#   prepin/frcmod, and all PDB models in a given range.  Each model is
#   solvated in an OPC water box, neutralised with Na⁺/Mg²⁺/Cl⁻ ions, and
#   saved as Amber prmtop/rst7 files.  Finally, ambpdb converts each rst7
#   to PDB format.
#
#   The script mirrors the directory structure of pdb_for_md/ and writes
#   all output into leap/CSS1_<LIGAND_UPPER>_constrained/.
#
# DEPENDENCIES
#   - AmberTools (tleap, ambpdb)
#   - Bash (Linux/WSL)
#
# USAGE
#   ./leap.sh [-l ligand] [-s start_model] [-e end_model] \
#             [-n num_na] [-m num_mg] [-c total_cl] [-h]
#
# ARGUMENTS
#   -l <name>    Ligand name (default: PLACEHOLDER_LIGAND).
#                Used for ff/${ligand}_resp.prepin/.frcmod and
#                pdb_for_md/CSS1_<UPPER>_constrained paths.
#   -s <int>     Starting model index, inclusive (default: 0)
#   -e <int>     Ending model index, inclusive (default: 24)
#   -n <int>     Number of Na⁺ ions to add (default: 65)
#   -m <int>     Number of Mg²⁺ ions to add (default: 1)
#   -c <int>     Total number of Cl⁻ ions (default: 22).
#                Automatically split between MgCl₂ and NaCl.
#   -h           Show this help message and exit.
#
# Cl⁻ DISTRIBUTION
#   The total Cl⁻ count passed with -c is automatically split:
#     - MgCl₂ requires  2 Cl⁻ per Mg²⁺  →  cl_mg = 2 × -m
#     - NaCl  receives the remainder     →  cl_na = total_cl - cl_mg
#
# INPUT FILES REQUIRED
#   ff/${ligand}_resp.prepin     — Amber prepin library (from ligand_prep.sh)
#   ff/${ligand}_resp.frcmod     — Force-field modifications (from ligand_prep.sh)
#   pdb_for_md/CSS1_<UPPER>_constrained/input_model_${i}.pdb  — PDB models
#
# OUTPUT STRUCTURE
#   leap/CSS1_<LIGAND_UPPER>_constrained/
#     tleap.in                    — Combined tLEAP input script
#     cplx_${i}.prmtop            — Amber topology (one per model)
#     cplx_${i}.rst7              — Amber restart/coords (one per model)
#     cplx_${i}.pdb               — PDB conversion (one per model)
#
# WORKFLOW (full pipeline order)
#   1. QM optimisation + RESP2 charges: orca_steps_wsl.sh → multiwfn_steps_windows.ps1
#   2. Ligand parameterisation:          ligand_prep.sh
#   3. System solvation + ionisation:    leap.sh (this script)
#   4. MD production:                    (user-defined)
#
# EXAMPLES
#   # Run with defaults (ligand=PLACEHOLDER_LIGAND, models 0–24, 65 Na⁺, 1 Mg²⁺, 22 Cl⁻)
#   ./leap.sh
#
#   # Custom ligand and model range
#   ./leap.sh -l hcy -s 5 -e 18
#
#   # Custom ion counts (24 Cl⁻ total → 2 for MgCl₂ + 22 for NaCl)
#   ./leap.sh -n 60 -m 1 -c 24
#
#   # Show help
#   ./leap.sh -h
#
# TROUBLESHOOTING
#   - Ensure the ligand prepin and frcmod files exist in ff/ before running.
#   - Check that PDB files exist in pdb_for_md/CSS1_<UPPER>_constrained/.
#   - If tleap reports missing parameters, review the frcmod file from
#     ligand_prep.sh for ATTN warnings.
#   - The OPC water model requires leaprc.water.opc (AmberTools ≥ 18).
#===============================================================================

# Default values
ligand="PLACEHOLDER_LIGAND"
start_model=0
end_model=24
num_na=65
num_mg=1
total_cl=22

# Parse command line arguments
while getopts "l:s:e:n:m:c:h" opt; do
  case $opt in
    l) ligand="$OPTARG" ;;
    s) start_model="$OPTARG" ;;
    e) end_model="$OPTARG" ;;
    n) num_na="$OPTARG" ;;
    m) num_mg="$OPTARG" ;;
    c) total_cl="$OPTARG" ;;
    h)
      echo "Usage: $0 [-l ligand] [-s start_model] [-e end_model] [-n num_na] [-m num_mg] [-c total_cl]"
      echo ""
      echo "Arguments:"
      echo "  -l  Ligand name (default: PLACEHOLDER_LIGAND)"
      echo "  -s  Starting model index (default: 0)"
      echo "  -e  Ending model index (default: 24)"
      echo "  -n  Number of Na+ ions (default: 65)"
      echo "  -m  Number of Mg2+ ions (default: 1)"
      echo "  -c  Total number of Cl- ions (default: 22)"
      echo "      (Cl- are divided between NaCl and MgCl2 additions)"
      exit 0
      ;;
    *) echo "Usage: $0 [-l ligand] [-s start_model] [-e end_model] [-n num_na] [-m num_mg] [-c total_cl]"
       exit 1 ;;
  esac
done

# Convert ligand to uppercase for directory/file naming
ligand_upper=$(echo "$ligand" | tr '[:lower:]' '[:upper:]')

# Calculate Cl- for MgCl2 and remaining for NaCl
# Each Mg2+ requires 2 Cl- → MgCl2
cl_mg=$((2 * num_mg))
cl_na=$((total_cl - cl_mg))

# Create leap directory and output subdirectory if they don't exist
leap_subdir="leap/CSS1_${ligand_upper}_constrained"
mkdir -p "$leap_subdir"

# Create a single tleap input file
cat > "${leap_subdir}/tleap.in" << EOF
source leaprc.gaff2
source leaprc.DNA.OL21
source leaprc.water.opc
loadAmberPrep ff/${ligand}_resp.prepin
loadAmberParams ff/${ligand}_resp.frcmod
EOF

for i in $(seq "$start_model" "$end_model"); do
  cat >> "${leap_subdir}/tleap.in" << EOF
oligo${i} = loadPDB pdb_for_md/CSS1_${ligand_upper}_constrained/input_model_${i}.pdb
solvateOct oligo${i} OPCBOX 14.0
addIonsRand oligo${i} Na+ ${num_na} Cl- ${cl_na}
addIonsRand oligo${i} MG ${num_mg} Cl- ${cl_mg}
saveamberparm oligo${i} ${leap_subdir}/cplx_${i}.prmtop ${leap_subdir}/cplx_${i}.rst7
EOF
done

# add quit at the end
echo "quit" >> "${leap_subdir}/tleap.in"

# Run tleap once on the full input
tleap -f "${leap_subdir}/tleap.in"

for i in $(seq "$start_model" "$end_model"); do
  ambpdb -p "${leap_subdir}/cplx_${i}.prmtop" -c "${leap_subdir}/cplx_${i}.rst7" > "${leap_subdir}/cplx_${i}.pdb"
done