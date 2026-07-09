#!/bin/bash
#===============================================================================
# ligand_prep.sh — AmberTools GAFF2 parameterisation pipeline for ligands
#
# DESCRIPTION
#   Takes an optimised .mol2 file from the QM/ directory and runs the full
#   AmberTools parameterisation pipeline:
#     1. Assign GAFF2 atom types and BCC charges with antechamber.
#     2. Replace BCC charges with RESP2 charges from QM (opt.chg file).
#     3. Generate an Amber prepin library file with prepgen.
#     4. Check missing force-field parameters with parmchk2.
#
#   Output files are organised into antechamber/, ff/, ff/prepin/, and
#   ff/parmchk2/ subdirectories.
#
# DEPENDENCIES
#   - AmberTools (antechamber, prepgen, parmchk2)
#   - Bash (Linux/WSL)
#
# USAGE
#   ./ligand_prep.sh
#
#   Before running, edit the per-project settings at the top of this file:
#     MOL  — residue/ligand name (default: "PLACEHOLDER_LIGAND")
#     NC   — net charge (default: 0)
#     S    — antechamber verbosity (0=silent, 2=verbose)
#
#   The script expects:
#     QM/${MOL}.mol2    — optimised ligand structure (output from QM ORCA/Multiwfn)
#     QM/${MOL}_opt.chg — RESP2 charges (output from multiwfn_steps_windows.ps1)
#
# OUTPUT STRUCTURE
#   antechamber/
#     ${MOL}.ac           — GAFF2 atom types with BCC charges
#     ${MOL}_resp.ac      — same, with RESP2 charges substituted
#     (auxiliary files from antechamber)
#   ff/
#     ${MOL}_resp.pdb     — PDB-like file for complex preparation in tLeap
#     ${MOL}_resp.prepin  — Amber prepin library file
#     ${MOL}_resp.frcmod  — force-field modification parameters
#   ff/prepin/            — PREP.INF and NEWPDB.PDB (auxiliary)
#   ff/parmchk2/          — ANTECHAMBER.FRCMOD and ATOMTYPE.INF (auxiliary)
#
# WORKFLOW (full pipeline order)
#   1. QM optimisation + RESP2 charges: orca_steps_wsl.sh → multiwfn_steps_windows.ps1
#   2. Ligand parameterisation:          ligand_prep.sh (this script)
#   3. System solvation + ionisation:    leap.sh
#   4. MD production:                    (user-defined)
#
# TROUBLESHOOTING
#   - If the .mol2 or .chg files are missing, the script exits with an error.
#   - Check the output of parmchk2 for ATTN warnings and non-zero penalty terms.
#   - The RESP2 charge substitution (Step 2) relies on the same atom ordering
#     in both the .mol2 and the .chg file.  Verify alignment if charges look
#     incorrect.
#
# EXAMPLE
#   # Prepare ligand "hcy" (default), neutral, with verbose output
#   ./ligand_prep.sh
#
#   # Prepare a different ligand
#   # Edit MOL="my_ligand" and NC=0 in the User settings section, then run.
#===============================================================================

# ─── User settings (EDIT FOR YOUR LIGAND) ─────────────────────
# MOL must match the residue/ligand name used in the complex structure file (PDB/CIF).
# It also determines the input file names: QM/${MOL}.mol2 and QM/${MOL}_opt.chg.
MOL="PLACEHOLDER_LIGAND"
S=2       # antechamber verbosity (0=silent, 2=verbose)
NC=0      # net charge
# ────────────────────────────────────────────────────────────────

# ─── Formatting helpers ──────────────────────────────────────
BOLD='\033[1m'
DIM='\033[2m'
RED='\033[1;31m'
ORANGE='\033[1;33m'
YELLOW='\033[1;93m'
GREEN='\033[1;32m'
BLUE='\033[1;34m'
MAGENTA='\033[1;35m'
CYAN='\033[1;36m'
RST='\033[0m' # No Color

header() {
  local msg="$1"
  local color="${2:-$BLUE}"
  local line
  printf -v line '%*s' 80 ''
  line="${line// /═}"
  echo -e "\n${color}${line}${RST}"
  echo -e "${color}  ${msg}${RST}"
  echo -e "${color}${line}${RST}\n"
}

# ─────────────────────────────────────────────────────────────
mkdir -p antechamber ff ff/prepin ff/parmchk2

header "Ligand parameterisation pipeline" $MAGENTA
echo -e "  ${BOLD}Ligand:${RST}  ${MOL}"
echo -e "  ${BOLD}Verbose:${RST} ${S}"
echo -e "  ${BOLD}Charge:${RST}   ${NC}\n"

header "Step 1/4: Assigning GAFF2 atom types and BCC charges with antechamber" $RED
cd antechamber

mol2_file=$(find ../QM -maxdepth 1 -iname "${MOL}.mol2" | head -1)
if [ -z "$mol2_file" ]; then
    echo -e "\033[1;31mERROR:${RST} ${MOL}.mol2 not found in QM/"
    exit 1
fi

antechamber \
-fi mol2 \
-i "$mol2_file" \
-fo ac \
-o ${MOL}.ac \
-at gaff2 \
-c bcc \
-s $S \
-nc $NC

cd ../
mv -f ANTECHAMBER.FRCMOD ATOMTYPE.INF antechamber/ 2>/dev/null
echo -e "  ${DIM}• antechamber/${MOL}.ac${RST}"
echo -e "  ${DIM}• antechamber/ (auxiliary files)${RST}"

header "Step 2/4: Replacing BCC charges with RESP2 charges from QM" $ORANGE

awk '{gsub(/\r/,"")} NR==FNR {charge[NR]=$5; next} /^ATOM/ {sub(/[^ ]+ +[^ ]+$/, charge[++i] " " $NF)} 1' \
  QM/${MOL}_opt.chg antechamber/${MOL}.ac > antechamber/${MOL}_resp.ac
echo -e "  ${DIM}• antechamber/${MOL}_resp.ac${RST}"
cp antechamber/${MOL}_resp.ac ff/${MOL}_resp.pdb
echo -e "  ${DIM}• ff/${MOL}_resp.pdb (for complex preparation in tLeap)${RST}"

header "Step 3/4: Generating Amber prepin library file" $GREEN

prepgen -i antechamber/${MOL}_resp.ac -o ff/${MOL}_resp.prepin -rn ${MOL^^}
mv -f PREP.INF NEWPDB.PDB ff/prepin/ 2>/dev/null
echo -e "  ${DIM}• ff/${MOL}_resp.prepin${RST}"
echo -e "  ${DIM}• ff/prepin/ (auxiliary files)${RST}"

header "Step 4/4: Checking missing force field parameters with parmchk2" $BLUE

parmchk2 -i ff/${MOL}_resp.prepin -f prepi -o ff/${MOL}_resp.frcmod -a Y -s gaff2
mv -f ANTECHAMBER.FRCMOD ATOMTYPE.INF ff/parmchk2/ 2>/dev/null
echo -e "  ${DIM}• ff/${MOL}_resp.frcmod${RST}"
echo -e "  ${DIM}• ff/parmchk2/ (auxiliary files)${RST}"

header "frcmod issues (ATTN or non-zero penalty)" $MAGENTA

awk '/ATTN/ || (/penalty/ && $NF !~ /^0\.0+$/)' ff/${MOL}_resp.frcmod

header "Done! Output files in antechamber/ and ff/" $CYAN