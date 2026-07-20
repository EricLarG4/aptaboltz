#!/bin/bash
#===============================================================================
# ligand_prep.sh — AmberTools GAFF2 parameterisation pipeline for ligands
#
# DESCRIPTION
#   Takes an optimised .mol2 file (or .xyz as fallback) from the QM/ directory
#   and runs the full AmberTools parameterisation pipeline:
#     0. (Fallback) Convert .xyz to .mol2 via antechamber if .mol2 not found.
#     1. Assign GAFF2 atom types and BCC charges with antechamber.
#     2. Replace BCC charges with RESP2 charges from QM (_opt.chg file).
#     3. Generate an Amber prepin library file with prepgen.
#     4. Check missing force-field parameters with parmchk2.
#
#   Output files are organised into antechamber/, ff/, ff/prepin/, and
#   ff/parmchk2/ subdirectories.
#
# DEPENDENCIES
#   - AmberTools (antechamber, prepgen, parmchk2)
#   - OpenBabel (optional, only if converting .xyz → .mol2 via obabel)
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
#     QM/${MOL}.mol2       — optimised ligand structure (preferred, see Step 0)
#     QM/${MOL}_opt.chg    — RESP2 charges (from multiwfn_steps*)
#
#   If QM/${MOL}.mol2 is NOT found, the script attempts to fall back to
#     QM/${MOL}.xyz        — optimised ligand in XYZ format (from QM package)
#   and converts it to mol2 using either antechamber or OpenBabel.
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
#   0. QM geometry optimisation (external): produces QM/${MOL}.xyz
#   1. ORCA single-point:                   orca_steps_wsl.sh
#   2. RESP2 charge fitting:                multiwfn_steps_windows.ps1 or multiwfn_steps_linux.sh
#   3. Ligand parameterisation:             ligand_prep.sh (this script)
#   4. Model preparation for MD:            model_prep.py
#   5. System solvation + ionisation:       leap.sh
#   6. MD production:                       (user-defined)
#
# TROUBLESHOOTING
#   - If both .mol2 and .xyz are missing, the script exits with an error.
#   - Check the output of parmchk2 for ATTN warnings and non-zero penalty terms.
#   - The RESP2 charge substitution (Step 2) relies on the same atom ordering
#     in both the .mol2 and the .chg file.  Verify alignment if charges look
#     incorrect.
#   - The .xyz geometry must be the same optimised geometry used for the
#     ORCA single-point calculations (orca_steps_wsl.sh).
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
    echo -e "  ${YELLOW}${MOL}.mol2 not found in QM/ — trying .xyz fallback...${RST}"
    xyz_file=$(find ../QM -maxdepth 1 -iname "${MOL}.xyz" | head -1)
    if [ -z "$xyz_file" ]; then
        echo -e "\033[1;31mERROR:${RST} Neither ${MOL}.mol2 nor ${MOL}.xyz found in QM/"
        echo -e "  ${DIM}Provide an optimised structure as QM/${MOL}.mol2 or QM/${MOL}.xyz${RST}"
        exit 1
    fi
    echo -e "  ${DIM}Found ${xyz_file} — converting to mol2...${RST}"
    # Try antechamber first (AmberTools); fall back to OpenBabel
    if command -v antechamber &>/dev/null; then
        antechamber -fi xyz -i "$xyz_file" -fo mol2 -o "${MOL}.mol2" -at gaff2 -c none -s 0 -nc $NC
        mol2_file="${MOL}.mol2"
    elif command -v obabel &>/dev/null; then
        obabel "$xyz_file" -O "${MOL}.mol2"
        mol2_file="${MOL}.mol2"
    else
        echo -e "\033[1;31mERROR:${RST} Need antechamber or obabel to convert .xyz to .mol2"
        exit 1
    fi
    echo -e "  ${GREEN}Converted ${xyz_file} → ${mol2_file}${RST}"
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
python templates/MD/python/fix_pdb_elements.py ff/${MOL}_resp.pdb
echo -e "  ${DIM}• ff/${MOL}_resp.pdb (for complex preparation in tLeap)${RST}"

header "Step 3/4: Generating Amber prepin library file" $GREEN

prepgen -i antechamber/${MOL}_resp.ac -o ff/${MOL}_resp.prepin -rn ${MOL^^}
mv -f PREP.INF NEWPDB.PDB ff/prepin/ 2>/dev/null
python templates/MD/python/fix_pdb_elements.py ff/prepin/NEWPDB.PDB
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