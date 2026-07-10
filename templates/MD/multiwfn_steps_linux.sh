#!/bin/bash
#===============================================================================
# multiwfn_steps_linux.sh — RESP2 charge calculation with Multiwfn (Linux)
#
# DESCRIPTION
#   Linux/Bash equivalent of multiwfn_steps_windows.ps1.  Reads gas-phase and
#   solvent-phase Molden files (produced by orca_steps_wsl.sh) and runs
#   Multiwfn to compute RESP charges for each phase.  It then combines the
#   two charge sets into RESP2 charges using the mixing parameter delta.
#
#   NOTE: Using Multiwfn on Windows (multiwfn_steps_windows.ps1) is
#   recommended for the latest developments.  This Linux script is provided
#   as a convenience for single-OS setups.
#
# DEPENDENCIES
#   - Multiwfn (Linux binary, e.g. Multiwfn_Linux_noGUI)
#   - Bash (Linux/WSL)
#
# USAGE
#   ./multiwfn_steps_linux.sh <input_name> [delta]
#
# ARGUMENTS
#   input_name  — Base name of the input molecule (e.g. "PLACEHOLDER_LIGAND").
#                 Used to name the output .chg file (e.g. PLACEHOLDER_LIGAND.chg).
#                 This should match the name used in ligand_prep.sh.
#   delta       — RESP2 mixing parameter (default: 0.5).
#                 RESP2_charge = (1 - delta) * gas_charge + delta * solv_charge.
#                 delta = 0.5 gives equal weighting (RESP2-0.5).
#                 Must match the delta used in orca_steps_wsl.sh.
#
# INPUT FILES REQUIRED (paths are editable — see CONFIGURATION below)
#   <MolDir>/SP_gas.molden     — Gas-phase Molden file (from orca_steps_wsl.sh)
#   <MolDir>/SP_solv.molden    — Solvent-phase Molden file (from orca_steps_wsl.sh)
#
# OUTPUT
#   <MolDir>/gas.chg            — Gas-phase RESP charges
#   <MolDir>/solv.chg           — Solvent-phase RESP charges
#   <MolDir>/<input_name>.chg   — Combined RESP2 charges
#                                 (last column = RESP2(delta) charges)
#
# WORKFLOW (full pipeline order)
#   1. QM geometry optimisation (external):          produces QM/${MOL}.xyz
#   2. ORCA single-point calculations:                orca_steps_wsl.sh
#   3. RESP2 charge fitting:                          multiwfn_steps_linux.sh (this script)
#   4. Ligand parameterisation:                       ligand_prep.sh
#   5. Model preparation for MD:                      model_prep.py
#   6. System solvation + ionisation:                 leap.sh
#   7. MD production:                                 (user-defined)
#
# STEP BY STEP
#   1. Run orca_steps_wsl.sh on WSL/Linux to produce SP_gas.molden and SP_solv.molden.
#   2. Edit the paths below if Molden files are in a different directory.
#   3. Run this script:
#        ./multiwfn_steps_linux.sh PLACEHOLDER_LIGAND
#   4. Use <input_name>.chg as input for ligand_prep.sh.
#
# TROUBLESHOOTING
#   - Ensure the Molden files exist at the specified paths.
#   - Multiwfn must be on $PATH or update the $Multiwfn path below.
#   - The RESP2 charge substitution in ligand_prep.sh assumes identical atom
#     ordering between .mol2 and .chg files.  Verify alignment if charges
#     appear incorrect.
#   - On Linux, Multiwfn accepts input commands via stdin piping; no batch
#     file wrapper is needed (unlike the Windows version).
#
# EXAMPLE
#   # Default delta=0.5
#   ./multiwfn_steps_linux.sh hcy
#
#   # Custom delta
#   ./multiwfn_steps_linux.sh hcy 0.6
#===============================================================================

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — Edit for your system                                    ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
Multiwfn="/path/to/Multiwfn_Linux_noGUI"   # ← Path to Multiwfn executable
MolDir="QM"                                # ← Directory with .molden files
OutDir="QM"                                # ← Output directory for .chg files
# ═════════════════════════════════════════════════════════════════════════════

input_name="${1:?Usage: $0 <input_name> [delta]}"
delta="${2:-0.5}"

gasMolden="${MolDir}/SP_gas.molden"
solvMolden="${MolDir}/SP_solv.molden"

if [ ! -f "$gasMolden" ]; then
    echo "ERROR: SP_gas.molden not found in ${MolDir}. Run orca_steps_wsl.sh first."
    exit 1
fi
if [ ! -f "$solvMolden" ]; then
    echo "ERROR: SP_solv.molden not found in ${MolDir}. Run orca_steps_wsl.sh first."
    exit 1
fi

echo "delta parameter is $delta"

run_multiwfn() {
    local molden="$1"
    local label="$2"
    echo "Running Multiwfn for $label..."

    # Multiwfn input:
    #   (empty line for initial prompt)
    #   7  → RESP charge fitting
    #   18 → Read from Molden file
    #   1  → CHELPG grid
    #   y  → Symmetry
    #   0  → Default
    #   0  → Quit current function
    #   q  → Quit Multiwfn
    echo -e "\n7\n18\n1\ny\n0\n0\nq" | "$Multiwfn" "$molden" -ispecial 1 > /dev/null 2>&1

    local expected="${molden%.*}.chg"
    if [ ! -f "$expected" ]; then
        echo "ERROR: Multiwfn did not produce $expected for $label phase."
        exit 1
    fi
    echo "$expected"
}

### Gas phase
gasChg=$(run_multiwfn "$gasMolden" "gas")
mv -f "$gasChg" "${OutDir}/gas.chg"
echo "Gas-phase RESP charges written to ${OutDir}/gas.chg"

### Solvent phase
solvChg=$(run_multiwfn "$solvMolden" "solvent")
mv -f "$solvChg" "${OutDir}/solv.chg"
echo "Solvent-phase RESP charges written to ${OutDir}/solv.chg"

### Combine into RESP2 charges
output="${OutDir}/${input_name}.chg"
awk -v delta="$delta" '
    NR == FNR {gas_charge[FNR] = $5; next}
    {
        chg = (1.0 - delta) * gas_charge[FNR] + delta * $5
        printf "%-3s %12.6f %12.6f %12.6f %15.10f\n", $1, $2, $3, $4, chg
    }
' "${OutDir}/gas.chg" "${OutDir}/solv.chg" > "$output"

echo ""
echo "Finished! RESP2 charges written to $output"
echo "The last column contains the RESP2($delta) charges."

### Optionally clean up intermediate files
# rm -f "${OutDir}/gas.chg" "${OutDir}/solv.chg"
