#!/bin/bash
#===============================================================================
# orca_steps_wsl.sh — ORCA single-point calculations (gas + solvent CPCM-SMD)
#
# DESCRIPTION
#   Performs two consecutive single-point energy calculations with ORCA:
#     1. Gas phase (B3LYP-D3/def2-TZVP with RIJCOSX)
#     2. Solvent phase (same level, CPCM-SMD continuum solvation)
#
#   Both calculations write Molden-format files with an Nval valence-electron
#   header (for def2 ECP basis sets), ready for RESP charge fitting with
#   Multiwfn (see multiwfn_steps_windows.ps1).
#
#   Adapted from RESP2_ORCA.sh by Tian Lu (sobereva@sina.com).
#
# DEPENDENCIES
#   - ORCA (tested with v6.x)
#   - orca_2mkl utility (included with ORCA)
#   - Bash (WSL/Linux)
#
# USAGE
#   ./orca_steps_wsl.sh input.xyz [charge] [multiplicity] [solvent]
#
# ARGUMENTS
#   input.xyz     — Molecular structure file (XYZ format).
#   charge        — Net charge of the molecule (default: 0).
#   multiplicity  — Spin multiplicity (default: 1).
#   solvent       — Solvent name for SMD (default: "Water").
#                    Names must match ORCA's SMDsolvent keyword.
#
# OUTPUT
#   SP_gas.molden   — Gas-phase Molden file with Nval header (for Multiwfn).
#   SP_solv.molden  — Solvent-phase Molden file with Nval header.
#
# WORKFLOW
#   1. Copy the optimised ligand XYZ from QM/ directory (e.g. QM/hcy_opt.xyz).
#   2. Run this script on WSL/Linux.
#   3. Copy SP_gas.molden and SP_solv.molden to the Windows QM/ directory.
#   4. Run multiwfn_steps_windows.ps1 to compute RESP2 charges.
#
# NOTES
#   - The script uses B3LYP-D3/def2-TZVP with RIJCOSX by default.  Adjust
#     the $keyword_SP variable for other levels of theory.
#   - The Nval.txt header is required for correct valence-electron counts
#     when using def2 ECP basis sets with Multiwfn.
#   - Solvent-phase SP reuses the gas-phase converged WF as initial guess.
#   - If a single-point task fails, the script exits immediately.
#
# EXAMPLE
#   # Default: neutral singlet in water
#   ./orca_steps_wsl.sh hcy_opt.xyz
#
#   # Charged, different solvent
#   ./orca_steps_wsl.sh ligand_opt.xyz -1 1 "Ethanol"
#
# TROUBLESHOOTING
#   - Check SP.out for ORCA errors if the script exits with a failure message.
#   - Ensure ORCA binaries are on $PATH or update the $ORCA variable below.
#   - For large systems, increase $nprocs and $maxcore.
#===============================================================================

# --- HARDCODED PATHS & SETTINGS (EDIT FOR YOUR SYSTEM) ---
ORCA="/home/PLACEHOLDER_USER/orca_PLACEHOLDER_VERSION/orca"
orca_2mkl="/home/PLACEHOLDER_USER/orca_PLACEHOLDER_VERSION/orca_2mkl"
nprocs=16
maxcore=32768
keyword_SP="! B3LYP/G D3 def2-TZVP def2/J RIJCOSX"
# ---------------------------------------------------------

# --- RESP2 mixing coefficient (delta). Used in multiwfn_steps_windows.ps1.
#     This value is passed to the PowerShell script, not used directly here,
#     but recorded in this template for traceability.
delta=0.5

export inname=$1
filename=${inname%.*}
suffix=${inname##*.}

if [ $2 ];then
    echo "Net charge = $2"
    chg=$2
else
    echo "Net charge was not defined. Default to 0"
    chg=0
fi

if [ $3 ];then
    echo "Spin multiplicity = $3"
    multi=$3
else
    echo "Spin multiplicity was not defined. Default to 1"
    multi=1
fi

if [ "$4" ];then
    echo Solvent is $4
    solvent=$4
else
    echo "Solvent name was not defined. Default to water"
    solvent="Water"
fi

echo delta parameter is $delta

### Nval.txt: valence electron count for def2 ECP basis sets
cat << NVALEOF > Nval.txt
[Nval]
Rb  9
Sr 10
Y  11
Zr 12
Nb 13
Mo 14
Tc 15
Ru 16
Rh 17
Pd 18
Ag 19
Cd 20
In 21
Sn 22
Sb 23
Te 24
I  25
Xe 26
Cs  9
Ba 10
La 11
Ce 30
Pr 31
Nd 32
Pm 33
Sm 34
Eu 35
Gd 36
Tb 37
Dy 38
Ho 39
Er 40
Tm 41
Yb 42
Lu 43
Hf 12
Ta 13
W  14
Re 15
Os 16
Ir 17
Pt 18
Au 19
Hg 20
Tl 21
Pb 22
Bi 23
Po 24
At 25
Rn 26
NVALEOF


### Run gas-phase single point
cat << SPEOF > SP.inp
$keyword_SP
%maxcore $maxcore
%pal nprocs $nprocs end
SPEOF
echo "* xyz $chg $multi" >> SP.inp
awk '{if (NR>2) print }' $1 >> SP.inp
echo "*" >> SP.inp

echo
echo Running single point task in gas via ORCA...
$ORCA SP.inp > SP.out

if grep -Fq "ORCA TERMINATED NORMALLY" SP.out
then
    echo Done!
else
    echo The single point task has failed! Please check content of SP.out to find reason
    echo The script is terminated
    mv SP.out tmp.out
    rm SP.* SP_*
    mv tmp.out SP.out
    exit 1
fi

echo Running orca_2mkl...
$orca_2mkl SP -molden > /dev/null
cat Nval.txt SP.molden.input > SP_gas.molden
echo "Gas-phase molden file written: SP_gas.molden"
rm SP.molden.input


### Run solvent-phase single point (uses gas converged WF as guess)
cat << SPEOF > SP.inp
$keyword_SP
%maxcore $maxcore
%pal nprocs $nprocs end
%cpcm
smd true
SMDsolvent "$solvent"
end
SPEOF
echo "* xyz $chg $multi" >> SP.inp
awk '{if (NR>2) print }' $1 >> SP.inp
echo "*" >> SP.inp

echo
echo Running single point task in solvent via ORCA...
$ORCA SP.inp > SP.out

if grep -Fq "ORCA TERMINATED NORMALLY" SP.out
then
    echo Done!
else
    echo The single point task has failed! Please check content of SP.out to find reason
    echo The script is terminated
    mv SP.out tmp.out
    rm SP.* SP_*
    mv tmp.out SP.out
    exit 1
fi

echo Running orca_2mkl...
$orca_2mkl SP -molden > /dev/null
cat Nval.txt SP.molden.input > SP_solv.molden
echo "Solvent-phase molden file written: SP_solv.molden"
rm SP.molden.input


echo
echo "=== ORCA steps complete ==="
echo "Output files:"
echo "  SP_gas.molden   (gas phase, with Nval header)"
echo "  SP_solv.molden  (solvent phase, with Nval header)"
echo
echo "Next: copy these .molden files to Windows and run multiwfn_steps_windows.ps1"