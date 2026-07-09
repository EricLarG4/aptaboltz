#===============================================================================
# multiwfn_steps_windows.ps1 — RESP2 charge calculation with Multiwfn
#
# DESCRIPTION
#   Reads gas-phase and solvent-phase Molden files (produced by
#   orca_steps_wsl.sh on WSL/Linux) and runs Multiwfn to compute
#   RESP charges for each phase.  It then combines the two charge
#   sets into RESP2 charges using the mixing parameter delta.
#
#   Adapted from RESP2_ORCA.sh by Tian Lu (sobereva@sina.com).
#
# DEPENDENCIES
#   - Multiwfn (tested with v2026.7.5)
#   - PowerShell (Windows)
#
# USAGE
#   .\multiwfn_steps_windows.ps1 input_name [-delta <value>]
#
# PARAMETERS
#   input_name  — Base name of the input molecule (e.g. "PLACEHOLDER_LIGAND").
#                 Used to name the output .chg file (e.g. PLACEHOLDER_LIGAND.chg).
#                 This should match the name used in ligand_prep.sh.
#   -delta      — RESP2 mixing parameter (default: 0.5).
#                 RESP2_charge = (1 - delta) * gas_charge + delta * solv_charge.
#                 delta = 0.5 gives equal weighting (RESP2-0.5).
#
# INPUT FILES REQUIRED (paths are hardcoded — edit for your system)
#   <MoldenDir>/SP_gas.molden    — Gas-phase Molden file (from WSL)
#   <MoldenDir>/SP_solv.molden   — Solvent-phase Molden file (from WSL)
#
# OUTPUT
#   <OutDir>/gas.chg              — Gas-phase RESP charges
#   <OutDir>/solv.chg             — Solvent-phase RESP charges
#   <OutDir>/<input_name>.chg     — Combined RESP2 charges
#                                   (last column = RESP2(delta) charges)
#
# WORKFLOW (full pipeline order)
#   1. QM optimisation + RESP2 charges:
#      orca_steps_wsl.sh       (WSL)
#      → multiwfn_steps_windows.ps1 (Windows — this script)
#   2. Ligand parameterisation:
#      ligand_prep.sh          (WSL)
#   3. System solvation + ionisation:
#      leap.sh                 (WSL)
#   4. MD production:
#      (user-defined)
#
# STEP BY STEP
#   1. Run orca_steps_wsl.sh on WSL to produce SP_gas.molden and SP_solv.molden.
#   2. Copy the two .molden files from WSL to <MoldenDir> on Windows.
#   3. Edit the hardcoded paths below for your system.
#   4. Run this script:
#      .\multiwfn_steps_windows.ps1 hcy
#   5. Copy <input_name>.chg back to QM/ on WSL for ligand_prep.sh.
#
# TROUBLESHOOTING
#   - Ensure the Molden files exist at the specified paths.
#   - Multiwfn must be able to run from cmd.exe via the batch wrapper.
#   - The script writes a temporary batch file (_run_multiwfn.bat) and input
#     file (_multiwfn_input.txt) which are deleted after execution.
#   - If Multiwfn exits with an error, check the host output for details.
#   - The RESP2 charge substitution assumes identical atom ordering between
#     the gas and solvent Molden files.  Verify alignment if charges appear
#     incorrect.
#
# EXAMPLE
#   # Default delta=0.5
#   .\multiwfn_steps_windows.ps1 hcy
#
#   # Custom delta
#   .\multiwfn_steps_windows.ps1 hcy -delta 0.6
#===============================================================================

param(
    [Parameter(Mandatory=$true)]
    [string]$input_name,
    [double]$delta = 0.5
)

# --- HARDCODED PATHS & SETTINGS (EDIT FOR YOUR SYSTEM) ---
$Multiwfn   = "E:\sDrive\Projects\boltz\CSS\MD\Multiwfn_PLACEHOLDER_VERSION_bin_Win64\Multiwfn.exe"
$MoldenDir  = "E:\sDrive\Projects\boltz\CSS\MD\QM"
$OutDir     = "E:\sDrive\Projects\boltz\CSS\MD\QM"
# ---------------------------------------------------------

$gasMolden  = Join-Path $MoldenDir "SP_gas.molden"
$solvMolden = Join-Path $MoldenDir "SP_solv.molden"

if (-not (Test-Path $gasMolden)) {
    Write-Error "SP_gas.molden not found in $MoldenDir. Copy from WSL first."
    exit 1
}
if (-not (Test-Path $solvMolden)) {
    Write-Error "SP_solv.molden not found in $MoldenDir. Copy from WSL first."
    exit 1
}

Write-Host "delta parameter is $delta"

function Run-Multiwfn {
    param([string]$molden, [string]$label)
    Write-Host "Running Multiwfn for $label..."

    # Write input commands to a temp file (CRLF for Windows cmd.exe)
    $inputFile = Join-Path $OutDir "_multiwfn_input.txt"
    $inputContent = "`r`n7`r`n18`r`n1`r`ny`r`n0`r`n0`r`nq`r`n"
    Set-Content -Path $inputFile -Value $inputContent -Encoding Ascii

    # Create a temporary batch file using here-string
    $batFile = Join-Path $OutDir "_run_multiwfn.bat"
    $batContent = @"
@echo off
cd /d "$OutDir"
"$Multiwfn" "$molden" -ispecial 1 < "$inputFile"
"@
    Set-Content -Path $batFile -Value $batContent -Encoding Ascii

    # Execute the batch file via cmd.exe, send output to host (not output stream)
    cmd.exe /c "$batFile" | Out-Host

    Remove-Item $inputFile -Force
    Remove-Item $batFile -Force -ErrorAction SilentlyContinue

    # Multiwfn names output after the input file (e.g. SP_gas.molden → SP_gas.chg)
    $expected = [System.IO.Path]::ChangeExtension([System.IO.Path]::GetFullPath($molden), ".chg")
    if (-not (Test-Path $expected)) {
        Write-Error "Multiwfn did not produce $expected for $label phase."
        exit 1
    }
    return $expected
}

### Gas phase
$gasChg = Run-Multiwfn -molden $gasMolden -label "gas"
$gasDest = Join-Path $OutDir "gas.chg"
Move-Item $gasChg $gasDest -Force
Write-Host "Gas-phase RESP charges written to $gasDest"

### Solvent phase
$solvChg = Run-Multiwfn -molden $solvMolden -label "solvent"
$solvDest = Join-Path $OutDir "solv.chg"
Move-Item $solvChg $solvDest -Force
Write-Host "Solvent-phase RESP charges written to $solvDest"

### Combine into RESP2 charges
$gasLines  = Get-Content $gasDest
$solvLines = Get-Content $solvDest
$outFile   = Join-Path $OutDir "$input_name.chg"

$lines = for ($i = 0; $i -lt $gasLines.Count; $i++) {
    $g = $gasLines[$i] -split '\s+', 5
    $s = $solvLines[$i] -split '\s+', 5
    if ($g.Count -ge 5 -and $s.Count -ge 5) {
        $chg = (1.0 - $delta) * [double]$g[4] + $delta * [double]$s[4]
        "{0,-3} {1,12:F6} {2,12:F6} {3,12:F6} {4,15:F10}" -f $g[0], [double]$g[1], [double]$g[2], [double]$g[3], $chg
    }
}
$lines | Set-Content $outFile -Encoding Ascii

Write-Host ""
Write-Host "Finished! RESP2 charges written to $outFile"
Write-Host "The last column contains the RESP2($delta) charges."

### Optionally clean up intermediate files
# Remove-Item $gasDest, $solvDest