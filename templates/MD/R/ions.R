#===============================================================================
# ions.R — Ion concentration calculator for leap.sh
#
# DESCRIPTION
#   Calculates the number of Na+, Mg2+, and Cl- ions needed to reach the
#   target ionic concentrations in a solvated MD system, given the water-box
#   volume and net system charge.  The results are printed to the console in
#   a format that tells you exactly what flags to pass to leap.sh.
#
#   Two methods are computed:
#     - SPLIT:  N0 = Nw * conc / 55.5  (classic approximation)
#     - SLTCAP: self-consistent solution (Machado + SLTCAP paper)
#               https://doi.org/10.1021/acs.jctc.7b01254
#
#   The SLTCAP method is preferred; it does not require N0 >> |Q|.
#
# DEPENDENCIES
#   Base R only — no external packages required.
#
# USAGE
#   1. Edit the four variables below (C_NaCl, C_MgCl2, Nw, Q) to match
#      your system.  Obtain Nw and box.vol.angst from the leap.log output
#      after running leap.sh once with approximate ion counts.
#   2. Run in R:
#        source("ions.R")
#   3. Use the printed leap.sh command in your leap.sh invocation.
#
# INPUT VARIABLES (EDIT THESE)
#   C_NaCl        — Target NaCl concentration (Molar). Default: 0.14
#   C_MgCl2       — Target MgCl2 concentration (Molar). Default: 0.01
#   Nw            — Number of water molecules in the solvated box.
#                   Obtain from leap.log: "Added X water molecules."
#   Q             — Net charge of the system (from phosphate backbone, etc.)
#                   e.g. -45 for a 46-mer DNA aptamer.
#
# METHOD
#   N0  = Nw * C_NaCl / 55.5
#   Na⁺ = N0 * sqrt(1 + (Q / (2 * N0))^2) - Q / 2    (SLTCAP)
#   Cl⁻ = N0 * sqrt(1 + (Q / (2 * N0))^2) + Q / 2    (SLTCAP)
#
#   Mg²⁺ and its accompanying Cl⁻ are then scaled proportionally from
#   the NaCl count by the concentration ratio C_MgCl2 / C_NaCl.
#
#   Final Cl⁻ = NaCl_Cl⁻ + 2 * Mg²⁺
#
# REFERENCES
#   - Matías Machado method: http://archive.ambermd.org/202002/0194.html
#   - SLTCAP: https://doi.org/10.1021/acs.jctc.7b01254
#===============================================================================

# ─── EDIT THESE VALUES FOR YOUR SYSTEM ──────────────────────────────
C_NaCl  <- 0.14  # Molar concentration of NaCl
C_MgCl2 <- 0.01  # Molar concentration of MgCl2
Nw      <- 14206 # Number of water molecules (from leap.log)
Q       <- -45   # Net system charge
# ─────────────────────────────────────────────────────────────────────

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  Ion Concentration Calculator                           \n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("  Target NaCl  concentration:  %.2f M\n", C_NaCl))
cat(sprintf("  Target MgCl2 concentration:  %.2f M\n", C_MgCl2))
cat(sprintf("  Water molecules (Nw):        %d\n", Nw))
cat(sprintf("  Net charge (Q):              %d\n", Q))
cat("───────────────────────────────────────────────────────────────\n")

# N0 = Nw * C_NaCl / 55.5 (water molarity at STP)
N0 <- Nw * C_NaCl / 55.5

# SLTCAP method (preferred)
sqrt_term <- sqrt(1 + (Q / (2 * N0))^2)
na_sltcap <- round(N0 * sqrt_term - Q / 2, 0)
cl_sltcap <- round(N0 * sqrt_term + Q / 2, 0)

# SPLIT method (for comparison)
na_split <- round(N0 - Q / 2, 0)
cl_split <- round(N0 + Q / 2, 0)

# Mg2+ scaled from NaCl by concentration ratio
n_mg <- round(C_MgCl2 * cl_sltcap / C_NaCl, 0)
n_cl_mg <- 2 * n_mg
n_cl_total <- cl_sltcap + n_cl_mg

cat("\n")
cat("  ── Results (SLTCAP method) ──\n")
cat(sprintf("  Na+  (NaCl):   %d\n", na_sltcap))
cat(sprintf("  Mg2+ (MgCl2):  %d\n", n_mg))
cat(sprintf("  Cl-  (total):  %d  (%d from NaCl + %d from MgCl2)\n",
    n_cl_total, cl_sltcap, n_cl_mg))

# Verify charge balance
charge_check <- na_sltcap * 1 + n_mg * 2 + n_cl_total * (-1) + Q
cat(sprintf("  Charge check:  %d  (should be 0)\n", charge_check))

cat("\n")
cat("  ── leap.sh command ──\n")
cat(sprintf("  leap.sh -n %d -m %d -c %d\n", na_sltcap, n_mg, n_cl_total))
cat("\n")
cat("  ── SPLIT method (for reference) ──\n")
cat(sprintf("  Na+  (NaCl):   %d\n", na_split))
cat(sprintf("  Cl-  (NaCl):   %d\n", cl_split))
cat("═══════════════════════════════════════════════════════════\n")
cat("\n")
