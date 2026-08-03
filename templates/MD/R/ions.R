#===============================================================================
# ions.R — Ion concentration calculator for leap.sh
#
# DESCRIPTION
#   Calculates the number of Na+, Mg2+, and Cl- ions needed to reach the
#   target ionic concentrations in a solvated MD system, given the water-box
#   volume and net system charge. Results are displayed as DT tables.
#
#   Two methods are computed:
#     - SPLIT:  N0 = Nw * conc / 55.5  (classic approximation)
#     - SLTCAP: self-consistent solution (Machado + SLTCAP paper)
#               https://doi.org/10.1021/acs.jctc.7b01254
#   The SLTCAP method is preferred; it does not require N0 >> |Q|.
#
# DEPENDENCIES
#   data.table, DT
#
# USAGE
#   1. Edit the variables below (C_NaCl, C_MgCl2, box.vol.angst, Nw, Q) to
#      match your system.  Obtain Nw and box.vol.angst from the leap.log
#      output after running leap.sh once with approximate ion counts.
#   2. Run in R:
#        source("ions.R")
#   3. Use the resulting n_Na / n_Mg / n_Cl counts in leap.sh:
#        ./leap.sh -n <n_Na> -m <n_Mg> -c <n_Cl>
#
# REFERENCES
#   - Matías Machado method: http://archive.ambermd.org/202002/0194.html
#   - SLTCAP: https://doi.org/10.1021/acs.jctc.7b01254
#===============================================================================

library(data.table)
library(DT)

# ─── EDIT THESE VALUES FOR YOUR SYSTEM ─────────────────────────────────
C_NaCl <- 0.14      # Target NaCl concentration (Molar)
C_MgCl2 <- 0.01     # Target MgCl2 concentration (Molar)
box.vol.angst <- 492689.905  # Volume of the box (cubic Angstroms, from leap.log)
Nw <- 14206         # Number of water molecules (from leap.log)
Q <- -45            # Net system charge (e.g. -45 for a 46-mer DNA aptamer)
# ────────────────────────────────────────────────────────────────────────


salt_dt <- data.table(
  conc = C_NaCl,
  box.vol.angst = box.vol.angst,
  Nw = Nw,
  Q = Q
) |>
  _[, `:=`(
    # "avogadro's method"
    box.vol.liter = box.vol.angst * 1e-27, # convert to liters
    # box.vol.reduced = box.vol.angst/(3.16655 * 3.16655 * 3.16655), # Lennard-Jones sigma parameter for OPC water model: https://pubs.acs.org/doi/10.1021/jz501780a
    # method of Matias Machado: http://archive.ambermd.org/202002/0194.html
    N0 = Nw * conc / 55.5 # 56 in original method
  )] |>
  _[, `:=`(
    test.SPLIT = (N0 > 10 * abs(Q)), # test hypothesis
    sodium.SPLIT = round(N0 - Q / 2, 0), # number of sodium ions
    chloride.SPLIT = round(N0 + Q / 2, 0) # number of chloride ions
  )] |>
  _[, `:=`(
    verif.SPLIT = (sodium.SPLIT * 1 + chloride.SPLIT * -1 + Q == 0), # N0 >> Q
    # SLTCAP does not require N0 >> Q: https://doi.org/10.1021/acs.jctc.7b01254
    sodium.SLTCAP = round(N0 * sqrt(1 + (Q / (2 * N0))^2) - Q / 2, 0),
    chloride.SLTCAP = round(N0 * sqrt(1 + (Q / (2 * N0))^2) + Q / 2, 0)
  )] |>
  _[, `:=`(
    verif.SLTCAP = (sodium.SLTCAP * 1 + chloride.SLTCAP * -1 + Q == 0),
    eq.SPLIT.SLTCAP = (sodium.SPLIT == sodium.SLTCAP &
      chloride.SPLIT == chloride.SLTCAP)
  )]

datatable(salt_dt)

# Number of Na and Cl corresponding to the NaCl concentration:
n_NaCl <- salt_dt$chloride.SLTCAP

# By proportion, the number of Mg and Cl2 corresponding to the MgCl2
# concentration:
n_Mg <- round(C_MgCl2 * n_NaCl / C_NaCl, 0)
n_Cl <- 2 * n_Mg

# Total number of chloride ions:
n_Cl <- n_Cl + n_NaCl

datatable(
  data.table(
    n_Na = salt_dt$sodium.SLTCAP,
    n_Mg = n_Mg,
    n_Cl = n_Cl
  )
)
