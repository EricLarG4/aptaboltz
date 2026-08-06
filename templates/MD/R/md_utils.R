#===============================================================================
# md_utils.R — MD trajectory analysis utilities
#
# DESCRIPTION
#   Utility functions for processing and visualizing molecular dynamics
#   simulation outputs. Currently includes:
#     - Reading 2D RMSD matrices from cpptraj's 2drms command
#     - Reading per-residue atomic fluctuation data from cpptraj's
#       atomicfluct command
#     - Custom ggplot2 themes for publication-quality figures
#     - Plotting 2D RMSD heatmaps and atomic fluctuation profiles
#     - Tracking ligand-DNA H-bond contacts over time (presence heatmaps,
#       summary tables)
#     - Tracking ligand-DNA pi-stacking interactions over time (presence
#       heatmaps, summary tables), driven by the external
#       templates/MD/python/pi_stacking.py precomputation script
#     - Reading per-model constraint verification CSVs for minimised
#       structures (written by process_minimized_constraint_verification
#       in templates/MD/python/pymol_utils.py)
#     - Extracting density time series from Amber production .out files
#       and assessing density equilibration via a single-exponential
#       plateau fit (Bogetti et al., J. Chem. Phys. 153, 054123 (2020))
#
# DEPENDENCIES
#   data.table, stringr, ggplot2, ggrepel, patchwork, ggtext, ggpattern
#   (jsonlite is used by pdb_pi_stacking() via the :: qualifier)
#
# FUNCTIONS
#   clean_dna_residue(res)
#     Strips the one-letter DNA prefix (D) and underscores from
#     cpptraj-style residue labels: "DG_21" -> "G21".
#
#   element_span(atom)
#     Wraps an atom label in a coloured HTML span (CPK-like element
#     colours) for DT tables and ggtext markdown axis labels.
#
#   ring_label(ring)
#     Maps ligand ring/unit names (Q0-Q1 fused-ring units by default, or
#     R0-R3 with --no-simplify) to descriptive names (e.g. "Q1"). The
#     mapping in ring_labels is ligand-specific — edit it for your ligand.
#     Unmapped names are returned unchanged.
#
#   read_rmsd(file, simulation_time, extract_frequency, time_step)
#     Reads a 2D RMSD .dat file into a long-format data.table with
#     time coordinates (t1, t2) and RMSD values. Optionally subsamples
#     to a coarser time resolution via time_step.
#
#   read_atomicfluct(file)
#     Reads an atomic fluctuation .dat file into a data.table with
#     residue numbers, fluctuation values, and metadata.
#
#   theme_custom(scaling)
#     Returns a ggplot2 theme with configurable font sizes via a
#     scaling factor. Removes grid lines and uses a white background
#     for clean publication output.
#
#   plot_rmsd(rmsd, max_t, time_step, scale)
#     Generates a 2D RMSD heatmap (geom_raster) faceted by experiment
#     and model, using viridis color scale. time_step controls the
#     time resolution (ns) for plotting.
#
#   plot_atomfluct(atomfluct, res_cutoff, label_cutoff, scale)
#     Plots per-residue atomic fluctuation as a line/point chart with
#     automatic labeling of high-fluctuation regions via ggrepel.
#
#   read_contacts_avg(file, ligand_residue)
#     Reads cpptraj average-geometry output (hbond avgout) for a run,
#     keeping only intermolecular ligand-DNA contacts.
#
#   read_contacts_matrix(file, ligand_residue, simulation_time,
#                        extract_frequency)
#     Reads a per-frame ligand-DNA contact count matrix into a
#     long-format data.table with time coordinates.
#
#   plot_contacts(dt_avg, dt_ts, max_t, scale, min_frac)
#     Two-panel contact plot: occupancy (fraction of frames) and a
#     time series of the ligand-DNA contact count.
#
#   read_contacts_long(file, ligand_residue, simulation_time,
#                      extract_frequency)
#     Reads a per-frame H-bond contact matrix (cpptraj contactseries)
#     into a long-format data.table, one row per contact per frame.
#
#   pdb_hbond_contacts(pdb_path)
#     Evaluates the H-bond set of a single structure (e.g. the final
#     minimised structure) in situ with the same geometric criteria as
#     cpptraj.
#
#   plot_contacts_tracked(dt_long, ref, max_t, scale, heat_bin_ns)
#     A per-contact presence heatmap for the reference H-bond sets, one
#     panel per model. Tile colour encodes presence and reference-set
#     membership (absent, initial-only, final-only, both). Rows are labelled
#     "<DNA residue>@<DNA atom>\u00b7<ligand atom>" with the ligand atom
#     coloured by element. Returns a ggplot object.
#
#   contacts_summary(ref, avg)
#     Summary data.table of the reference contacts with set membership,
#     occupancy and mean H-bond geometry, for display (e.g. DT).
#     Model column is labelled "model N".
#
#   read_pi_stacking(file, stride, simulation_time)
#     Reads an external pi-stacking CSV (written by
#     templates/MD/python/pi_stacking.py, copied to the project) into a
#     long-format data.table, one row per pair per analysed frame, with
#     absent frames zero-padded.
#
#   pdb_pi_stacking(pdb, ring_defs_file, ligand_resname, ...)
#     Evaluates the pi-stacking set of a single structure (e.g. the final
#     minimised structure) in situ from a PDB and ring_defs.json, using the
#     same geometric criteria as the Python script.
#
#   plot_pi_stacking_tracked(dt_long, ref, max_t, scale, heat_bin_ns)
#     Per-pair presence heatmap for the reference pi-stacking sets,
#     mirroring plot_contacts_tracked(). Returns a ggplot object.
#
#   pi_stacking_summary(ref, dt_long)
#     Summary data.table of the reference pi-stacking pairs with set
#     membership, occupancy, mean geometry and dominant stacking mode.
#
#   read_minimized_constraints(project)
#     Reads all "*_minimized_constraints.csv" files under
#     "{project}/MD/pmemd/out/" into a single data.table with sequence,
#     condition and job metadata, keeping only the most recent job per
#     experiment.
#
#   read_density(file)
#     Extracts the instantaneous-density time series from an Amber pmemd
#     production .out file (one point per NSTEP block, stopping before
#     the block-average section) into a data.table with t_ns, t_ps and
#     seq/ligand/experiment/job/model metadata.
#
#   density_plateau(dt, slope_cut, df_cut, chi2_cut)
#     Fits the density series to D(t) = D_i + (D_f - D_i)(1 - exp(-k t))
#     per group (Bogetti et al. 2020 / cpptraj evalplateau) and tests the
#     plateau criteria (final slope, |D_f - second-half mean|, reduced
#     chi-squared). Returns one summary row per group.
#
#   plot_density(dt, plateau, max_t, scale)
#     Line chart of the density series faceted by model, with the fitted
#     plateau exponential overlaid.
#
# USAGE
#   source("md_utils.R")
#
#   # Then use in a processing script:
#   rmsd_dt    <- lapply(rmsd_files, read_rmsd) |> rbindlist()
#   fluct_dt   <- lapply(fluct_files, read_atomicfluct) |> rbindlist()
#   plot_rmsd(rmsd_dt, max_t = 100, time_step = 1, scale = 0.7)
#   plot_atomfluct(fluct_dt, res_cutoff = 46, label_cutoff = 4)
#
# FILE STRUCTURE ASSUMPTIONS
#   The file path passed to read_rmsd() and read_atomicfluct() is expected
#   to follow the pattern:
#     <project>/<category>/<method>/out/<experiment>_<replicate>/<model>/<file>.dat
#   where:
#     - path component [5]  = experiment identifier (e.g. "CSS1_rep1")
#     - path component [7]  = model name (e.g. "wildtype")
#   The experiment name is split on '_' and joined with a middle dot (·).
#
# NOTES
#   - Frame-to-time conversion assumes uniform extraction frequency.
#   - simulation_time and extract_frequency default to global variables
#     sim_time and extract_freq if not provided explicitly.
#   - The plot functions use rleid() to assign subgroup IDs for faceting,
#     so models within each experiment are stacked vertically.
#===============================================================================

# ─── DEPENDENCIES ──────────────────────────────────────────────────────
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggtext)
library(ggpattern)


# ─── LABEL HELPERS ─────────────────────────────────────────────────────
# clean_dna_residue() strips the one-letter DNA prefix (D) and any
# underscores from cpptraj-style residue labels: "DG_21" -> "G21",
# "DC3_46" -> "C346". Handles both the H-bond ("DG_21") and the
# pi-stacking ("DA10") label conventions.
clean_dna_residue <- function(res) {
  gsub("_", "", sub("^D", "", res))
}

# pad_dna_residue() appends a non-breaking space (U+00A0) to single-digit
# residue numbers so that stacked "Q1\u2225C9" labels align vertically with
# multi-digit ones ("Q1\u2225G21"). A regular space would be stripped by the
# ggtext markdown parser; NBSP renders identically and survives.
pad_dna_residue <- function(res) {
  ifelse(grepl("^[A-Z][0-9]$", res), paste0(res, "\u00A0"), res)
}

# CPK-like element colours for colouring ligand atoms by element (N blue,
# O red, ...). Used for HTML spans in DT tables and ggtext axis labels.
element_colour <- c(H = "#FFFFFF", C = "#909090", N = "#3050F8",
                    O = "#FF0D0D", P = "#FF8000", S = "#FFFF30",
                    Cl = "#1FF01F", BR = "#A62929", I = "#940094")

# element_span(atom) wraps an atom label in a coloured HTML span
# (e.g. "O2" -> "<span style='color:#FF0D0D'>O2</span>"). Valid for both
# DT (escape = FALSE) and ggtext element_markdown().
element_span <- function(atom) {
  el <- toupper(sub("[0-9]+$", "", atom))
  col <- unname(element_colour[el])
  col[is.na(col)] <- "#909090"
  paste0("<span style='color:", col, ";'>", atom, "</span>")
}

# Descriptive labels for the piperaquine aromatic systems. With the default
# --simplify output of pi_stacking.py, fused aromatic rings are merged into
# single units Q0, Q1, ... (the two quinoline units map to Q1/Q2); the
# per-ring names R0-R3 (--no-simplify) are also mapped, each quinoline
# contributing a pyridine ring (pyr) and a benzene ring (benz): R0/R1 =
# quinoline Q1, R2/R3 = quinoline Q2.
ring_labels <- c(Q0 = "Q1", Q1 = "Q2",
                 R0 = "Q1\u00b7pyr", R1 = "Q1\u00b7benz",
                 R2 = "Q2\u00b7pyr", R3 = "Q2\u00b7benz")

ring_label <- function(ring) {
  unname(ifelse(ring %in% names(ring_labels), ring_labels[ring], ring))
}

# Separator for pi-stacking pair labels (deliberately distinct from the \u00b7
# used for H-bond contacts): two parallel bars evoke a stacked pair, e.g.
# "Q1\u2225A10".
stack_sep <- "\u2225"


# ═══════════════════════════════════════════════════════════════════════════
# I/O FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

# ─── read_rmsd ──────────────────────────────────────────────────────────
# Reads a 2D RMSD matrix file (from cpptraj's 2drms command) and converts
# it from a wide-format triangular matrix into a long-format data.table.
#
# The raw file has frames as both rows and columns, with the first column
# being the row frame index. After melting, each row becomes a pairwise
# comparison (frame_1 vs frame_2) with its RMSD value.
#
# The 2drms command in cpptraj computes pairwise RMSD between all frames
# in a trajectory. The atom mask used (!:Na+,Cl-,MG) excludes ions and
# the output format is a space-separated matrix with "#Frame" as the
# header for the first column.
#
# Frame indices are converted to time (ns) using:
#   t = frame * simulation_time / total_frames
# where total_frames = number_of_rows * extract_frequency.
#
# Args:
#   file              — Path to the 2D RMSD .dat file.
#   simulation_time   — Total simulation length in ns (default: sim_time).
#   extract_frequency — Frame extraction interval (default: extract_freq).
#   time_step         — Subsample to this time resolution (ns). The
#                       nearest available frames are kept. NA = no
#                       filtering (default).
#
# Returns:
#   A data.table with columns: t1, t2, rmsd, model, experiment.
read_rmsd <- function(
  file,
  simulation_time = sim_time,
  extract_frequency = extract_freq,
  time_step = NA
) {
  # Extract experiment name from the 5th path component
  # e.g. "CSS/MD/pmemd/out/CSS1_rep1/..." → "CSS1_rep1" → "CSS1·rep1"
  experiment <- strsplit(file, split = "[/\\\\]")[[1]][5]
  sequence <- strsplit(experiment, '_')[[1]][1]
  ligand <- strsplit(experiment, '_')[[1]][2]
  experiment <- paste0(
    sequence,
    "·",
    ligand
  )

  # Extract model name from the 7th path component
  # e.g. ".../CSS1_rep1/wildtype/2drmsd.dat" → "wildtype"
  model <- strsplit(file, split = "[/\\\\]")[[1]][7]

  cat("Reading", model, "for experiment", experiment, "\n")

  # Read the wide-format matrix, melt to long format, and rename columns
  rmsd_dt <- fread(file) |>
    melt(
      id.vars = "#Frame"
    ) |>
    setNames(c("frame_1", "frame_2", "rmsd")) |>
    _[, frame_2 := as.integer(frame_2)] |>
    _[, `:=`(
      frame_1 = (frame_1 - 1) * extract_frequency + 1,
      frame_2 = (frame_2 - 1) * extract_frequency + 1
    )]

  # Total number of frames in the trajectory (for time conversion)
  total_frames <- length(unique(rmsd_dt$frame_1)) * extract_frequency

  # Convert frame numbers to time points and attach metadata
  rmsd_dt[,
    `:=`(
      t1 = frame_1 * simulation_time / total_frames,
      t2 = frame_2 * simulation_time / total_frames,
      model = model,
      experiment = experiment,
      seq = sequence,
      ligand = ligand
    )
  ]

  # Optionally subsample to a coarser time resolution
  if (!is.na(time_step)) {
    unique_t <- sort(unique(rmsd_dt$t1))
    dt <- median(diff(unique_t))
    n <- max(1, round(time_step / dt))
    keep_t <- unique_t[seq(1, length(unique_t), by = n)]
    cat(
      "Time step filter: requested",
      time_step,
      "ns,",
      "data resolution",
      round(dt, 4),
      "ns,",
      "keeping every",
      n,
      "th frame (effective:",
      round(n * dt, 4),
      "ns)\n"
    )
    rmsd_dt <- rmsd_dt[t1 %in% keep_t & t2 %in% keep_t]
  }

  rmsd_dt
}


# ─── read_atomicfluct ───────────────────────────────────────────────────
# Reads an atomic fluctuation file (from cpptraj's atomicfluct command)
# and returns a tidy data.table.
#
# The cpptraj atomicfluct command with the "byres" keyword computes
# per-residue atomic fluctuation (RMSF) values. The atom mask used
# (!:Na+,Cl-,MG) excludes ions from the calculation.
#
# The raw file has columns: #Res (residue number), AtomicFlx (fluctuation
# in Angstroms). Metadata (experiment, model) is extracted from the file
# path using the same convention as read_rmsd().
#
# Args:
#   file — Path to the atomicfluct .dat file.
#
# Returns:
#   A data.table with columns: res, atomfluct, experiment, model.
read_atomicfluct <- function(file) {
  # Extract experiment name from path component [5]
  experiment <- strsplit(file, split = "[/\\\\]")[[1]][5]
  sequence <- strsplit(experiment, '_')[[1]][1]
  ligand <- strsplit(experiment, '_')[[1]][2]

  experiment <- paste0(
    sequence,
    "·",
    ligand
  )

  # Extract model name from path component [7]
  model <- strsplit(file, split = "[/\\\\]")[[1]][7]

  cat("Reading", model, "for experiment", experiment, "\n")

  # Read and attach metadata, rename columns for clarity
  fread(file) |>
    _[, experiment := experiment] |>
    _[, model := model] |>
    _[, seq := sequence] |>
    _[, ligand := ligand] |>
    _[, .(res = `#Res`, atomfluct = AtomicFlx, experiment, model, seq, ligand)]
}


# ═══════════════════════════════════════════════════════════════════════════
# PLOTTING UTILITIES
# ═══════════════════════════════════════════════════════════════════════════

# ─── theme_custom ───────────────────────────────────────────────────────
# A reusable ggplot2 theme optimized for publication figures.
#
# Strips all background/grid elements and applies consistent typography
# with a global scaling factor. Designed for faceted plots with
# experiment on the x-axis and model subgroups on the y-axis.
#
# Args:
#   scaling — Multiplier for all font sizes and line widths.
#             Default 1 = base sizes; 0.7 for smaller figures, etc.
#
# Returns:
#   A ggplot2 theme object.
theme_custom <- function(scaling = 1) {
  theme(
    axis.line = element_line(linewidth = 0.75 * scaling, color = 'black'),
    axis.ticks = element_line(linewidth = 0.75 * scaling, color = 'black'),
    axis.title = element_text(size = 20 * scaling, face = 'bold'),
    axis.text = element_text(size = 16 * scaling),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18 * scaling, face = 'bold'),
    strip.text.y = element_text(
      size = 18 * scaling,
      face = 'bold',
      angle = -90
    ),
    legend.title = element_text(
      size = 18 * scaling,
      face = 'bold',
      angle = -90
    ),
    legend.text = element_text(size = 16 * scaling),
    legend.key.width = unit(10, 'pt'),
    legend.key.height = unit(25, 'pt'),
    legend.background = element_blank(),
    legend.position = 'right',
    plot.background = element_rect(fill = 'white', colour = NA),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )
}


# ─── plot_rmsd ──────────────────────────────────────────────────────────
# Generates a 2D RMSD heatmap using geom_raster().
#
# The plot shows time-vs-time with RMSD as the fill color (viridis scale).
# Models are faceted vertically within each experiment using rleid() to
# assign subgroup IDs based on the ordering of model names.
#
# Args:
#   rmsd       — A data.table from read_rmsd() (must have t1, t2, rmsd,
#                 model, experiment columns).
#   max_t      — Maximum time (ns) for the x-axis. NA = auto (default).
#   time_step  — Plot at this time resolution (ns). The nearest available
#                frames are kept. NA = use all frames (default).
#   scale      — Scaling factor passed to theme_custom().
#
# Returns:
#   A ggplot object.
plot_rmsd <- function(rmsd, max_t = NA, time_step = NA, scale = 1,
                      clean_labels = TRUE) {
  cat(
    "Plotting RMSD for sequence(s):",
    paste(unique(rmsd$seq), collapse = ", "),
    "\n"
  )

  # Optionally subsample to a coarser time resolution
  if (!is.na(time_step)) {
    unique_t <- sort(unique(rmsd$t1))
    dt <- median(diff(unique_t))
    n <- max(1, round(time_step / dt))
    keep_t <- unique_t[seq(1, length(unique_t), by = n)]
    cat(
      "Time step filter: requested",
      time_step,
      "ns,",
      "data resolution",
      round(dt, 4),
      "ns,",
      "keeping every",
      n,
      "th frame (effective:",
      round(n * dt, 4),
      "ns)\n"
    )
    rmsd <- rmsd[t1 %in% keep_t & t2 %in% keep_t]
  }

  # Clean experiment labels for facet headers
  if (isTRUE(clean_labels)) {
    rmsd[, experiment := gsub(".*·", "", experiment)]
    rmsd[experiment == "free", experiment := "Free"]
  }

  # Assign subgroup IDs for vertical faceting within each experiment
  rmsd[
    order(experiment, model),
    subgroup_id := rleid(model),
    by = experiment
  ] |>
    ggplot(aes(t1, t2, fill = rmsd)) +
    geom_raster() +
    theme_custom(scale) +
    scale_fill_viridis_c("RMSD (Å)", guide = guide_colorbar(title.vjust = 0)) +
    facet_grid(subgroup_id ~ experiment) +
    scale_x_continuous("t (ns)", expand = c(0, 0), limits = c(1, max_t)) +
    scale_y_continuous("t (ns)", expand = c(0, 0))
}


# ─── plot_atomfluct ─────────────────────────────────────────────────────
# Plots per-residue atomic fluctuation as a line chart with points.
#
# Residues with fluctuation above label_cutoff are annotated with their
# residue number using ggrepel to avoid overlapping labels. High-fluctuation
# regions are identified by grouping consecutive residues above the
# cutoff (using rleid on the difference between residue index and position).
#
# Only the peak residue (max fluctuation) in each high-fluctuation
# tract is labeled to reduce clutter.
#
# Args:
#   atomfluct      — A data.table from read_atomicfluct() with columns
#                    res, atomfluct, model, experiment.
#   res_cutoff     — Maximum residue number to display (default: 46).
#   label_cutoff   — Fluctuation threshold (Å) above which residues are
#                    labeled (default: 4).
#   scale          — Scaling factor passed to theme_custom().
#
# Returns:
#   A ggplot object.
plot_atomfluct <- function(
  atomfluct = atomicfluct_dt,
  res_cutoff = 46,
  label_cutoff = 4,
  scale = 1
) {
  # Assign subgroup IDs and filter to the residue range of interest
  atomfluct <- atomfluct[
    order(experiment, model),
    subgroup_id := rleid(model),
    by = experiment
  ] |>
    _[res <= res_cutoff]

  # Identify high-fluctuation residues and group consecutive ones into tracts
  high_fluct <- atomfluct[
    atomfluct > label_cutoff,
    .(res, atomfluct),
    by = .(model, experiment)
  ][
    order(model, experiment, res)
  ][,
    # Group consecutive residues above the cutoff into "tracts"
    # tract_id changes whenever the gap between residue numbers > 1
    tract_id := rleid(res - seq_len(.N)),
    by = .(model, experiment)
  ][,
    # Keep only the peak residue (max fluctuation) in each tract
    .SD[which.max(atomfluct)],
    by = .(model, experiment, tract_id)
  ][,
    subgroup_id := rleid(model),
    by = experiment
  ]

  # Build the plot: line + points, with repel labels for high-fluct regions
  atomfluct |>
    ggplot(
      aes(x = res, y = atomfluct)
    ) +
    geom_line(linewidth = 0.75, color = 'grey') +
    geom_point(aes(color = atomfluct), size = 4, show.legend = FALSE) +
    geom_text_repel(
      data = high_fluct,
      inherit.aes = FALSE,
      aes(label = res, x = res, y = atomfluct),
      color = "#525252",
      size = 6,
      fontface = 'bold',
      force = 50,
      direction = 'y',
      ylim = label_cutoff + 3
    ) +
    scale_x_continuous("") +
    scale_y_continuous("Atomic fluctuation (Å)") +
    scale_color_viridis_c() +
    theme_custom(scaling = scale) +
    facet_grid(subgroup_id ~ experiment)
}


# ─── read_contacts_avg ──────────────────────────────────────────────────
read_contacts_avg <- function(file, ligand_residue = "HCY_47") {
  experiment <- strsplit(file, split = "[/\\\\]")[[1]][5]
  sequence <- strsplit(experiment, "_")[[1]][1]
  lgd <- strsplit(experiment, "_")[[1]][2]
  model <- strsplit(file, split = "[/\\\\]")[[1]][7]

  cat("Reading", model, "for experiment",
      paste0(sequence, "\u00B7", lgd), "\n")

  dt <- fread(file, na.strings = "", strip.white = TRUE)
  setnames(dt, gsub("^#", "", names(dt)))

  # Keep only ligand-DNA (intermolecular), exclude ligand-ligand (intra)
  dt <- dt[(grepl(ligand_residue, Acceptor) |
            grepl(ligand_residue, Donor)) &
           !(grepl(ligand_residue, Acceptor) &
             grepl(ligand_residue, Donor))]
  if (nrow(dt) == 0) return(NULL)

  dt[, dna_residue := fifelse(
    grepl(ligand_residue, Acceptor),
    sub("@.*", "", Donor),
    sub("@.*", "", Acceptor))]

  dt[, dna_clean := gsub("_", "", dna_residue)]

  dt[, dna_atom := fifelse(
    grepl(ligand_residue, Acceptor),
    sub(".*@", "", Donor),
    sub(".*@", "", Acceptor))]

  dt[, bond_label := paste0(dna_clean, " ", dna_atom)]

  dt[, `:=`(
    model = model,
    experiment = paste0(sequence, "\u00B7", lgd),
    seq = sequence,
    ligand = lgd
  )]

  # Contact label in cpptraj convention: Acceptor-Donor-DonorH
  # (DonorH in the avg file carries its residue: HCY_47@H16 -> H16)
  dt[, contact := paste0(Acceptor, "-", Donor, "-",
                        sub(".*@", "", DonorH))]
  dt[, ligand_atom := fifelse(
    grepl(ligand_residue, Acceptor),
    sub(".*@", "", Acceptor),
    sub(".*@", "", Donor))]

  dt[]
}


# ─── read_contacts_matrix ───────────────────────────────────────────────
read_contacts_matrix <- function(file, ligand_residue = "HCY_47",
                                 simulation_time = sim_time,
                                 extract_frequency = extract_freq) {
  experiment <- strsplit(file, split = "[/\\\\]")[[1]][5]
  sequence <- strsplit(experiment, "_")[[1]][1]
  lgd <- strsplit(experiment, "_")[[1]][2]
  model <- strsplit(file, split = "[/\\\\]")[[1]][7]

  cat("Reading matrix", model, "for experiment",
      paste0(sequence, "\u00B7", lgd), "\n")

  dt <- fread(file, na.strings = "", strip.white = TRUE, header = TRUE)

  lgd_cols <- grep(ligand_residue, names(dt), value = TRUE)
  if (length(lgd_cols) == 0) return(NULL)

  n <- nrow(dt)
  mat <- as.matrix(dt[, ..lgd_cols])

  traj_frame <- (seq_len(n) - 1) * extract_frequency + 1
  data.table(
    frame = 1:n,
    time = traj_frame * simulation_time / (n * extract_frequency),
    n_contacts = rowSums(mat, na.rm = TRUE),
    model = model,
    experiment = paste0(sequence, "\u00B7", lgd),
    seq = sequence,
    ligand = lgd
  )
}


# ─── plot_contacts ──────────────────────────────────────────────────────
plot_contacts <- function(dt_avg, dt_ts, max_t = 100, scale = 1,
                          min_frac = 0.05) {
  dt_avg <- dt_avg[Frac >= min_frac]
  occ_p <- ggplot(dt_avg, aes(x = Frac, y = reorder(bond_label, Frac))) +
    geom_segment(aes(xend = 0, yend = bond_label),
                 color = "grey70", linewidth = 0.4) +
    geom_point(aes(color = dna_residue), size = 2) +
    facet_grid(model ~ ., scales = "free_y", space = "free_y") +
    theme_custom(scale) +
    scale_x_continuous("Fraction of frames",
                       limits = c(0, 1), expand = c(0.02, 0)) +
    scale_color_brewer(palette = "Paired") +
    labs(y = NULL, color = "DNA residue") +
    theme(axis.text.y = element_text(size = 8),
          strip.text.y = element_text(size = 7))

  ts_p <- ggplot(dt_ts, aes(x = time, y = n_contacts)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE,
                linewidth = 0.6, color = "#D55E00") +
    facet_grid(model ~ .) +
    theme_custom(scale) +
    scale_x_continuous("Time (ns)",
                       limits = c(0, max_t), expand = c(0, 0)) +
    scale_y_continuous("Ligand-DNA contacts",
                       expand = c(0.1, 0)) +
    theme(strip.text.y = element_text(size = 7))

  occ_p / ts_p + plot_layout(heights = c(2, 1))
}


# ═══════════════════════════════════════════════════════════════════════════
# LIGAND-DNA CONTACT TRACKING
#
# These functions support tracking a predefined set of ligand-DNA H-bond
# contacts over time, in order to detect (and date) changes in the ligand
# binding mode. The workflow is:
#   1. read_contacts_long()  — per-frame presence of every ligand-DNA H-bond
#                              (from cpptraj hbond series matrices).
#   2. pdb_hbond_contacts()  — H-bond set of a single structure (e.g. the
#                              final minimised structure) evaluated in situ,
#                              using the same geometric criteria as cpptraj.
#   3. Reference sets are built by the caller: "initial" contacts (present
#      at the start of the production run, within a tunable time window) and
#      "final" contacts (present in the final minimised structure).
#   4. plot_contacts_tracked()  — presence heatmap for the reference sets.
#   5. contacts_summary()       — table of reference contacts with occupancy
#                                 and mean geometry.
# ═══════════════════════════════════════════════════════════════════════════


# ─── read_contacts_long ─────────────────────────────────────────────────
# Reads a per-frame H-bond contact matrix (cpptraj `hbond` + `writedata
# contactseries`, e.g. step10_hbond_contacts_*.csv) into a long-format
# data.table. One row per contact per frame.
#
# The matrix header is "#Frame, contacts[UU], <Acceptor>-<Donor>-<DonorH>...",
# with 0/1 entries giving H-bond presence per stored frame. Only
# intermolecular ligand-DNA contacts are kept (ligand-ligand H-bonds such as
# HCY_47@O4-HCY_47@O2-H29 are excluded). The ligand residue token (e.g.
# "HCY_47", "PQ_48") is inferred from the column labels if not supplied.
#
# Args:
#   file              — Path to the contact matrix CSV.
#   ligand_residue    — Ligand residue token (e.g. "HCY_47"). NULL = infer.
#   simulation_time   — Total simulation length in ns (default: sim_time).
#   extract_frequency — Trajectory frame extraction interval (default:
#                       extract_freq).
#
# Returns:
#   A data.table with columns: contact, present, frame, time, ligand_atom,
#   dna_residue, dna_atom, model, experiment, seq, ligand.
read_contacts_long <- function(file, ligand_residue = NULL,
                               simulation_time = sim_time,
                               extract_frequency = extract_freq) {
  experiment <- strsplit(file, split = "[/\\\\]")[[1]][5]
  sequence <- strsplit(experiment, "_")[[1]][1]
  lgd <- strsplit(experiment, "_")[[1]][2]
  model <- strsplit(file, split = "[/\\\\]")[[1]][7]

  cat("Reading contact series", model, "for experiment",
      paste0(sequence, "\u00B7", lgd), "\n")

  dt <- fread(file, na.strings = "", strip.white = TRUE, header = TRUE)
  setnames(dt, gsub("^#", "", names(dt)))

  contact_cols <- setdiff(names(dt), c("Frame", "contacts[UU]"))
  if (length(contact_cols) == 0) return(NULL)

  # Infer the ligand residue token from the column labels if not supplied
  if (is.null(ligand_residue)) {
    tokens <- unique(unlist(strsplit(contact_cols, "-", fixed = TRUE)))
    tokens <- tokens[grepl("@", tokens)]
    res <- sub("@.*", "", tokens)
    lig <- unique(res[!grepl("^(DC|DA|DG|DT)", res)])
    if (length(lig) == 0) return(NULL)
    ligand_residue <- lig
  }

  long <- melt(dt, id.vars = "Frame", measure.vars = contact_cols,
               variable.name = "contact", value.name = "present")
  long[, contact := as.character(contact)]
  long[, c("acceptor", "donor", "donorH") :=
         tstrsplit(contact, "-", fixed = TRUE)]

  # Keep only intermolecular ligand-DNA contacts (XOR on ligand presence)
  long <- long[grepl(ligand_residue, acceptor) != grepl(ligand_residue, donor)]

  long[, `:=`(
    ligand_atom = fifelse(grepl(ligand_residue, acceptor),
                          sub(".*@", "", acceptor),
                          sub(".*@", "", donor)),
    dna_part = fifelse(grepl(ligand_residue, acceptor), donor, acceptor)
  )]
  long[, `:=`(
    dna_residue = sub("@.*", "", dna_part),
    dna_atom = sub(".*@", "", dna_part)
  )]
  long[, dna_part := NULL]

  # Time conversion (same convention as read_rmsd / read_contacts_matrix)
  n <- nrow(dt)
  long[, time := ((as.integer(Frame) - 1) * extract_frequency + 1) *
         simulation_time / (n * extract_frequency)]
  long[, frame := as.integer(Frame)]

  long[, `:=`(
    model = model,
    experiment = paste0(sequence, "\u00B7", lgd),
    seq = sequence,
    ligand = lgd
  )]
  long[, c("Frame", "acceptor", "donor", "donorH") := NULL]
  long[]
}


# ─── pdb_hbond_contacts ─────────────────────────────────────────────────
# Detects H-bonds in a single structure (PDB file) using the same geometric
# criteria as the cpptraj `hbond` analysis used for the trajectories
# (see CSS/MD/hbond_analysis.cpptraj and PQ4/MD/hbond_analysis.cpptraj):
#   - donor-acceptor distance <= dist_cut (default 3.2 A)
#   - D-H-A angle >= angle_cut (default 135 degrees)
# Hydrogen-to-donor bonding is assigned to the nearest heavy atom within
# 1.2 A. Contacts within a single residue are excluded (this removes
# ligand-ligand H-bonds).
#
# Args:
#   pdb              — Path to a PDB file (e.g. the final minimised
#                      structure, *_stripped.pdb).
#   dist_cut         — Donor-acceptor distance cutoff in A (default: 3.2).
#   angle_cut        — D-H-A angle cutoff in degrees (default: 135).
#   ligand_residue   — If supplied, keep only contacts in which exactly one
#                      side is this residue (e.g. "HCY_47", "PQ_48").
#
# Returns:
#   A data.table with columns: acceptor, donor, donorH, contact, ligand_atom,
#   dna_residue, dna_atom.
pdb_hbond_contacts <- function(pdb, dist_cut = 3.2, angle_cut = 135,
                               ligand_residue = NULL) {
  lines <- readLines(pdb, warn = FALSE)
  lines <- lines[grepl("^(ATOM  |HETATM)", lines)]
  if (length(lines) == 0) return(data.table(contact = character()))

  resnum <- as.integer(substr(lines, 23, 26))
  resname <- gsub(" ", "", substr(lines, 18, 20))
  atom <- gsub(" ", "", substr(lines, 13, 16))
  x <- as.numeric(substr(lines, 31, 38))
  y <- as.numeric(substr(lines, 39, 46))
  z <- as.numeric(substr(lines, 47, 54))
  elem <- gsub(" ", "", substr(lines, 77, 78))
  missing <- elem == ""
  elem[missing] <- substr(atom[missing], 1, 1)

  atoms <- data.table(atom = atom, resname = resname, resnum = resnum,
                      x = x, y = y, z = z, elem = elem)
  atoms[, resid := paste0(resname, "_", resnum)]

  heavy <- atoms[elem %in% c("C", "N", "O", "P", "S")]
  hydr <- atoms[elem == "H"]
  if (nrow(heavy) == 0 || nrow(hydr) == 0) {
    return(data.table(contact = character()))
  }

  # Assign each hydrogen to its nearest heavy atom (bonding cutoff 1.2 A)
  hx <- hydr$x; hy <- hydr$y; hz <- hydr$z
  gx <- heavy$x; gy <- heavy$y; gz <- heavy$z
  h_donor <- vector("integer", nrow(hydr))
  for (i in seq_len(nrow(hydr))) {
    dx <- hx[i] - gx; dy <- hy[i] - gy; dz <- hz[i] - gz
    j <- which.min(dx^2 + dy^2 + dz^2)
    if (sqrt(dx[j]^2 + dy[j]^2 + dz[j]^2) <= 1.2) h_donor[i] <- j
  }
  keep <- h_donor > 0L
  dh <- data.table(donor_idx = h_donor[keep], hydr_row = which(keep))
  if (nrow(dh) == 0) return(data.table(contact = character()))

  # Donors are N/O heavy atoms bearing at least one hydrogen
  dh[, donor_elem := heavy$elem[donor_idx]]
  dh <- dh[donor_elem %in% c("N", "O")]
  if (nrow(dh) == 0) return(data.table(contact = character()))

  # Acceptors are N/O heavy atoms
  acc <- heavy[elem %in% c("N", "O")]

  cos_cut <- cos(angle_cut * pi / 180)
  out <- vector("list", nrow(dh))
  for (i in seq_len(nrow(dh))) {
    di <- dh$donor_idx[i]
    hi <- dh$hydr_row[i]
    dx <- heavy$x[di] - hydr$x[hi]; dy <- heavy$y[di] - hydr$y[hi]
    dz <- heavy$z[di] - hydr$z[hi]
    ax <- acc$x - hydr$x[hi]; ay <- acc$y - hydr$y[hi]; az <- acc$z - hydr$z[hi]
    adx <- acc$x - heavy$x[di]; ady <- acc$y - heavy$y[di]; adz <- acc$z - heavy$z[di]
    d_h <- sqrt(dx^2 + dy^2 + dz^2)
    a_h <- sqrt(ax^2 + ay^2 + az^2)
    d_a <- sqrt(adx^2 + ady^2 + adz^2)
    cosv <- (dx * ax + dy * ay + dz * az) / (d_h * a_h)
    ok <- which(a_h > 0 & d_a <= dist_cut & cosv <= cos_cut &
                  acc$resid != heavy$resid[di])
    if (length(ok)) {
      out[[i]] <- data.table(
        acceptor = paste0(acc$resid[ok], "@", acc$atom[ok]),
        donor    = paste0(heavy$resid[di], "@", heavy$atom[di]),
        donorH   = hydr$atom[hi]
      )
    }
  }
  res <- rbindlist(out)
  if (is.null(res) || nrow(res) == 0) {
    return(data.table(contact = character()))
  }
  res[, contact := paste0(acceptor, "-", donor, "-", donorH)]

  if (!is.null(ligand_residue)) {
    res <- res[grepl(ligand_residue, acceptor) != grepl(ligand_residue, donor)]
  }
  if (nrow(res) == 0) return(data.table(contact = character()))

  res[, `:=`(
    ligand_atom = fifelse(grepl(ligand_residue, acceptor),
                          sub(".*@", "", acceptor),
                          sub(".*@", "", donor)),
    dna_part = fifelse(grepl(ligand_residue, acceptor), donor, acceptor)
  )]
  res[, `:=`(
    dna_residue = sub("@.*", "", dna_part),
    dna_atom = sub(".*@", "", dna_part)
  )]
  res[, dna_part := NULL]
  res[]
}


# ─── plot_contacts_tracked ──────────────────────────────────────────────
# Plots the time evolution of a predefined set of ligand-DNA contacts.
#
# A single per-contact presence heatmap is produced, with one panel per
# sequence model (all models are shown even if they have no reference
# contacts). Tile colour combines presence with reference-set membership:
# absent, present and initial-only, present and final-only, present in both
# sets. Rows are ordered final-set first, then initial-only.
#
# Row labels use "<DNA residue>@<DNA atom>\u00b7<ligand atom>" (e.g.
# "G39@O6\u00b7O2"), where the DNA residue is cleaned (DG_21 -> G21) and the
# ligand atom is coloured by element via ggtext markdown.
#
# Models are displayed as "model N" (the raw "task_N" names are kept
# internally in model_levels / row_key, so ordering and merges are
# unaffected).
#
# Args:
#   dt_long      — Long contact data from read_contacts_long().
#   ref          — Reference-set membership: one row per contact with columns
#                  contact, model, seq, experiment, in_initial (0/1),
#                  in_final (0/1), plus dna_residue/dna_atom/ligand_atom for
#                  labelling.
#   max_t        — Maximum time (ns) for the x-axis.
#   scale        — Scaling factor passed to theme_custom().
#   heat_bin_ns  — Time resolution (ns) used to bin the heatmap.
#
# Returns:
#   A ggplot object.
plot_contacts_tracked <- function(dt_long, ref, max_t = 100, scale = 1,
                                  heat_bin_ns = 0.1) {
  ref_sub <- ref[, .(contact, model, seq, experiment, in_initial, in_final,
                     dna_residue, dna_atom, ligand_atom)]
  # drop dna_residue/dna_atom/ligand_atom from the long data to avoid merge
  # suffix clash
  dt <- copy(dt_long)[, c("dna_residue", "dna_atom", "ligand_atom") := NULL]
  dt <- merge(dt, ref_sub, by = c("contact", "model", "seq", "experiment"))
  dt <- dt[in_initial == 1 | in_final == 1]
  if (nrow(dt) == 0) {
    msg <- paste0("No ligand-DNA H-bonds (3.2 A / 135 deg) were present in the ",
                  "initial state or the final minimised structure for these ",
                  "models.")
    empty_p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
      annotate("text", 0, 0, label = msg, size = 3, lineheight = 0.9) +
      theme_void() +
      theme(panel.border = element_rect(fill = NA, colour = "grey70"))
    return(empty_p)
  }
  dt[, bin := floor(time / heat_bin_ns) * heat_bin_ns]

  # Facet over every model of the sequence, even those without reference
  # contacts (shown as empty panels)
  model_levels <- unique(dt_long[, .(model, seq)])
  model_levels <- model_levels[order(as.integer(sub("task_", "", model)))]$model
  model_labels <- paste0("model ", sub("task_", "", model_levels))

  # ---- Presence heatmap ----
  # Tile colour combines presence with reference-set membership of the contact
  heat_agg <- dt[, .(present_frac = mean(present),
                     in_initial = max(in_initial),
                     in_final = max(in_final)),
                 by = .(contact, model, seq, bin)]
  heat_agg[, fill_key := fcase(
    present_frac >= 0.5 & in_initial == 1 & in_final == 1, "both",
    present_frac >= 0.5 & in_initial == 1, "initial",
    present_frac >= 0.5 & in_final == 1, "final",
    default = "absent")]

  # Row ordering: per model, final-set contacts first, then initial-only
  order_dt <- unique(dt[, .(contact, model, seq, in_initial, in_final,
                            dna_residue, dna_atom, ligand_atom)])
  order_dt[,
    row_key := paste0(model, " | ", contact)
  ][order(model, seq, -in_final, -in_initial, dna_residue, contact),
    row_level := seq_len(.N)
  ]

  heat_agg <- merge(heat_agg,
                    order_dt[, .(contact, model, seq, row_key)],
                    by = c("contact", "model", "seq"))

  # Row labels: "<cleaned DNA residue>@<DNA atom>\u00b7<ligand atom>" with the
  # ligand atom coloured by element (ggtext markdown)
  row_lab <- setNames(
    paste0(clean_dna_residue(order_dt$dna_residue), "@", order_dt$dna_atom,
           "\u00b7", element_span(order_dt$ligand_atom)),
    order_dt$row_key)

  # placeholder rows so models without reference contacts still get a panel
  missing_models <- setdiff(model_levels, unique(as.character(heat_agg$model)))
  dummy_keys <- character()
  if (length(missing_models) > 0) {
    dummy_keys <- paste0(missing_models, " | \u2014")
    heat_agg <- rbind(heat_agg, data.table(
      bin = NA_real_,
      model = missing_models,
      seq = unique(dt_long$seq),
      contact = NA_character_,
      present_frac = NA_real_,
      in_initial = 0L,
      in_final = 0L,
      fill_key = NA_character_,
      row_key = dummy_keys), fill = TRUE)
    row_lab <- c(row_lab, setNames(rep("", length(dummy_keys)), dummy_keys))
  }

  # ---- Minimised-structure membership column ----
  # One tile per row to the right of the trajectory showing whether the
  # contact is present in the final minimised structure (step-11 PDB),
  # independent of the trajectory occupancy used for the time bins.
  marker_x <- max_t + 1.5
  marker_w <- 0.8
  marker_dt <- order_dt[, .(model, seq, row_key, in_initial, in_final)]
  marker_dt[, `:=`(
    bin = marker_x,
    present_frac = NA_real_,
    fill_key = fifelse(in_final == 1, "min_yes", "min_no")
  )]
  heat_agg <- rbind(heat_agg, marker_dt, fill = TRUE)

  heat_agg[, row_key := factor(row_key,
                               levels = c(order_dt$row_key[order(order_dt$row_level)],
                                          dummy_keys))]
  heat_agg[, model := factor(model, levels = model_levels, labels = model_labels)]

  time_agg <- heat_agg[fill_key %in% c("min_yes", "min_no") == FALSE]
  marker_agg <- heat_agg[fill_key %in% c("min_yes", "min_no")]

  heat_p <- ggplot(heat_agg, aes(bin, row_key, fill = fill_key,
                                 pattern = fill_key)) +
    geom_tile_pattern(
      data = time_agg,
      width = heat_bin_ns * 1.02, height = 0.9, na.rm = TRUE,
      pattern_type = "stripe",
      pattern_colour = NA,
      pattern_fill = "#FFA6D9",
      pattern_fill2 = "#C9303E",
      pattern_density = 0.5,
      pattern_angle = 45
    ) +
    geom_tile_pattern(
      data = marker_agg,
      width = marker_w, height = 0.9, na.rm = TRUE,
      pattern_type = "none",
      pattern_colour = NA
    ) +
    geom_vline(xintercept = max_t, colour = "grey60",
               linewidth = 0.4, na.rm = TRUE) +
    facet_grid(model ~ ., scales = "free_y", space = "free_y", drop = FALSE) +
    theme_custom(scale) +
    scale_fill_manual(
      values = c(absent = "#F2F2F2", initial = "#FFA6D9",
                 final = "#C9303E", both = "#C9303E",
                 min_yes = "#C9303E", min_no = "#FFFFFF"),
      breaks = c("absent", "initial", "final", "both"),
      labels = c("absent", "present, initial", "present, minimized",
                 "present, initial & minimized"),
      name = NULL
    ) +
    scale_pattern_manual(
      values = c(absent = "none", initial = "none", final = "none",
                 both = "stripe", min_yes = "none", min_no = "none"),
      breaks = c("absent", "initial", "final", "both"),
      labels = c("absent", "present, initial", "present, minimized",
                 "present, initial & minimized"),
      name = NULL
    ) +
    scale_x_continuous("Time (ns)", breaks = seq(0, max_t, by = 25),
                       expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_t + 3.2)) +
    scale_y_discrete(labels = function(x) row_lab[x], name = NULL) +
    guides(fill = guide_legend(nrow = 2)) +
    theme(legend.position = "top",
          legend.key.height = unit(30, "pt"),
          legend.key.spacing.y = unit(0.4, "cm"),
          panel.spacing.y = unit(0.35, "lines"),
          axis.text.y = element_markdown(size = 14 * scale, face = "bold"),
          strip.text.y = element_text(size = 14 * scale, face = "bold"),
          strip.background = element_rect(fill = "grey92", colour = NA))

  heat_p
}


# ─── contacts_summary ───────────────────────────────────────────────────
# Builds a summary data.table of the reference contacts for display (e.g.
# with DT::datatable). One row per contact of the union of the initial and
# final reference sets, with set membership, occupancy over the full run and
# mean H-bond geometry (from the cpptraj average output).
#
# Args:
#   ref  — Reference-set membership (see plot_contacts_tracked).
#   avg  — Output of read_contacts_avg() for the same runs (may be NULL).
#
# Returns:
#   A data.table with columns: Model (labelled "model N"), Contact,
#   DNA residue (cleaned, e.g. "G39@O6"), Ligand atom (coloured by element
#   via an HTML span), Set, Occupancy, Mean D-A (A), Mean angle (degrees).
contacts_summary <- function(ref, avg = NULL) {
  summ <- copy(ref)
  summ <- summ[in_initial == 1 | in_final == 1]
  summ[, set := fifelse(in_initial == 1 & in_final == 1, "both",
                        fifelse(in_initial == 1, "initial", "minimized"))]

  if (!is.null(avg) && nrow(avg) > 0) {
    avg_sub <- avg[, .(contact, model,
                       occupancy = Frac, avg_dist = AvgDist, avg_ang = AvgAng)]
    summ <- merge(summ, avg_sub, by = c("contact", "model"), all.x = TRUE)
  } else {
    summ[, `:=`(occupancy = NA_real_, avg_dist = NA_real_, avg_ang = NA_real_)]
  }
  summ[is.na(occupancy), occupancy := 0]

  summ[, `:=`(
    Model = paste0("model ", sub("task_", "", model)),
    Contact = contact,
    `DNA residue` = paste0(clean_dna_residue(dna_residue), "@", dna_atom),
    `Ligand atom` = element_span(ligand_atom),
    Set = set,
    `Occupancy` = round(occupancy, 3),
    `Mean D-A (A)` = round(avg_dist, 2),
    `Mean angle (deg)` = round(avg_ang, 1)
  )]
  summ <- summ[order(model, -in_final, -in_initial, dna_residue)]
  summ[, .(Model, Contact, `DNA residue`, `Ligand atom`, Set,
           `Occupancy`, `Mean D-A (A)`, `Mean angle (deg)`)]
}


# ═══════════════════════════════════════════════════════════════════════════
# LIGAND-DNA PI-STACKING TRACKING
#
# These functions mirror the ligand-DNA H-bond tracking above, but for
# pi-stacking interactions between the aromatic rings of the ligand (only
# piperaquine, PQ, has aromatic rings in this project) and the aromatic base
# rings of the DNA. The trajectory data are precomputed externally by
# PQ4/MD/python/pi_stacking.py (MDAnalysis + RDKit; see the re-run note in
# the report), which writes one CSV per model with per-frame geometry for
# every ligand-ring x base-ring pair that comes within a buffer distance.
# The workflow is:
#   1. read_pi_stacking()        — per-frame presence/geometry of every pair
#                                  from the external CSVs, zero-padded to the
#                                  full analysed frame range.
#   2. pdb_pi_stacking()         — pi-stacking set of a single structure
#                                  (e.g. the final minimised structure)
#                                  evaluated in situ from a PDB and
#                                  ring_defs.json, using the same geometric
#                                  criteria as the Python script.
#   3. Reference sets are built by the caller exactly like the H-bonds:
#      "initial" pairs (present within a tunable window) and "final" pairs
#      (present in the final minimised structure).
#   4. plot_pi_stacking_tracked() — presence heatmap for the reference sets.
#   5. pi_stacking_summary()      — table of reference pairs with occupancy,
#                                  mean geometry and dominant stacking mode.
#
# Geometric criteria (matching the Python script):
#   parallel   : centroid distance <= 5.5 A, acute interplanar angle <= 30 deg,
#                lateral offset <= 2.0 A
#   T-shaped   : centroid distance <= 5.5 A, angle 60-90 deg
# A pair is "present" at a frame if it satisfies either mode.
# ═══════════════════════════════════════════════════════════════════════════


# ─── read_pi_stacking ───────────────────────────────────────────────────
# Reads an external pi-stacking CSV (written by PQ4/MD/python/pi_stacking.py)
# into a long-format data.table, one row per pair per analysed frame, with
# absent frames zero-padded so `present` is defined for every frame.
#
# The CSV files live in <project>/MD/pi_stacking/ and are named
# <experiment>_task_<n>.csv (e.g. "PQ4_PQ_constrained_task_0.csv"); the
# experiment and model (task_N) are parsed from the file name. Only pairs
# that come within the buffer distance get any rows, so padding is required
# to define absence over the full frame range.
#
# Args:
#   file              — Path to the pi-stacking CSV.
#   stride            — Analysis frame stride used by the Python script
#                       (default: 5).
#   simulation_time   — Total simulation length in ns (default: sim_time).
#
# Returns:
#   A data.table with columns: frame, time, pair, lig_ring, base_residue,
#   present, mode, d, angle, offset, model, experiment, experiment_raw, seq,
#   ligand.
read_pi_stacking <- function(file, stride = 5, simulation_time = sim_time) {
  base <- basename(file)
  experiment_raw <- sub("_task_.*\\.csv$", "", base)
  experiment <- experiment_raw
  sequence <- strsplit(experiment, "_")[[1]][1]
  lgd <- strsplit(experiment, "_")[[1]][2]
  model <- sub("^.*_task_", "task_", sub("\\.csv$", "", base))

  cat("Reading pi-stacking", model, "for experiment",
      paste0(sequence, "\u00B7", lgd), "\n")

  dt <- fread(file, na.strings = "", strip.white = TRUE)
  setnames(dt, gsub("^#", "", names(dt)))

  n_frames <- max(dt$Frame, na.rm = TRUE)

  # Zero-pad absent frames so `present` is defined over the whole run
  grid <- CJ(frame = seq_len(n_frames), pair = unique(dt$pair))
  grid <- merge(grid, unique(dt[, .(pair, lig_ring, base_residue)]),
                by = "pair")
  dt[, frame := as.integer(Frame)]
  dt[, Frame := NULL]
  dt <- merge(grid, dt[, .(frame, pair, present, mode, d, angle, offset)],
              by = c("frame", "pair"), all.x = TRUE)
  dt[is.na(present), present := 0L]
  dt[is.na(mode), mode := ""]
  dt[is.na(d), d := NA_real_]
  dt[is.na(angle), angle := NA_real_]
  dt[is.na(offset), offset := NA_real_]

  dt[, time := ((frame - 1) * stride + 1) *
       simulation_time / (n_frames * stride)]
  dt[, `:=`(
    model = model,
    experiment = paste0(sequence, "\u00B7", lgd),
    experiment_raw = experiment_raw,
    seq = sequence,
    ligand = lgd
  )]
  dt[]
}


# ─── pdb_pi_stacking ────────────────────────────────────────────────────
# Evaluates the pi-stacking set of a single structure (PDB file, e.g. the
# final minimised structure) in situ, using the same geometric criteria as
# the Python script (PQ4/MD/python/pi_stacking.py):
#   - aromatic ligand rings from ring_defs.json (written by the script)
#   - aromatic base rings from the standard DNA base ring atom sets
#     (same ring_defs.json)
#   - parallel : centroid d <= 5.5 A, acute angle <= 30 deg, offset <= 2 A
#   - T-shaped : centroid d <= 5.5 A, angle 60-90 deg
#
# Args:
#   pdb             — Path to a PDB file (e.g. *_stripped.pdb).
#   ring_defs_file  — Path to ring_defs.json (PQ4/MD/pi_stacking/ring_defs.json).
#   ligand_resname  — Ligand residue name (default: "PQ").
#   dist_cut        — Centroid distance cutoff in A (default: 5.5).
#   angle_cut       — Parallel interplanar angle cutoff in degrees (default: 30).
#   offset_cut      — Lateral offset cutoff in A (default: 2.0).
#   buffer          — Distance prefilter in A (default: 6.5).
#
# Returns:
#   A data.table with columns: pair, lig_ring, base_residue, present, mode,
#   d, angle, offset.
pdb_pi_stacking <- function(pdb, ring_defs_file, ligand_resname = "PQ",
                            dist_cut = 5.5, angle_cut = 30, offset_cut = 2.0,
                            buffer = 6.5) {
  defs <- jsonlite::fromJSON(ring_defs_file)
  lig_rings <- defs$ligand[[ligand_resname]]
  # jsonlite auto-simplifies single-ring bases to vectors; normalise to a
  # list of ring atom-name vectors (matrices -> list of rows) for uniformity
  base_sets <- lapply(defs$bases, function(g) {
    if (is.matrix(g)) g <- lapply(seq_len(nrow(g)), function(i) g[i, ])
    if (is.list(g)) g else list(g)
  })

  lines <- readLines(pdb, warn = FALSE)
  lines <- lines[grepl("^(ATOM  |HETATM)", lines)]
  if (length(lines) == 0) return(data.table(pair = character()))

  resnum <- as.integer(substr(lines, 23, 26))
  resname <- gsub(" ", "", substr(lines, 18, 20))
  atom <- gsub(" ", "", substr(lines, 13, 16))
  x <- as.numeric(substr(lines, 31, 38))
  y <- as.numeric(substr(lines, 39, 46))
  z <- as.numeric(substr(lines, 47, 54))
  elem <- gsub(" ", "", substr(lines, 77, 78))
  missing <- elem == ""
  elem[missing] <- substr(atom[missing], 1, 1)

  atoms <- data.table(atom = atom, resname = resname, resnum = resnum,
                      x = x, y = y, z = z, elem = elem)

  # Ring geometry: centroid + normal (smallest singular vector)
  ring_geom <- function(tbl) {
    ctr <- colMeans(tbl[, .(x, y, z)])
    svd_t <- svd(scale(tbl[, .(x, y, z)], center = TRUE, scale = FALSE))
    nml <- svd_t$v[, 3]
    list(ctr = ctr, nml = nml)
  }

  rings <- list()
  for (rn in names(lig_rings)) {
    sub <- atoms[resname == ligand_resname & atom %in% lig_rings[[rn]]]
    if (nrow(sub) == length(lig_rings[[rn]])) {
      rings[[rn]] <- ring_geom(sub)
    }
  }

  base_rings <- list()
  for (bs in names(base_sets)) {
    rng <- atoms[resname == bs]
    for (ring_atoms in base_sets[[bs]]) {
      for (res in unique(rng$resnum)) {
        sub <- rng[resnum == res & atom %in% ring_atoms]
        if (nrow(sub) == length(ring_atoms)) {
          base_rings[[paste0(bs, res)]] <- ring_geom(sub)
        }
      }
    }
  }

  out <- vector("list", length(rings) * length(base_rings))
  k <- 0L
  for (rn in names(rings)) {
    c1 <- rings[[rn]]$ctr; n1 <- rings[[rn]]$nml
    for (key in names(base_rings)) {
      c2 <- base_rings[[key]]$ctr; n2 <- base_rings[[key]]$nml
      dvec <- c2 - c1
      d <- sqrt(sum(dvec^2))
      if (d > buffer) next
      cosv <- abs(sum(n1 * n2)) / (sqrt(sum(n1^2)) * sqrt(sum(n2^2)))
      ang <- acos(pmin(1, cosv)) * 180 / pi
      h <- abs(sum(dvec * n1))
      off <- sqrt(max(0, d^2 - h^2))
      if (d <= dist_cut && ang <= angle_cut && off <= offset_cut) {
        mode <- "parallel"; present <- 1L
      } else if (d <= dist_cut && ang >= 60 && ang <= 90) {
        mode <- "T-shaped"; present <- 1L
      } else {
        mode <- ""; present <- 0L
      }
      k <- k + 1L
      out[[k]] <- data.table(
        pair = paste0(rn, "-", key),
        lig_ring = rn,
        base_residue = key,
        present = present, mode = mode,
        d = round(d, 3), angle = round(ang, 2), offset = round(off, 3)
      )
    }
  }
  res <- rbindlist(out[seq_len(k)])
  if (is.null(res) || nrow(res) == 0) return(data.table(pair = character()))
  res[]
}


# ─── plot_pi_stacking_tracked ───────────────────────────────────────────
# Plots the time evolution of a predefined set of ligand-DNA pi-stacking
# pairs. Mirrors plot_contacts_tracked() (same single-panel heatmap layout
# and colour scheme), with rows being ligand-ring/unit x base-ring pairs
# instead of H-bonds. Row labels use "<descriptive ligand ring>\u2225<cleaned
# DNA base>" (e.g. "Q1\u2225A10"), ordered by ring first so each ring's rows
# group together and reading top-to-bottom traces the ring across nucleotides:
# by default the pair is a fused quinoline unit (Q0/Q1 -> Q1/Q2); with
# --no-simplify it is a single ring (R0-R3 -> "Q1\u00b7pyr" etc.), where the
# ligand ring is named by its quinoline system (Q1/Q2) and ring type
# (pyr/benz).
#
# Args:
#   dt_long      — Long stacking data from read_pi_stacking().
#   ref          — Reference-set membership: one row per pair with columns
#                  pair, model, seq, experiment, in_initial (0/1),
#                  in_final (0/1), plus base_residue/lig_ring for labelling.
#   max_t        — Maximum time (ns) for the x-axis.
#   scale        — Scaling factor passed to theme_custom().
#   heat_bin_ns  — Time resolution (ns) used to bin the heatmap.
#
# Returns:
#   A ggplot object.
plot_pi_stacking_tracked <- function(dt_long, ref, max_t = 100, scale = 1,
                                     heat_bin_ns = 0.1) {
  ref_sub <- ref[, .(pair, model, seq, experiment, in_initial, in_final,
                     base_residue, lig_ring)]
  dt <- copy(dt_long)[, c("base_residue", "lig_ring") := NULL]
  dt <- merge(dt, ref_sub, by = c("pair", "model", "seq", "experiment"))
  dt <- dt[in_initial == 1 | in_final == 1]
  if (nrow(dt) == 0) {
    msg <- paste0("No ligand-DNA pi-stacking (5.5 A centroid distance, ",
                  "30 deg parallel / 60-90 deg T-shaped) was present in the ",
                  "initial state or the final minimised structure for these ",
                  "models.")
    empty_p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
      annotate("text", 0, 0, label = msg, size = 3, lineheight = 0.9) +
      theme_void() +
      theme(panel.border = element_rect(fill = NA, colour = "grey70"))
    return(empty_p)
  }
  dt[, bin := floor(time / heat_bin_ns) * heat_bin_ns]

  model_levels <- unique(dt_long[, .(model, seq)])
  model_levels <- model_levels[order(as.integer(sub("task_", "", model)))]$model
  model_labels <- paste0("model ", sub("task_", "", model_levels))

  # ---- Presence heatmap ----
  heat_agg <- dt[, .(present_frac = mean(present),
                     in_initial = max(in_initial),
                     in_final = max(in_final)),
                 by = .(pair, model, seq, bin)]
  heat_agg[, fill_key := fcase(
    present_frac >= 0.5 & in_initial == 1 & in_final == 1, "both",
    present_frac >= 0.5 & in_initial == 1, "initial",
    present_frac >= 0.5 & in_final == 1, "final",
    default = "absent")]

  order_dt <- unique(dt[, .(pair, model, seq, in_initial, in_final,
                            base_residue, lig_ring)])
  order_dt[, base_num := as.integer(gsub("[^0-9]", "", base_residue))]
  order_dt[,
    row_key := paste0(model, " | ", pair)
  ][order(model, seq, lig_ring, base_num, base_residue),
    row_level := seq_len(.N)
  ]

  heat_agg <- merge(heat_agg,
                    order_dt[, .(pair, model, seq, row_key)],
                    by = c("pair", "model", "seq"))

  # Row labels: "<descriptive ligand ring>\u2225<cleaned DNA base>"
  row_lab <- setNames(
    paste0(ring_label(order_dt$lig_ring), stack_sep,
           pad_dna_residue(clean_dna_residue(order_dt$base_residue))),
    order_dt$row_key)

  missing_models <- setdiff(model_levels, unique(as.character(heat_agg$model)))
  dummy_keys <- character()
  if (length(missing_models) > 0) {
    dummy_keys <- paste0(missing_models, " | \u2014")
    heat_agg <- rbind(heat_agg, data.table(
      bin = NA_real_,
      model = missing_models,
      seq = unique(dt_long$seq),
      pair = NA_character_,
      present_frac = NA_real_,
      in_initial = 0L,
      in_final = 0L,
      fill_key = NA_character_,
      row_key = dummy_keys), fill = TRUE)
    row_lab <- c(row_lab, setNames(rep("", length(dummy_keys)), dummy_keys))
  }

  # ---- Minimised-structure membership column ----
  # One tile per row to the right of the trajectory showing whether the
  # pair is present in the final minimised structure (step-11 PDB),
  # independent of the trajectory occupancy used for the time bins.
  marker_x <- max_t + 1.5
  marker_w <- 0.8
  marker_dt <- order_dt[, .(model, seq, row_key, in_initial, in_final)]
  marker_dt[, `:=`(
    bin = marker_x,
    present_frac = NA_real_,
    fill_key = fifelse(in_final == 1, "min_yes", "min_no")
  )]
  heat_agg <- rbind(heat_agg, marker_dt, fill = TRUE)

  # Q1 series on top, Q2 below (row order is inverted relative to the sort)
  heat_agg[, row_key := factor(row_key,
                               levels = c(rev(order_dt$row_key[order(order_dt$row_level)]),
                                          dummy_keys))]
  heat_agg[, model := factor(model, levels = model_levels, labels = model_labels)]

  time_agg <- heat_agg[fill_key %in% c("min_yes", "min_no") == FALSE]
  marker_agg <- heat_agg[fill_key %in% c("min_yes", "min_no")]

  heat_p <- ggplot(heat_agg, aes(bin, row_key, fill = fill_key,
                                 pattern = fill_key)) +
    geom_tile_pattern(
      data = time_agg,
      width = heat_bin_ns * 1.02, height = 0.9, na.rm = TRUE,
      pattern_type = "stripe",
      pattern_colour = NA,
      pattern_fill = "#A6E6DB",
      pattern_fill2 = "#0D2B52",
      pattern_density = 0.5,
      pattern_angle = 45
    ) +
    geom_tile_pattern(
      data = marker_agg,
      width = marker_w, height = 0.9, na.rm = TRUE,
      pattern_type = "none",
      pattern_colour = NA
    ) +
    geom_vline(xintercept = max_t, colour = "grey60",
               linewidth = 0.4, na.rm = TRUE) +
    facet_grid(model ~ ., scales = "free_y", space = "free_y", drop = FALSE) +
    theme_custom(scale) +
    scale_fill_manual(
      values = c(absent = "#F2F2F2", initial = "#A6E6DB",
                 final = "#0D2B52", both = "#0D2B52",
                 min_yes = "#0D2B52", min_no = "#FFFFFF"),
      breaks = c("absent", "initial", "final", "both"),
      labels = c("absent", "present, initial", "present, minimized",
                 "present, initial & minimized"),
      name = NULL
    ) +
    scale_pattern_manual(
      values = c(absent = "none", initial = "none", final = "none",
                 both = "stripe", min_yes = "none", min_no = "none"),
      breaks = c("absent", "initial", "final", "both"),
      labels = c("absent", "present, initial", "present, minimized",
                 "present, initial & minimized"),
      name = NULL
    ) +
    scale_x_continuous("Time (ns)", breaks = seq(0, max_t, by = 25),
                       expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_t + 3.2)) +
    scale_y_discrete(labels = function(x) row_lab[x], name = NULL) +
    guides(fill = guide_legend(nrow = 2)) +
    theme(legend.position = "top",
          legend.key.height = unit(30, "pt"),
          legend.key.spacing.y = unit(0.4, "cm"),
          panel.spacing.y = unit(0.35, "lines"),
          axis.text.y = element_markdown(size = 14 * scale, face = "bold"),
          strip.text.y = element_text(size = 14 * scale, face = "bold"),
          strip.background = element_rect(fill = "grey92", colour = NA))

  heat_p
}


# ─── pi_stacking_summary ────────────────────────────────────────────────
# Builds a summary data.table of the reference pi-stacking pairs for display
# (e.g. with DT::datatable). One row per pair of the union of the initial and
# final reference sets, with set membership, occupancy over the full run and
# mean geometry / dominant stacking mode (over the frames where the pair is
# present).
#
# Args:
#   ref      — Reference-set membership (see plot_pi_stacking_tracked).
#   dt_long  — Long stacking data from read_pi_stacking().
#
# Returns:
#   A data.table with columns: Model (labelled "model N"), Pair (as
#   "<ligand ring>\u2225<DNA base>", e.g. "Q1\u2225A10"), DNA residue (cleaned),
#   Ligand ring (descriptive), Set, Occupancy, Mean d (A), Mean angle (deg),
#   Mode.
pi_stacking_summary <- function(ref, dt_long) {
  summ <- copy(ref)
  summ <- summ[in_initial == 1 | in_final == 1]
  summ[, set := fifelse(in_initial == 1 & in_final == 1, "both",
                        fifelse(in_initial == 1, "initial", "minimized"))]

  occ <- dt_long[, .(occupancy = mean(present)), by = .(model, pair)]
  geo <- dt_long[present == 1, .(
    mean_d = mean(d, na.rm = TRUE),
    mean_ang = mean(angle, na.rm = TRUE),
    mode = names(which.max(table(mode)))
  ), by = .(model, pair)]

  summ <- merge(summ, occ, by = c("model", "pair"), all.x = TRUE)
  summ <- merge(summ, geo, by = c("model", "pair"), all.x = TRUE)
  summ[is.na(occupancy), `:=`(occupancy = 0, mean_d = NA_real_,
                              mean_ang = NA_real_, mode = NA_character_)]
  summ[is.na(mode), mode := ""]

  summ[, `:=`(
    Model = paste0("model ", sub("task_", "", model)),
    Pair = paste0(ring_label(lig_ring), stack_sep,
                  pad_dna_residue(clean_dna_residue(base_residue))),
    `DNA residue` = clean_dna_residue(base_residue),
    `Ligand ring` = ring_label(lig_ring),
    Set = set,
    `Occupancy` = round(occupancy, 3),
    `Mean d (A)` = round(mean_d, 2),
    `Mean angle (deg)` = round(mean_ang, 1),
    Mode = mode
  )]
  summ[, base_num := as.integer(gsub("[^0-9]", "", base_residue))]
  summ <- summ[order(model, lig_ring, base_num, base_residue)]
  summ[, .(Model, Pair, `DNA residue`, `Ligand ring`, Set,
           `Occupancy`, `Mean d (A)`, `Mean angle (deg)`, Mode)]
}


#===============================================================================
# Constraint verification for minimised structures
#===============================================================================

# read_minimized_constraints(project)
#   Reads all "*_minimized_constraints.csv" files under
#   "{project}/MD/pmemd/out/" (one per experiment/job, written by
#   process_minimized_constraint_verification in templates/MD/python/pymol_utils.py)
#   into a single data.table. Derives 'sequence' (e.g. "CSS1") and 'condition'
#   (Free / HCY / PQ) labels from the experiment name and adds the job id.
#
# Returns:
#   A data.table with columns sequence, condition, experiment, job, model,
#   n_constraints, n_satisfied, all_verified, failed_constraints,
#   worst_min_distance, ordered by sequence -> condition -> model.
read_minimized_constraints <- function(project) {
  constraint_files <- list.files(
    path = file.path(project, "MD", "pmemd", "out"),
    pattern = "_minimized_constraints\\.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(constraint_files) == 0) {
    return(data.table(
      sequence = character(), condition = character(), experiment = character(),
      job = character(), model = integer(), n_constraints = integer(),
      n_satisfied = integer(), all_verified = integer(),
      failed_constraints = character(), worst_min_distance = numeric()
    ))
  }

  constraint_dt <- lapply(constraint_files, function(x) {
    dt <- fread(x)
    exp_dir <- dirname(x)
    dt[, c("experiment", "job") := .(
      basename(dirname(exp_dir)),
      basename(exp_dir)
    )]
    dt
  }) |>
    rbindlist(fill = TRUE)

  constraint_dt[, sequence := sub("_.*$", "", experiment)]
  model_base <- sub("_constrained$", "", constraint_dt$experiment)
  cond_raw <- sub("^[^_]+_", "", model_base)
  cond_map <- c(
    free = "Free",
    hcy = "HCY",
    pq = "PQ",
    piperaquine = "PQ"
  )
  cond_key <- tolower(cond_raw)
  constraint_dt[, condition := ifelse(
    cond_key %in% names(cond_map),
    unname(cond_map[cond_key]),
    cond_raw
  )]

  # Keep only the most recent (highest-numbered) job per experiment, so that a
  # superseded MD run (e.g. a relaunched production simulation) does not add a
  # duplicate set of rows next to the latest one.  All model rows of that job
  # are retained (one row per task replicate).
  constraint_dt[, job_num := as.numeric(sub("^J", "", job))]
  constraint_dt <- constraint_dt[
    constraint_dt[, .I[job_num == max(job_num)], by = experiment]$V1
  ]
  constraint_dt[, job_num := NULL]

  setcolorder(
    constraint_dt,
    c(
      "sequence", "condition", "experiment", "job", "model",
      "n_constraints", "n_satisfied", "all_verified",
      "failed_constraints", "worst_min_distance"
    )
  )
  setorder(constraint_dt, sequence, condition, experiment, model)
  constraint_dt
}


# ─── read_density ───────────────────────────────────────────────────────
# Extracts the time series of instantaneous density from an Amber pmemd
# .out file (the production run) and returns a tidy data.table.
#
# Amber prints one data block per "NSTEP" line in the "4. RESULTS" section:
# the block starting with "NSTEP = ..." reports instantaneous values
# (including "Density") for that step. Parsing stops at the
# "A V E R A G E S   O V E R  100000 S T E P S" marker that precedes the
# final block-average section, so that the run-mean density and its RMS
# fluctuation (which would otherwise appear as an artefactual jump in the
# series) are not included.
#
# Metadata (seq, ligand, experiment, job, model) is extracted from the file
# path using the same convention as read_rmsd()/read_contacts_long(): path
# components 5 (experiment), 6 (job) and 7 (model, e.g. "task_0").
#
# Args:
#   file — Path to a step10_*.out file.
#
# Returns:
#   A data.table with columns: seq, ligand, experiment, job, model,
#   t_ns, t_ps, density.
read_density <- function(file) {
  comps <- strsplit(file, split = "[/\\\\]")[[1]]
  experiment <- comps[5]
  sequence <- strsplit(experiment, "_")[[1]][1]
  ligand <- strsplit(experiment, "_")[[1]][2]
  job <- comps[6]
  model <- comps[7]

  cat("Reading", model, "for experiment",
      paste0(sequence, "\u00B7", ligand), "\n")

  lines <- readLines(file)

  # Stop at the block-average section: the average NSTEP block that follows
  # it carries the run-mean density / RMS fluctuation, not trajectory data.
  marker <- grep("A V E R A G E S   O V E R", lines, fixed = TRUE)
  marker <- if (length(marker) > 0) marker[1] else length(lines)

  nstep_idx <- grep("NSTEP\\s*=\\s*\\d+\\s+TIME\\(PS\\)", lines)
  nstep_idx <- nstep_idx[nstep_idx < marker]
  if (length(nstep_idx) == 0) {
    stop("no NSTEP blocks found before the block-average marker in: ", file)
  }

  time_ps <- as.numeric(sub(".*TIME\\(PS\\)\\s*=\\s*([0-9.]+).*", "\\1",
                            lines[nstep_idx]))

  # Density is printed on the 5th line after the NSTEP line. Fall back to a
  # short forward scan if the offset is ever off (e.g. extra output lines).
  density_idx <- nstep_idx + 5L
  has_density <- grepl("Density", lines[density_idx])
  density <- rep(NA_real_, length(nstep_idx))
  density[has_density] <- as.numeric(sub(".*Density\\s*=\\s*([0-9.]+).*", "\\1",
                                         lines[density_idx][has_density]))
  off <- which(!has_density)
  for (m in off) {
    for (j in (nstep_idx[m] + 1L):(nstep_idx[m] + 6L)) {
      if (grepl("Density", lines[j])) {
        density[m] <- as.numeric(sub(".*Density\\s*=\\s*([0-9.]+).*", "\\1",
                                     lines[j]))
        break
      }
    }
  }

  data.table(
    seq = sequence,
    ligand = ligand,
    experiment = paste0(sequence, "\u00B7", ligand),
    job = job,
    model = model,
    t_ps = time_ps,
    t_ns = time_ps / 1000,
    density = density
  )
}


# ─── density_plateau ────────────────────────────────────────────────────
# Assesses whether the density time series of a production run has reached
# a plateau, following Bogetti et al., J. Chem. Phys. 153, 054123 (2020)
# (as implemented in cpptraj's "evalplateau" command).
#
# The density series D(t) is fitted to a single exponential:
#   D(t) = D_initial + (D_final - D_initial) * (1 - exp(-k * t))
# with time shifted so that the first point is at t = 0, and D_initial
# seeded from the mean of the first 1% of the data. A plateau is reached
# when all three criteria hold:
#   1. the final slope of the fitted curve is < slope_cut (g cm-3 ps-1);
#   2. |D_final - mean(second half of the data)| < df_cut (g cm-3);
#   3. the (reduced) chi-squared of the fit is < chi2_cut.
#
# The chi-squared is computed as the mean squared residual (SSR / N), which
# matches the paper's statement that the cutoff of 0.5 corresponds to a
# total deviation of about 0.71 g cm-3 (sqrt(0.5)), and is scale-invariant
# with respect to the number of points (cpptraj's raw SSR is calibrated for
# its short, ~1 ns datasets and would spuriously fail a 100,000-point
# series).
#
# Args:
#   dt        — A data.table from read_density() with columns t_ns and
#               density, plus grouping columns (seq, experiment, model, ...).
#   slope_cut — Final-slope cutoff in g cm-3 ps-1 (default 1e-6).
#   df_cut    - |D_final - second-half mean| cutoff in g cm-3 (default 0.02).
#   chi2_cut  - Reduced chi-squared cutoff (default 0.5).
#
# Returns:
#   A data.table with one row per group: n, d_initial, d_final, k, chi2,
#   final_slope, df_diff, plateau_ns, reached, converged.
density_plateau <- function(dt, slope_cut = 1e-6, df_cut = 0.02,
                            chi2_cut = 0.5) {
  group_cols <- setdiff(names(dt), c("t_ns", "t_ps", "density"))
  fit <- function(t_ns, density) {
    ord <- order(t_ns)
    t0_ns <- t_ns[ord][1]
    t <- t_ns[ord] - t0_ns
    D <- density[ord]
    n <- length(t)

    # Least-squares fit of D(t) = DF - (DF - DI) * exp(-k * t). For a fixed
    # k the parameters DF and (DF - DI) are linear in the model, so k is
    # profiled over a log-spaced grid and the linear problem solved for each
    # candidate. This is equivalent to the non-linear least squares of
    # Bogetti et al. (2020) but is stable on near-constant density series,
    # where the exponential amplitude is within the noise and Gauss-Newton
    # type solvers diverge. The series is decimated for the fit (the trend
    # is smooth); the reduced chi-squared is evaluated on the full series.
    #
    # k is bounded below by 1/t_max so that exp(-k*t) decays meaningfully
    # within the observed window (exp(-k*t_max) <= exp(-1)). Without this
    # bound, essentially-flat series return degenerate fits (k -> 0, with
    # D_final unconstrained), which spurious fail the |D_final - second-half
    # mean| criterion.
    idx <- seq.int(1L, n, by = 100)
    t_max <- t[n]
    k_min <- max(1e-6, 1 / t_max)
    k_grid <- 10^seq(log10(k_min), 2, length.out = 300)
    best <- NULL
    for (k in k_grid) {
      X <- cbind(1, -exp(-k * t[idx]))
      Q <- qr(X)
      if (Q$rank < 2) next
      b <- tryCatch(qr.coef(Q, D[idx]), error = function(e) NULL)
      if (is.null(b) || anyNA(b)) next
      resid <- D[idx] - drop(X %*% b)
      sse <- sum(resid^2)
      if (is.null(best) || sse < best$sse) {
        best <- list(k = k, DF = b[1], DI = b[1] - b[2], sse = sse)
      }
    }

    converged <- !is.null(best)
    if (converged) {
      k <- best$k
      DI <- best$DI
      DF <- best$DF
      pred <- DI + (DF - DI) * (1 - exp(-k * t))
    } else {
      k <- NA_real_
      DI <- mean(D)
      DF <- mean(D)
      pred <- rep(DF, n)
    }
    chi2 <- mean((D - pred)^2)
    second_half <- D[seq.int(n %/% 2 + 1L, n)]
    df_diff <- abs(DF - mean(second_half))

    if (converged && is.finite(k) && k > 0) {
      amp <- k * abs(DF - DI)                 # slope magnitude at t = 0 (ns-1)
      final_slope <- amp * exp(-k * t_max) / 1000  # g cm-3 ps-1
      if (amp > slope_cut * 1000) {
        plateau_ns <- max(0, log(amp / (slope_cut * 1000)) / k)
      } else {
        plateau_ns <- 0
      }
    } else {
      final_slope <- 0
      plateau_ns <- 0
    }

    reached <- isTRUE(final_slope < slope_cut) &&
      isTRUE(df_diff < df_cut) &&
      isTRUE(chi2 < chi2_cut)

    data.table(
      n = n,
      t0_ns = t0_ns,
      d_initial = DI,
      d_final = DF,
      k = k,
      chi2 = chi2,
      final_slope = final_slope,
      df_diff = df_diff,
      plateau_ns = plateau_ns,
      reached = as.integer(reached),
      converged = as.integer(converged)
    )
  }

  if (length(group_cols) == 0) {
    fit(dt$t_ns, dt$density)
  } else {
    dt[, fit(t_ns, density), by = c(group_cols)]
  }
}


# ─── plot_density ───────────────────────────────────────────────────────
# Plots the density time series of a production run as a line chart,
# faceted by model, with the fitted plateau exponential overlaid.
#
# Args:
#   dt       — A data.table from read_density() (t_ns, density, model,
#              experiment columns).
#   plateau  — A data.table from density_plateau() with the same grouping
#              columns plus d_initial, d_final, k and the first time point
#              (t0_ns). NULL = no fitted curve (default).
#   max_t    — Maximum time (ns) for the x-axis (default: 100).
#   scale    — Scaling factor passed to theme_custom().
#
# Returns:
#   A ggplot object.
plot_density <- function(dt, plateau = NULL, max_t = 100, scale = 1) {
  dt[
    order(experiment, model),
    subgroup_id := rleid(model),
    by = experiment
  ]
  if (!is.null(plateau) && nrow(plateau) > 0) {
    plateau[
      order(experiment, model),
      subgroup_id := rleid(model),
      by = experiment
    ]
  }
  max_plot <- max(c(max_t, max(dt$t_ns)))

  p <- dt |>
    ggplot(aes(t_ns, density)) +
    geom_point(size = 0.3, color = "grey45", alpha = 0.5, stroke = 0) +
    facet_grid(subgroup_id ~ experiment) +
    scale_x_continuous("t (ns)", expand = c(0, 0), limits = c(0, max_plot)) +
    scale_y_continuous("Density (g cm\u207B\u00B3)") +
    theme_custom(scaling = scale)

  if (!is.null(plateau) && nrow(plateau) > 0) {
    curve_dt <- plateau[, {
      tt <- seq(t0_ns, max_plot, length.out = 400)
      k_eff <- ifelse(is.na(k), 0, pmax(k, 0))
      .(
        t_ns = tt,
        d_fit = d_initial + (d_final - d_initial) *
          (1 - exp(-k_eff * (tt - t0_ns)))
      )
    }, by = .(experiment, subgroup_id, model)]
    p <- p +
      geom_line(
        data = curve_dt,
        aes(t_ns, d_fit, group = model),
        color = "firebrick",
        linewidth = 0.75,
        inherit.aes = FALSE
      )
  }

  p
}
