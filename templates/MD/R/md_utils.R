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
#
# DEPENDENCIES
#   data.table, stringr, ggplot2, ggrepel
#
# FUNCTIONS
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
#     scaling factor. Removes grid lines and background for clean
#     publication output.
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
  experiment <- strsplit(file, split = "/")[[1]][5]
  sequence <- strsplit(experiment, '_')[[1]][1]
  ligand <- strsplit(experiment, '_')[[1]][2]
  experiment <- paste0(
    sequence,
    "·",
    ligand
  )

  # Extract model name from the 7th path component
  # e.g. ".../CSS1_rep1/wildtype/2drmsd.dat" → "wildtype"
  model <- strsplit(file, split = "/")[[1]][7]

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
      "Time step filter: requested", time_step, "ns,",
      "data resolution", round(dt, 4), "ns,",
      "keeping every", n, "th frame (effective:",
      round(n * dt, 4), "ns)\n"
    )
    rmsd_dt <- rmsd_dt[t1 %in% keep_t & t2 %in% keep_t]
  }
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
# in Angstroms). Metadata (experiment, model, seq, ligand) is extracted
# from the file path using the same convention as read_rmsd().
#
# Args:
#   file — Path to the atomicfluct .dat file.
#
# Returns:
#   A data.table with columns: res, atomfluct, experiment, model, seq,
#   ligand.
read_atomicfluct <- function(file) {
  # Extract experiment name from path component [5]
  experiment <- strsplit(file, split = "/")[[1]][5]
  sequence <- strsplit(experiment, '_')[[1]][1]
  ligand <- strsplit(experiment, '_')[[1]][2]

  experiment <- paste0(
    sequence,
    "·",
    ligand
  )

  # Extract model name from path component [7]
  model <- strsplit(file, split = "/")[[1]][7]

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
    plot.background = element_blank(),
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
plot_rmsd <- function(rmsd, max_t = NA, time_step = NA, scale = 1) {
  cat("Plotting RMSD for sequence(s):", paste(unique(rmsd$seq), collapse = ", "), "\n")

  # Optionally subsample to a coarser time resolution
  if (!is.na(time_step)) {
    unique_t <- sort(unique(rmsd$t1))
    dt <- median(diff(unique_t))
    n <- max(1, round(time_step / dt))
    keep_t <- unique_t[seq(1, length(unique_t), by = n)]
    cat(
      "Time step filter: requested", time_step, "ns,",
      "data resolution", round(dt, 4), "ns,",
      "keeping every", n, "th frame (effective:",
      round(n * dt, 4), "ns)\n"
    )
    rmsd <- rmsd[t1 %in% keep_t & t2 %in% keep_t]
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
