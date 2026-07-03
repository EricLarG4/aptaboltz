# processor.R
# ===========
# Boltz-2 result visualisation in R.
#
# Reads the aggregated CSV files produced by process_boltz_results.py and
# produces:
#   - An interactive HTML confidence-score table (DT + bslib theming)
#   - PAE / PDE heatmap faceted plots (ggplot2)
#   - pLDDT per-residue line plots (ggplot2)
#
# Entry-point functions (called from a project-specific process.R):
#   process_confidence(project)  -> data.table
#   table_confidence(dt)         -> DT::datatable (side effect: saves HTML)
#   process_pxe(project, type)   -> named list of data.tables
#   plot_pxe(pxe, type, ...)     -> list of ggplot2 objects (side effect: saves PNGs)
#   plot_plddt(project)          -> list of ggplot2 objects (side effect: saves PNGs)


library(data.table)
library(DT)
library(ggplot2)
library(khroma)
library(jsonlite)
library(yaml)
library(bslib)


# =====================================================================
#  Confidence table: read aggregated CSVs, parse dict columns, reorder
# =====================================================================

process_confidence <- function(project) {
  """
  Read all *confidence.csv files under *project* into a single data.table.

  Columns containing JSON-like dict strings (e.g. `chains_ptm`,
  `pair_chains_iptm`) are parsed and expanded into one column per chain
  (e.g. `chains_ptm_chain0`, `chains_ptm_chain1`).  The original dict
  columns are removed.

  Returns a data.table with columns ordered:
      base (model, experiment, id) -> metrics -> chain-specific columns.
  """

  # 1. List all confidence CSV files recursively
  confidence_files <- list.files(
    path = project,
    pattern = "confidence.csv",
    recursive = TRUE,
    full.names = TRUE
  )

  # 2. Read each CSV, add model and experiment info from the directory path
  confidence_dt <- lapply(
    confidence_files,
    function(x) {
      dt <- fread(x)
      dt[,
        c("model", "experiment") := tstrsplit(
          gsub(paste0(project, "/"), "", dirname(x)),
          "/",
          fixed = TRUE
        )
      ]
      dt[, id := .I]
      return(dt)
    }
  ) |>
    rbindlist()

  # 3. Parse dict columns and expand into chain-specific columns
  #    These columns contain JSON-like strings: {'0': value, '1': value, ...}
  #    We need to parse them because Python's pd.DataFrame.to_csv() serialises
  #    nested dict cells as Python-repr strings, not valid JSON.
  dict_cols <- c(
    "chains_ptm",
    "pair_chains_iptm",
    "chains_pae",
    "pair_chains_pae"
  )
  for (col in intersect(dict_cols, names(confidence_dt))) {
    # Replace single quotes with double quotes for JSON parsing
    parsed <- lapply(confidence_dt[[col]], function(x) {
      fromJSON(gsub("'", '"', x))
    })

    # Determine all chain keys present (e.g., "0", "1", ...)
    chains <- unique(unlist(lapply(parsed, names)))
    chains <- sort(chains) # ensure consistent ordering
    # Create one column per chain
    for (ch in chains) {
      values <- sapply(parsed, function(x) as.numeric(x[[ch]]))
      confidence_dt[, paste0(col, "_chain", ch) := values]
    }
    # Remove the original dict column
    confidence_dt[, (col) := NULL]
  }

  # 4. Reorder columns so base identifiers come first,
  #    then metric columns, then chain-specific expansions
  all_cols <- names(confidence_dt)
  base_cols <- c("model", "experiment", "id", "confidence_score", "ptm")
  base_cols <- intersect(base_cols, all_cols)

  chain_cols <- grep("_chain", all_cols, value = TRUE)
  metric_cols <- setdiff(
    all_cols,
    c(
      "model",
      "experiment",
      "id",
      "ptm",
      grep("_chain", all_cols, value = TRUE)
    )
  )

  # Build ordered column list: base -> metrics -> chain columns
  ordered_cols <- c(
    intersect(base_cols, all_cols),
    setdiff(metric_cols, base_cols),
    chain_cols
  )
  confidence_dt <- confidence_dt[, ..ordered_cols]

  return(confidence_dt)
}


# =====================================================================
#  Render interactive confidence table (DT) and save as HTML
# =====================================================================

table_confidence <- function(confidence_dt) {
  """
  Build an interactive DT::datatable with colour-bar styling, wrap it in
  a Bootstrap 5 theme (bslib), and save the result as an HTML file in the
  project directory.

  Styling rules:
    - `Smaller is better` metrics (PAE, PDE, IPDE) get a horizontal bar
      with angle = 270 (i.e. bars grow right-to-left so smaller values
      appear fuller).
    - All other metrics get angle = 90 (bars grow left-to-right).
    - Chain-specific columns are assigned cycling hues derived from their
      parent dict-column type.
    - Ligand_iptm / protein_iptm columns are removed if all values are zero.
  """

  # 1. Define colour palette for each metric column
  color_list <- list(
    confidence_score = "lightblue",
    ptm = "thistle",
    iptm = "#e2cfe2",
    ligand_iptm = "#ecdfec",
    protein_iptm = "#f6eff6",
    complex_plddt = "#9fc69f",
    complex_iplddt = "#bfd9bf",
    complex_pde = "#D3D3D3",
    complex_ipde = "#e0e0e0",
    complex_pae = "#b0c4de",
    complex_ipae = "#c2d1e6"
  )

  # Dynamically add colours for chain-specific columns (e.g. "chains_ptm_chain0")
  chain_hues <- list(
    chains_ptm = c("#FFB6C1", "#ffd0d7"),
    pair_chains_iptm = c("#ffcdc8", "#FFE4E1"),
    chains_pae = c("#40E0D0", "#40E0D0"),
    pair_chains_pae = c("#40E0D0", "#40E0D0")
  )

  for (col_name in names(confidence_dt)) {
    matched <- names(chain_hues)[sapply(names(chain_hues), function(prefix) {
      grepl(paste0("^", prefix, "_chain"), col_name)
    })]
    if (length(matched) == 1) {
      chain_idx <- gsub("^.*_chain", "", col_name)
      chain_num <- as.numeric(chain_idx) + 1  # 0 -> 1 for 1-based cycling
      hues <- chain_hues[[matched]]
      color_list[[col_name]] <- hues[(chain_num - 1) %% length(hues) + 1]
    }
  }

  # 2. Remove zero-valued iptm columns (no ligand was present)
  if (all(confidence_dt$ligand_iptm == 0, na.rm = TRUE)) {
    confidence_dt[, ligand_iptm := NULL]
    confidence_dt[, protein_iptm := NULL]
    color_list <- color_list[
      !names(color_list) %in% c("ligand_iptm", "protein_iptm")
    ]
  }

  # 3. Bootstrap 5 theme via bslib
  bs_theme <- bslib::bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = bslib::font_google("Inter")
  )

  # 4. Build the datatable with copy / CSV / Excel / print buttons
  display_table <- DT::datatable(
    confidence_dt,
    extensions = 'Buttons',
    filter = 'top',
    class = 'display cell-border stripe hover compact',
    options = list(
      pageLength = max(confidence_dt$id),
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'print'),
      lengthMenu = list(c(10, 25, 50, -1), c('10', '25', '50', 'All'))
    ),
    rownames = FALSE
  ) |>
    # Round all metric columns to 3 significant digits
    formatRound(columns = 4:ncol(confidence_dt), digits = 3) |>
    # Apply colour-bar background to each metric column
    (\(x) {
      Reduce(
        function(dt, col) {
          # "Smaller is better" metrics get reverse-direction colour bar
          smaller_better <- grepl("pae|pde|ipde", col, ignore.case = TRUE)
          bar_angle <- if (smaller_better) 270 else 90
          formatStyle(
            dt,
            columns = col,
            background = styleColorBar(
              range(confidence_dt[[col]], na.rm = TRUE),
              color_list[[col]],
              angle = bar_angle
            ),
            backgroundSize = '100% 90%',
            backgroundRepeat = 'no-repeat',
            backgroundPosition = 'center'
          )
        },
        names(color_list),
        init = x
      )
    })()

  # 5. Wrap in themed page and save as HTML
  themed_page <- bslib::page_fluid(
    title = "Confidence Scores",
    theme = bs_theme,
    display_table
  )
  htmltools::save_html(themed_page, file = paste0(project, "/confidence_table.html"))

  display_table
}


# =====================================================================
#  PAE / PDE: read CSVs and melt into long format for heatmap plotting
# =====================================================================

process_pxe <- function(project, type = c("pae", "pde")) {
  """
  Read all *_{type}.csv files under *project* into a named list of
  data.tables, where names are the relative paths (model/experiment).

  Parameters
  ----------
  project : str
  type : `pae` or `pde`
  """
  pxe_files <- list.files(
    path = project,
    pattern = paste0(type, ".csv"),
    recursive = TRUE,
    full.names = TRUE
  )
  pxe_dt <- lapply(pxe_files, fread) |>
    setNames(gsub(paste0(project, "/"), "", dirname(pxe_files)))
}

melt_pxe <- function(pxe) {
  """
  Convert a wide PAE/PDE matrix (residue × residue per model) into long
  format suitable for ggplot2 geom_raster().

  Each model's square matrix is melted into (res_1, res_2, pxe) triples,
  and a 'model' column is added (1-indexed).
  """
  melted_pxe <- lapply(
    unique(pxe$model),
    function(m) {
      pxe[model == m] |>
        _[, model := NULL] |>
        _[, res_1 := .I] |>
        melt(
          id.vars = "res_1",
          measure.vars = 1:(ncol(pxe) - 1),
          variable.name = "res_2",
          value.name = "pxe"
        ) |>
        _[, res_2 := as.integer(res_2)] |>
        _[, model := m + 1]
    }
  ) |>
    rbindlist()
}


# =====================================================================
#  Plot PAE / PDE heatmap faceted by model
# =====================================================================

plot_pxe <- function(pxe, type = c("pae", "pde"), ligand_number = NULL) {
  """
  Generate and save faceted PAE or PDE heatmaps, one PNG per region.

  If *ligand_number* is provided, the matrix is split into three regions:
    - `DNA`    : residues 1..ligand_number-1 x 1..ligand_number-1
    - `Ligand` : residues ligand_number..N x ligand_number..N
    - `other`  : cross-region entries (not plotted)
  When *ligand_number* is NULL, the full matrix is plotted as a single
  "all" region.

  PAE uses the khroma `nuuk` scale; PDE uses reversed `batlowW` scale.

  Returns a list of ggplot objects (invisible).  Saves PNGs as a side effect.
  """
  title <- names(pxe)

  # Extract model / experiment from the list name
  model <- gsub("/.*", "", title)
  experiment <- gsub(".*/", "", title)
  title <- paste0(model, " [", experiment, "]")
  cat(paste0("Plotting ", toupper(type), " matrices for ", title, "...\n"))

  data_to_plot <- melt_pxe(pxe[[1]])

  # Split into regions if a ligand_number boundary is given
  if (!is.null(ligand_number)) {
    data_to_plot <- data_to_plot[, region := "all"] |>
      _[,
        region := fcase(
          res_1 < ligand_number & res_2 < ligand_number , "DNA"    ,
          res_1 >= ligand_number & res_2 >= ligand_number , "Ligand" ,
          default = "other"
        )
      ]
  } else {
    data_to_plot[, region := "all"]
  }

  # Plot each region separately, excluding "other"
  lapply(
    unique(data_to_plot$region)[unique(data_to_plot$region) != "other"],
    function(r) {
      pxe_plot <- ggplot(
        data = data_to_plot[region == r],
        aes(x = res_1, y = res_2, fill = pxe)
      ) +
        facet_wrap(~model, ncol = 5) +
        {
          if (type == "pae") {
            scale_fill_nuuk(name = "Predicted aligned error (Å)")
          } else {
            scale_fill_batlowW(
              name = "Predicted distance error (Å)",
              transform = 'reverse'
            )
          }
        } +
        scale_x_continuous(
          name = if (r == 'DNA') "Scored residue" else "Scored ligand atom",
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          name = if (r == 'DNA') "Aligned residue" else "Aligned ligand atom",
          expand = c(0, 0),
          transform = 'reverse'
        ) +
        geom_raster() +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          plot.subtitle = element_text(size = 16, face = "italic"),
          panel.grid = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold"),
          strip.text = element_text(size = 18, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          legend.position = 'bottom',
          legend.ticks = element_line(colour = 'white', linewidth = 1)
        ) +
        guides(
          fill = guide_colorbar(
            barwidth = 20,
            barheight = 1,
            title.position = "top",
            title.hjust = 0.5
          )
        ) +
        labs(
          title = title,
          subtitle = paste0("Region: ", r)
        )

      # 4. Save to PNG
      plot_file_path <- paste0(project, "/", model, "/", experiment, "/")
      plot_file_name <- paste0(
        gsub("/", "_", title), "_", r, "_", type, ".png"
      )
      cat(paste0(
        "Saving ", toupper(type), " plot for region ", r,
        " as ", plot_file_name, " in ", plot_file_path, "\n"
      ))
      ggsave(
        plot = pxe_plot,
        filename = paste0(plot_file_path, plot_file_name),
        width = 12,
        height = 10
      )
      pxe_plot
    }
  )
}


# =====================================================================
#  pLDDT per-residue line plots
# =====================================================================

plot_plddt <- function(project) {
  """
  Generate pLDDT per-residue line plots for each experiment, faceted by
  model (5 columns).

  For each pLDDT CSV file the function:
    1. Reads the corresponding YAML input to extract the DNA sequence.
    2. Maps numeric residue indices to labelled residues (e.g. `G1`, `G2`).
    3. Plots a line + colour-point trace of pLDDT vs. residue per model.

  Returns a list of ggplot objects (invisible).  Saves PNGs as a side effect.
  """

  # 1. Find all pLDDT CSV files
  plddt_files <- list.files(
    path = paste0(project),
    pattern = "plddt.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  cat(paste0(
    "Found ", length(plddt_files), " pLDDT CSV files in ", project, "/...\n"
  ))

  # 2. Read all YAML files to extract DNA sequences
  yaml_files <- list.files(
    paste0(project, '/yaml/'),
    pattern = '\\.yaml$',
    full.names = TRUE
  )
  seq_names <- lapply(
    yaml_files,
    function(f) {
      basename(f) |>
        strsplit('_') |>
        unlist() |>
        head(1)
    }
  ) |>
    unique()
  seq_list <- lapply(yaml_files, function(f) {
    read_yaml(f)$sequences[[1]]$dna$sequence
  }) |>
    unique() |>
    setNames(seq_names)
  cat("Found ", length(seq_list), " unique sequences in YAML files:\n")
  cat(paste0("  - ", names(seq_list), ": ", seq_list, "\n", collapse = ""))

  # 3. Generate one plot per pLDDT CSV file
  plots <- lapply(
    plddt_files,
    function(f) {
      cat(paste0("Processing ", f, "...\n"))

      # Parse file path for metadata
      seq_name <- strsplit(strsplit(f, '/')[[1]][2], '_')[[1]][1]
      lgd_name <- strsplit(strsplit(f, '/')[[1]][2], '_')[[1]][2]
      rst <- strsplit(strsplit(f, '/')[[1]][2], '_')[[1]][3]
      job <- strsplit(f, "_")[[1]] |>
        grep("^J\\d+$", x = _, value = TRUE)

      # Get the DNA sequence and create labelled residue vector (e.g. "G1", "G2")
      seq <- seq_list[[seq_name]]
      seq_vec <- strsplit(seq, "")[[1]]
      names(seq_vec) <- seq_along(seq_vec)
      seq_vector <- paste0(seq_vec, names(seq_vec), sep = "")
      cat("The sequence is: ", paste(seq_vector, collapse = ""), "\n")
      seq_length <- length(seq_vec)
      cat(paste0("The sequence length is: ", seq_length, "\n"))

      # Read pLDDT data, filter to sequence residues only, add labels
      fread(f) |>
        _[, res := seq_len(.N), by = model] |>
        _[res <= seq_length] |>
        _[, residue := seq_vector[res]] |>
        _[, residue := factor(residue, levels = seq_vector)] |>
        _[, model := factor(model + 1)] |>
        ggplot(aes(x = residue, y = pLDDT, group = model)) +
        geom_line() +
        geom_point(aes(color = pLDDT), size = 3) +
        facet_wrap(~model, ncol = 5) +
        scale_x_discrete(
          name = "Residue",
          breaks = seq_vector[seq(1, length(seq_vector), by = 5)],
          labels = seq_vector[seq(1, length(seq_vector), by = 5)],
          expand = c(0.05, 0.05)
        ) +
        scale_color_gradientn(
          name = "pLDDT",
          colors = c("#8B0000", "#DAA520", "#87CEEB", "#191970"),
          values = scales::rescale(
            c(0.5, 0.6, 0.7, 0.9),
            to = c(0, 1),
            from = c(0.5, 0.9)
          ),
          limits = c(0.5, 0.9),
          oob = scales::squish
        ) +
        theme_minimal() +
        theme(
          panel.grid.major = element_line(
            colour = 'grey', linewidth = 0.5, linetype = 'dashed'
          ),
          panel.grid.minor = element_line(
            colour = 'grey', linewidth = 0.5, linetype = 'dashed'
          ),
          axis.line = element_line(colour = 'black', linewidth = 1),
          axis.text.x = element_text(
            size = 18, angle = -90, hjust = 0, vjust = 0
          ),
          axis.text.y = element_text(size = 18),
          axis.ticks = element_line(colour = 'black', linewidth = 1),
          axis.title = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 20, face = "bold"),
          plot.subtitle = element_text(size = 16, face = "italic"),
          strip.text = element_text(size = 20, face = "bold"),
          legend.position = 'none',
          legend.ticks = element_line(colour = 'white', linewidth = 1)
        ) +
        labs(
          title = {
            title <- paste0("pLDDT per residue for ", seq_name)
            if (lgd_name != 'free' & lgd_name != 'Free') {
              title <- paste0(title, " + ", lgd_name)
            } else {
              title <- paste0(title, " (unbound)")
            }
            if (!is.na(rst)) {
              title <- paste0(title, " [", rst, "]")
            }
            title
          },
          subtitle = job,
          x = "Residue",
          y = "pLDDT"
        )

      # Handle optional constraint suffix
      if (is.na(rst)) {
        rst <- ""
      } else {
        rst <- paste0("_", rst)
      }
      cat(paste0(
        "Saving pLDDT plot for ", seq_name, " + ", lgd_name, rst,
        " as ", seq_name, "_", lgd_name, rst, "_plddt.png in ",
        project, "/", seq_name, "_", lgd_name, rst, "/...\n"
      ))

      # Save the plot
      ggsave(
        filename = paste0(
          project, "/", seq_name, "_", lgd_name, rst, "/", job, "/",
          seq_name, "_", lgd_name, "_", rst, "_plddt.png"
        ),
        width = 12,
        height = 10
      )
    }
  )
  plots
}