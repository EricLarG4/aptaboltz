# install_packages.R
# ------------------
# One-time setup script for the boltz project.
# Installs required R packages and creates the 'boltz' conda environment
# with all Python dependencies (defined in environment.yml).
#
# Usage:
#   source("install_packages.R")
#
# This is idempotent — safe to run multiple times.

required_r_pkgs <- c(
  "data.table",
  "ggplot2",
  "DT",
  "khroma",
  "jsonlite",
  "yaml",
  "bslib",
  "reticulate"
)

# Install any missing R packages from CRAN
missing_r <- required_r_pkgs[!sapply(required_r_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_r) > 0) {
  message("Installing missing R packages: ", paste(missing_r, collapse = ", "))
  install.packages(missing_r, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}

# Create / update the 'boltz' conda environment from environment.yml
if (!"boltz" %in% reticulate::conda_list()$name) {
  message("Creating 'boltz' conda environment from environment.yml ...")
  reticulate::conda_create(
    envname = "boltz",
    environment = "environment.yml"
  )
  message("'boltz' conda environment created.")
} else {
  message("'boltz' conda environment already exists.")
}

# Tell reticulate to use the boltz environment
reticulate::use_condaenv("boltz", required = TRUE)
message("reticulate is now using the 'boltz' conda environment.")