# dependencies.R
# Required packages for all R analysis scripts in this directory
#
# This file centralizes package management for reproducibility
# Run this script before executing any analysis scripts

# Required packages for all analysis scripts
required_packages <- c(
  # Core tidyverse and data manipulation
  "tidyverse",
  "dplyr",
  "readr",
  "tibble",
  "tidyr",
  "stringr",
  "rlang",

  # Plotting and visualization
  "ggplot2",
  "patchwork",
  "gplots",
  "ggsignif",
  "ggrepel",
  "ggtext",
  "scales",

  # Statistical modeling
  "mgcv",
  "broom",
  "rstatix",
  "DHARMa",

  # Utilities
  "forcats",
  "progress",
  "gt",

  # Path management
  "here"
)

# Check which packages are not installed
missing_packages <- setdiff(required_packages, rownames(installed.packages()))

# Install missing packages
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
} else {
  cat("All required packages are already installed.\n")
}

# Load all packages
cat("Loading packages...\n")
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat("  -", pkg, "\n")
}))

cat("\nAll packages loaded successfully!\n")
cat("R version:", R.version.string, "\n")
