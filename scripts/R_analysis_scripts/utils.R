# utils.R
# Shared utility functions for R analysis scripts
#
# This file contains helper functions used across multiple analysis scripts
# to reduce code duplication and improve maintainability

library(here)
library(readr)

#' Read metadata CSV file from the assemblies directory
#'
#' @param filename Name of the CSV file (e.g., "metadata_updated_ecoli.csv")
#' @param show_col_types Logical, whether to show column types (default: FALSE)
#' @return A tibble containing the metadata
#' @export
read_metadata <- function(filename, show_col_types = FALSE) {
  file_path <- here("assemblies", filename)
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  read_csv(file_path, show_col_types = show_col_types)
}

#' Create output directory if it doesn't exist
#'
#' @param dir_path Directory path (can be relative or absolute)
#' @param recursive Logical, create parent directories if needed (default: TRUE)
#' @return NULL (invisibly)
#' @export
ensure_dir <- function(dir_path, recursive = TRUE) {
  dir.create(dir_path, recursive = recursive, showWarnings = FALSE)
  invisible(NULL)
}

#' Convert p-values to significance stars
#'
#' @param p P-value(s)
#' @return Character vector with significance stars
#' @export
star <- function(p) {
  sapply(p, function(x) {
    if (x < 0.001) {
      "***"
    } else if (x < 0.01) {
      "**"
    } else if (x < 0.05) {
      "*"
    } else {
      "ns"
    }
  })
}
