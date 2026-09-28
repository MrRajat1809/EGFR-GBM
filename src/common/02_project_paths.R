# .egfr_src is located by the calling script before this file is sourced.
PROJECT_ROOT <- normalizePath(file.path(.egfr_src, ".."), mustWork = TRUE)

project_path <- function(relative_path) {
  file.path(PROJECT_ROOT, relative_path)
}

use_working_directory <- function(relative_path) {
  directory <- project_path(relative_path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  setwd(directory)
}
