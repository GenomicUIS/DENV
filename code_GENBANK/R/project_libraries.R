# Per-project R library management. Packages live under .Rlib/ so runs are
# self-contained and do not touch the user's global R library.

pipeline_packages <- function() {
  c("httr2", "xml2", "openxlsx", "ggplot2", "dplyr", "tidyr", "patchwork", "sf", "rnaturalearth", "rnaturalearthdata", "ggridges", "ggalluvial")
}

project_library_paths <- function(project_root) {
  c(file.path(project_root, ".Rlib"),
    file.path(project_root, "r_libs"))
}

activate_project_libraries <- function(project_root) {
  paths <- project_library_paths(project_root)
  paths <- paths[dir.exists(paths)]
  if (length(paths) == 0L) return(invisible(.libPaths()))
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
  .libPaths(unique(c(paths, .libPaths())))
  invisible(.libPaths())
}

#' Install any missing packages into the project library.
#' Returns the (possibly empty) character vector of newly installed packages.
ensure_packages <- function(packages,
                            project_root,
                            repos = Sys.getenv("CRAN_REPO",
                                               "https://cloud.r-project.org"),
                            lib = file.path(project_root, ".Rlib")) {
  packages <- unique(packages[nzchar(packages)])
  if (length(packages) == 0L) return(invisible(character(0)))
  
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)
  .libPaths(unique(c(lib, .libPaths())))
  
  missing <- packages[!vapply(packages, requireNamespace, logical(1L),
                              quietly = TRUE)]
  if (length(missing) == 0L) return(invisible(character(0)))
  
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(as.numeric(old_timeout), 600))
  
  for (pkg in missing) {
    message(sprintf("Installing '%s' into %s", pkg, lib))
    utils::install.packages(pkg, lib = lib, repos = repos)
  }
  invisible(missing)
}

#' Ensure the pipeline packages are installed and available.
load_pipeline_packages <- function(project_root,
                                   packages = pipeline_packages()) {
  activate_project_libraries(project_root)
  ensure_packages(packages, project_root)
  invisible(TRUE)
}