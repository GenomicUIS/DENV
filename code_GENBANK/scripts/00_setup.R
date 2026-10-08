#!/usr/bin/env Rscript
# 00_setup.R — install the pipeline's R packages into the project library.

project_root <- local({
  env <- Sys.getenv("DENGUE_SURAMERICA_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]),
                                         winslash = "/", mustWork = FALSE))))
  }
  probe <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(probe, "DESCRIPTION"))) {
    parent <- dirname(probe)
    if (identical(parent, probe)) stop("Cannot find project root.")
    probe <- parent
  }
  probe
})
source(file.path(project_root, "R", "utils.R"))
load_project_helpers(project_root)

installed <- ensure_packages(pipeline_packages(), project_root)

cat("Project root:", project_root, "\n")
cat("Project library:", file.path(project_root, ".Rlib"), "\n")
cat("Library paths:\n")
cat(paste0("  - ", .libPaths(), collapse = "\n"), "\n\n")

if (length(installed) == 0L) {
  cat("All pipeline packages already installed.\n")
} else {
  cat("Installed:", paste(installed, collapse = ", "), "\n")
}
for (pkg in pipeline_packages()) {
  ok  <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(pkg)) else "NA"
  cat(sprintf("  %-10s available=%-5s version=%s\n", pkg, ok, ver))
}