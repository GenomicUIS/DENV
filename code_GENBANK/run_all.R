#!/usr/bin/env Rscript
# run_all.R - run the pipeline end to end (00 -> 03) as separate processes.

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

steps <- c("scripts/00_setup.R", "scripts/01_download.R",
           "scripts/02_qc_tiers.R", "scripts/03_report.R")
rscript <- file.path(R.home("bin"),
                     if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# Child processes inherit the environment, so setting the root here is enough.
Sys.setenv(DENGUE_SURAMERICA_ROOT = project_root)

for (step in steps) {
  cat(sprintf("\n=== %s ===\n", step))
  status <- system2(rscript, shQuote(file.path(project_root, step)))
  if (status != 0L) {
    stop(sprintf("Step %s failed (exit code %d)", step, status))
  }
}
cat("\nAll steps completed successfully.\n")