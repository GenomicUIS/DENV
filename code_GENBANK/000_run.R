# 000_run.R
# Calls the main scripts to execute steps 00 to 04 in order.

# Run the complete pipeline as a clean system process
source("scripts/00_setup.R")
Sys.setenv(DENV_DEDUPE_SEQUENCE = "false")
source("scripts/01_download.R")
source("scripts/02_qc_tiers.R")
source("scripts/03_report.R")
source("scripts/04_plot_results.R")
source("scripts/05_advanced_plots.R")
