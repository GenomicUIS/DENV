# 000_run.R
# Calls the main scripts to execute steps 00 to 04 in order.

# Before running the full workflow with these simple commands, 
# verify that you are in the correct directory (folder); 
# you can check this using `getwd()`.

# Important: You need a GenBank API key (log in to GenBank and obtain one first);
# enter your email and API key into the "ncbi_client.R" file to ensure faster download speeds.

# Run the complete pipeline as a clean system process
source("scripts/00_setup.R")
Sys.setenv(DENV_DEDUPE_SEQUENCE = "false")
source("scripts/01_download.R")
source("scripts/02_qc_tiers.R")
source("scripts/03_report.R")
source("scripts/04_plot_results.R")
source("scripts/05_advanced_plots.R")
