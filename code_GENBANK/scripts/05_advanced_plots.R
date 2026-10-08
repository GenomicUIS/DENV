#!/usr/bin/env Rscript
# 05_advanced_plots.R — Generates advanced analytical visualizations for DENV genomes

# ---- bootstrap ------------------------------------------------------------
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

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ggalluvial)

# ---- configuration --------------------------------------------------------
summary_dir <- file.path(project_root, "results", "summary")
qc_dir      <- file.path(project_root, "results", "sequence_curation", "summary")
plots_dir   <- file.path(project_root, "results", "plots", "advanced")
ensure_dir(plots_dir)

# Load data
metadata_path <- file.path(summary_dir, "metadata_DENV_south_america_up_to_2026.csv")
qc_path       <- file.path(qc_dir, "qc_coding_region_sequences.csv")

if (!file.exists(metadata_path) || !file.exists(qc_path)) {
  stop("Required CSV files not found. Ensure previous scripts ran successfully.")
}

df_meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
df_qc   <- read.csv(qc_path, stringsAsFactors = FALSE)

# Color palettes
denv_palette <- c("DENV1" = "#E63946", "DENV2" = "#457B9D", 
                  "DENV3" = "#2A9D8F", "DENV4" = "#F4A261")

theme_manuscript <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 15)),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", margin = margin(b = 15)),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text = element_text(face = "bold", size = 12)
    )
}

# ---- Plot 1: Spatiotemporal Surveillance Heatmap ---------------------------
# Visualizes historical sequencing gaps and outbreak peaks by country
p_heat <- df_meta %>%
  filter(!is.na(year)) %>%
  count(year, country) %>%
  ggplot(aes(x = year, y = reorder(country, n, sum), fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Genomes", trans = "log1p",
                       breaks = c(1, 10, 50, 200)) +
  scale_x_continuous(breaks = seq(min(df_meta$year, na.rm=T), max(df_meta$year, na.rm=T), by = 5)) +
  theme_manuscript() +
  labs(title = "Spatiotemporal Density of DENV Genomic Surveillance",
       x = "Collection Year", y = "Country") +
  theme(legend.position = "right", panel.grid.major = element_blank())

ggsave(file.path(plots_dir, "05_heatmap_surveillance.tiff"), plot = p_heat,
       device = "tiff", width = 11, height = 7, dpi = 300, compression = "lzw")

# ---- Plot 2: Ridgeline Plot for Genome Lengths ----------------------------
# Justifies the 9000nt cutoff by showing the length distribution
p_ridge <- df_qc %>%
  filter(!is.na(coding_length), coding_length > 0) %>%
  ggplot(aes(x = coding_length, y = serotype, fill = serotype)) +
  geom_density_ridges(alpha = 0.8, color = "black", scale = 1.2, rel_min_height = 0.01) +
  geom_vline(xintercept = 9000, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = denv_palette) +
  theme_manuscript() +
  labs(title = "Distribution of Extracted CDS Lengths by Serotype",
       subtitle = "Red dashed line indicates the 9000 nt automated exclusion threshold",
       x = "Coding Sequence Length (Nucleotides)", y = "Serotype") +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, "06_ridgeline_genome_length.tiff"), plot = p_ridge,
       device = "tiff", width = 9, height = 6, dpi = 300, compression = "lzw")

# ---- Plot 3: Sankey/Alluvial Diagram for QC Methodology -------------------
# Tracks the fate of sequences through the curation tiers
qc_flow <- df_qc %>%
  mutate(
    Stage1 = "L1: Raw CDS",
    Stage2 = ifelse(level2_pass, "L2: No Ns/Gaps", "Dropped"),
    Stage3 = ifelse(level3_pass, "L3: Pure ACGT", ifelse(level2_pass, "Dropped", "Dropped")),
    Stage4 = ifelse(level3_pass & !auto_flag_required, "L4: Retained", 
                    ifelse(level3_pass & auto_flag_required, "Auto-Flagged", Stage3))
  ) %>%
  count(Stage1, Stage2, Stage3, Stage4, serotype)

p_sankey <- ggplot(qc_flow,
                   aes(y = n, axis1 = Stage1, axis2 = Stage2, axis3 = Stage3, axis4 = Stage4)) +
  geom_alluvium(aes(fill = serotype), width = 1/12, alpha = 0.7, color = "white") +
  geom_stratum(width = 1/4, fill = "gray20", color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "white", size = 3.5, fontface = "bold") +
  scale_fill_manual(values = denv_palette) +
  scale_x_discrete(limits = c("Level 1\nExtraction", "Level 2\nAmbiguities", "Level 3\nCharacters", "Level 4\nFinal Filter"), expand = c(0.05, 0.05)) +
  theme_minimal(base_size = 14) +
  labs(title = "Sequence Curation Pipeline Flow",
       y = "Number of Sequences", fill = "Serotype") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank())

ggsave(file.path(plots_dir, "07_alluvial_qc_flow.tiff"), plot = p_sankey,
       device = "tiff", width = 11, height = 7, dpi = 300, compression = "lzw")

# ---- Plot 4: Seasonality in Polar Coordinates -----------------------------
# Evaluates seasonal bias in sample collection
polar_df <- df_meta %>%
  mutate(month = case_when(
    grepl("^[0-9]{4}-[0-9]{2}", collection_date) ~ as.integer(substr(collection_date, 6, 7)),
    TRUE ~ NA_integer_
  )) %>%
  filter(!is.na(month), month >= 1, month <= 12) %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb)) %>%
  count(month_name, serotype)

p_polar <- ggplot(polar_df, aes(x = month_name, y = n, fill = serotype)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = denv_palette) +
  coord_polar(start = -pi/12) +
  theme_manuscript() +
  labs(title = "Seasonal Distribution of Sample Collections",
       subtitle = "Aggregated monthly counts across all years",
       x = NULL, y = NULL, fill = "Serotype") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"))

ggsave(file.path(plots_dir, "08_polar_seasonality.tiff"), plot = p_polar,
       device = "tiff", width = 8, height = 8, dpi = 300, compression = "lzw")

# ---- Plot 5: Host Distribution (Waffle/Stacked Bar) -----------------------
# Cleans the messy host metadata to compare Clinical vs Vector sequences
host_df <- df_meta %>%
  mutate(host_clean = case_when(
    grepl("homo|human|paciente|patient", tolower(host)) ~ "Human (Clinical)",
    grepl("aedes|mosquito|vector", tolower(host)) ~ "Mosquito (Vector)",
    grepl("monkey|macaque|nhp", tolower(host)) ~ "Non-Human Primate",
    TRUE ~ "Unknown / Not Specified"
  )) %>%
  count(country, host_clean) %>%
  group_by(country) %>%
  mutate(prop = n / sum(n))

p_host <- ggplot(host_df, aes(x = reorder(country, -n, sum), y = prop, fill = host_clean)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_manuscript() +
  labs(title = "Proportion of Sequences by Host Origin",
       x = "Country", y = "Percentage of Total Sequences", fill = "Host Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(plots_dir, "09_stacked_host_distribution.tiff"), plot = p_host,
       device = "tiff", width = 10, height = 6, dpi = 300, compression = "lzw")

cat("Advanced analytical plots successfully generated in:", plots_dir, "\n")