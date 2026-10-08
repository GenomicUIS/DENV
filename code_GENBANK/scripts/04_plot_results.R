#!/usr/bin/env Rscript
# 04_plot_results.R — Generates high-resolution TIFF plots of DENV genomes

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
library(patchwork)
library(sf)
library(rnaturalearth)

# ---- configuration --------------------------------------------------------
summary_dir <- file.path(project_root, "results", "summary")
plots_dir   <- file.path(project_root, "results", "plots")
ensure_dir(plots_dir)

metadata_path <- file.path(summary_dir, "metadata_DENV_south_america_up_to_2026.csv")
if (!file.exists(metadata_path)) stop("Metadata not found.")
df <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Custom highly visible palette for serotypes
denv_palette <- c("DENV1" = "#E63946", # Bold Red
                  "DENV2" = "#457B9D", # Steel Blue
                  "DENV3" = "#2A9D8F", # Teal/Green
                  "DENV4" = "#F4A261") # Sandy Orange

# Custom theme for clean aesthetics and centered titles
theme_manuscript <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold")
    )
}

# ---- Plot 1: Bars by country and serotype -----------------------------------
p1 <- df %>%
  count(country, serotype) %>%
  ggplot(aes(x = reorder(country, -n), y = n, fill = serotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = denv_palette) +
  theme_manuscript() +
  labs(title = "DENV Sequences by Country and Serotype",
       x = "Country", y = "Number of Sequences", fill = "Serotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(plots_dir, "01_bars_country_serotype.tiff"), plot = p1,
       device = "tiff", width = 10, height = 6, dpi = 300, compression = "lzw")

# ---- Plot 2: Time series (Faceted for clarity) ---------------------------
p2 <- df %>%
  filter(!is.na(year)) %>%
  count(year, serotype) %>%
  ggplot(aes(x = year, y = n, color = serotype, group = serotype)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = denv_palette) +
  facet_wrap(~ serotype, ncol = 1, scales = "free_y") +
  theme_manuscript() +
  labs(title = "Temporal Evolution of Sequences",
       x = "Collection Year", y = "Number of Sequences", color = "Serotype") +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"))

ggsave(file.path(plots_dir, "02_time_series_serotype.tiff"), plot = p2,
       device = "tiff", width = 8, height = 10, dpi = 300, compression = "lzw")

# ---- Plot 3: Multipanel Integration (Patchwork) ---------------------------
p_combined <- (p1 / p2) + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = 'A')
ggsave(file.path(plots_dir, "03_multipanel_statistics.tiff"), plot = p_combined,
       device = "tiff", width = 12, height = 15, dpi = 300, compression = "lzw")

# ---- Plot 4: South America Maps by Serotype -----------------------------
sa_map <- ne_countries(continent = "south america", returnclass = "sf")

map_data <- df %>%
  count(country, serotype) %>%
  complete(country, serotype = c("DENV1", "DENV2", "DENV3", "DENV4"), fill = list(n = 0))

map_joined <- sa_map %>%
  left_join(map_data, by = c("name" = "country")) %>%
  filter(!is.na(serotype))

p_map <- ggplot(map_joined) +
  geom_sf(aes(fill = n), color = "gray30", linewidth = 0.3) +
  geom_sf_text(aes(label = ifelse(n > 0, n, "")), size = 3, color = "black", fontface = "bold") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Sequences", 
                       trans = "log1p", breaks = c(0, 10, 50, 100, 500)) +
  facet_wrap(~ serotype, ncol = 2) +
  theme_void(base_size = 14) +
  labs(title = "Spatial Distribution of Dengue Serotypes in South America") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 15)),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "#f8f9fa", color = NA),
        strip.text = element_text(face = "bold", size = 14, margin = margin(b = 10)))

ggsave(file.path(plots_dir, "04_maps_serotype.tiff"), plot = p_map,
       device = "tiff", width = 10, height = 10, dpi = 300, compression = "lzw")

cat("Plots successfully generated in:", plots_dir, "\n")