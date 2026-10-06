source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## COMPARATIVE INTEGRATION: DENGUE, SEVERE DENGUE AND MORTALITY
## Purpose: Generate comparative figures by consolidating social data.
## If the labeled columns do not exist in the RDS files, they are translated 
## from the original codes.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

results_dir <- file.path(SIVIGILA_PATH, "resultados")
integrated_figures_dir <- file.path(results_dir, "integracion_eventos", "figuras")
dir.create(integrated_figures_dir, recursive = TRUE, showWarnings = FALSE)

# Verify inputs from script 14
.required_paths <- file.path(
  results_dir, "dengue", "reportes", "tablas_figuras", "08_perfil_social.csv"
)
.missing <- .required_paths[!file.exists(.required_paths)]
if (length(.missing)) {
  stop("Missing inputs. Execute R/12_generate_dengue_graphics.R first:\n",
       paste(.missing, collapse = "\n"), call. = FALSE)
}

events <- c("dengue", "dengue_grave", "mortalidad_dengue")
labels <- c(dengue = "Dengue (210)", dengue_grave = "Severe Dengue (220)", mortalidad_dengue = "Mortality (580)")
palette <- c(dengue = "#0072B2", dengue_grave = "#D55E00", mortalidad_dengue = "#8E4585")

integrated_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), colour = "#17202A"),
      plot.subtitle = element_text(colour = "#52616F", margin = margin(b = 10)),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E2E8F0", linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

## -----------------------------------------------------------------------------
## 1. Failsafe extraction for Social Profile
## -----------------------------------------------------------------------------
calculate_social_from_rds <- function(rds_path, event_name) {
  if (!file.exists(rds_path)) return(NULL)
  datos <- readRDS(rds_path)
  
  # Emergency translator if labeled columns are missing
  if (!"area_etiqueta" %in% names(datos) && "area" %in% names(datos)) {
    datos$area_etiqueta <- factor(datos$area, levels = 1:3, labels = c("Cabecera municipal", "Centro poblado", "Rural disperso"))
  }
  if (!"regimen_salud_etiqueta" %in% names(datos) && "tip_ss" %in% names(datos)) {
    datos$regimen_salud_etiqueta <- factor(datos$tip_ss, levels = c("C","S","P","E","N","I"), labels = c("Contributivo", "Subsidiado", "Excepción", "Especial", "No asegurado", "Indeterminado/Pendiente"))
  }
  if (!"pertenencia_etnica_etiqueta" %in% names(datos) && "per_etn" %in% names(datos)) {
    datos$pertenencia_etnica_etiqueta <- factor(datos$per_etn, levels = 1:6, labels = c("Indígena", "ROM/Gitano", "Raizal", "Palenquero", "Negro, mulato o afrocolombiano", "Otro"))
  }
  
  vars <- c(area_etiqueta = "Área de procedencia", 
            regimen_salud_etiqueta = "Régimen de salud", 
            pertenencia_etnica_etiqueta = "Pertenencia étnica")
  
  do.call(rbind, lapply(names(vars), function(v) {
    if (!v %in% names(datos)) return(NULL)
    cat <- as.character(datos[[v]])
    cat[is.na(cat) | cat == ""] <- "Sin información"
    if (v == "pertenencia_etnica_etiqueta") {
      cat[cat %in% c("ROM/Gitano", "Raizal", "Palenquero")] <- "Otros grupos étnicos minoritarios"
    }
    t <- as.data.frame(table(cat), stringsAsFactors = FALSE)
    names(t) <- c("categoria", "n")
    t$variable <- unname(vars[[v]])
    t$evento <- factor(event_name, levels = events)
    t
  }))
}

# Consolidate the data
df_dengue <- read.csv(file.path(results_dir, "dengue", "reportes", "tablas_figuras", "08_perfil_social.csv"), stringsAsFactors = FALSE)
df_dengue$evento <- factor("dengue", levels = events)

severe_path <- file.path(results_dir, "dengue_grave", "datos_procesados", "dengue_grave_analitica_sin_duplicados.rds")
df_severe <- calculate_social_from_rds(severe_path, "dengue_grave")

mortality_path <- file.path(results_dir, "mortalidad_dengue", "datos_procesados", "mortalidad_dengue_canonica_sin_duplicados.rds")

if (file.exists(mortality_path)) {
  mortality_data <- readRDS(mortality_path)
  mortality_data <- mortality_data[!is.na(mortality_data$con_fin) & mortality_data$con_fin == "2", ]
  temp_path <- tempfile(fileext = ".rds")
  saveRDS(mortality_data, temp_path)
  df_mortality <- calculate_social_from_rds(temp_path, "mortalidad_dengue")
  unlink(temp_path)
} else {
  df_mortality <- NULL
}

social_data <- bind_rows(df_dengue, df_severe, df_mortality)

## -----------------------------------------------------------------------------
## 2. Generation of Figure A: Comparative Social Profile
## -----------------------------------------------------------------------------
global_social <- social_data %>%
  group_by(evento, variable, categoria) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(evento, variable) %>%
  mutate(porcentaje = n / sum(n)) %>%
  ungroup()

plot_social_variable <- function(var_name) {
  df <- global_social %>% filter(variable == var_name)
  ggplot(df, aes(x = reorder(categoria, porcentaje), y = porcentaje, fill = evento)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = percent(porcentaje, accuracy = 0.1, decimal.mark = ",")),
              position = position_dodge(width = 0.8), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_manual(values = palette, labels = labels) +
    scale_y_continuous(labels = percent_format(decimal.mark = ","), limits = c(0, max(df$porcentaje) * 1.2)) +
    labs(title = var_name, x = NULL, y = "Percentage") +
    integrated_theme()
}

p_social <- (plot_social_variable("Área de procedencia") / 
               plot_social_variable("Régimen de salud") / 
               plot_social_variable("Pertenencia étnica")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(integrated_figures_dir, "01_perfil_social_comparativo.tiff"), p_social, 
       width = 10, height = 12, dpi = 600, compression = "lzw", bg = "white")

cat("Social profile integration completed. Figure saved in:", integrated_figures_dir, "\n")
