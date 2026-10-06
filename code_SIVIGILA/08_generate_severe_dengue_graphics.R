source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## GENERATION OF GRAPHICS AND MAPS FOR SEVERE DENGUE
## Purpose: Create the complete series of visualizations (trends, maps, 
## pyramids, timeliness of care, and quality) using the deduplicated base.
## -----------------------------------------------------------------------------

## 1. Verify that the necessary packages to process and plot are installed.
packages <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "sf", "lwgeom")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    ". Install them with install.packages()."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(sf)
  library(lwgeom)
})

## Load support scripts for figure formatting and departmental maps.
source("figures_quality.R", encoding = "UTF-8")
source("departmental_maps.R", encoding = "UTF-8")

## 2. Load the analytical database without duplicates.
data_path <- file.path(
  SIVIGILA_PATH, "resultados", "dengue_grave", "datos_procesados",
  "dengue_grave_analitica_sin_duplicados.rds"
)
if (!file.exists(data_path)) {
  stop("The deduplicated base does not exist. Execute the deduplication script first.")
}
datos <- readRDS(data_path)
stopifnot(
  nrow(datos) == valor_deduplicacion("dengue_grave", INDICADORES_DEDUP[["base_deduplicada"]])
)

## 3. Create the organized folders to classify each type of figure.
figures_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue_grave", "figuras")
directories <- file.path(
  figures_dir,
  c(
    "01_temporales", "02_demograficas", "03_geograficas",
    "04_atencion", "05_sociales", "06_desenlace", "07_calidad",
    "_previsualizaciones"
  )
)
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))

## 4. Define official color palette for the thesis (legible in print).
blue <- "#0072B2"
light_blue <- "#56B4E9"
green <- "#009E73"
orange <- "#D55E00"
yellow <- "#E69F00"
purple <- "#CC79A7"
gray <- "#6B7280"

## Unified graphic theme for all figures.
thesis_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(face = "bold", size = rel(1.25), colour = "#17202A"),
      plot.subtitle = element_text(colour = "#475569", margin = margin(b = 8)),
      plot.caption = element_text(colour = "#64748B", hjust = 0, size = rel(0.78)),
      axis.title = element_text(colour = "#334155"),
      axis.text = element_text(colour = "#334155"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E2E8F0", linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", colour = "#334155"),
      plot.margin = margin(12, 16, 12, 12)
    )
}

wrap_text <- function(x, width = 28) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

## 5. Auxiliary function to save figures in TIFF format (high resolution),
## PNG (for quick views) and register metadata in a general index.
manifest <- list()
save_figure <- function(plot, subfolder, file_name, width, height, title_text, sample_size, note) {
  path <- file.path(figures_dir, subfolder, file_name)
  ggsave(
    filename = path, plot = plot, device = "tiff",
    width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white", limitsize = FALSE
  )
  png_path <- file.path(
    figures_dir, "_previsualizaciones",
    sub("\\.tiff$", ".png", file_name, ignore.case = TRUE)
  )
  ggsave(
    filename = png_path, plot = plot, device = "png",
    width = width, height = height, units = "in", dpi = 150,
    bg = "white", limitsize = FALSE
  )
  manifest[[length(manifest) + 1L]] <<- data.frame(
    orden = length(manifest) + 1L,
    categoria = subfolder,
    archivo_tiff = normalizePath(path, winslash = "/", mustWork = TRUE),
    titulo = title_text,
    muestra = sample_size,
    nota = note,
    ancho_pulgadas = width,
    alto_pulgadas = height,
    dpi = 600L,
    compresion = "LZW",
    stringsAsFactors = FALSE
  )
  guardar_evidencia_figura(plot, path, width, height)
  invisible(path)
}

caption_base <- paste0(
  "Source: SIVIGILA, severe dengue (event 220), 2007-2025. ",
  "Deduplicated dataset; n = ", format(nrow(datos), big.mark = ",", decimal.mark = "."), "."
)

## -----------------------------------------------------------------------------
## 6. CONSTRUCTION OF THE 11 THESIS FIGURES
## -----------------------------------------------------------------------------

# Figure 1. Annual series of reported cases
annual_series <- datos %>% count(ano, name = "casos") %>% arrange(ano)

p_series <- ggplot(annual_series, aes(ano, casos)) +
  geom_line(colour = blue, linewidth = 1) +
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.8, stroke = 0.9) +
  geom_label(label.padding = grid::unit(0.07, "lines"), linewidth = 0, fill = "white",
             aes(label = label_number(big.mark = ",", decimal.mark = ".")(casos)),
             vjust = -0.75, size = 3, colour = "#334155"
  ) +
  scale_x_continuous(breaks = annual_series$ano) +
  scale_y_continuous(
    labels = label_number(big.mark = ",", decimal.mark = "."),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  labs(
    title = "Severe dengue reported per year",
    subtitle = "Number of event 220 records; corresponds to counts, not population rates",
    x = "Epidemiological year", y = "Reported severe cases", caption = caption_base
  ) +
  thesis_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(
  p_series, "01_temporales", "01_serie_anual_dengue_grave.tiff", 11, 6.5,
  "Severe dengue reported per year", nrow(datos),
  "Counts; denominators are required for incidence."
)

# Figure 2. Seasonality heatmap (Week × Year)
heatmap_data <- datos %>%
  count(ano, semana, name = "casos") %>%
  complete(ano = SIVIGILA_ANIOS, semana = 1:53, fill = list(casos = 0L)) %>%
  mutate(ano_factor = factor(ano, levels = rev(SIVIGILA_ANIOS)))

p_heatmap <- ggplot(heatmap_data, aes(semana, ano_factor, fill = casos)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", breaks = pretty_breaks(5)) +
  scale_x_continuous(breaks = c(1, seq(4, 52, 4), 53), expand = c(0, 0)) +
  labs(
    title = "Seasonality of severe dengue",
    subtitle = "Reported cases by epidemiological week and year",
    x = "Epidemiological week", y = "Year", fill = "Cases", caption = caption_base
  ) +
  thesis_theme() +
  theme(panel.grid = element_blank(), legend.key.width = grid::unit(1.3, "cm"))

save_figure(
  p_heatmap, "01_temporales", "02_mapa_calor_semana_anio.tiff", 11, 7.5,
  "Seasonality per epidemiological week", nrow(datos),
  "Color scale with square root to preserve contrast."
)

# Figure 3. Demographic profile by age group and sex
age_sex <- datos %>%
  filter(!is.na(grupo_edad_tesis), sexo_etiqueta %in% c("Femenino", "Masculino")) %>%
  count(grupo_edad_tesis, sexo_etiqueta, name = "casos") %>%
  mutate(
    valor = if_else(sexo_etiqueta == "Femenino", -casos, casos),
    grupo_edad_tesis = factor(
      grupo_edad_tesis,
      levels = c("<1", "1-4", "5-14", "15-24", "25-44", "45-64", "65+")
    )
  )

p_age_sex <- ggplot(age_sex, aes(grupo_edad_tesis, valor, fill = sexo_etiqueta)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = label_number(big.mark = ",", decimal.mark = ".")(casos), hjust = if_else(valor < 0, 1.08, -0.08)),
    size = 3, colour = "#334155"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("Femenino" = purple, "Masculino" = blue)) +
  scale_y_continuous(
    labels = function(x) label_number(big.mark = ",", decimal.mark = ".")(abs(x)),
    expand = expansion(mult = c(0.16, 0.16))
  ) +
  labs(
    title = "Severe cases by age group and sex",
    subtitle = "Age was converted to years using both EDAD and UNI_MED collectively",
    x = "Age group (years)", y = "Reported severe cases",
    fill = "Sex", caption = caption_base
  ) +
  thesis_theme()

save_figure(
  p_age_sex, "02_demograficas", "03_edad_sexo.tiff", 10, 7,
  "Severe dengue by age and sex", sum(age_sex$casos),
  "Includes female or male sex and interpretable age."
)

# Figures 4 and 5. Departmental maps of residence and occurrence
dept_map <- cargar_cartografia_departamental()

create_map <- function(variable_code, title_text, subtitle_text) {
  counts <- datos %>%
    filter(!is.na(.data[[variable_code]]), .data[[variable_code]] != "00") %>%
    count(codigo_dpto = .data[[variable_code]], name = "casos")
  map_data <- dept_map %>%
    left_join(counts, by = "codigo_dpto") %>%
    mutate(casos = replace_na(casos, 0L))
  unmapped <- sum(
    is.na(datos[[variable_code]]) |
      datos[[variable_code]] == "00" |
      !datos[[variable_code]] %in% dept_map$codigo_dpto
  )
  plot <- dibujar_mapa_departamental(map_data, "casos", title_text, subtitle_text,
                                     wrap_text(paste0(caption_base, " Unmapped records: ", unmapped,
                                                      ". Cartography: DANE, MGN 2025."), 105), "Cases", nrow(datos))
  list(grafico = plot, no_mapeados = unmapped)
}

map_residence <- create_map(
  "cod_dpto_r",
  "Severe dengue according to department of residence",
  "Reported counts; resident population by department and year is required for rates"
)
save_figure(
  map_residence$grafico, "03_geograficas", "04_mapa_residencia_departamento.tiff",
  8.5, 9, "Map by department of residence", nrow(datos),
  paste("Counts, not rates. Unmapped:", map_residence$no_mapeados)
)

map_occurrence <- create_map(
  "cod_dpto_o",
  "Severe dengue according to department of occurrence",
  "Occurrence should not be confused with residence or notification location"
)
save_figure(
  map_occurrence$grafico, "03_geograficas", "05_mapa_ocurrencia_departamento.tiff",
  8.5, 9, "Map by department of occurrence", nrow(datos),
  paste("Counts, not rates. Unmapped:", map_occurrence$no_mapeados)
)

# Figure 6. Timeliness of medical care
period_data <- datos %>%
  mutate(periodo = cut(
    ano,
    breaks = c(2006, 2012, 2018, 2025),
    labels = c("2007-2012", "2013-2018", "2019-2025")
  ))

interval_labels <- c(
  dias_inicio_consulta = "Symptoms → consultation",
  dias_inicio_hospitalizacion = "Symptoms → hospitalization",
  dias_inicio_notificacion = "Symptoms → notification"
)

times_data <- period_data %>%
  select(periodo, all_of(names(interval_labels))) %>%
  pivot_longer(-periodo, names_to = "intervalo", values_to = "dias") %>%
  filter(!is.na(dias), dias >= 0) %>%
  group_by(periodo, intervalo) %>%
  summarise(
    n = n(), mediana = median(dias),
    p25 = quantile(dias, 0.25), p75 = quantile(dias, 0.75),
    .groups = "drop"
  ) %>%
  mutate(intervalo = factor(
    interval_labels[intervalo], levels = unname(interval_labels)
  ))

p_times <- ggplot(times_data, aes(periodo, mediana, colour = periodo)) +
  geom_errorbar(aes(ymin = p25, ymax = p75), width = 0.14, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0("Median: ", label_number()(mediana))), vjust = -0.9, size = 3) +
  facet_wrap(~intervalo, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(
    "2007-2012" = blue, "2013-2018" = green, "2019-2025" = orange
  )) +
  labs(
    title = "Timeliness of severe dengue care",
    subtitle = "Median and interquartile range by period",
    x = "Period", y = "Days", colour = "Period",
    caption = paste0(
      caption_base,
      " Missing or negative values are excluded only from the corresponding interval; ",
      "vertical scales vary between panels."
    )
  ) +
  thesis_theme() +
  theme(legend.position = "none")

save_figure(
  p_times, "04_atencion", "06_oportunidad_atencion.tiff", 13, 5.8,
  "Timeliness of care", nrow(datos),
  "Medians and IQR; missing and negative values are excluded per interval."
)

# Figure 7. Social profile, insurance, and ethnic belonging
distribution_plot <- function(variable, title_text, color) {
  table_data <- datos %>%
    mutate(categoria = replace_na(as.character(.data[[variable]]), "Sin información")) %>%
    count(categoria, name = "n") %>%
    mutate(porcentaje = n / sum(n)) %>%
    arrange(porcentaje) %>%
    mutate(categoria = factor(wrap_text(categoria, 23), levels = wrap_text(categoria, 23)))
  ggplot(table_data, aes(categoria, porcentaje)) +
    geom_col(fill = color, width = 0.72) +
    geom_text(aes(label = porcentaje_es(porcentaje, accuracy = 0.1)), hjust = -0.08, size = 3) +
    coord_flip() +
    scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.18))) +
    labs(title = title_text, x = NULL, y = "Percentage") +
    thesis_theme(10) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

datos$etnia_grafica <- ifelse(
  datos$pertenencia_etnica_etiqueta %in% c("ROM/Gitano", "Raizal", "Palenquero"),
  "Otros grupos étnicos minoritarios",
  as.character(datos$pertenencia_etnica_etiqueta)
)
p_area <- distribution_plot("area_etiqueta", "Area of origin", blue)
p_regimen <- distribution_plot("regimen_salud_etiqueta", "Health insurance regime", green)
p_ethnic <- distribution_plot("etnia_grafica", "Ethnic belonging", orange)
p_social <- (p_area | p_regimen | p_ethnic) +
  plot_annotation(
    title = "Territorial, insurance, and ethnic belonging profile",
    subtitle = "Percentage distribution of reported severe cases",
    caption = caption_base,
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"),
      plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"),
      plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0)
    )
  )

save_figure(
  p_social, "05_sociales", "07_perfil_social.tiff", 15, 7.5,
  "Social profile of severe dengue", nrow(datos),
  "Descriptive percentages; they do not represent relative risks."
)

# Figure 8. Hospitalization and temporal evolution
annual_hospitalization <- datos %>%
  group_by(ano) %>%
  summarise(
    n = n(), hospitalizados = sum(pac_hos == "1", na.rm = TRUE),
    porcentaje = hospitalizados / n, .groups = "drop"
  )

p_hosp_annual <- ggplot(annual_hospitalization, aes(ano, porcentaje)) +
  geom_line(colour = blue, linewidth = 0.9) +
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.5, stroke = 0.8) +
  scale_x_continuous(breaks = seq(2007, 2025, 2)) +
  scale_y_continuous(labels = porcentaje_es, limits = c(0, 1)) +
  labs(title = "Hospitalization per year", x = "Year", y = "Hospitalized cases") +
  thesis_theme()

hosp_delay <- datos %>%
  filter(pac_hos == "1", !is.na(dias_inicio_hospitalizacion), dias_inicio_hospitalizacion >= 0) %>%
  mutate(
    demora = case_when(
      dias_inicio_hospitalizacion <= 2 ~ "0-2 días",
      dias_inicio_hospitalizacion <= 5 ~ "3-5 días",
      dias_inicio_hospitalizacion <= 7 ~ "6-7 días",
      TRUE ~ "8 días o más"
    ),
    demora = factor(demora, levels = c("0-2 días", "3-5 días", "6-7 días", "8 días o más"))
  ) %>%
  count(demora, name = "n") %>%
  mutate(porcentaje = n / sum(n))

p_hosp_delay <- ggplot(hosp_delay, aes(demora, porcentaje)) +
  geom_col(fill = green, width = 0.68) +
  geom_text(
    aes(label = paste0(label_number(big.mark = ",", decimal.mark = ".")(n), " (", porcentaje_es(porcentaje, accuracy = 0.1), ")")),
    vjust = -0.5, size = 3.1
  ) +
  scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.17))) +
  labs(title = "Symptoms until hospitalization", x = NULL, y = "Percentage") +
  thesis_theme() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

p_hospitalization <- (p_hosp_annual | p_hosp_delay) +
  plot_annotation(
    title = "Hospitalization of severe dengue cases",
    subtitle = "Annual trend and time since symptom onset",
    caption = caption_base,
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"),
      plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"),
      plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0)
    )
  )

save_figure(
  p_hospitalization, "04_atencion", "08_hospitalizacion_evolucion.tiff", 13, 7,
  "Hospitalization and evolution time", nrow(datos),
  paste(
    "Valid delay for",
    format(sum(hosp_delay$n), big.mark = ",", decimal.mark = "."),
    "hospitalized patients."
  )
)

# Figure 9. Annual fatal outcome within event 220
annual_outcome <- datos %>%
  group_by(ano) %>%
  summarise(
    casos = n(), fallecidos = sum(con_fin == "2", na.rm = TRUE),
    proporcion = fallecidos / casos, .groups = "drop"
  )

intervals <- t(vapply(seq_len(nrow(annual_outcome)), function(i) {
  test_result <- prop.test(
    annual_outcome$fallecidos[[i]], annual_outcome$casos[[i]], correct = FALSE
  )
  as.numeric(test_result$conf.int)
}, numeric(2)))
annual_outcome$li <- intervals[, 1]
annual_outcome$ls <- intervals[, 2]

p_outcome <- ggplot(annual_outcome, aes(ano, proporcion)) +
  geom_ribbon(aes(ymin = li, ymax = ls), fill = orange, alpha = 0.16) +
  geom_line(colour = orange, linewidth = 1) +
  geom_point(colour = orange, fill = "white", shape = 21, size = 2.7, stroke = 0.9) +
  scale_x_continuous(breaks = SIVIGILA_ANIOS) +
  scale_y_continuous(labels = porcentaje_formato_es(accuracy = 0.1), expand = expansion(mult = c(0.025, 0.08))) +
  labs(
    title = "Fatal outcome among severe dengue notifications",
    subtitle = "Annual percentage and 95% confidence interval",
    x = "Epidemiological year", y = "Deceased among severe cases",
    caption = wrap_text(paste0(
      caption_base, " This proportion does not equal the lethality of all dengue, ",
      "because the denominator contains only event 220 notifications.",
      if (any(datos$ano == 2025) &&
          all(datos$con_fin[datos$ano == 2025] == "1", na.rm = TRUE))
        " In 2025 all records in this file have CON_FIN=1; event 580 deaths are analyzed separately." else ""
    ), 115)
  ) +
  thesis_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(
  p_outcome, "06_desenlace", "09_desenlace_mortal_anual.tiff", 11, 6.8,
  "Fatal outcome within event 220", nrow(datos),
  paste("Deceased:", sum(datos$con_fin == "2", na.rm = TRUE),
        "; percentage over severe cases, not over all dengue.")
)

# Figure 10. Registered vulnerability groups
vulnerability_vars <- c(
  gp_discapa = "Discapacidad",
  gp_desplaz = "Población desplazada",
  gp_migrant = "Población migrante",
  gp_carcela = "Privación de la libertad",
  gp_gestan = "Gestante",
  gp_indigen = "Habitante de calle",
  gp_pobicfb = "Población infantil ICBF",
  gp_mad_com = "Madre comunitaria",
  gp_desmovi = "Población desmovilizada",
  gp_psiquia = "Centro psiquiátrico",
  gp_vic_vio = "Víctima de violencia armada"
)

vulnerability <- do.call(rbind, lapply(names(vulnerability_vars), function(variable) {
  x <- as.character(datos[[variable]])
  valid_responses <- x %in% c("1", "2")
  data.frame(
    grupo = unname(vulnerability_vars[[variable]]),
    n_si = sum(x == "1", na.rm = TRUE), n_validos = sum(valid_responses),
    porcentaje = if (sum(valid_responses)) sum(x == "1", na.rm = TRUE) / sum(valid_responses) else NA_real_,
    stringsAsFactors = FALSE
  )
})) %>%
  filter(n_si >= 5) %>%
  arrange(porcentaje) %>%
  mutate(grupo = factor(grupo, levels = grupo))

p_vulnerability <- ggplot(vulnerability, aes(grupo, porcentaje)) +
  geom_col(fill = purple, width = 0.7) +
  geom_text(
    aes(label = paste0(label_number(big.mark = ",", decimal.mark = ".")(n_si), " (", porcentaje_es(porcentaje, accuracy = 0.01), ")")),
    hjust = -0.08, size = 3.2
  ) +
  coord_flip() +
  scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.23))) +
  labs(
    title = "Registered vulnerability groups",
    subtitle = "Percentage among records with valid response for each indicator",
    x = NULL, y = "Percentage with Yes response",
    caption = paste0(
      caption_base,
      " For confidentiality reasons, indicators with less than five affirmative responses are not shown."
    )
  ) +
  thesis_theme()

save_figure(
  p_vulnerability, "05_sociales", "10_grupos_vulnerabilidad.tiff", 10, 7,
  "Vulnerability groups", nrow(datos),
  "Categories with fewer than five affirmative responses are omitted."
)

# Figure 11. Data quality issues by year
quality_flags <- c(
  flag_edad_fecha_nacimiento_discordante = "Discordant age/birth date",
  flag_anio_fec_not_distinto = "FEC_NOT year different from ANO",
  flag_ajuste_anterior_notificacion = "Adjustment prior to notification",
  flag_cronologia_clinica_inconsistente = "Negative clinical chronology",
  flag_hospitalizacion_inconsistente = "Inconsistent hospitalization/date",
  flag_muerto_sin_fecha_defuncion = "Deceased without date of death"
)

quality_data <- datos %>%
  select(ano, all_of(names(quality_flags))) %>%
  pivot_longer(-ano, names_to = "bandera", values_to = "marcado") %>%
  group_by(ano, bandera) %>%
  summarise(porcentaje = mean(marcado, na.rm = TRUE), .groups = "drop") %>%
  mutate(indicador = factor(
    quality_flags[bandera], levels = rev(unname(quality_flags))
  ))

p_quality <- ggplot(quality_data, aes(ano, indicador, fill = porcentaje)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "B", labels = porcentaje_es, breaks = pretty_breaks(5)) +
  scale_x_continuous(breaks = SIVIGILA_ANIOS) +
  labs(
    title = "Data quality issues by year",
    subtitle = "Percentage of records with each flag; they are not clinical outcomes",
    x = "Epidemiological year", y = NULL, fill = "Records", caption = caption_base
  ) +
  thesis_theme() +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = grid::unit(1.3, "cm")
  )

save_figure(
  p_quality, "07_calidad", "11_calidad_datos_por_anio.tiff", 12, 6.8,
  "Data quality by year", nrow(datos),
  "Quality control; do not interpret as a clinical pattern."
)

## 7. Export indexes, package versions, and README documentation file.
figures_index <- do.call(rbind, manifest)
write.csv(
  figures_index, file.path(figures_dir, "INDICE_FIGURAS.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

package_versions <- data.frame(
  paquete = packages,
  version = vapply(packages, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
write.csv(
  package_versions, file.path(figures_dir, "VERSIONES_PAQUETES.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

writeLines(c(
  "# Severe Dengue Figures",
  "",
  paste("All figures use the deduplicated base of", nrow(datos), "records."),
  "Exclusions due to matching content are documented in reportes/deduplicacion/resumen_deduplicacion.csv.",
  "TIFF files were saved at 600 dpi with LZW compression, accompanied by vector PDFs and PNGs for quick review.",
  "",
  "## Folders",
  "",
  "- 01_temporales: annual trend and seasonality.",
  "- 02_demograficas: age and sex.",
  "- 03_geograficas: residence and occurrence, always separated.",
  "- 04_atencion: timeliness and hospitalization.",
  "- 05_sociales: social profile and vulnerability.",
  "- 06_desenlace: fatal outcome within event 220.",
  "- 07_calidad: quality flags per year.",
  "- _previsualizaciones: low resolution PNG copies.",
  "",
  "Check INDICE_FIGURAS.csv for titles, sizes, samples, and methodological notes."
), file.path(figures_dir, "README.md"), useBytes = TRUE)

tiff_files <- list.files(
  figures_dir, pattern = "\\.tiff$", recursive = TRUE, full.names = TRUE
)
stopifnot(
  nrow(figures_index) == 11L,
  length(tiff_files) == 11L,
  all(file.info(tiff_files)$size > 10000)
)

cat("Severe dengue figures generated:", length(tiff_files), "\n")
cat("Directory:", normalizePath(figures_dir, winslash = "/"), "\n")
