## 05_generate_dengue_mortality_graphics.R

source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

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
source("figures_quality.R", encoding = "UTF-8")
# Direct connection to your cartography script
source("departmental_maps.R", encoding = "UTF-8")

data_path <- file.path(
  SIVIGILA_PATH, "resultados", "mortalidad_dengue", "datos_procesados",
  "mortalidad_dengue_muestra_graficas.rds"
)
if (!file.exists(data_path)) {
  stop("The deduplicated sample does not exist. Execute R/04_deduplicate_dengue_mortality.R")
}
datos <- readRDS(data_path)
# stopifnot(nrow(datos) == valor_deduplicacion("mortalidad_dengue", "Registros incluidos en las gráficas de mortalidad"), all(datos$con_fin == "2"))

figures_dir <- file.path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "figuras")
directories <- file.path(
  figures_dir,
  c(
    "01_temporales", "02_demograficas", "03_geograficas",
    "04_atencion", "05_sociales", "06_causas", "07_calidad",
    "_previsualizaciones"
  )
)
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))

blue <- "#0072B2"
light_blue <- "#56B4E9"
green <- "#009E73"
orange <- "#D55E00"
yellow <- "#E69F00"
purple <- "#CC79A7"
gray <- "#6B7280"

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

manifest <- list()
save_figure <- function(plot, subfolder, file_name, width, height, title, sample_size, note) {
  path <- file.path(figures_dir, subfolder, file_name)
  ggsave(
    filename = path,
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  png_path <- file.path(
    figures_dir, "_previsualizaciones",
    sub("\\.tiff$", ".png", file_name, ignore.case = TRUE)
  )
  ggsave(
    filename = png_path,
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 150,
    bg = "white",
    limitsize = FALSE
  )
  manifest[[length(manifest) + 1L]] <<- data.frame(
    orden = length(manifest) + 1L,
    categoria = subfolder,
    archivo_tiff = normalizePath(path, winslash = "/", mustWork = TRUE),
    titulo = title,
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
  "Source: SIVIGILA, dengue mortality, 2007–2025. ",
  "Deduplicated dataset; n = ", format(nrow(datos), big.mark = ",", decimal.mark = "."), "."
)

# 1. Annual series -------------------------------------------------------------
annual_series <- datos %>%
  count(ano, name = "muertes") %>%
  arrange(ano)

p_series <- ggplot(annual_series, aes(ano, muertes)) +
  geom_line(colour = blue, linewidth = 1) +
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.8, stroke = 0.9) +
  geom_label(label.padding = grid::unit(0.07, "lines"), linewidth = 0, fill = "white", aes(label = muertes), vjust = -0.8, size = 3.1, colour = "#334155") +
  scale_x_continuous(breaks = annual_series$ano) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.13)), labels = label_number()) +
  labs(
    title = "Reported dengue mortality per year",
    subtitle = "Counts of records classified as deceased; they do not correspond to population rates",
    x = "Epidemiological year",
    y = "Reported deaths",
    caption = caption_base
  ) +
  thesis_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(
  p_series, "01_temporales", "01_serie_anual_mortalidad.tiff", 11, 6.5,
  "Reported dengue mortality per year", nrow(datos),
  "Counts; population denominators are required for rates."
)

# 2. Seasonality per week -----------------------------------------------
heatmap_data <- datos %>%
  count(ano, semana, name = "muertes") %>%
  complete(ano = SIVIGILA_ANIOS, semana = 1:53, fill = list(muertes = 0L)) %>%
  mutate(ano_factor = factor(ano, levels = rev(SIVIGILA_ANIOS)))

p_heatmap <- ggplot(heatmap_data, aes(semana, ano_factor, fill = muertes)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", breaks = pretty_breaks(5)) +
  scale_x_continuous(breaks = c(1, seq(4, 52, 4), 53), expand = c(0, 0)) +
  labs(
    title = "Seasonality of dengue mortality",
    subtitle = "Number of deaths reported by epidemiological week and year",
    x = "Epidemiological week",
    y = "Year",
    fill = "Deaths",
    caption = caption_base
  ) +
  thesis_theme() +
  theme(panel.grid = element_blank(), legend.key.width = grid::unit(1.3, "cm"))

save_figure(
  p_heatmap, "01_temporales", "02_mapa_calor_semana_anio.tiff", 11, 7.5,
  "Seasonality per epidemiological week", nrow(datos),
  "The color scale uses a square root transformation to preserve contrast."
)
names(datos)
datos <- datos %>%
  mutate(sexo_etiqueta = case_when(
    sexo == "F" ~ "Femenino",
    sexo == "M" ~ "Masculino",
    TRUE ~ "Sin información"
  ))

# 3. Profile by age and sex ---------------------------------------------------
age_sex <- datos %>%
  filter(!is.na(grupo_edad_tesis), sexo_etiqueta %in% c("Femenino", "Masculino")) %>%
  count(grupo_edad_tesis, sexo_etiqueta, name = "muertes") %>%
  mutate(
    valor = if_else(sexo_etiqueta == "Femenino", -muertes, muertes),
    grupo_edad_tesis = factor(
      grupo_edad_tesis,
      levels = c("<1", "1-4", "5-14", "15-24", "25-44", "45-64", "65+")
    )
  )

p_age_sex <- ggplot(age_sex, aes(grupo_edad_tesis, valor, fill = sexo_etiqueta)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = muertes, hjust = if_else(valor < 0, 1.1, -0.1)),
    size = 3.2,
    colour = "#334155"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("Femenino" = purple, "Masculino" = blue)) +
  scale_y_continuous(
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  labs(
    title = "Reported deaths by age group and sex",
    subtitle = "Age was converted to years using both EDAD and UNI_MED collectively",
    x = "Age group (years)",
    y = "Number of deaths",
    fill = "Sex",
    caption = caption_base
  ) +
  thesis_theme()

save_figure(
  p_age_sex, "02_demograficas", "03_edad_sexo.tiff", 10, 7,
  "Deaths by age and sex", sum(age_sex$muertes),
  "Includes only female or male sex and interpretable age."
)

# 4. Departmental cartography ------------------------------------------------
dept_map <- cargar_cartografia_departamental()

create_map <- function(variable_code, title_text, subtitle_text) {
  counts <- datos %>%
    filter(!is.na(.data[[variable_code]]), .data[[variable_code]] != "00") %>%
    count(codigo_dpto = .data[[variable_code]], name = "muertes")
  map_data <- dept_map %>%
    left_join(counts, by = "codigo_dpto") %>%
    mutate(muertes = replace_na(muertes, 0L))
  unmapped <- sum(
    is.na(datos[[variable_code]]) |
      datos[[variable_code]] == "00" |
      !datos[[variable_code]] %in% dept_map$codigo_dpto
  )
  plot <- dibujar_mapa_departamental(map_data, "muertes", title_text, subtitle_text,
                                     wrap_text(paste0(caption_base, " Unmapped records: ", unmapped,
                                                      ". Cartography: DANE, MGN 2025."), 105), "Deaths", nrow(datos))
  list(grafico = plot, no_mapeados = unmapped)
}

map_residence <- create_map(
  "cod_dpto_r",
  "Dengue mortality according to department of residence",
  "Reported counts; resident population by department and year is required for rates"
)
save_figure(
  map_residence$grafico, "03_geograficas", "04_mapa_residencia_departamento.tiff", 8.5, 9,
  "Map by department of residence", nrow(datos),
  paste("Counts, not rates. Unmapped:", map_residence$no_mapeados)
)

map_occurrence <- create_map(
  "cod_dpto_o",
  "Dengue mortality according to department of occurrence",
  "Occurrence/origin should not be confused with residence or notification"
)
save_figure(
  map_occurrence$grafico, "03_geograficas", "05_mapa_ocurrencia_departamento.tiff", 8.5, 9,
  "Map by department of occurrence", nrow(datos),
  paste("Counts, not rates. Unmapped:", map_occurrence$no_mapeados)
)

# 5. Care times and clinical evolution -----------------------------------------
period_data <- datos %>%
  filter(ano >= 2013) %>%
  mutate(periodo = cut(
    ano,
    breaks = c(2012, 2016, 2020, 2025),
    labels = c("2013–2016", "2017–2020", "2021–2025")
  ))

interval_labels <- c(
  dias_inicio_consulta = "Symptoms → consultation",
  dias_inicio_defuncion = "Symptoms → death",
  dias_inicio_notificacion = "Symptoms → notification"
)

times_data <- period_data %>%
  select(periodo, all_of(names(interval_labels))) %>%
  pivot_longer(-periodo, names_to = "intervalo", values_to = "dias") %>%
  filter(!is.na(dias), dias >= 0) %>%
  group_by(periodo, intervalo) %>%
  summarise(
    n = n(),
    mediana = median(dias),
    p25 = quantile(dias, 0.25),
    p75 = quantile(dias, 0.75),
    .groups = "drop"
  ) %>%
  mutate(intervalo = factor(
    interval_labels[intervalo],
    levels = unname(interval_labels)
  ))

p_times <- ggplot(times_data, aes(periodo, mediana, colour = periodo)) +
  geom_errorbar(aes(ymin = p25, ymax = p75), width = 0.14, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0("Median: ", label_number()(mediana))), vjust = -0.9, size = 3) +
  facet_wrap(~intervalo, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("2013–2016" = blue, "2017–2020" = green, "2021–2025" = orange)) +
  labs(
    title = "Timeliness of care and clinical evolution",
    subtitle = "Period 2013–2025, with greater temporal consistency; median and interquartile range",
    x = "Period",
    y = "Days",
    colour = "Period",
    caption = paste0(
      "Source: SIVIGILA, dengue mortality, 2013–2025. Deduplicated dataset; n = ",
      format(nrow(period_data), big.mark = ",", decimal.mark = "."), ". ",
      " Missing or negative values are excluded only from each calculation; vertical scales vary between panels."
    )
  ) +
  thesis_theme() +
  theme(legend.position = "none")

save_figure(
  p_times, "04_atencion", "06_tiempos_atencion_evolucion.tiff", 13, 5.8,
  "Care times and evolution", nrow(period_data),
  "Period 2013–2025. Medians and IQR; missing and negative values are excluded per interval."
)

# 6. Area, insurance, and ethnic belonging --------------------------------
distribution_plot <- function(variable, title_text, color) {
  table_data <- datos %>%
    mutate(categoria = as.character(.data[[variable]])) %>%
    mutate(categoria = replace_na(categoria, "Sin información")) %>%
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
names(datos)
datos <- datos %>%
  mutate(
    pertenencia_etnica_etiqueta = case_when(
      per_etn == "1" ~ "Indígena",
      per_etn == "2" ~ "ROM/Gitano",
      per_etn == "3" ~ "Raizal",
      per_etn == "4" ~ "Palenquero",
      per_etn == "5" ~ "Negro/Mulato/Afrocolombiano",
      TRUE ~ "Otro/Sin información"
    ),
    area_etiqueta = case_when(
      area == "1" ~ "Cabecera municipal",
      area == "2" ~ "Centro poblado",
      area == "3" ~ "Rural disperso",
      TRUE ~ "Sin información"
    ),
    regimen_salud_etiqueta = case_when(
      tip_ss == "C" ~ "Contributivo",
      tip_ss == "S" ~ "Subsidiado",
      tip_ss == "P" ~ "Excepción",
      tip_ss == "E" ~ "Especial",
      tip_ss == "N" ~ "No asegurado",
      TRUE ~ "Sin información"
    ),
    hospitalizado_etiqueta = case_when(
      pac_hos == "1" ~ "Sí",
      pac_hos == "2" ~ "No",
      TRUE ~ "Sin información"
    )
  )

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
    subtitle = "Percentage distribution of reported deaths",
    caption = caption_base,
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"),
      plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"),
      plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0)
    )
  )

save_figure(
  p_social, "05_sociales", "07_perfil_social.tiff", 15, 7.5,
  "Social profile of reported deaths", nrow(datos),
  "Descriptive percentages; they do not represent relative risks."
)

# 7. Hospitalization and delay -------------------------------------------------
hospitalization <- datos %>%
  mutate(categoria = replace_na(as.character(hospitalizado_etiqueta), "Sin información")) %>%
  count(categoria, name = "n") %>%
  mutate(porcentaje = n / sum(n))

p_hosp_prop <- ggplot(hospitalization, aes(categoria, porcentaje, fill = categoria)) +
  geom_col(width = 0.66, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", porcentaje_es(porcentaje, accuracy = 0.1), ")")), vjust = -0.5, size = 3.4) +
  scale_fill_manual(values = c("Sí" = blue, "No" = gray, "Sin información" = yellow)) +
  scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Hospitalized proportion", x = NULL, y = "Percentage") +
  thesis_theme()

hosp_delay <- datos %>%
  filter(
    hospitalizado_etiqueta == "Sí",
    !is.na(dias_consulta_hospitalizacion),
    dias_consulta_hospitalizacion >= 0
  ) %>%
  mutate(
    demora = case_when(
      dias_consulta_hospitalizacion == 0 ~ "Mismo día",
      dias_consulta_hospitalizacion == 1 ~ "1 día",
      dias_consulta_hospitalizacion <= 3 ~ "2–3 días",
      TRUE ~ "4 días o más"
    ),
    demora = factor(demora, levels = c("Mismo día", "1 día", "2–3 días", "4 días o más"))
  ) %>%
  count(demora, name = "n") %>%
  mutate(porcentaje = n / sum(n))

p_hosp_delay <- ggplot(hosp_delay, aes(demora, porcentaje)) +
  geom_col(fill = green, width = 0.68) +
  geom_text(aes(label = paste0(n, " (", porcentaje_es(porcentaje, accuracy = 0.1), ")")), vjust = -0.5, size = 3.2) +
  scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.17))) +
  labs(title = "Time between consultation and hospitalization", x = NULL, y = "Percentage") +
  thesis_theme() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_hospitalization <- (p_hosp_prop | p_hosp_delay) +
  plot_annotation(
    title = "Hospitalization of deceased cases",
    subtitle = "The delay is calculated only for hospitalized patients with valid and non-negative dates",
    caption = caption_base,
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"),
      plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"),
      plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0)
    )
  )

save_figure(
  p_hospitalization, "04_atencion", "08_hospitalizacion_demora.tiff", 12, 7,
  "Hospitalization and delay", nrow(datos),
  paste("Valid delay for", sum(hosp_delay$n), "hospitalized patients.")
)

# 8. Basic cause of death ---------------------------------------------------
causes <- datos %>%
  mutate(codigo = toupper(trimws(as.character(cbmte)))) %>%
  filter(!is.na(codigo), grepl("^[A-Z][0-9]{2}[A-Z0-9]?$", codigo)) %>%
  count(codigo, sort = TRUE, name = "n") %>%
  slice_head(n = 12) %>%
  mutate(codigo = factor(codigo, levels = rev(codigo)))

p_causes <- ggplot(causes, aes(codigo, n)) +
  geom_col(fill = orange, width = 0.72) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = label_number()) +
  labs(
    title = "Most frequent basic cause of death codes",
    subtitle = "Codes with ICD-10 structure are shown; they must be validated and decoded prior to clinical interpretation",
    x = "Code registered in CBMTE",
    y = "Reported deaths",
    caption = paste0(
      caption_base, " Missing CBMTE: ", sum(is.na(datos$cbmte)),
      "; numerical codes 8888/9999 are not presented as ICD-10 causes."
    )
  ) +
  thesis_theme()

save_figure(
  p_causes, "06_causas", "09_causas_basicas_frecuentes.tiff", 10, 7.5,
  "Most frequent basic causes", nrow(datos),
  "Exact codes; require validation and ICD-10 grouping."
)

# 9. Vulnerability groups -------------------------------------------------
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
    n_si = sum(x == "1", na.rm = TRUE),
    n_validos = sum(valid_responses),
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
    aes(label = paste0(n_si, " (", porcentaje_es(porcentaje, accuracy = 0.1), ")")),
    hjust = -0.08,
    size = 3.3
  ) +
  coord_flip() +
  scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Registered vulnerability groups",
    subtitle = "Percentage among records with valid response for each indicator",
    x = NULL,
    y = "Percentage with Yes response",
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

# 10. Quality per year ---------------------------------------------------------
quality_flags <- c(
  flag_condicion_final_no_muerto = "Final condition is not dead",
  flag_fecha_defuncion_faltante = "Missing date of death",
  flag_cronologia_inconsistente = "Negative clinical chronology",
  flag_edad_fuera_rango = "Age out of range"
)

quality_data <- datos %>%
  select(ano, all_of(names(quality_flags))) %>%
  pivot_longer(-ano, names_to = "bandera", values_to = "marcado") %>%
  group_by(ano, bandera) %>%
  summarise(porcentaje = mean(marcado, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    indicador = factor(
      quality_flags[bandera],
      levels = rev(unname(quality_flags))
    )
  )

p_quality <- ggplot(quality_data, aes(ano, indicador, fill = porcentaje)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "B", labels = porcentaje_es, breaks = pretty_breaks(5)) +
  scale_x_continuous(breaks = SIVIGILA_ANIOS) +
  labs(
    title = "Data quality issues by year",
    subtitle = "Percentage of records with each flag; these indicators are not clinical outcomes",
    x = "Epidemiological year",
    y = NULL,
    fill = "Records",
    caption = caption_base
  ) +
  thesis_theme() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = grid::unit(1.3, "cm")
  )

save_figure(
  p_quality, "07_calidad", "11_calidad_datos_por_anio.tiff", 12, 6.8,
  "Data quality by year", nrow(datos),
  "Quality control; do not interpret as a clinical pattern."
)

figures_index <- do.call(rbind, manifest)
write.csv(
  figures_index,
  file.path(figures_dir, "INDICE_FIGURAS.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

package_versions <- data.frame(
  paquete = packages,
  version = vapply(packages, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
write.csv(
  package_versions,
  file.path(figures_dir, "VERSIONES_PAQUETES.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

tiff_files <- list.files(figures_dir, pattern = "\\.tiff$", recursive = TRUE, full.names = TRUE)
stopifnot(
  nrow(figures_index) == 11L,
  length(tiff_files) == 11L,
  all(file.info(tiff_files)$size > 10000)
)

cat("Generated figures:", length(tiff_files), "\n")
cat("Directory:", normalizePath(figures_dir, winslash = "/"), "\n")
