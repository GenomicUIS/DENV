source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## DENGUE GRAPHICS GENERATION (EVENT 210)
## Purpose: Create visualizations by reading deduplicated annual bases 
## iteratively to optimize RAM usage.
## -----------------------------------------------------------------------------

packages <- c("ggplot2", "dplyr", "tidyr", "scales", "patchwork", "sf", "lwgeom")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), ". Install them with install.packages().")
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
source("departmental_maps.R", encoding = "UTF-8")

results_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue")
data_dir <- file.path(results_dir, "datos_procesados", "analitica_sin_duplicados_por_anio")
figures_dir <- file.path(results_dir, "figuras")
tables_dir <- file.path(results_dir, "reportes", "tablas_figuras")
years <- SIVIGILA_ANIOS

paths <- file.path(data_dir, sprintf("dengue_analitica_sin_duplicados_%d.rds", years))
if (any(!file.exists(paths))) {
  stop("Missing deduplicated bases. Execute the deduplication script first.")
}

subfolders <- c(
  "01_temporales", "02_demograficas", "03_geograficas", "04_clasificacion",
  "05_atencion", "06_sociales", "07_integracion", "08_calidad",
  "_previsualizaciones"
)
invisible(lapply(c(file.path(figures_dir, subfolders), tables_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

## Palette and theme
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

wrap_text <- function(x, width = 30) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

manifest <- list()
save_figure <- function(plot, subfolder, file_name, width, height, title_text, sample_size, note) {
  path <- file.path(figures_dir, subfolder, file_name)
  ggsave(path, plot, device = "tiff", width = width, height = height, units = "in", dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  png_path <- file.path(figures_dir, "_previsualizaciones", sub("\\.tiff$", ".png", file_name, ignore.case = TRUE))
  ggsave(png_path, plot, device = "png", width = width, height = height, units = "in", dpi = 150, bg = "white", limitsize = FALSE)
  
  manifest[[length(manifest) + 1L]] <<- data.frame(
    orden = length(manifest) + 1L, categoria = subfolder,
    archivo_tiff = normalizePath(path, winslash = "/", mustWork = TRUE),
    titulo = title_text, muestra = sample_size, nota = note,
    ancho_pulgadas = width, alto_pulgadas = height,
    dpi = 600L, compresion = "LZW", stringsAsFactors = FALSE
  )
  guardar_evidencia_figura(plot, path, width, height)
  invisible(path)
}

write_table <- function(table_data, name) {
  write.csv(table_data, file.path(tables_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

social_variables <- c(area_etiqueta = "Área de procedencia", regimen_salud_etiqueta = "Régimen de salud", pertenencia_etnica_etiqueta = "Pertenencia étnica")
vulnerability_variables <- c(gp_discapa = "Discapacidad", gp_desplaz = "Población desplazada", gp_migrant = "Población migrante", gp_carcela = "Privación de la libertad", gp_gestan = "Gestante", gp_indigen = "Habitante de calle", gp_pobicfb = "Población infantil ICBF", gp_mad_com = "Madre comunitaria", gp_desmovi = "Población desmovilizada", gp_psiquia = "Centro psiquiátrico", gp_vic_vio = "Víctima de violencia armada")
quality_flags <- c(flag_edad_fecha_nacimiento_discordante = "Discordant age/birth date", flag_anio_fec_not_distinto = "FEC_NOT year different from ANO", flag_ajuste_anterior_notificacion = "Adjustment prior to notification", flag_cronologia_clinica_inconsistente = "Negative clinical chronology", flag_hospitalizacion_inconsistente = "Inconsistent hospitalization/date", flag_muerto_sin_fecha_defuncion = "Deceased without date of death")
care_intervals <- c(dias_inicio_consulta = "Symptoms → consultation", dias_inicio_hospitalizacion = "Symptoms → hospitalization", dias_inicio_notificacion = "Symptoms → notification")

## 1. Iterative data aggregation per year
aggregates <- lapply(seq_along(paths), function(i) {
  datos <- readRDS(paths[[i]])
  year <- years[[i]]
  
  times <- do.call(rbind, lapply(names(care_intervals), function(variable) {
    x <- datos[[variable]]
    x <- x[!is.na(x) & x >= 0]
    data.frame(
      ano = year, intervalo = unname(care_intervals[[variable]]), n = length(x),
      mediana = if (length(x)) median(x) else NA_real_,
      p25 = if (length(x)) unname(quantile(x, 0.25)) else NA_real_,
      p75 = if (length(x)) unname(quantile(x, 0.75)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  
  social <- do.call(rbind, lapply(names(social_variables), function(variable) {
    category <- as.character(datos[[variable]])
    category[is.na(category) | category == ""] <- "Sin información"
    if (variable == "pertenencia_etnica_etiqueta") {
      category[category %in% c("ROM/Gitano", "Raizal", "Palenquero")] <- "Otros grupos étnicos minoritarios"
    }
    table_data <- as.data.frame(table(category), stringsAsFactors = FALSE)
    names(table_data) <- c("categoria", "n")
    table_data$ano <- year
    table_data$variable <- unname(social_variables[[variable]])
    table_data[, c("ano", "variable", "categoria", "n")]
  }))
  
  vulnerability <- do.call(rbind, lapply(names(vulnerability_variables), function(variable) {
    x <- as.character(datos[[variable]])
    data.frame(ano = year, grupo = unname(vulnerability_variables[[variable]]), n_si = sum(x == "1", na.rm = TRUE), n_validos = sum(x %in% c("1", "2")), stringsAsFactors = FALSE)
  }))
  
  quality <- do.call(rbind, lapply(names(quality_flags), function(variable) {
    x <- as.logical(datos[[variable]])
    data.frame(ano = year, indicador = unname(quality_flags[[variable]]), n_marcados = sum(x, na.rm = TRUE), n_evaluables = sum(!is.na(x)), stringsAsFactors = FALSE)
  }))
  
  list(
    anual = data.frame(ano = year, casos = nrow(datos), hospitalizados = sum(datos$pac_hos == "1", na.rm = TRUE), stringsAsFactors = FALSE),
    semanal = count(datos, ano, semana, name = "casos"),
    edad_sexo = count(datos, grupo_edad_tesis, sexo_etiqueta, name = "casos"),
    residencia = count(datos, cod_dpto_r, name = "casos"),
    ocurrencia = count(datos, cod_dpto_o, name = "casos"),
    estado_final = count(datos, ano, nom_est_f_caso, name = "casos"),
    tiempos = times, sociales = social, vulnerabilidad = vulnerability, calidad = quality
  )
})

merge_data <- function(name) bind_rows(lapply(aggregates, `[[`, name))

annual_series <- merge_data("anual") %>% mutate(p_hospitalizados = hospitalizados / casos)
heatmap <- merge_data("semanal") %>% group_by(ano, semana) %>% summarise(casos = sum(casos), .groups = "drop") %>% complete(ano = SIVIGILA_ANIOS, semana = 1:53, fill = list(casos = 0L))
age_sex <- merge_data("edad_sexo") %>% group_by(grupo_edad_tesis, sexo_etiqueta) %>% summarise(casos = sum(casos), .groups = "drop")
residence <- merge_data("residencia") %>% group_by(cod_dpto_r) %>% summarise(casos = sum(casos), .groups = "drop")
occurrence <- merge_data("ocurrencia") %>% group_by(cod_dpto_o) %>% summarise(casos = sum(casos), .groups = "drop")
final_status <- merge_data("estado_final") %>% group_by(ano, nom_est_f_caso) %>% summarise(casos = sum(casos), .groups = "drop") %>% group_by(ano) %>% mutate(porcentaje = casos / sum(casos)) %>% ungroup()
times <- merge_data("tiempos")
social <- merge_data("sociales") %>% mutate(periodo = cut(ano, breaks = c(2006, 2012, 2018, 2025), labels = c("2007-2012", "2013-2018", "2019-2025"))) %>% group_by(periodo, variable, categoria) %>% summarise(n = sum(n), .groups = "drop") %>% group_by(periodo, variable) %>% mutate(porcentaje = n / sum(n)) %>% ungroup()
vulnerability <- merge_data("vulnerabilidad") %>% group_by(grupo) %>% summarise(n_si = sum(n_si), n_validos = sum(n_validos), .groups = "drop") %>% mutate(porcentaje = n_si / n_validos)
quality <- merge_data("calidad") %>% mutate(porcentaje = n_marcados / n_evaluables)

total_n <- sum(annual_series$casos)
caption_base <- paste0("Source: SIVIGILA, dengue (event 210), 2007-2025. Deduplicated base; n = ", format(total_n, big.mark = ".", decimal.mark = ","), ".")

write_table(annual_series, "01_serie_anual.csv")
write_table(heatmap, "02_semana_anio.csv")
write_table(age_sex, "03_edad_sexo.csv")
write_table(residence, "04_residencia_departamento.csv")
write_table(occurrence, "05_ocurrencia_departamento.csv")
write_table(final_status, "06_estado_final.csv")
write_table(times, "07_oportunidad_atencion.csv")
write_table(social, "08_perfil_social.csv")
write_table(vulnerability, "09_vulnerabilidad.csv")
write_table(quality, "12_calidad_datos.csv")

## 2. Figures Construction
# 2.1 Annual trend
p_series <- ggplot(annual_series, aes(ano, casos)) +
  geom_line(colour = blue, linewidth = 1) +
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.8, stroke = 0.9) +
  geom_label(label.padding = grid::unit(0.07, "lines"), linewidth = 0, fill = "white", aes(label = label_number(scale_cut = escala_miles_es(), decimal.mark = ",")(casos)), vjust = -0.75, size = 3, colour = "#334155") +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ","), expand = expansion(mult = c(0.02, 0.14))) +
  labs(title = "Reported dengue per year", subtitle = "Number of event 210 records; corresponds to counts, not population rates", x = "Epidemiological year", y = "Reported cases", caption = caption_base) +
  thesis_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p_series, "01_temporales", "01_serie_anual_dengue.tiff", 11, 6.5, "Reported dengue per year", total_n, "Counts; denominators are required for incidence.")

# 2.2 Seasonality
heatmap$ano_factor <- factor(heatmap$ano, levels = rev(years))
p_heatmap <- ggplot(heatmap, aes(semana, ano_factor, fill = casos)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", breaks = pretty_breaks(5)) +
  scale_x_continuous(breaks = c(1, seq(4, 52, 4), 53), expand = c(0, 0)) +
  labs(title = "Seasonality of dengue", subtitle = "Reported cases by epidemiological week and year", x = "Epidemiological week", y = "Year", fill = "Cases", caption = caption_base) +
  thesis_theme() + theme(panel.grid = element_blank(), legend.position = "right") + guides(fill = guide_colourbar(barheight = grid::unit(4.2, "cm")))
save_figure(p_heatmap, "01_temporales", "02_mapa_calor_semana_anio.tiff", 11, 7.5, "Seasonality per epidemiological week", total_n, "Color scale with square root to preserve contrast.")

# 2.3 Age and sex
age_sex <- age_sex %>% filter(!is.na(grupo_edad_tesis), sexo_etiqueta %in% c("Femenino", "Masculino")) %>% mutate(valor = if_else(sexo_etiqueta == "Femenino", -casos, casos), grupo_edad_tesis = factor(grupo_edad_tesis, levels = c("<1", "1-4", "5-14", "15-24", "25-44", "45-64", "65+")))
p_age_sex <- ggplot(age_sex, aes(grupo_edad_tesis, valor, fill = sexo_etiqueta)) +
  geom_col(width = 0.78) +
  geom_text(aes(label = label_number(scale_cut = escala_miles_es(), decimal.mark = ",")(casos), hjust = if_else(valor < 0, 1.08, -0.08)), size = 3, colour = "#334155") +
  coord_flip() +
  scale_fill_manual(values = c("Femenino" = purple, "Masculino" = blue)) +
  scale_y_continuous(labels = function(x) label_number(big.mark = ".", decimal.mark = ",")(abs(x)), expand = expansion(mult = c(0.18, 0.18))) +
  labs(title = "Dengue cases by age group and sex", subtitle = "Age was converted to years using both EDAD and UNI_MED collectively", x = "Age group (years)", y = "Reported cases", fill = "Sex", caption = caption_base) + thesis_theme()
save_figure(p_age_sex, "02_demograficas", "03_edad_sexo.tiff", 10, 7, "Dengue by age and sex", sum(age_sex$casos), "Includes female or male sex and interpretable age.")

# 2.4 and 2.5 Geographical distribution
dept_map <- cargar_cartografia_departamental()
create_map <- function(counts, code_var, title_text, subtitle_text) {
  counts <- counts %>% rename(codigo_dpto = all_of(code_var))
  valid_codes <- dept_map$codigo_dpto
  unmapped <- sum(counts$casos[is.na(counts$codigo_dpto) | counts$codigo_dpto == "00" | !counts$codigo_dpto %in% valid_codes])
  counts <- counts %>% filter(!is.na(codigo_dpto), codigo_dpto != "00")
  map_data <- dept_map %>% left_join(counts, by = "codigo_dpto") %>% mutate(casos = replace_na(casos, 0L))
  plot <- dibujar_mapa_departamental(map_data, "casos", title_text, subtitle_text, wrap_text(paste0(caption_base, " Records with unmappable code: ", format(unmapped, big.mark=".", decimal.mark=","), ". Cartography: DANE, MGN 2025."), 105), "Cases", total_n)
  list(grafico = plot, no_mapeados = unmapped)
}

map_residence <- create_map(residence, "cod_dpto_r", "Dengue according to department of residence", "Reported counts; resident population by department and year is required for rates")
save_figure(map_residence$grafico, "03_geograficas", "04_mapa_residencia_departamento.tiff", 8.5, 9, "Map by department of residence", total_n, paste("Counts, not rates. Unmapped:", map_residence$no_mapeados))

map_occurrence <- create_map(occurrence, "cod_dpto_o", "Dengue by department of occurrence", "Occurrence should not be confused with residence or notification location")
save_figure(map_occurrence$grafico, "03_geograficas", "05_mapa_ocurrencia_departamento.tiff", 8.5, 9, "Map by department of occurrence", total_n, paste("Counts, not rates. Unmapped:", map_occurrence$no_mapeados))

# 2.6 Final classification
status_levels <- c("Probable", "Confirmado por Nexo Epidemiológico", "Confirmado por laboratorio")
final_status$nom_est_f_caso <- factor(final_status$nom_est_f_caso, levels = status_levels)
p_status <- ggplot(final_status, aes(ano, porcentaje, fill = nom_est_f_caso)) +
  geom_col(width = 0.84) +
  scale_fill_manual(values = c("Probable" = gray, "Confirmado por Nexo Epidemiológico" = yellow, "Confirmado por laboratorio" = blue), drop = FALSE) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(labels = porcentaje_es) +
  labs(title = "Final case classification by year", subtitle = "Percentage composition; final status is used, not initial classification", x = "Epidemiological year", y = "Percentage of cases", fill = "Final status", caption = caption_base) + thesis_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p_status, "04_clasificacion", "06_estado_final_por_anio.tiff", 12, 7, "Final classification by year", total_n, "100% bars; describes composition, not risk.")

# 2.7 Hospitalization and timeliness
p_hospital <- ggplot(annual_series, aes(ano, p_hospitalizados)) +
  geom_line(colour = blue, linewidth = 0.95) +
  geom_point(colour = blue, fill = "white", shape = 21, size = 2.5, stroke = 0.8) +
  scale_x_continuous(breaks = seq(2007, 2025, 2)) +
  scale_y_continuous(labels = porcentaje_es, limits = c(0, max(annual_series$p_hospitalizados) * 1.12)) +
  labs(title = "A. Hospitalized cases", x = "Year", y = "Percentage") + thesis_theme()

times$intervalo <- factor(times$intervalo, levels = unname(care_intervals))
p_times <- ggplot(times, aes(ano, mediana)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = light_blue, alpha = 0.24) +
  geom_line(colour = blue, linewidth = 0.9) +
  geom_point(colour = blue, size = 2.1) +
  geom_label(data = data.frame(ano = 2013, mediana = 138, intervalo = factor("Symptoms → notification", levels = unname(care_intervals))), aes(ano, mediana, label = "Quality anomaly\nin 2013"), inherit.aes = FALSE, size = 2.7, colour = "#7C2D12", fill = "#FFF7ED", linewidth = 0.2) +
  facet_wrap(~intervalo, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(2007, 2025, 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(title = "B. Median and interquartile range of timeliness", x = "Year", y = "Days") + thesis_theme(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_care <- (p_hospital / p_times) + plot_layout(heights = c(0.9, 1.25)) + plot_annotation(title = "Hospitalization and timeliness of care", subtitle = "Annual evolution of care registered in SIVIGILA", caption = wrap_text(paste0(caption_base, " For medians, only missing or negative values of the corresponding interval are excluded; the 2013 value is kept and flagged as a quality anomaly."), 140), theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"), plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"), plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0)))
save_figure(p_care, "05_atencion", "07_hospitalizacion_oportunidad.tiff", 15, 10, "Hospitalization and timeliness", total_n, paste("Hospitalized percentage, medians, and annual IQR; negative intervals excluded by indicator. The 2013 reporting anomaly is kept and flagged."))

# 2.8 Social profile by period
period_palette <- c("2007-2012" = blue, "2013-2018" = green, "2019-2025" = orange)
social_plot <- function(variable_name) {
  table_data <- social %>% filter(variable == variable_name) %>% mutate(categoria = factor(wrap_text(categoria, 22), levels = rev(unique(wrap_text(categoria, 22)))))
  ggplot(table_data, aes(categoria, porcentaje, fill = periodo)) + geom_col(position = position_dodge(width = 0.78), width = 0.7) + coord_flip() + scale_fill_manual(values = period_palette) + scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.12))) + labs(title = variable_name, x = NULL, y = "Percentage", fill = "Period") + thesis_theme(9.5) + theme(plot.title = element_text(size = 11, face = "bold"))
}
p_social <- (social_plot("Área de procedencia") | social_plot("Régimen de salud") | social_plot("Pertenencia étnica")) + plot_layout(guides = "collect") + plot_annotation(title = "Social profile of reported cases", subtitle = "Descriptive comparison between three periods", caption = wrap_text(paste0(caption_base, " Percentages describe notification composition and not relative risks."), 140), theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 15, colour = "#17202A"), plot.subtitle = element_text(family = "sans", size = 11, colour = "#475569"), plot.caption = element_text(family = "sans", size = 8.5, colour = "#64748B", hjust = 0))) & theme(legend.position = "bottom")
save_figure(p_social, "06_sociales", "08_perfil_social_por_periodo.tiff", 16, 8, "Social profile by period", total_n, "Descriptive percentages by period; includes the category without information.")

# 2.9 Vulnerability groups
vulnerability_plot_data <- vulnerability %>% filter(n_si >= 5) %>% arrange(porcentaje) %>% mutate(grupo = factor(grupo, levels = grupo))
p_vulnerability <- ggplot(vulnerability_plot_data, aes(grupo, porcentaje)) + geom_col(fill = purple, width = 0.7) + geom_text(aes(label = paste0(label_number(big.mark = ".", decimal.mark = ",")(n_si), " (", porcentaje_es(porcentaje, accuracy = 0.01), ")")), hjust = -0.08, size = 3.1) + coord_flip() + scale_y_continuous(labels = porcentaje_es, expand = expansion(mult = c(0, 0.25))) + labs(title = "Registered vulnerability groups", subtitle = "Percentage among records with valid response for each indicator", x = NULL, y = "Percentage with Yes response", caption = wrap_text(paste0(caption_base, " For confidentiality, indicators with fewer than five affirmative responses are not shown."), 120)) + thesis_theme()
save_figure(p_vulnerability, "06_sociales", "09_grupos_vulnerabilidad.tiff", 10.5, 7.5, "Vulnerability groups", total_n, "Categories with fewer than five affirmative responses are omitted.")

# 2.10 and 2.11 Events integration
severe_path <- file.path(SIVIGILA_PATH, "resultados", "dengue_grave", "datos_procesados", "dengue_grave_analitica_sin_duplicados.rds")
mortality_path <- file.path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "datos_procesados", "mortalidad_dengue_muestra_graficas.rds")

if (!file.exists(severe_path) || !file.exists(mortality_path)) {
  stop("Missing deduplicated bases for severe dengue or dengue mortality.")
}
severe <- readRDS(severe_path) %>% count(ano, name = "dengue_grave")
mortality <- readRDS(mortality_path) %>% count(ano, name = "muertes")
integration <- annual_series %>% select(ano, dengue = casos) %>% left_join(severe, by = "ano") %>% left_join(mortality, by = "ano") %>% mutate(dengue_grave = replace_na(dengue_grave, 0L), muertes = replace_na(muertes, 0L), proporcion_grave = dengue_grave / (dengue + dengue_grave))
write_table(integration, "10_integracion_eventos.csv")

p_severe_proportion <- ggplot(integration, aes(ano, proporcion_grave)) + geom_line(colour = orange, linewidth = 1) + geom_point(colour = orange, fill = "white", shape = 21, size = 2.8, stroke = 0.9) + scale_x_continuous(breaks = years) + scale_y_continuous(labels = porcentaje_formato_es(accuracy = 0.1), expand = expansion(mult = c(0.03, 0.1))) + labs(title = "Proportion of notifications classified as severe dengue", subtitle = "Event 220 divided by the sum of events 210 and 220 in each year", x = "Epidemiological year", y = "Severe dengue among events 210 + 220", caption = wrap_text(paste0("Source: SIVIGILA, events 210 and 220, 2007-2025; deduplicated bases. It is a proportion among notifications, not a population rate."), 115)) + thesis_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(p_severe_proportion, "07_integracion", "10_proporcion_dengue_grave.tiff", 11, 6.5, "Proportion of severe dengue", sum(integration$dengue + integration$dengue_grave), "Annual denominator: event 210 + event 220; it is not incidence.")

three_events <- integration %>% select(ano, `Dengue (event 210)` = dengue, `Severe dengue (event 220)` = dengue_grave, `Dengue mortality (event 580)` = muertes) %>% pivot_longer(-ano, names_to = "evento", values_to = "registros")
p_three_events <- ggplot(three_events, aes(ano, registros, colour = evento)) + geom_line(linewidth = 0.9) + geom_point(size = 2.2) + facet_wrap(~evento, ncol = 1, scales = "free_y") + scale_colour_manual(values = c(blue, purple, orange)) + scale_x_continuous(breaks = seq(2007, 2025, 2)) + scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",")) + labs(title = "Comparative evolution of dengue, severe dengue, and mortality", subtitle = "Each panel has its own vertical scale to make the three series visible", x = "Epidemiological year", y = "Notified records", colour = "Event", caption = wrap_text(paste0("Source: SIVIGILA, events 210, 220, and 580, 2007-2025; deduplicated bases. Counts are not rates and vertical scales differ between panels."), 120)) + thesis_theme() + theme(legend.position = "none")
save_figure(p_three_events, "07_integracion", "11_tendencia_tres_eventos.tiff", 11, 9.5, "Comparative trend of three events", sum(three_events$registros), "Panels with free scale; compare temporal patterns, not heights between panels.")

# 2.12 Quality
quality$indicador <- factor(quality$indicador, levels = rev(unname(quality_flags)))
p_quality <- ggplot(quality, aes(ano, indicador, fill = porcentaje)) + geom_tile(colour = "white", linewidth = 0.25) + scale_fill_viridis_c(option = "B", labels = porcentaje_es, breaks = pretty_breaks(5)) + scale_x_continuous(breaks = years) + labs(title = "Data quality issues by year", subtitle = "Percentage of records with each flag after deduplication", x = "Epidemiological year", y = NULL, fill = "Records", caption = caption_base) + thesis_theme() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") + guides(fill = guide_colourbar(barheight = grid::unit(4.2, "cm")))
save_figure(p_quality, "08_calidad", "12_calidad_datos_por_anio.tiff", 12.5, 7, "Data quality by year", total_n, "Flags are consistency controls and not clinical outcomes.")

## 3. Closure and metadata
figures_index <- do.call(rbind, manifest)
write.csv(figures_index, file.path(figures_dir, "INDICE_FIGURAS.csv"), row.names = FALSE, fileEncoding = "UTF-8")

package_versions <- data.frame(paquete = packages, version = vapply(packages, function(x) as.character(packageVersion(x)), character(1)), stringsAsFactors = FALSE)
write.csv(package_versions, file.path(figures_dir, "VERSIONES_PAQUETES.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("DENGUE figures successfully generated. Directory:", normalizePath(figures_dir, winslash = "/"), "\n")
