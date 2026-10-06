source("00_config.R", encoding = "UTF-8")

## -----------------------------------------------------------------------------
## DENGUE DOCUMENTATION AND DICTIONARY (EVENT 210)
## Consolidates annual profiles, generates the variable dictionary and creates
## quality and duplicates audit reports.
## -----------------------------------------------------------------------------

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## 1. Working paths
results_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue")
data_dir      <- file.path(results_dir, "datos_procesados")
reports_dir   <- file.path(results_dir, "reportes")
control_dir    <- file.path(data_dir, "control_interno")

## -----------------------------------------------------------------------------
## 2. Consolidation of clinical intervals
## 09_prepare_dengue.R saves annual fragments of the intervals in
## control_interno/intervalos_YYYY.rds. Here they are joined to view the global
## distribution of care times.
## -----------------------------------------------------------------------------
intervals_path     <- file.path(reports_dir, "resumen_intervalos.csv")
intervals_files <- file.path(control_dir,
                             paste0("intervalos_", SIVIGILA_ANIOS, ".rds"))

if (all(file.exists(intervals_files))) {
  interval_vars <- names(readRDS(intervals_files[[1L]]))
  
  intervals_summary <- do.call(rbind, lapply(interval_vars, function(variable) {
    x <- unlist(lapply(intervals_files,
                       function(ruta) readRDS(ruta)[[variable]]),
                use.names = FALSE)
    valid_vals <- x[!is.na(x)]
    
    q <- if (length(valid_vals) > 0L) {
      stats::quantile(valid_vals, c(0, .25, .5, .75, .95, 1), names = FALSE)
    } else {
      rep(NA_real_, 6L)
    }
    
    data.frame(
      intervalo     = variable,
      n_disponibles = length(valid_vals),
      n_faltantes   = sum(is.na(x)),
      n_negativos   = sum(valid_vals < 0),
      minimo        = q[1L],
      p25           = q[2L],
      mediana       = q[3L],
      p75           = q[4L],
      p95           = q[5L],
      maximo        = q[6L],
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(intervals_summary, intervals_path,
            row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

## -----------------------------------------------------------------------------
## 3. Verify inputs
## Assumes 09_prepare_dengue.R has already been executed.
## -----------------------------------------------------------------------------
required_files <- c(
  "perfil_variables_por_anio.csv",
  "calidad_por_anio.csv",
  "resumen_duplicados.csv",
  "disponibilidad_variables_por_anio.csv",
  "resumen_intervalos.csv"
)

missing_files <- required_files[!file.exists(file.path(reports_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("Missing results from 09_prepare_dengue.R: ",
       paste(missing_files, collapse = ", "), call. = FALSE)
}

read_report <- function(name, check.names = FALSE) {
  read.csv(file.path(reports_dir, name),
           stringsAsFactors = FALSE,
           check.names = check.names,
           encoding = "UTF-8")
}

annual_profile   <- read_report("perfil_variables_por_anio.csv", check.names = TRUE)
annual_quality  <- read_report("calidad_por_anio.csv", check.names = TRUE)
availability <- read_report("disponibilidad_variables_por_anio.csv")
duplicates     <- read_report("resumen_duplicados.csv")

variables        <- unique(annual_profile$variable)
source_vars <- unique(c("fuente_evento", "fuente_archivo", "fuente_anio",
                        "fuente_fila_excel", availability$variable))

## -----------------------------------------------------------------------------
## 4. Global profile per variable
## -----------------------------------------------------------------------------
global_profile <- do.call(rbind, lapply(split(annual_profile, annual_profile$variable),
                                        function(x) {
                                          data.frame(
                                            variable               = x$variable[[1L]],
                                            tipo_r                 = x$tipo_r[[1L]],
                                            n_registros            = sum(x$n_registros),
                                            n_faltantes            = sum(x$n_faltantes),
                                            pct_faltantes          = 100 * sum(x$n_faltantes) / sum(x$n_registros),
                                            n_anios_con_algun_dato = sum(x$n_faltantes < x$n_registros),
                                            n_unicos_anuales_max   = max(x$n_unicos_anuales, na.rm = TRUE),
                                            origen                 = if (x$variable[[1L]] %in% source_vars) {
                                              "Harmonized original"
                                            } else {
                                              "Derived"
                                            },
                                            stringsAsFactors = FALSE
                                          )
                                        }))
rownames(global_profile) <- NULL

## -----------------------------------------------------------------------------
## 5. Thematic groups and analytical priorities
## -----------------------------------------------------------------------------
def_groups <- list(
  "Traceability and system" = c("fuente_evento", "fuente_archivo", "fuente_anio",
                                "fuente_fila_excel", "separacion", "particion",
                                "consecutive", "consecutive_origen", "cod_eve",
                                "cod_eve__1", "nombre_evento", "confirmados", "va_sispro"),
  "Time" = c("fec_not", "semana", "ano", "fec_con", "ini_sin", "fec_hos",
             "fec_def", "fecha_nto", "fec_aju", "fec_arc_xl",
             grep("^dias_", variables, value = TRUE)),
  "Demographics" = c("edad", "uni_med", "edad_anios", "grupo_edad_tesis",
                     "sexo", "sexo_etiqueta", "nacionalidad",
                     "nombre_nacionalidad", "ocupacion"),
  "Geography" = c(grep("^cod_pais_|^cod_dpto_|^cod_mun_|^divipola_", variables, value = TRUE),
                  "pais_ocurrencia", "departamento_ocurrencia", "municipio_ocurrencia",
                  "pais_residencia", "departamento_residencia", "municipio_residencia",
                  "departamento_notificacion", "municipio_notificacion",
                  "area", "area_etiqueta"),
  "Health services" = c("cod_pre", "cod_sub", "llave_upgd", "nom_upgd",
                        "tip_ss", "regimen_salud_etiqueta", "cod_ase",
                        "pac_hos", "hospitalizado_etiqueta"),
  "Ethnicity and vulnerability" = c("per_etn", "pertenencia_etnica_etiqueta",
                                    "gru_pob", "nom_grupo", "estrato",
                                    grep("^gp_", variables, value = TRUE),
                                    "sem_ges", "fuente"),
  "Classification and outcome" = c("tip_cas", "tipo_caso_inicial_etiqueta",
                                   "con_fin", "condicion_final_etiqueta",
                                   "ajuste", "ajuste_etiqueta",
                                   "estado_final_de_caso", "nom_est_f_caso",
                                   "cer_def", "cbmte"),
  "Military forces" = c("fm_fuerza", "fm_unidad", "fm_grado"),
  "Quality control" = grep("^flag_|^grupo_", variables, value = TRUE)
)

group <- setNames(rep("Other", length(variables)), variables)
for (name in names(def_groups)) {
  matches <- intersect(def_groups[[name]], variables)
  if (length(matches) > 0L) group[matches] <- name
}

essentials <- c(
  "ano", "semana", "fec_not", "ini_sin", "fec_con", "fec_hos",
  "edad", "uni_med", "edad_anios", "grupo_edad_tesis", "sexo", "sexo_etiqueta",
  "area", "area_etiqueta", "tip_ss", "regimen_salud_etiqueta", "per_etn",
  "pertenencia_etnica_etiqueta", "pac_hos", "hospitalizado_etiqueta",
  "ajuste", "ajuste_etiqueta", "estado_final_de_caso", "nom_est_f_caso",
  "divipola_ocurrencia", "departamento_ocurrencia", "municipio_ocurrencia",
  "divipola_residencia", "departamento_residencia", "municipio_residencia",
  "dias_inicio_consulta", "dias_inicio_hospitalizacion", "dias_inicio_notificacion"
)

complementaries <- c(
  "ocupacion", "cod_ase", "nacionalidad", "nombre_nacionalidad", "gru_pob",
  "nom_grupo", "estrato", grep("^gp_", variables, value = TRUE), "sem_ges",
  "cod_pre", "cod_sub", "llave_upgd", "nom_upgd", "tip_cas",
  "tipo_caso_inicial_etiqueta", "con_fin", "condicion_final_etiqueta",
  "fec_def", "cer_def", "cbmte", "fecha_nto",
  "cod_pais_o", "cod_pais_r", "pais_ocurrencia", "pais_residencia",
  "cod_dpto_o", "cod_mun_o", "cod_dpto_r", "cod_mun_r",
  "cod_dpto_n", "cod_mun_n", "divipola_notificacion",
  "departamento_notificacion", "municipio_notificacion",
  "dias_consulta_hospitalizacion", "dias_notificacion_ajuste", "dias_inicio_defuncion"
)

control_vars <- c(
  "fuente_evento", "fuente_archivo", "fuente_anio", "fuente_fila_excel",
  "separacion", "particion", "consecutive", "consecutive_origen",
  "cod_eve", "cod_eve__1", "nombre_evento", "fec_arc_xl", "fec_aju",
  "confirmados", "va_sispro", grep("^flag_|^grupo_", variables, value = TRUE)
)

not_recommended <- c("separacion", "particion", "cod_eve__1",
                     "fm_fuerza", "fm_unidad", "fm_grado",
                     "gp_otros", "fuente")

priority <- setNames(rep("Complementary", length(variables)), variables)
priority[intersect(essentials, variables)]      <- "Essential"
priority[intersect(complementaries, variables)] <- "Complementary"
priority[intersect(control_vars, variables)]         <- "Control/traceability"
priority[intersect(not_recommended, variables)] <- "Not recommended for main analysis"

## -----------------------------------------------------------------------------
## 6. Graphical suggestions, cautions, and descriptions
## -----------------------------------------------------------------------------
graphical_use <- setNames(
  rep("Only if it answers a specific question and there is comparable data.",
      length(variables)),
  variables
)
graphical_use[intersect(c("ano", "semana", "fec_not"), variables)] <- "Time series or year x week heatmap."
graphical_use[intersect(c("edad_anios", "grupo_edad_tesis", "sexo_etiqueta"), variables)] <- "Histogram, age/sex bars, or pyramid."
graphical_use[intersect(c("departamento_ocurrencia", "municipio_ocurrencia", "divipola_ocurrencia"), variables)] <- "Occurrence map; use rates with compatible denominator."
graphical_use[intersect(c("departamento_residencia", "municipio_residencia", "divipola_residencia"), variables)] <- "Residence map; align with resident population."
graphical_use[intersect(c("departamento_notificacion", "municipio_notificacion", "divipola_notificacion"), variables)] <- "Surveillance flow; do not interpret as transmission location."
graphical_use[intersect(c("area_etiqueta", "regimen_salud_etiqueta", "pertenencia_etnica_etiqueta"), variables)] <- "Composition bars and trends by period."
graphical_use[intersect(c("hospitalizado_etiqueta", grep("^dias_", variables, value = TRUE)), variables)] <- "Proportions, boxplots, or care opportunity intervals."
graphical_use[intersect(c("nom_est_f_caso", "estado_final_de_caso", "ajuste_etiqueta"), variables)] <- "Stacked bars and evolution of final classification."
graphical_use[intersect(c("ocupacion", "cod_ase", "nom_upgd", "cbmte"), variables)] <- "Group or decode first; then ordered bars or Pareto."
graphical_use[intersect(control_vars, variables)]         <- "Do not plot as outcome; use for auditing and filters."
graphical_use[intersect(not_recommended, variables)] <- "Do not plot in the main analysis."

caution <- setNames(
  rep("Document categories, missing values, and historical changes.", length(variables)),
  variables
)
caution[intersect(c("edad", "uni_med"), variables)] <- "Interpret together; EDAD is not always expressed in years."
caution[intersect(c("tip_cas", "tipo_caso_inicial_etiqueta"), variables)] <- "Initial classification; does not replace final status."
caution[intersect(c("con_fin", "condicion_final_etiqueta", "fec_def", "cer_def", "cbmte"), variables)] <- "Very few deceased in event 210; analyze mortality with event 580."
caution[intersect(c("departamento_ocurrencia", "municipio_ocurrencia", "divipola_ocurrencia"), variables)] <- "Do not mix occurrence with residence or notification."
caution[intersect(c("departamento_residencia", "municipio_residencia", "divipola_residencia"), variables)] <- "For rates, use resident population denominators."
caution[intersect(c("departamento_notificacion", "municipio_notificacion", "divipola_notificacion"), variables)] <- "Represents where it was notified, not where it was transmitted."
caution[intersect(c("nacionalidad", "nombre_nacionalidad", "cod_pais_r", "pais_residencia", "gru_pob", "nom_grupo", "estrato"), variables)] <- "High historical absence; restrict to comparable periods or exclude."
caution[intersect(c("ocupacion", "cod_ase", "cbmte"), variables)] <- "Code that requires a reference table or grouping."
caution[intersect("cod_eve__1", variables)] <- "Technical duplicate exclusive to 2016; matches COD_EVE and should not be used as a second variable."
caution[intersect(c("fm_fuerza", "fm_unidad", "fm_grado"), variables)] <- "More than 99% missing; insufficient for main analysis."
caution[intersect(c("gp_otros", "fuente"), variables)] <- "Insufficient utility or semantics for central analysis."
caution[grep("^flag_|^grupo_", variables, value = TRUE)] <- "Auditing variable; not an original clinical feature."

description <- setNames(paste0("Variable ", variables, "."), variables)
description[c("ano", "semana", "fec_not")] <- c("Epidemiological year.", "Epidemiological week (1 to 53).", "Notification date.")
description[c("edad", "uni_med", "edad_anios", "grupo_edad_tesis")] <- c("Original age value.", "Original age unit.", "Age converted to years.", "Derived age group for the thesis.")
description[c("tip_cas", "ajuste", "estado_final_de_caso", "nom_est_f_caso")] <- c("Initial case classification.", "Last recorded adjustment.", "Consolidated final status code.", "Consolidated final status label.")
description[c("con_fin", "fec_def", "cer_def", "cbmte")] <- c("Final condition: alive or dead.", "Date of death, when applicable.", "Death certificate, when applicable.", "Basic cause of death, when applicable.")
description[c("cod_eve", "cod_eve__1")] <- c("Official event code; for DENGUE it must be 210.", "Second COD_EVE header present only in 2016; identical content to the main one.")
description[grep("^dias_", variables, value = TRUE)]  <- "Interval in days derived from clinical dates."
description[grep("^flag_", variables, value = TRUE)]  <- "Logical quality control flag."
description[grep("^grupo_duplicado|^grupo_posible", variables, value = TRUE)] <- "Internal identifier of the matching group."

dict_url <- "https://www.ins.gov.co/BibliotecaDigital/diccionario-datos-sivigila-2026-VF.pdf"
manual_url      <- "https://www.ins.gov.co/BibliotecaDigital/manual-del-usuario-sivigila-4-0.pdf"
ficha_url       <- "https://www.ins.gov.co/Direcciones/Vigilancia/sivigila/FichasdeNotificacion/DENGUE%20F210-220-580.pdf"

doc_source <- setNames(rep(dict_url, length(variables)), variables)
doc_source[intersect(control_vars, variables)] <- manual_url
doc_source[setdiff(variables, source_vars)] <-
  "Calculation documented in R/09_prepare_dengue.R"
doc_source[intersect(c("cod_eve", "nombre_evento"), variables)] <- ficha_url

## -----------------------------------------------------------------------------
## 7. Final dictionary
## -----------------------------------------------------------------------------
dictionary <- cbind(
  global_profile,
  data.frame(
    grupo             = unname(group[global_profile$variable]),
    descripcion       = unname(description[global_profile$variable]),
    prioridad_tesis   = unname(priority[global_profile$variable]),
    grafico_sugerido  = unname(graphical_use[global_profile$variable]),
    cautela           = unname(caution[global_profile$variable]),
    fuente_documental = unname(doc_source[global_profile$variable]),
    stringsAsFactors  = FALSE
  )
)

priority_order <- c("Essential", "Complementary",
                    "Control/traceability",
                    "Not recommended for main analysis")
dictionary <- dictionary[order(match(dictionary$prioridad_tesis, priority_order),
                               dictionary$grupo,
                               dictionary$variable), ]
rownames(dictionary) <- NULL

write.csv(dictionary,
          file.path(reports_dir, "diccionario_variables_dengue.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## -----------------------------------------------------------------------------
## 8. Consolidated quality report
## -----------------------------------------------------------------------------
flags_labels <- c(
  flag_cod_evento_no_210                      = "Event code different from 210",
  flag_cod_evento_duplicado_inconsistente     = "Second COD_EVE from 2016 inconsistent",
  flag_unidad_edad_invalida                   = "Invalid age unit",
  flag_edad_no_calculable                     = "Age not calculable in years",
  flag_edad_fuera_rango                       = "Age out of 0 to 130 years range",
  flag_edad_fecha_nacimiento_discordante      = "Discordant age with birth date (>2 years)",
  flag_inicio_sintomas_faltante               = "Symptom onset missing",
  flag_fecha_consulta_faltante                = "Consultation date missing",
  flag_hospitalizacion_inconsistente          = "Hospitalization inconsistent with FEC_HOS",
  flag_cronologia_clinica_inconsistente       = "Clinical chronology with negative interval",
  flag_ajuste_anterior_notificacion           = "Adjustment prior to notification",
  flag_muerto_sin_fecha_defuncion             = "Dead without date of death",
  flag_vivo_con_fecha_defuncion               = "Alive with date of death",
  flag_muerto_sin_certificado                 = "Dead without death certificate",
  flag_muerto_sin_causa_basica                = "Dead without basic cause",
  flag_gestante_sin_semana_gestacion          = "Pregnant without gestation week",
  flag_semana_gestacion_sin_gestante          = "Gestation week reported without marking pregnant",
  flag_anio_fec_not_distinto                  = "FEC_NOT year different from ANO",
  flag_confirmados_estado_final_inconsistente = "CONFIRMADOS inconsistent with final status",
  flag_duplicado_consecutive                  = "Repeated consecutive identifier",
  flag_duplicado_contenido_analitico          = "Matching analytical content",
  flag_posible_duplicado_clave_estricta       = "Matching strict epidemiological key",
  flag_cert_def_sentinel                      = "Death certificate with textual sentinel",
  flag_residencia_municipio_cero              = "Municipality of residence at sentinel 000",
  flag_municipio_texto_desconocido            = "Municipality name marked as unknown"
)

quality_flags <- grep("^flag_", names(annual_quality), value = TRUE)
flags_totals <- vapply(quality_flags,
                       function(f) sum(annual_quality[[f]], na.rm = TRUE),
                       numeric(1))

quality_summary <- data.frame(
  indicador  = c("Kept records", "Annual files", "Period",
                 unname(flags_labels[quality_flags])),
  valor      = c(sum(annual_quality$n_registros), nrow(annual_quality),
                 paste0(min(SIVIGILA_ANIOS), "-", max(SIVIGILA_ANIOS)),
                 flags_totals),
  porcentaje = c(NA, NA, NA,
                 100 * flags_totals / sum(annual_quality$n_registros)),
  stringsAsFactors = FALSE
)

write.csv(quality_summary,
          file.path(reports_dir, "resumen_calidad.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## -----------------------------------------------------------------------------
## 9. Global distributions
## -----------------------------------------------------------------------------
annual_dist <- read_report("distribuciones_clave_por_anio.csv")
global_dist <- aggregate(n ~ variable + valor, data = annual_dist, sum)
global_dist <- global_dist[order(global_dist$variable, -global_dist$n), ]
rownames(global_dist) <- NULL

write.csv(global_dist,
          file.path(reports_dir, "distribuciones_clave_globales.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## -----------------------------------------------------------------------------
## 10. Graph proposals and guide
## -----------------------------------------------------------------------------
proposals <- data.frame(
  prioridad = c("High", "High", "High", "High", "High", "High",
                "High", "High", "Medium", "Medium", "Medium", "Control"),
  pregunta = c(
    "How did reported dengue change between 2007 and 2025?",
    "Is there seasonality by epidemiological week?",
    "What ages and sexes concentrate the notifications?",
    "Where is it concentrated according to occurrence and residence?",
    "How did territorial incidence change when incorporating population?",
    "How did the final probable/confirmed classification evolve?",
    "What proportion was hospitalized and how timely was care?",
    "What proportion of notifications corresponded to severe dengue?",
    "How does it vary by area, insurance, and ethnic belonging?",
    "What vulnerable populations are registered?",
    "How do dengue, severe dengue, and mortality relate over the period?",
    "What quality and duplication problems change by year?"
  ),
  variables = c(
    "ano", "ano + semana", "grupo_edad_tesis + sexo_etiqueta",
    "divipola_ocurrencia + divipola_residencia",
    "casos + population by territory/year", "nom_est_f_caso + ano",
    "hospitalizado_etiqueta + dias_inicio_consulta + dias_inicio_hospitalizacion",
    "event 210 + event 220",
    "area_etiqueta + regimen_salud_etiqueta + pertenencia_etnica_etiqueta",
    "gp_discapa ... gp_vic_vio", "events 210 + 220 + 580", "ano + flag_*"
  ),
  grafico = c(
    "Annual line", "Year × week heatmap",
    "Pyramid or bars by age and sex", "Separate maps",
    "Maps/rate trends", "100% stacked bars per year",
    "Hospitalization trend and medians/IQR for delays",
    "Annual line of severe proportion", "Composition bars by period",
    "Descriptive prevalence bars", "Comparative panel of the three events",
    "Flags heatmap"
  ),
  cautela = c(
    "Counts, not incidence.", "Compare low-volume years with caution.",
    "Use edad_anios, not EDAD alone.",
    "Do not mix occurrence, residence, and notification.",
    "The denominator must match residence and period.",
    "TIP_CAS is initial; prefer final status.",
    "Exclude negative intervals only from the corresponding calculation and report them.",
    "Keep events 210 and 220 separate and define the denominator.",
    "Monitor missing values and historical changes.", "Group small cells.",
    "Cross-referencing events requires an explicit linking and deduplication strategy.",
    "Flags are controls, not clinical outcomes."
  ),
  stringsAsFactors = FALSE
)
write.csv(proposals,
          file.path(reports_dir, "propuestas_graficas_tesis.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

originals        <- dictionary[dictionary$origen == "Harmonized original", ]
priority_counts <- table(originals$prioridad_tesis)

guide <- c(
  "# Dengue Analysis Guide - Master in Biology UIS",
  paste("Period:", min(SIVIGILA_ANIOS), "to", max(SIVIGILA_ANIOS)),
  paste("Original records:", sum(annual_quality$n_registros)),
  "Canonical bases keep all original records.",
  "For analysis, use the analytical version without duplicates after the full pipeline.",
  "Check the variable dictionary, quality reports, and deduplication audit.",
  "Maps are counts; there are no population denominators.",
  "Documentation is regenerated by running R/10_document_dengue.R again."
)
writeLines(guide, file.path(reports_dir, "GUIA_ANALISIS_TESIS.md"), useBytes = TRUE)

## -----------------------------------------------------------------------------
## 11. Final validations
## -----------------------------------------------------------------------------
# The control_interno/ folder contains:
#   - intervalos_YYYY.rds (read in section 2)
#   - claves/contenido_YYYY.rds, estricta_YYYY.rds, id_YYYY.rds (used by 13)
# That's why it is NOT deleted here. If at some point you want to free up space,
# do it manually when you are sure that 13 has already run and will not be
# run again on the same data.
message("Intermediate artifacts preserved in: ", control_dir)

stopifnot(
  setequal(dictionary$variable, unique(annual_profile$variable))
)

cat("DENGUE documentation completed.\n")
print(as.data.frame(priority_counts), row.names = FALSE)