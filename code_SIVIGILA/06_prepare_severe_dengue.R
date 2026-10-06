#### 06_prepare_severe_dengue.R ####

source("01_load_databases.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## PREPARATION OF THE HISTORICAL SEVERE DENGUE DATABASE (EVENT 220)
## Unifies the data, calculates derived variables (age in years, care times)
## and marks errors or duplicates without deleting any records.
## Auxiliary functions live in 99_utils.R.
## -----------------------------------------------------------------------------

## 1. Output folders
results_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue_grave")
data_dir    <- file.path(results_dir, "datos_procesados")
reports_dir <- file.path(results_dir, "reportes")
dir.create(data_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

## 2. Loading and consolidating the raw database
raw_tables <- cargar_evento("dengue_grave")

# If it is already a single table, we use it; if it is a list, we bind it.
if (is.data.frame(raw_tables)) {
  text_data   <- raw_tables
  text_tables <- split(text_data, as.character(text_data$fuente_anio))
} else {
  text_data   <- dplyr::bind_rows(raw_tables)
  text_tables <- raw_tables
}

provenance_columns <- c("fuente_evento", "fuente_archivo",
                        "fuente_anio", "fuente_fila_excel")
source_columns <- setdiff(names(text_data), provenance_columns)

# Strict validation against the inventory
stopifnot(
  nrow(text_data) == filas_originales_esperadas("dengue_grave")
)

## 3. Schema reports by year
schemas_by_year <- do.call(rbind, lapply(names(text_tables), function(anio) {
  tabla <- text_tables[[anio]]
  data.frame(
    ano = as.integer(anio),
    registros = nrow(tabla),
    columnas_fuente = length(setdiff(names(tabla), provenance_columns)),
    variables_exclusivas_o_ausentes = paste(setdiff(source_columns, names(tabla)),
                                            collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
rownames(schemas_by_year) <- NULL
write.csv(schemas_by_year,
          file.path(reports_dir, "esquemas_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

availability_by_year <- do.call(rbind, lapply(names(text_tables), function(anio) {
  tabla <- text_tables[[anio]]
  data.frame(
    ano = as.integer(anio),
    variable = source_columns,
    presente_en_archivo = source_columns %in% names(tabla),
    stringsAsFactors = FALSE
  )
}))
rownames(availability_by_year) <- NULL
write.csv(availability_by_year,
          file.path(reports_dir, "disponibilidad_variables_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## Copy of harmonized raw data
saveRDS(text_data,
        file.path(data_dir, "dengue_grave_original_armonizada.rds"),
        compress = "xz")

## 4. Type transformation
datos <- text_data
date_columns <- intersect(c(
  "fec_not", "fec_con", "ini_sin", "fec_hos", "fec_def", "fecha_nto",
  "fec_arc_xl", "fec_aju"
), names(datos))
datos[date_columns] <- lapply(datos[date_columns], parsear_fecha)

for (variable in intersect(c("ano", "semana", "edad", "uni_med", "sem_ges"), names(datos))) {
  datos[[variable]] <- suppressWarnings(as.integer(datos[[variable]]))
}

code_widths <- c(
  cod_pais_o = 3L, cod_dpto_o = 2L, cod_mun_o = 3L,
  cod_pais_r = 3L, cod_dpto_r = 2L, cod_mun_r = 3L,
  cod_dpto_n = 2L, cod_pre = 10L, cod_sub = 2L
)
for (variable in intersect(names(code_widths), names(datos))) {
  datos[[variable]] <- normalizar_codigo(datos[[variable]], code_widths[[variable]])
}

## 5. Derived variables (DIVIPOLA, age, times)
datos$divipola_ocurrencia <- ifelse(
  is.na(datos$cod_dpto_o) | is.na(datos$cod_mun_o),
  NA_character_,
  paste0(datos$cod_dpto_o, datos$cod_mun_o)
)
datos$divipola_residencia <- ifelse(
  is.na(datos$cod_dpto_r) | is.na(datos$cod_mun_r),
  NA_character_,
  paste0(datos$cod_dpto_r, datos$cod_mun_r)
)
datos$divipola_notificacion <- unir_divipola_notificacion(datos$cod_dpto_n, datos$cod_mun_n)
datos$llave_upgd <- ifelse(
  is.na(datos$cod_pre) | is.na(datos$cod_sub),
  NA_character_,
  paste0(datos$cod_pre, "-", datos$cod_sub)
)

datos$edad_anios       <- calcular_edad_anios(datos$edad, datos$uni_med)
datos$grupo_edad_tesis <- grupo_edad_tesis(datos$edad_anios)

## 6. Plain text labels
datos$sexo_etiqueta               <- etiquetar_codigo(datos$sexo,   c(M = "Masculino", F = "Femenino", I = "Indeterminado"))
datos$area_etiqueta               <- etiquetar_codigo(datos$area,   c(`1` = "Cabecera municipal", `2` = "Centro poblado", `3` = "Rural disperso"))
datos$regimen_salud_etiqueta      <- etiquetar_codigo(datos$tip_ss, c(C = "Contributivo", S = "Subsidiado", P = "Excepción", E = "Especial", N = "No asegurado", I = "Indeterminado/Pendiente"))
datos$pertenencia_etnica_etiqueta <- etiquetar_codigo(datos$per_etn, c(`1` = "Indígena", `2` = "ROM/Gitano", `3` = "Raizal", `4` = "Palenquero", `5` = "Negro, mulato o afrocolombiano", `6` = "Otro"))
datos$tipo_caso_inicial_etiqueta  <- etiquetar_codigo(datos$tip_cas, c(`1` = "Sospechoso", `2` = "Probable", `3` = "Confirmado por laboratorio", `4` = "Confirmado por clínica", `5` = "Confirmado por nexo epidemiológico"))
datos$hospitalizado_etiqueta      <- etiquetar_codigo(datos$pac_hos, c(`1` = "Sí", `2` = "No"))
datos$condicion_final_etiqueta    <- etiquetar_codigo(datos$con_fin, c(`0` = "No sabe/No responde", `1` = "Vivo", `2` = "Muerto"))
datos$ajuste_etiqueta             <- etiquetar_codigo(datos$ajuste,  c(`0` = "Sin ajuste/Primera vez", `3` = "Confirmado por laboratorio", `4` = "Confirmado por clínica", `5` = "Confirmado por nexo epidemiológico", `6` = "Descartado", `7` = "Otro ajuste", D = "Descarte por error de digitación"))

## 7. Care times in days
datos$dias_inicio_consulta          <- as.integer(datos$fec_con - datos$ini_sin)
datos$dias_consulta_hospitalizacion <- as.integer(datos$fec_hos - datos$fec_con)
datos$dias_inicio_hospitalizacion   <- as.integer(datos$fec_hos - datos$ini_sin)
datos$dias_inicio_notificacion      <- as.integer(datos$fec_not - datos$ini_sin)
datos$dias_notificacion_ajuste      <- as.integer(datos$fec_aju - datos$fec_not)
datos$dias_inicio_defuncion         <- as.integer(datos$fec_def - datos$ini_sin)

## 8. Duplicates marking (they are marked, not deleted)
content_exclusion_vars <- c(
  "fuente_evento", "fuente_archivo", "fuente_anio", "fuente_fila_excel",
  "separacion", "particion", "consecutive", "consecutive_origen",
  "fec_arc_xl", "fec_aju", "confirmados", "va_sispro",
  "estado_final_de_caso", "nom_est_f_caso", "nombre_evento", "cod_eve"
)
content_vars <- setdiff(names(text_data), content_exclusion_vars)
content_groups <- marcar_grupos_repetidos(crear_clave(text_data, content_vars))

id <- trimws(as.character(text_data$consecutive))
id[es_vacio(id)] <- paste0("<FALTANTE_", which(es_vacio(id)), ">")
id_groups <- marcar_grupos_repetidos(id)

strict_key_vars <- c(
  "cod_eve", "ano", "semana", "cod_pre", "cod_sub", "edad", "uni_med",
  "sexo", "cod_dpto_o", "cod_mun_o", "fec_not", "ini_sin"
)
key_groups <- marcar_grupos_repetidos(crear_clave(text_data, strict_key_vars))

datos$flag_duplicado_consecutive               <- id_groups$marca
datos$grupo_duplicado_consecutive              <- id_groups$grupo
datos$flag_duplicado_contenido_analitico       <- content_groups$marca
datos$grupo_duplicado_contenido_analitico      <- content_groups$grupo
datos$flag_posible_duplicado_clave_estricta    <- key_groups$marca
datos$grupo_posible_duplicado_clave_estricta   <- key_groups$grupo

## 9. Quality flags
datos$flag_cod_evento_no_220       <- is.na(datos$cod_eve) | datos$cod_eve != "220"
datos$flag_unidad_edad_invalida    <- is.na(datos$uni_med) | !(datos$uni_med %in% 1:5)
datos$flag_edad_no_calculable      <- is.na(datos$edad_anios)
datos$flag_edad_fuera_rango        <- !is.na(datos$edad_anios) & (datos$edad_anios < 0 | datos$edad_anios > 130)

age_by_date <- as.numeric(datos$fec_not - datos$fecha_nto) / 365.25
datos$flag_edad_fecha_nacimiento_discordante <- !is.na(age_by_date) &
  !is.na(datos$edad_anios) & abs(age_by_date - datos$edad_anios) > 2

datos$flag_inicio_sintomas_faltante  <- is.na(datos$ini_sin)
datos$flag_fecha_consulta_faltante   <- is.na(datos$fec_con)
datos$flag_hospitalizacion_inconsistente <- (datos$pac_hos == "1" & is.na(datos$fec_hos)) |
  (datos$pac_hos == "2" & !is.na(datos$fec_hos))
datos$flag_cronologia_clinica_inconsistente <-
  (!is.na(datos$dias_inicio_consulta)          & datos$dias_inicio_consulta < 0L) |
  (!is.na(datos$dias_consulta_hospitalizacion) & datos$dias_consulta_hospitalizacion < 0L) |
  (!is.na(datos$dias_inicio_hospitalizacion)   & datos$dias_inicio_hospitalizacion < 0L) |
  (!is.na(datos$dias_inicio_notificacion)      & datos$dias_inicio_notificacion < 0L) |
  (!is.na(datos$dias_inicio_defuncion)         & datos$dias_inicio_defuncion < 0L)

datos$flag_ajuste_anterior_notificacion    <- !is.na(datos$dias_notificacion_ajuste) &
  datos$dias_notificacion_ajuste < 0L
datos$flag_muerto_sin_fecha_defuncion      <- datos$con_fin == "2" & is.na(datos$fec_def)
datos$flag_vivo_con_fecha_defuncion        <- datos$con_fin == "1" & !is.na(datos$fec_def)
datos$flag_muerto_sin_certificado          <- datos$con_fin == "2" & es_vacio(datos$cer_def)
datos$flag_muerto_sin_causa_basica         <- datos$con_fin == "2" & es_vacio(datos$cbmte)
datos$flag_gestante_sin_semana_gestacion   <- datos$gp_gestan == "1" & is.na(datos$sem_ges)
datos$flag_semana_gestacion_sin_gestante   <- datos$gp_gestan != "1" & !is.na(datos$sem_ges)
datos$flag_anio_fec_not_distinto           <- !is.na(datos$fec_not) &
  as.integer(format(datos$fec_not, "%Y")) != datos$ano

# Flags derived from raw files review
datos$flag_cert_def_sentinel           <- datos$con_fin == "2" &
  !is.na(datos$cer_def) & es_vacio(datos$cer_def)
datos$flag_residencia_municipio_cero   <- !is.na(datos$cod_mun_r) & datos$cod_mun_r == "000"
datos$flag_municipio_texto_desconocido <- !is.na(datos$municipio_residencia) &
  grepl("DESCONOCIDO|SIN DATO|^\\*", datos$municipio_residencia)

## 10. Canonical base
saveRDS(datos, file.path(data_dir, "dengue_grave_canonica.rds"), compress = "xz")
write.csv(datos, file.path(data_dir, "dengue_grave_canonica.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## 11. Analytical base
analytical_vars <- unique(c(
  "fuente_archivo", "fuente_anio", "fuente_fila_excel", "consecutive",
  "ano", "semana", "fec_not", "ini_sin", "fec_con", "fec_hos", "fec_def", "fec_aju",
  "edad", "uni_med", "edad_anios", "grupo_edad_tesis", "sexo", "sexo_etiqueta",
  "area", "area_etiqueta", "ocupacion", "tip_ss", "regimen_salud_etiqueta",
  "cod_ase", "per_etn", "pertenencia_etnica_etiqueta", "gru_pob", "nom_grupo",
  "estrato", grep("^gp_", names(datos), value = TRUE), "sem_ges",
  "pais_ocurrencia", "cod_dpto_o", "cod_mun_o", "divipola_ocurrencia",
  "departamento_ocurrencia", "municipio_ocurrencia", "cod_pais_r", "pais_residencia",
  "cod_dpto_r", "cod_mun_r", "divipola_residencia", "departamento_residencia",
  "municipio_residencia", "cod_dpto_n", "cod_mun_n", "divipola_notificacion",
  "departamento_notificacion", "municipio_notificacion", "cod_pre", "cod_sub",
  "llave_upgd", "nom_upgd", "tip_cas", "tipo_caso_inicial_etiqueta",
  "pac_hos", "hospitalizado_etiqueta", "con_fin", "condicion_final_etiqueta",
  "ajuste", "ajuste_etiqueta", "estado_final_de_caso", "nom_est_f_caso",
  "cer_def", "cbmte", grep("^dias_", names(datos), value = TRUE),
  grep("^flag_", names(datos), value = TRUE),
  grep("^grupo_duplicado|^grupo_posible", names(datos), value = TRUE)
))
analytical_vars <- intersect(analytical_vars, names(datos))
analytical_data <- datos[, analytical_vars, drop = FALSE]

saveRDS(analytical_data, file.path(data_dir, "dengue_grave_analitica.rds"), compress = "xz")
write.csv(analytical_data, file.path(data_dir, "dengue_grave_analitica.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## 12. Quality and duplicate reports
duplicates_summary <- data.frame(
  control = c(
    "Repeated consecutive identifier",
    "Matching analytical content, without identifiers or technical fields",
    "Match in strict epidemiological key"
  ),
  filas_marcadas = c(
    sum(datos$flag_duplicado_consecutive),
    sum(datos$flag_duplicado_contenido_analitico),
    sum(datos$flag_posible_duplicado_clave_estricta)
  ),
  grupos = c(
    length(unique(na.omit(datos$grupo_duplicado_consecutive))),
    length(unique(na.omit(datos$grupo_duplicado_contenido_analitico))),
    length(unique(na.omit(datos$grupo_posible_duplicado_clave_estricta)))
  ),
  decision = c(
    "Keep and review provenance: CONSECUTIVE can be reused between years; it is not enough to exclude rows.",
    "Authorized redundancy to exclude from analysis using R/08; keep canonical database intact.",
    "Could be a coincidence or a legitimate adjustment; keep in the analyses."
  ),
  stringsAsFactors = FALSE
)
write.csv(duplicates_summary, file.path(reports_dir, "resumen_duplicados.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

possible_dups <- datos[datos$flag_duplicado_contenido_analitico |
                         datos$flag_posible_duplicado_clave_estricta, , drop = FALSE]
if (nrow(possible_dups)) {
  possible_dups$tipo_coincidencia <- ifelse(
    possible_dups$flag_duplicado_contenido_analitico,
    "Matching analytical content",
    "Only strict epidemiological key matches"
  )
  group_ids <- sort(unique(na.omit(possible_dups$grupo_posible_duplicado_clave_estricta)))
  differences <- setNames(vapply(group_ids, function(g) {
    x <- possible_dups[possible_dups$grupo_posible_duplicado_clave_estricta == g, , drop = FALSE]
    candidate_vars <- setdiff(names(text_data),
                              c("fuente_evento", "fuente_archivo",
                                "fuente_anio", "fuente_fila_excel"))
    differ <- candidate_vars[vapply(candidate_vars, function(variable) {
      val <- trimws(as.character(x[[variable]]))
      val[es_vacio(val)] <- "<FALTANTE>"
      length(unique(val)) > 1L
    }, logical(1))]
    paste(differ, collapse = " | ")
  }, character(1)), as.character(group_ids))
  possible_dups$variables_que_difieren <- differences[as.character(possible_dups$grupo_posible_duplicado_clave_estricta)]
  possible_cols <- unique(c(
    "grupo_duplicado_contenido_analitico", "grupo_posible_duplicado_clave_estricta",
    "tipo_coincidencia", "consecutive", "consecutive_origen", "fuente_archivo",
    "fuente_anio", "fuente_fila_excel", strict_key_vars, "area", "pac_hos",
    "con_fin", "ajuste", "fecha_nto", "fec_arc_xl", "fec_aju", "variables_que_difieren"
  ))
  possible_dups <- possible_dups[, intersect(possible_cols, names(possible_dups)), drop = FALSE]
}
write.csv(possible_dups, file.path(reports_dir, "posibles_duplicados.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

quality_indicators <- c(
  flag_cod_evento_no_220              = "Event code different from 220",
  flag_unidad_edad_invalida           = "Invalid age unit",
  flag_edad_no_calculable             = "Age not calculable in years",
  flag_edad_fuera_rango               = "Age out of 0 to 130 years range",
  flag_edad_fecha_nacimiento_discordante = "Discordant age with birth date (>2 years)",
  flag_inicio_sintomas_faltante       = "Symptom onset missing",
  flag_fecha_consulta_faltante        = "Consultation date missing",
  flag_hospitalizacion_inconsistente  = "Hospitalization inconsistent with FEC_HOS",
  flag_cronologia_clinica_inconsistente = "Clinical chronology with negative interval",
  flag_ajuste_anterior_notificacion   = "Adjustment prior to notification",
  flag_muerto_sin_fecha_defuncion     = "Dead without date of death",
  flag_vivo_con_fecha_defuncion       = "Alive with date of death",
  flag_muerto_sin_certificado         = "Dead without death certificate",
  flag_muerto_sin_causa_basica        = "Dead without basic cause",
  flag_gestante_sin_semana_gestacion  = "Pregnant without gestation week",
  flag_semana_gestacion_sin_gestante  = "Gestation week reported without being marked pregnant",
  flag_anio_fec_not_distinto          = "FEC_NOT year different from ANO",
  flag_cert_def_sentinel              = "Death certificate with textual sentinel",
  flag_residencia_municipio_cero      = "Municipality of residence at sentinel 000",
  flag_municipio_texto_desconocido    = "Municipality name marked as unknown"
)
quality_summary <- data.frame(
  indicador = c("Kept records", "Annual files", "Period",
                unname(quality_indicators)),
  valor = c(
    nrow(datos),
    length(unique(datos$fuente_archivo)),
    paste0(min(datos$ano), "-", max(datos$ano)),
    vapply(names(quality_indicators), function(v) sum(datos[[v]], na.rm = TRUE), numeric(1))
  ),
  stringsAsFactors = FALSE
)
quality_summary$accion_recomendada <- c(
  "Keep all cases; apply specific filters according to each analysis.",
  "Maintain provenance per file.",
  "Analyze as a historical series and monitor availability changes per year.",
  rep("Review the flag before using the related variable; do not delete the entire case.",
      length(quality_indicators))
)
write.csv(quality_summary, file.path(reports_dir, "resumen_calidad.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

flag_vars <- grep("^flag_", names(datos), value = TRUE)
quality_by_year <- do.call(rbind, lapply(split(datos, datos$ano), function(x) {
  counts <- vapply(x[, flag_vars, drop = FALSE],
                   function(z) sum(z, na.rm = TRUE), numeric(1))
  data.frame(
    ano = unique(x$ano)[[1L]],
    n_registros = nrow(x),
    as.list(counts),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))
rownames(quality_by_year) <- NULL
write.csv(quality_by_year, file.path(reports_dir, "calidad_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

interval_vars <- c(
  "dias_inicio_consulta", "dias_consulta_hospitalizacion",
  "dias_inicio_hospitalizacion", "dias_inicio_notificacion",
  "dias_notificacion_ajuste", "dias_inicio_defuncion"
)
intervals_summary <- do.call(rbind, lapply(interval_vars, function(variable) {
  x <- datos[[variable]]
  valid_vals <- x[!is.na(x)]
  q <- if (length(valid_vals)) {
    stats::quantile(valid_vals, c(0, .25, .5, .75, .95, 1), names = FALSE)
  } else rep(NA_real_, 6L)
  data.frame(
    intervalo = variable,
    n_disponibles = length(valid_vals),
    n_faltantes = sum(is.na(x)),
    n_negativos = sum(valid_vals < 0),
    minimo = q[1], p25 = q[2], mediana = q[3], p75 = q[4],
    p95 = q[5], maximo = q[6], stringsAsFactors = FALSE
  )
}))
rownames(intervals_summary) <- NULL
write.csv(intervals_summary, file.path(reports_dir, "resumen_intervalos.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## 13. Data dictionary
variables <- names(datos)
group <- setNames(rep("Other", length(variables)), variables)
assign_group <- function(name, vars) group[intersect(vars, variables)] <<- name

assign_group("Traceability and system", c("fuente_evento", "fuente_archivo", "fuente_anio", "fuente_fila_excel", "separacion", "particion", "consecutive", "consecutive_origen", "cod_eve", "nombre_evento", "confirmados", "va_sispro", "estado_final_de_caso", "nom_est_f_caso"))
assign_group("Time", c("fec_not", "semana", "ano", "fec_con", "ini_sin", "fec_hos", "fec_def", "fecha_nto", "fec_aju", "fec_arc_xl", grep("^dias_", variables, value = TRUE)))
assign_group("Demographics", c("edad", "uni_med", "edad_anios", "grupo_edad_tesis", "sexo", "sexo_etiqueta", "nacionalidad", "nombre_nacionalidad", "ocupacion"))
assign_group("Geography", c(grep("^cod_pais_|^cod_dpto_|^cod_mun_|^divipola_", variables, value = TRUE), "pais_ocurrencia", "departamento_ocurrencia", "municipio_ocurrencia", "pais_residencia", "departamento_residencia", "municipio_residencia", "departamento_notificacion", "municipio_notificacion", "area", "area_etiqueta"))
assign_group("Health services", c("cod_pre", "cod_sub", "llave_upgd", "nom_upgd", "tip_ss", "regimen_salud_etiqueta", "cod_ase", "pac_hos", "hospitalizado_etiqueta"))
assign_group("Ethnicity and vulnerability", c("per_etn", "pertenencia_etnica_etiqueta", "gru_pob", "nom_grupo", "estrato", grep("^gp_", variables, value = TRUE), "sem_ges", "fuente"))
assign_group("Classification and outcome", c("tip_cas", "tipo_caso_inicial_etiqueta", "con_fin", "condicion_final_etiqueta", "ajuste", "ajuste_etiqueta", "cer_def", "cbmte"))
assign_group("Military forces", c("fm_fuerza", "fm_unidad", "fm_grado"))
assign_group("Quality control", grep("^flag_|^grupo_", variables, value = TRUE))

essentials <- c(
  "ano", "semana", "fec_not", "ini_sin", "fec_con", "fec_hos", "edad", "uni_med",
  "edad_anios", "grupo_edad_tesis", "sexo", "sexo_etiqueta", "area", "area_etiqueta",
  "tip_ss", "regimen_salud_etiqueta", "per_etn", "pertenencia_etnica_etiqueta",
  "pac_hos", "hospitalizado_etiqueta", "con_fin", "condicion_final_etiqueta",
  "ajuste", "ajuste_etiqueta", "estado_final_de_caso", "nom_est_f_caso",
  "divipola_ocurrencia", "departamento_ocurrencia", "municipio_ocurrencia",
  "divipola_residencia", "departamento_residencia", "municipio_residencia",
  "dias_inicio_consulta", "dias_inicio_hospitalizacion", "dias_inicio_notificacion"
)
complementaries <- c(
  "ocupacion", "cod_ase", "nacionalidad", "nombre_nacionalidad", "gru_pob",
  "nom_grupo", "estrato", grep("^gp_", variables, value = TRUE), "sem_ges",
  "cod_pre", "cod_sub", "llave_upgd", "nom_upgd", "tip_cas", "tipo_caso_inicial_etiqueta",
  "fec_def", "cer_def", "cbmte", "dias_consulta_hospitalizacion",
  "dias_notificacion_ajuste", "dias_inicio_defuncion"
)
control_vars <- c(
  "fuente_evento", "fuente_archivo", "fuente_anio", "fuente_fila_excel",
  "separacion", "particion", "consecutive", "consecutive_origen", "cod_eve",
  "nombre_evento", "fec_arc_xl", "fec_aju", "confirmados", "va_sispro",
  grep("^flag_|^grupo_", variables, value = TRUE)
)
not_recommended <- c("fm_fuerza", "fm_unidad", "fm_grado", "separacion",
                     "particion", "gp_otros", "fuente")

priority <- setNames(rep("Complementary", length(variables)), variables)
priority[intersect(essentials, variables)]        <- "Essential"
priority[intersect(complementaries, variables)]   <- "Complementary"
priority[intersect(control_vars, variables)]           <- "Control/traceability"
priority[intersect(not_recommended, variables)]   <- "Not recommended for main analysis"

graphical_use <- setNames(rep("Only if it answers a specific question and there is enough data.", length(variables)), variables)
graphical_use[intersect(c("ano", "semana", "fec_not"), variables)] <- "Time series or year x week heatmap."
graphical_use[intersect(c("edad_anios", "grupo_edad_tesis", "sexo_etiqueta"), variables)] <- "Histogram, age/sex bars, or population pyramid."
graphical_use[intersect(c("departamento_ocurrencia", "municipio_ocurrencia", "divipola_ocurrencia"), variables)] <- "Occurrence map; use rates when a compatible denominator exists."
graphical_use[intersect(c("departamento_residencia", "municipio_residencia", "divipola_residencia"), variables)] <- "Residence map; align with resident population."
graphical_use[intersect(c("departamento_notificacion", "municipio_notificacion", "divipola_notificacion"), variables)] <- "Notification flow; do not interpret as transmission location."
graphical_use[intersect(c("area_etiqueta", "regimen_salud_etiqueta", "pertenencia_etnica_etiqueta"), variables)] <- "Composition bars and trends by period."
graphical_use[intersect(c("hospitalizado_etiqueta", grep("^dias_", variables, value = TRUE)), variables)] <- "Proportions, boxplots, or violins for care and evolution."
graphical_use[intersect(c("condicion_final_etiqueta", "con_fin"), variables)] <- "Proportion of deceased among severe dengue notifications."
graphical_use[intersect(c("ocupacion", "cod_ase", "nom_upgd", "cbmte"), variables)] <- "Group/decode first; then ordered bars or Pareto."
graphical_use[intersect(control_vars, variables)]          <- "Do not plot as outcome; use for auditing and filters."
graphical_use[intersect(not_recommended, variables)]  <- "Do not plot in the main analysis."

caution <- setNames(rep("Document categories, missing values, and historical changes.", length(variables)), variables)
caution[intersect(c("edad", "uni_med"), variables)] <- "Interpret together; EDAD is not always expressed in years."
caution[intersect(c("tip_cas", "tipo_caso_inicial_etiqueta"), variables)] <- "Initial classification; does not replace the final status of the case."
caution[intersect(c("con_fin", "condicion_final_etiqueta"), variables)] <- "Allows describing the outcome within event 220; it does not equate to the lethality of all dengue by itself."
caution[intersect(c("fec_def", "cer_def", "cbmte", "dias_inicio_defuncion"), variables)] <- "Evaluate only among the deceased; global absence is expected in living persons."
caution[intersect(c("nacionalidad", "nombre_nacionalidad", "cod_pais_r", "pais_residencia", "gru_pob", "nom_grupo", "estrato"), variables)] <- "High historical absence; restrict years or leave out of the main analysis."
caution[intersect(c("ocupacion", "cod_ase", "cbmte"), variables)] <- "Code that requires a reference table or grouping before interpretation."
caution[intersect(c("departamento_ocurrencia", "municipio_ocurrencia", "divipola_ocurrencia"), variables)] <- "Do not mix occurrence with residence or notification."
caution[intersect(c("departamento_residencia", "municipio_residencia", "divipola_residencia"), variables)] <- "For rates, use resident population denominators."
caution[intersect(c("departamento_notificacion", "municipio_notificacion", "divipola_notificacion"), variables)] <- "Represents where it was notified, not where it was transmitted."
caution[grep("^flag_|^grupo_", variables, value = TRUE)] <- "Auditing variable; not an original clinical feature."
caution[intersect(c("fm_fuerza", "fm_unidad", "fm_grado"), variables)] <- "More than 99.9% missing; insufficient for main analysis."
caution[intersect(c("gp_otros", "fuente"), variables)] <- "Semantics/filling less informative in this series; do not use without additional official clarification."

description <- setNames(paste0("Field ", variables, "."), variables)
description[c("ano", "semana", "fec_not")] <- c("Epidemiological year of the record.", "Epidemiological week (1 to 53).", "Notification date.")
description[c("edad", "uni_med", "edad_anios", "grupo_edad_tesis")] <- c("Original age value.", "Original age unit.", "Age converted to years.", "Derived age group for the thesis.")
description[c("cod_dpto_o", "cod_mun_o", "divipola_ocurrencia")] <- c("Department code of occurrence.", "Municipality code of occurrence.", "Derived DIVIPOLA code of occurrence.")
description[c("cod_dpto_r", "cod_mun_r", "divipola_residencia")] <- c("Department code of residence.", "Municipality code of residence.", "Derived DIVIPOLA code of residence.")
description[c("cod_dpto_n", "cod_mun_n", "divipola_notificacion")] <- c("Notifying department code.", "Notification municipality code; in some years it already contains the complete DIVIPOLA.", "Normalized DIVIPOLA notification code without duplicating the department.")
description[c("tip_cas", "ajuste", "estado_final_de_caso", "nom_est_f_caso")] <- c("Initial case classification.", "Last recorded adjustment.", "Consolidated final status code.", "Consolidated final status label.")
description[c("con_fin", "fec_def", "cer_def", "cbmte")] <- c("Final condition: alive or dead.", "Date of death, when applicable.", "Death certificate number, when applicable.", "Coded basic cause of death, when applicable.")
description[grep("^dias_", variables, value = TRUE)]  <- "Interval in days derived from clinical dates."
description[grep("^flag_", variables, value = TRUE)]  <- "Logical quality control flag."
description[grep("^grupo_duplicado|^grupo_posible", variables, value = TRUE)] <- "Internal identifier of the matching group."

profile <- do.call(rbind, lapply(variables, function(variable) {
  x <- datos[[variable]]
  empty_val <- es_vacio(x)
  data.frame(
    variable = variable,
    origen = if (variable %in% names(text_data)) "Harmonized original" else "Derived",
    tipo_r = paste(class(x), collapse = "/"),
    n_registros = length(x),
    n_faltantes = sum(empty_val),
    pct_faltantes = 100 * mean(empty_val),
    n_unicos_no_vacios = length(unique(x[!empty_val])),
    stringsAsFactors = FALSE
  )
}))

dict_url <- "https://www.ins.gov.co/BibliotecaDigital/diccionario-datos-sivigila-2026-VF.pdf"
manual_url      <- "https://www.ins.gov.co/BibliotecaDigital/manual-del-usuario-sivigila-4-0.pdf"
ficha_url       <- "https://www.ins.gov.co/Direcciones/Vigilancia/sivigila/FichasdeNotificacion/DENGUE%20F210-220-580.pdf"
doc_source <- setNames(rep(dict_url, length(variables)), variables)
doc_source[intersect(control_vars, variables)] <- manual_url
doc_source[setdiff(variables, names(text_data))] <- "Calculation documented in R/06_prepare_severe_dengue.R"
doc_source[intersect(c("cod_eve", "nombre_evento"), variables)] <- ficha_url

dictionary <- cbind(
  profile,
  data.frame(
    grupo = unname(group[profile$variable]),
    descripcion = unname(description[profile$variable]),
    prioridad_tesis = unname(priority[profile$variable]),
    grafico_sugerido = unname(graphical_use[profile$variable]),
    cautela = unname(caution[profile$variable]),
    fuente_documental = unname(doc_source[profile$variable]),
    stringsAsFactors = FALSE
  )
)
dictionary <- dictionary[order(
  match(dictionary$prioridad_tesis,
        c("Essential", "Complementary", "Control/traceability",
          "Not recommended for main analysis")),
  dictionary$grupo, dictionary$variable
), ]
rownames(dictionary) <- NULL
write.csv(dictionary,
          file.path(reports_dir, "diccionario_variables_dengue_grave.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

graphical_proposals <- data.frame(
  prioridad = c("High", "High", "High", "High", "High", "High", "High", "Medium", "Medium", "Control"),
  pregunta = c(
    "How did the number of severe cases reported change between 2007 and 2025?",
    "Is there seasonality by epidemiological week?",
    "Which age and sex groups concentrate severe dengue?",
    "Where is it concentrated according to occurrence and residence?",
    "How timely were consultation, hospitalization, and notification?",
    "What proportion was hospitalized and how did it change over time?",
    "How is the alive/dead outcome distributed within event 220?",
    "How does it vary by area, insurance, and ethnic belonging?",
    "What vulnerable populations appear among severe cases?",
    "What quality and duplication problems change by year?"
  ),
  variables = c(
    "ano", "ano + semana", "grupo_edad_tesis + sexo_etiqueta",
    "divipola_ocurrencia + divipola_residencia",
    "dias_inicio_consulta + dias_inicio_hospitalizacion + dias_inicio_notificacion",
    "hospitalizado_etiqueta + ano", "condicion_final_etiqueta + ano/grupo/territorio",
    "area_etiqueta + regimen_salud_etiqueta + pertenencia_etnica_etiqueta",
    "gp_discapa ... gp_vic_vio", "ano + flag_*"
  ),
  grafico = c(
    "Annual line with points", "Year × week heatmap",
    "Pyramid or bars by age group and sex",
    "Separate maps; prefer rates with compatible denominators",
    "Boxplots or violins by year/period", "Proportion bars by year",
    "Proportion of deceased with confidence intervals",
    "100% stacked bars by period", "Descriptive prevalence bars",
    "Heatmap of quality flags"
  ),
  cautela = c(
    "These are counts; denominators are needed for incidence.",
    "Use caution when comparing low-volume years.",
    "Use edad_anios, not EDAD alone.",
    "Do not mix occurrence, residence, and notification.",
    "Exclude negative intervals only from the calculation and report them.",
    "Verify consistency between PAC_HOS and FEC_HOS.",
    "It is lethality among severe notifications, not among all dengue cases.",
    "Monitor missing values and capture changes by year.",
    "Group small cells to protect confidentiality.",
    "Flags are controls, not clinical outcomes."
  ),
  stringsAsFactors = FALSE
)
write.csv(graphical_proposals,
          file.path(reports_dir, "propuestas_graficas_tesis.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 14. Final validations
stopifnot(
  nrow(text_data) == nrow(datos),
  nrow(datos) == nrow(analytical_data),
  nrow(datos) == filas_originales_esperadas("dengue_grave"),
  # Fails only if there is a non-NA value different from "220".
  # NA values are allowed and already captured by flag_cod_evento_no_220.
  !any(datos$cod_eve != "220" & !is.na(datos$cod_eve)),
  !anyDuplicated(clave_registro(datos)),
  all(nchar(datos$divipola_ocurrencia[!is.na(datos$divipola_ocurrencia)]) == 5L),
  all(nchar(datos$divipola_residencia[!is.na(datos$divipola_residencia)]) == 5L),
  all(nchar(datos$divipola_notificacion[!is.na(datos$divipola_notificacion)]) == 5L)
)

cat("Severe dengue preparation completed.\n")
cat("Kept records:", nrow(datos), "\n")
cat("Annual files:", length(unique(datos$fuente_archivo)), "\n")
cat("Harmonized original columns:", ncol(text_data), "\n")
cat("Canonical columns, derived and flags:", ncol(datos), "\n")
print(duplicates_summary, row.names = FALSE)
