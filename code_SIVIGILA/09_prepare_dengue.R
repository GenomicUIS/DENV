source("01_load_databases.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## HISTORICAL PREPARATION OF DENGUE (EVENT 210), 2007-2025
## More than 1.5 million records. Processed year by year to avoid overloading RAM.
## Auxiliary functions live in 99_utils.R.
## -----------------------------------------------------------------------------

## 1. Folder structure
results_dir  <- file.path(SIVIGILA_PATH, "resultados", "dengue")
data_dir       <- file.path(results_dir, "datos_procesados")
canonical_dir    <- file.path(data_dir, "canonica_por_anio")
analytical_dir   <- file.path(data_dir, "analitica_por_anio")
reports_dir    <- file.path(results_dir, "reportes")
control_dir     <- file.path(data_dir, "control_interno")
provisional_dir <- file.path(control_dir, "provisional")
keys_dir      <- file.path(control_dir, "claves")

invisible(lapply(
  c(canonical_dir, analytical_dir, reports_dir, provisional_dir, keys_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

## 2. Inventory
files <- inventario_sivigila[tolower(inventario_sivigila$evento) == "dengue", ]
files <- files[order(files$anio_archivo), ]

stopifnot(
  nrow(files) == length(SIVIGILA_ANIOS),
  files$anio_archivo == SIVIGILA_ANIOS
)

## 3. Pre-analysis of headers
headers <- lapply(seq_len(nrow(files)), function(i) {
  path <- files$ruta[[i]]
  names_found <- names(readxl::read_excel(
    path, n_max = 0, .name_repair = "minimal", progress = FALSE
  ))
  limpiar_nombres_sivigila(names_found)
})
names(headers) <- as.character(files$anio_archivo)
source_columns    <- unique(unlist(headers, use.names = FALSE))
provenance_columns <- c("fuente_evento", "fuente_archivo",
                        "fuente_anio", "fuente_fila_excel")

schemas_by_year <- do.call(rbind, lapply(names(headers), function(anio) {
  vars <- headers[[anio]]
  data.frame(
    ano = as.integer(anio),
    archivo = files$archivo[files$anio_archivo == as.integer(anio)],
    registros_esperados = NA_integer_,
    columnas_fuente = length(vars),
    variables_ausentes = paste(setdiff(source_columns, vars), collapse = " | "),
    variables_especiales = paste(intersect(vars, c("separacion", "particion", "cod_eve__1")),
                                 collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
rownames(schemas_by_year) <- NULL

availability_by_year <- do.call(rbind, lapply(names(headers), function(anio) {
  data.frame(
    ano = as.integer(anio), variable = source_columns,
    presente_en_archivo = source_columns %in% headers[[anio]],
    stringsAsFactors = FALSE
  )
}))
rownames(availability_by_year) <- NULL

## Duplicates criteria
content_exclusion_vars <- c(
  "separacion", "particion",
  "consecutive", "consecutive2", "consecutive_12", "consecutive_origen",
  "fec_arc_xl", "fec_aju", "confirmados", "va_sispro",
  "estado_final_de_caso", "nom_est_f_caso", "nombre_evento",
  "cod_eve", "cod_eve__1"
)
content_vars <- setdiff(source_columns, content_exclusion_vars)
strict_key_vars <- c(
  "cod_eve", "ano", "semana", "cod_pre", "cod_sub", "edad", "uni_med",
  "sexo", "cod_dpto_o", "cod_mun_o", "fec_not", "ini_sin"
)

date_columns <- c(
  "fec_not", "fec_con", "ini_sin", "fec_hos", "fec_def", "fecha_nto",
  "fec_arc_xl", "fec_aju"
)

rows_per_year <- integer(length = nrow(files))
names(rows_per_year) <- as.character(files$anio_archivo)

## Reading an annual file
# In addition to cleaning the names, we standardize columns that SIVIGILA
# delivered with alternative names in some years (CONSECUTIVE2,
# CONSECUTIVE_12, NOMBRE_UPGD). The duplicated COD_EVE that appears 
# in the 2016 file is also discarded.
leer_archivo_sivigila <- function(ruta, evento, agregar_procedencia = TRUE) {
  datos <- readxl::read_excel(ruta, col_types = "text",
                              .name_repair = "minimal", progress = FALSE)
  names(datos) <- limpiar_nombres_sivigila(names(datos))
  
  # Historical column name corrections
  if ("consecutive2" %in% names(datos)) {
    names(datos)[names(datos) == "consecutive2"] <- "consecutive"
  }
  if ("consecutive_12" %in% names(datos)) {
    names(datos)[names(datos) == "consecutive_12"] <- "consecutive"
  }
  if ("nombre_upgd" %in% names(datos)) {
    names(datos)[names(datos) == "nombre_upgd"] <- "nom_upgd"
  }
  
  # Some years bring duplicated COD_EVE. The first appearance is kept
  # to avoid breaking the column selection by name.
  if (anyDuplicated(names(datos))) {
    datos <- datos[, !duplicated(names(datos)), drop = FALSE]
  }
  
  if (agregar_procedencia) {
    datos$fuente_evento     <- evento
    datos$fuente_archivo    <- basename(ruta)
    detected_year          <- regmatches(basename(ruta),
                                         regexpr("[0-9]{4}", basename(ruta)))
    datos$fuente_anio       <- as.integer(detected_year)
    datos$fuente_fila_excel <- seq_len(nrow(datos)) + 1L
  }
  as.data.frame(datos)
}

## -----------------------------------------------------------------------------
## 4. ANNUAL PROCESSING
## -----------------------------------------------------------------------------
for (i in seq_len(nrow(files))) {
  anio <- files$anio_archivo[[i]]
  cat(sprintf("[%d/%d] Processing DENGUE %d...\n", i, nrow(files), anio))
  
  text_data <- leer_archivo_sivigila(files$ruta[[i]], evento = "dengue",
                                     agregar_procedencia = TRUE)
  
  # Align columns with the general schema
  missing_cols <- setdiff(source_columns, names(text_data))
  for (variable in missing_cols) text_data[[variable]] <- NA_character_
  text_data <- text_data[, c(provenance_columns, source_columns), drop = FALSE]
  
  rows_per_year[as.character(anio)] <- nrow(text_data)
  
  # Save keys for global duplicates analysis
  content_key <- crear_clave(text_data, content_vars)
  strict_key  <- crear_clave(text_data, strict_key_vars)
  identifier   <- trimws(as.character(text_data$consecutive))
  identifier[es_vacio(identifier)] <- paste0(
    "<FALTANTE_", anio, "_", which(es_vacio(identifier)), ">"
  )
  saveRDS(content_key, file.path(keys_dir, paste0("contenido_", anio, ".rds")), compress = "gzip")
  saveRDS(strict_key,  file.path(keys_dir, paste0("estricta_",  anio, ".rds")), compress = "gzip")
  saveRDS(identifier,   file.path(keys_dir, paste0("id_",        anio, ".rds")), compress = "gzip")
  rm(content_key, strict_key, identifier)
  
  # Type conversion
  datos <- text_data
  present_dates <- intersect(date_columns, names(datos))
  datos[present_dates] <- lapply(datos[present_dates], parsear_fecha)
  
  # Integer columns. estado_final_de_caso is included here because in several
  # years it arrives with trailing spaces ("3          ") and its comparison with
  # c(3L, 5L) requires it to be a clean numeric.
  integer_vars <- c("ano", "semana", "edad", "uni_med", "sem_ges",
                    "estado_final_de_caso")
  for (variable in intersect(integer_vars, names(datos))) {
    datos[[variable]] <- suppressWarnings(
      as.integer(trimws(as.character(datos[[variable]])))
    )
  }
  
  widths <- c(
    cod_pais_o = 3L, cod_dpto_o = 2L, cod_mun_o = 3L,
    cod_pais_r = 3L, cod_dpto_r = 2L, cod_mun_r = 3L,
    cod_dpto_n = 2L, cod_pre = 10L, cod_sub = 2L
  )
  for (variable in intersect(names(widths), names(datos))) {
    datos[[variable]] <- normalizar_codigo(datos[[variable]], widths[[variable]])
  }
  
  # Derived variables
  datos$divipola_ocurrencia <- ifelse(
    is.na(datos$cod_dpto_o) | is.na(datos$cod_mun_o), NA_character_,
    paste0(datos$cod_dpto_o, datos$cod_mun_o)
  )
  datos$divipola_residencia <- ifelse(
    is.na(datos$cod_dpto_r) | is.na(datos$cod_mun_r), NA_character_,
    paste0(datos$cod_dpto_r, datos$cod_mun_r)
  )
  datos$divipola_notificacion <- unir_divipola_notificacion(datos$cod_dpto_n, datos$cod_mun_n)
  datos$llave_upgd <- ifelse(
    is.na(datos$cod_pre) | is.na(datos$cod_sub), NA_character_,
    paste0(datos$cod_pre, "-", datos$cod_sub)
  )
  
  datos$edad_anios       <- calcular_edad_anios(datos$edad, datos$uni_med)
  datos$grupo_edad_tesis <- grupo_edad_tesis(datos$edad_anios)
  
  datos$sexo_etiqueta               <- etiquetar_codigo(datos$sexo,   c(M = "Masculino", F = "Femenino", I = "Indeterminado"))
  datos$area_etiqueta               <- etiquetar_codigo(datos$area,   c(`1` = "Cabecera municipal", `2` = "Centro poblado", `3` = "Rural disperso"))
  datos$regimen_salud_etiqueta      <- etiquetar_codigo(datos$tip_ss, c(C = "Contributivo", S = "Subsidiado", P = "Excepción", E = "Especial", N = "No asegurado", I = "Indeterminado/Pendiente"))
  datos$pertenencia_etnica_etiqueta <- etiquetar_codigo(datos$per_etn, c(`1` = "Indígena", `2` = "ROM/Gitano", `3` = "Raizal", `4` = "Palenquero", `5` = "Negro, mulato o afrocolombiano", `6` = "Otro"))
  datos$tipo_caso_inicial_etiqueta  <- etiquetar_codigo(datos$tip_cas, c(`1` = "Sospechoso", `2` = "Probable", `3` = "Confirmado por laboratorio", `4` = "Confirmado por clínica", `5` = "Confirmado por nexo epidemiológico"))
  datos$hospitalizado_etiqueta      <- etiquetar_codigo(datos$pac_hos, c(`1` = "Sí", `2` = "No"))
  datos$condicion_final_etiqueta    <- etiquetar_codigo(datos$con_fin, c(`0` = "No sabe/No responde", `1` = "Vivo", `2` = "Muerto"))
  datos$ajuste_etiqueta             <- etiquetar_codigo(datos$ajuste,  c(`0` = "Sin ajuste/Primera vez", `3` = "Confirmado por laboratorio", `4` = "Confirmado por clínica", `5` = "Confirmado por nexo epidemiológico", `6` = "Descartado", `7` = "Otro ajuste", D = "Descarte por error de digitación"))
  
  datos$dias_inicio_consulta          <- as.integer(datos$fec_con - datos$ini_sin)
  datos$dias_consulta_hospitalizacion <- as.integer(datos$fec_hos - datos$fec_con)
  datos$dias_inicio_hospitalizacion   <- as.integer(datos$fec_hos - datos$ini_sin)
  datos$dias_inicio_notificacion      <- as.integer(datos$fec_not - datos$ini_sin)
  datos$dias_notificacion_ajuste      <- as.integer(datos$fec_aju - datos$fec_not)
  datos$dias_inicio_defuncion         <- as.integer(datos$fec_def - datos$ini_sin)
  
  # Quality flags
  datos$flag_cod_evento_no_210 <- is.na(datos$cod_eve) | datos$cod_eve != "210"
  
  if ("cod_eve__1" %in% names(datos)) {
    datos$flag_cod_evento_duplicado_inconsistente <- !es_vacio(datos$cod_eve__1) &
      datos$cod_eve__1 != datos$cod_eve
  } else {
    datos$flag_cod_evento_duplicado_inconsistente <- FALSE
  }
  
  datos$flag_unidad_edad_invalida <- is.na(datos$uni_med) | !(datos$uni_med %in% 1:5)
  datos$flag_edad_no_calculable   <- is.na(datos$edad_anios)
  datos$flag_edad_fuera_rango     <- !is.na(datos$edad_anios) &
    (datos$edad_anios < 0 | datos$edad_anios > 130)
  
  age_by_date <- as.numeric(datos$fec_not - datos$fecha_nto) / 365.25
  datos$flag_edad_fecha_nacimiento_discordante <- !is.na(age_by_date) &
    !is.na(datos$edad_anios) & abs(age_by_date - datos$edad_anios) > 2
  
  datos$flag_inicio_sintomas_faltante <- is.na(datos$ini_sin)
  datos$flag_fecha_consulta_faltante  <- is.na(datos$fec_con)
  datos$flag_hospitalizacion_inconsistente <- (datos$pac_hos == "1" & is.na(datos$fec_hos)) |
    (datos$pac_hos == "2" & !is.na(datos$fec_hos))
  
  datos$flag_cronologia_clinica_inconsistente <-
    (!is.na(datos$dias_inicio_consulta)          & datos$dias_inicio_consulta < 0L) |
    (!is.na(datos$dias_consulta_hospitalizacion) & datos$dias_consulta_hospitalizacion < 0L) |
    (!is.na(datos$dias_inicio_hospitalizacion)   & datos$dias_inicio_hospitalizacion < 0L) |
    (!is.na(datos$dias_inicio_notificacion)      & datos$dias_inicio_notificacion < 0L) |
    (!is.na(datos$dias_inicio_defuncion)         & datos$dias_inicio_defuncion < 0L)
  
  datos$flag_ajuste_anterior_notificacion  <- !is.na(datos$dias_notificacion_ajuste) &
    datos$dias_notificacion_ajuste < 0L
  datos$flag_muerto_sin_fecha_defuncion    <- datos$con_fin == "2" & is.na(datos$fec_def)
  datos$flag_vivo_con_fecha_defuncion      <- datos$con_fin == "1" & !is.na(datos$fec_def)
  datos$flag_muerto_sin_certificado        <- datos$con_fin == "2" & es_vacio(datos$cer_def)
  datos$flag_muerto_sin_causa_basica       <- datos$con_fin == "2" & es_vacio(datos$cbmte)
  datos$flag_gestante_sin_semana_gestacion <- datos$gp_gestan == "1" & is.na(datos$sem_ges)
  datos$flag_semana_gestacion_sin_gestante <- datos$gp_gestan != "1" & !is.na(datos$sem_ges)
  datos$flag_anio_fec_not_distinto         <- !is.na(datos$fec_not) &
    as.integer(format(datos$fec_not, "%Y")) != datos$ano
  
  # estado_final_de_caso was already converted to clean integer; it is compared
  # against 3 and 5 ("Confirmado por laboratorio" and "Confirmado por
  # nexo epidemiológico" codes).
  expected_confirmation <- ifelse(
    datos$estado_final_de_caso %in% c(3L, 5L), "1", "0"
  )
  datos$flag_confirmados_estado_final_inconsistente <- !es_vacio(datos$confirmados) &
    datos$confirmados != expected_confirmation
  
  # Flags derived from raw files review.
  # They go INSIDE the loop because `datos` is freed at the end of each iteration.
  datos$flag_cert_def_sentinel           <- datos$con_fin == "2" &
    !is.na(datos$cer_def) & es_vacio(datos$cer_def)
  datos$flag_residencia_municipio_cero   <- !is.na(datos$cod_mun_r) & datos$cod_mun_r == "000"
  datos$flag_municipio_texto_desconocido <- !is.na(datos$municipio_residencia) &
    grepl("DESCONOCIDO|SIN DATO|^\\*", datos$municipio_residencia)
  
  # Critical validations of the block.
  # The cod_eve check fails only if there is a non-NA value different from
  # "210". NA values are allowed and already captured by
  # flag_cod_evento_no_210.
  stopifnot(
    !any(datos$cod_eve != "210" & !is.na(datos$cod_eve)),
    all(nchar(datos$divipola_ocurrencia[!is.na(datos$divipola_ocurrencia)]) == 5L),
    all(nchar(datos$divipola_residencia[!is.na(datos$divipola_residencia)]) == 5L),
    all(nchar(datos$divipola_notificacion[!is.na(datos$divipola_notificacion)]) == 5L)
  )
  
  # Save annual provisional and free memory
  saveRDS(datos,
          file.path(provisional_dir, paste0("dengue_canonica_", anio, ".rds")),
          compress = "gzip")
  intervals <- datos[, grep("^dias_", names(datos), value = TRUE), drop = FALSE]
  saveRDS(intervals,
          file.path(control_dir, paste0("intervalos_", anio, ".rds")),
          compress = "gzip")
  
  rm(datos, text_data, intervals, age_by_date, expected_confirmation)
  invisible(gc())
}

schemas_by_year$registros_esperados <- unname(rows_per_year[as.character(schemas_by_year$ano)])
write.csv(schemas_by_year,
          file.path(reports_dir, "esquemas_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(availability_by_year,
          file.path(reports_dir, "disponibilidad_variables_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Strict validation against the initial diagnosis
stopifnot(sum(rows_per_year) == filas_originales_esperadas("dengue"))

## -----------------------------------------------------------------------------
## 5. GLOBAL DUPLICATES REVIEW
## -----------------------------------------------------------------------------
combine_keys <- function(prefix) {
  unlist(lapply(SIVIGILA_ANIOS, function(anio) {
    readRDS(file.path(keys_dir, paste0(prefix, "_", anio, ".rds")))
  }), use.names = FALSE)
}

cat("Evaluating repeated identifiers across the series...\n")
key <- combine_keys("id")
id_duplicates <- marcar_grupos_repetidos(key)
rm(key); invisible(gc())

cat("Evaluating analytical content matching globally...\n")
key <- combine_keys("contenido")
content_duplicates <- marcar_grupos_repetidos(key)
rm(key); invisible(gc())

cat("Evaluating strict epidemiological key matches...\n")
key <- combine_keys("estricta")
strict_duplicates <- marcar_grupos_repetidos(key)
rm(key); invisible(gc())

## -----------------------------------------------------------------------------
## 6. DEFINITIVE BASES BY YEAR (CANONICAL AND ANALYTICAL)
## -----------------------------------------------------------------------------
offsets <- c(0L, cumsum(rows_per_year))
analytical_vars_base <- c(
  "fuente_archivo", "fuente_anio", "fuente_fila_excel", "consecutive",
  "ano", "semana", "fec_not", "ini_sin", "fec_con", "fec_hos", "fec_def", "fec_aju",
  "edad", "uni_med", "edad_anios", "grupo_edad_tesis", "sexo", "sexo_etiqueta",
  "area", "area_etiqueta", "ocupacion", "tip_ss", "regimen_salud_etiqueta",
  "cod_ase", "per_etn", "pertenencia_etnica_etiqueta", "gru_pob", "nom_grupo",
  "estrato", "sem_ges", "cod_dpto_o", "cod_mun_o", "divipola_ocurrencia",
  "departamento_ocurrencia", "municipio_ocurrencia", "cod_dpto_r", "cod_mun_r",
  "divipola_residencia", "departamento_residencia", "municipio_residencia",
  "cod_dpto_n", "cod_mun_n", "divipola_notificacion", "departamento_notificacion",
  "municipio_notificacion", "cod_pre", "cod_sub", "llave_upgd", "nom_upgd",
  "tip_cas", "tipo_caso_inicial_etiqueta", "pac_hos", "hospitalizado_etiqueta",
  "con_fin", "condicion_final_etiqueta", "ajuste", "ajuste_etiqueta",
  "estado_final_de_caso", "nom_est_f_caso", "cer_def", "cbmte"
)

annual_profile <- list()
annual_quality <- list()
annual_distributions <- list()
candidates <- list()
distribution_vars <- c(
  "sexo_etiqueta", "area_etiqueta", "regimen_salud_etiqueta",
  "pertenencia_etnica_etiqueta", "tipo_caso_inicial_etiqueta",
  "hospitalizado_etiqueta", "condicion_final_etiqueta", "ajuste_etiqueta",
  "nom_est_f_caso", "confirmados"
)

for (i in seq_along(rows_per_year)) {
  anio <- as.integer(names(rows_per_year)[[i]])
  cat(sprintf("Adding global flags and definitive profile %d...\n", anio))
  
  datos <- readRDS(file.path(provisional_dir, paste0("dengue_canonica_", anio, ".rds")))
  idx <- (offsets[[i]] + 1L):offsets[[i + 1L]]
  
  datos$flag_duplicado_consecutive             <- id_duplicates$marca[idx]
  datos$grupo_duplicado_consecutive            <- id_duplicates$grupo[idx]
  datos$flag_duplicado_contenido_analitico     <- content_duplicates$marca[idx]
  datos$grupo_duplicado_contenido_analitico    <- content_duplicates$grupo[idx]
  datos$flag_posible_duplicado_clave_estricta  <- strict_duplicates$marca[idx]
  datos$grupo_posible_duplicado_clave_estricta <- strict_duplicates$grupo[idx]
  
  saveRDS(datos,
          file.path(canonical_dir, paste0("dengue_canonica_", anio, ".rds")),
          compress = "gzip")
  
  analytical_vars <- unique(c(
    analytical_vars_base,
    grep("^gp_", names(datos), value = TRUE),
    grep("^dias_", names(datos), value = TRUE),
    grep("^flag_", names(datos), value = TRUE),
    grep("^grupo_duplicado|^grupo_posible", names(datos), value = TRUE)
  ))
  analytical_vars <- intersect(analytical_vars, names(datos))
  
  saveRDS(datos[, analytical_vars, drop = FALSE],
          file.path(analytical_dir, paste0("dengue_analitica_", anio, ".rds")),
          compress = "gzip")
  
  annual_profile[[as.character(anio)]] <- do.call(rbind, lapply(names(datos), function(variable) {
    x <- datos[[variable]]
    vacio <- es_vacio(x)
    data.frame(
      ano = anio, variable = variable,
      tipo_r = paste(class(x), collapse = "/"),
      n_registros = length(x), n_faltantes = sum(vacio),
      n_unicos_anuales = length(unique(x[!vacio])),
      stringsAsFactors = FALSE
    )
  }))
  
  flags <- grep("^flag_", names(datos), value = TRUE)
  annual_quality[[as.character(anio)]] <- data.frame(
    ano = anio, n_registros = nrow(datos),
    as.list(vapply(datos[, flags, drop = FALSE],
                   function(x) sum(x, na.rm = TRUE), numeric(1))),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  annual_distributions[[as.character(anio)]] <- do.call(rbind, lapply(
    intersect(distribution_vars, names(datos)), function(variable) {
      x <- as.character(datos[[variable]])
      x[es_vacio(x)] <- "Sin información"
      tab <- table(x, useNA = "no")
      data.frame(ano = anio, variable = variable,
                 valor = names(tab), n = as.integer(tab))
    }
  ))
  
  candidate_mark <- datos$flag_duplicado_contenido_analitico |
    datos$flag_posible_duplicado_clave_estricta |
    datos$flag_duplicado_consecutive
  
  if (any(candidate_mark)) {
    candidates[[as.character(anio)]] <- datos[
      candidate_mark,
      unique(c(
        provenance_columns, source_columns,
        "grupo_duplicado_consecutive", "grupo_duplicado_contenido_analitico",
        "grupo_posible_duplicado_clave_estricta",
        "flag_duplicado_consecutive", "flag_duplicado_contenido_analitico",
        "flag_posible_duplicado_clave_estricta"
      )), drop = FALSE
    ]
  }
  rm(datos); invisible(gc())
}

## 7. Final reports export
annual_profile <- do.call(rbind, annual_profile)
rownames(annual_profile) <- NULL
write.csv(annual_profile,
          file.path(reports_dir, "perfil_variables_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

annual_quality <- do.call(rbind, annual_quality)
rownames(annual_quality) <- NULL
write.csv(annual_quality,
          file.path(reports_dir, "calidad_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

distributions <- do.call(rbind, annual_distributions)
rownames(distributions) <- NULL
write.csv(distributions,
          file.path(reports_dir, "distribuciones_clave_por_anio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

duplicates_summary <- data.frame(
  control = c(
    "Repeated consecutive identifier",
    "Matching analytical content, without identifiers or technical fields",
    "Match in strict epidemiological key"
  ),
  filas_marcadas = c(
    sum(id_duplicates$marca), sum(content_duplicates$marca),
    sum(strict_duplicates$marca)
  ),
  grupos = c(
    length(unique(na.omit(id_duplicates$grupo))),
    length(unique(na.omit(content_duplicates$grupo))),
    length(unique(na.omit(strict_duplicates$grupo)))
  ),
  decision = c(
    "Review provenance; do not automatically exclude.",
    "Content redundancy; exclude via R/11 and keep audit and original databases.",
    "Could be a coincidence or legitimate update; keep until review."
  ),
  stringsAsFactors = FALSE
)
write.csv(duplicates_summary,
          file.path(reports_dir, "resumen_duplicados.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  data.frame(
    ano = as.integer(names(rows_per_year)), registros = unname(rows_per_year),
    ruta_canonica = file.path(canonical_dir, paste0("dengue_canonica_",  names(rows_per_year), ".rds")),
    ruta_analitica = file.path(analytical_dir, paste0("dengue_analitica_", names(rows_per_year), ".rds")),
    stringsAsFactors = FALSE
  ),
  file.path(data_dir, "dengue_particiones.rds"), compress = "gzip"
)

cat("Preparation, quality, and duplicates phase completed.\n")
cat("Kept records:", sum(rows_per_year), "\n")
print(duplicates_summary, row.names = FALSE)