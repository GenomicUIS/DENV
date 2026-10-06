# 02_initial_diagnosis.R
# Initial audit of the SIVIGILA databases: rows per year, duplicates,
# variable completeness, and date consistency.
# Produces salidas/diagnostico_inicial/01_resumen_eventos_y_duplicados.csv,
# which in turn feeds filas_originales_esperadas() in 00_config.R.

source("00_config.R", encoding = "UTF-8")
source("01_load_databases.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(fs)
})

## 1. Output folder
dir_salidas <- fs::path(SIVIGILA_PATH, "salidas", "diagnostico_inicial")
fs::dir_create(dir_salidas)

## 2. Audit functions
auditar_evento <- function(datos, nombre_evento) {
  message("Auditing: ", nombre_evento)
  
  columnas_contenido <- setdiff(names(datos),
                                grep("^fuente_", names(datos), value = TRUE))
  
  # A. Summary by year.
  # filas_vacias uses es_vacio() so that empty strings and SIVIGILA
  # sentinels ("SIN DATO", "DESCONOCIDO", "*", ...) are counted.
  resumen_anual <- datos %>%
    group_by(fuente_anio) %>%
    summarise(
      n_filas            = n(),
      filas_total        = n(),
      archivos_leidos    = n_distinct(fuente_archivo),
      errores_lectura    = 0L,
      filas_vacias       = sum(if_all(all_of(columnas_contenido), ~ es_vacio(.x))),
      duplicados_exactos = sum(duplicated(pick(all_of(columnas_contenido)))),
      id_duplicados      = if ("consecutive" %in% names(datos)) {
        sum(duplicated(datos$consecutive, incomparables = NA))
      } else {
        NA_integer_
      },
      .groups = "drop"
    ) %>%
    mutate(evento = nombre_evento, .before = 1)
  
  # B. Completeness by variable.
  perfil_variables <- datos %>%
    pivot_longer(
      cols = -starts_with("fuente_"),
      names_to = "variable",
      values_to = "valor",
      values_transform = as.character
    ) %>%
    group_by(fuente_anio, variable) %>%
    summarise(
      n_total       = n(),
      n_faltantes   = sum(es_vacio(valor)),
      pct_faltantes = round(100 * mean(es_vacio(valor)), 2),
      .groups = "drop"
    ) %>%
    mutate(evento = nombre_evento)
  
  # C. Date audit.
  columnas_fecha <- names(datos)[grepl("fec_|fecha|ini_sin",
                                       names(datos), ignore.case = TRUE)]
  
  perfil_fechas <- tibble()
  if (length(columnas_fecha) > 0) {
    perfil_fechas <- datos %>%
      select(fuente_anio, all_of(columnas_fecha)) %>%
      pivot_longer(-fuente_anio, names_to = "variable_fecha",
                   values_to = "valor_original") %>%
      filter(!es_vacio(valor_original)) %>%
      mutate(
        fecha_limpia      = parsear_fecha(valor_original),
        es_fecha_invalida = is.na(fecha_limpia),
        fuera_rango       = !is.na(fecha_limpia) &
          (fecha_limpia < as.Date("1900-01-01") |
             fecha_limpia > Sys.Date() + 366)
      ) %>%
      group_by(fuente_anio, variable_fecha) %>%
      summarise(
        n_registros_no_vacios = n(),
        n_no_parseables       = sum(es_fecha_invalida),
        n_fuera_rango         = sum(fuera_rango, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(evento = nombre_evento)
  }
  
  list(resumen   = resumen_anual,
       variables = perfil_variables,
       fechas    = perfil_fechas)
}

## 3. Execution
# Warning: this loads dengue, severe dengue, and mortality into RAM
# simultaneously. For >1.5 M dengue rows it can consume several GB.
message("Loading all databases. This may take a few minutes...")
bases_completas <- cargar_todas_las_bases()

resultados_auditoria <- purrr::map2(bases_completas,
                                    names(bases_completas),
                                    auditar_evento)

tabla_resumen   <- purrr::map_dfr(resultados_auditoria, "resumen")
tabla_variables <- purrr::map_dfr(resultados_auditoria, "variables")
tabla_fechas    <- purrr::map_dfr(resultados_auditoria, "fechas")

## 4. Export
readr::write_csv(tabla_resumen,
                 fs::path(dir_salidas, "01_resumen_eventos_y_duplicados.csv"))
readr::write_csv(tabla_variables,
                 fs::path(dir_salidas, "02_perfil_completitud_variables.csv"))
readr::write_csv(tabla_fechas,
                 fs::path(dir_salidas, "03_auditoria_fechas.csv"))

message("\nInitial diagnosis successfully finished.")
message("Audit files in: ", dir_salidas)