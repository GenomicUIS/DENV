# 01_load_databases.R
# Import and unification of SIVIGILA databases (Dengue, Severe Dengue, Mortality)

source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(janitor)
})

## -----------------------------------------------------------------------------
## Rename helper: does nothing if the source name is absent; drops the source
## if the target name already exists, preventing duplicate column names.
## -----------------------------------------------------------------------------
renombrar_si_no_existe <- function(datos, viejo, nuevo) {
  if (!viejo %in% names(datos)) return(datos)
  if (nuevo %in% names(datos)) {
    datos[, !names(datos) %in% viejo, drop = FALSE]
  } else {
    names(datos)[names(datos) == viejo] <- nuevo
    datos
  }
}

## -----------------------------------------------------------------------------
## Main import function
## -----------------------------------------------------------------------------
cargar_evento <- function(evento_buscado) {
  evento_buscado <- validar_evento(evento_buscado)
  
  archivos_evento <- inventario_sivigila %>% filter(evento == evento_buscado)
  
  message("\nStarting the loading of the ", nrow(archivos_evento),
          " files for: ", evento_buscado)
  
  datos_unificados <- purrr::pmap_dfr(
    list(
      ruta    = archivos_evento$ruta,
      archivo = archivos_evento$archivo,
      anio    = archivos_evento$anio_archivo
    ),
    function(ruta, archivo, anio) {
      message("  -> Importing: ", archivo)
      
      # Everything as text: prevents type conflicts across years.
      datos <- readxl::read_excel(ruta, col_types = "text", .name_repair = "minimal")
      datos <- janitor::clean_names(datos)
      
      # Historical SIVIGILA corrections, safe against duplicate targets.
      datos <- renombrar_si_no_existe(datos, "consecutive2",   "consecutive")
      datos <- renombrar_si_no_existe(datos, "consecutive_12", "consecutive")
      datos <- renombrar_si_no_existe(datos, "nombre_upgd",    "nom_upgd")
      
      # Drop any remaining duplicates, keeping the first occurrence.
      if (anyDuplicated(names(datos))) {
        datos <- datos[, !duplicated(names(datos)), drop = FALSE]
      }
      
      datos %>%
        mutate(
          fuente_evento     = evento_buscado,
          fuente_archivo    = archivo,
          fuente_anio       = anio,
          fuente_fila_excel = row_number() + 1L   # +1 = header row in Excel
        ) %>%
        relocate(starts_with("fuente_"), .before = 1)
    }
  )
  
  message("Loading finished for ", evento_buscado,
          ". Total reports read: ",
          format(nrow(datos_unificados), big.mark = ","))
  
  return(datos_unificados)
}

## -----------------------------------------------------------------------------
## Optional: load every event into a named list
## -----------------------------------------------------------------------------
cargar_todas_las_bases <- function() {
  resultado <- purrr::map(names(SIVIGILA_DIRS), cargar_evento)
  names(resultado) <- names(SIVIGILA_DIRS)
  return(resultado)
}