# 00_config.R
# Initial configuration, paths, and validation of the SIVIGILA environment.

## -----------------------------------------------------------------------------
## 0. Verify working directory
## -----------------------------------------------------------------------------
if (!file.exists("01_load_databases.R")) {
  stop(
    "Execute this script from the project's R/ folder.\n",
    "Current directory: ", getwd(), "\n",
    "Use setwd() or open the corresponding .Rproj project.",
    call. = FALSE
  )
}

options(readxl.show_progress = FALSE)

## -----------------------------------------------------------------------------
## 1. Load common packages
## -----------------------------------------------------------------------------
.paquetes_requeridos <- c(
  "dplyr", "stringr", "fs", "purrr", "readxl", "janitor",
  "tidyr", "ggplot2", "scales", "patchwork", "sf", "lubridate",
  "readr"
)

.faltantes <- .paquetes_requeridos[!vapply(
  .paquetes_requeridos, requireNamespace, logical(1), quietly = TRUE
)]
if (length(.faltantes)) {
  stop(
    "Missing packages: ", paste(.faltantes, collapse = ", "),
    ". Install them with install.packages().",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(fs)
  library(purrr)
})

## -----------------------------------------------------------------------------
## 2. Root path and local library
## -----------------------------------------------------------------------------
# By default, the R/ folder is inside the project root; the event folders
# (DENGUE, DENGUE GRAVE, ...) live one level up. If SIVIGILA_PATH is set in
# the environment, that value wins.
resolver_raiz_sivigila <- function() {
  candidata <- Sys.getenv("SIVIGILA_PATH", unset = NA_character_)
  if (is.na(candidata) || !nzchar(candidata)) {
    if (file.exists("01_load_databases.R")) {
      candidata <- dirname(getwd())      # parent of R/
    } else {
      candidata <- getwd()
    }
  }
  candidata <- fs::path_norm(candidata)
  if (!fs::dir_exists(candidata)) {
    stop("The SIVIGILA_PATH path does not exist: ", candidata, call. = FALSE)
  }
  candidata
}

SIVIGILA_PATH <- resolver_raiz_sivigila()

version_r <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
SIVIGILA_R_LIB <- fs::path(SIVIGILA_PATH, "R_library", version_r)
if (fs::dir_exists(SIVIGILA_R_LIB)) {
  .libPaths(c(SIVIGILA_R_LIB, .libPaths()))
}

## -----------------------------------------------------------------------------
## 3. Event folders
## -----------------------------------------------------------------------------
SIVIGILA_CARPETAS <- c(
  dengue            = "DENGUE",
  dengue_grave      = "DENGUE GRAVE",
  mortalidad_dengue = "MORTALIDAD DENGUE"
)

SIVIGILA_DIRS <- fs::path(SIVIGILA_PATH, SIVIGILA_CARPETAS)
names(SIVIGILA_DIRS) <- names(SIVIGILA_CARPETAS)

.faltantes <- names(SIVIGILA_DIRS)[!fs::dir_exists(SIVIGILA_DIRS)]
if (length(.faltantes)) {
  stop("Missing expected folders: ", paste(.faltantes, collapse = ", "),
       call. = FALSE)
}

SIVIGILA_ANIOS <- 2007:2025

## -----------------------------------------------------------------------------
## 4. File inventory
## -----------------------------------------------------------------------------
inventariar_sivigila <- function() {
  purrr::imap_dfr(SIVIGILA_DIRS, function(ruta_dir, evento) {
    archivos <- fs::dir_ls(ruta_dir, regexp = "\\.(xls|xlsx)$", ignore.case = TRUE)
    info <- fs::file_info(archivos)
    tibble::tibble(
      evento        = evento,
      anio_archivo  = as.integer(str_extract(fs::path_file(archivos), "(19|20)\\d{2}")),
      archivo       = fs::path_file(archivos),
      extension     = fs::path_ext(archivos),
      tamano_bytes  = info$size,
      modificado    = info$modification_time,
      ruta          = as.character(archivos)
    )
  }) %>% arrange(anio_archivo, archivo)
}

inventario_sivigila <- inventariar_sivigila()

## -----------------------------------------------------------------------------
## 5. Validations
## -----------------------------------------------------------------------------
validar_evento <- function(evento) {
  if (length(evento) != 1 || !evento %in% names(SIVIGILA_DIRS)) {
    stop("Invalid event. Use one of: ",
         paste(names(SIVIGILA_DIRS), collapse = ", "), call. = FALSE)
  }
  evento
}

purrr::walk(names(SIVIGILA_DIRS), function(evento) {
  anios_presentes <- inventario_sivigila %>%
    filter(evento == !!evento) %>%
    pull(anio_archivo)
  if (!identical(sort(anios_presentes), SIVIGILA_ANIOS)) {
    stop("Incomplete inventory or out of period for ", evento,
         ". One file per year from 2007 to 2025 is expected.", call. = FALSE)
  }
})

## -----------------------------------------------------------------------------
## 6. Reading external reports
## -----------------------------------------------------------------------------
INDICADORES_DEDUP <- c(
  base_canonica         = "Records in canonical base",
  grupos_redundantes    = "Groups with matching analytical content",
  pertenecen_grupos     = "Records belonging to those groups",
  excluidos             = "Redundant records excluded from analyses",
  base_deduplicada      = "Records in deduplicated base",
  representantes        = "Kept group representatives",
  grupos_clave_estricta = "Strict key groups with kept differences"
)

filas_originales_esperadas <- function(evento) {
  ruta <- fs::path(SIVIGILA_PATH, "salidas", "diagnostico_inicial",
                   "01_resumen_eventos_y_duplicados.csv")
  if (!fs::file_exists(ruta)) {
    stop("First execute R/02_initial_diagnosis.R", call. = FALSE)
  }
  x <- read.csv(ruta, fileEncoding = "UTF-8")
  col_filas <- intersect(c("filas_total", "n_filas"), names(x))[1]
  if (is.na(col_filas)) {
    stop("The diagnostic report does not contain a total rows column.",
         call. = FALSE)
  }
  x_ev <- x[x$evento == evento, , drop = FALSE]
  if (nrow(x_ev) == 0L) {
    stop("The diagnostic report does not contain rows for ", evento,
         call. = FALSE)
  }
  stopifnot(all(x_ev$errores_lectura == 0))
  sum(as.numeric(x_ev[[col_filas]]))
}

valor_deduplicacion <- function(evento, indicador_buscado) {
  ruta <- fs::path(SIVIGILA_PATH, "resultados", evento, "reportes",
                   "deduplicacion", "resumen_deduplicacion.csv")
  if (!fs::file_exists(ruta)) {
    stop("The deduplication summary for ", evento, " does not exist.",
         call. = FALSE)
  }
  x <- read.csv(ruta, fileEncoding = "UTF-8")
  x_ev <- x[x$indicador == indicador_buscado, , drop = FALSE]
  if (nrow(x_ev) != 1L) {
    stop(
      "Expected exactly 1 row for indicator '", indicador_buscado,
      "' in ", basename(ruta), ", found ", nrow(x_ev), ".\n",
      "Available indicators: ",
      paste(unique(x$indicador), collapse = " | "),
      call. = FALSE
    )
  }
  if (is.na(x_ev$valor)) {
    stop("Indicator '", indicador_buscado, "' has an NA value.", call. = FALSE)
  }
  as.numeric(x_ev$valor)
}

## -----------------------------------------------------------------------------
## 7. Unique record key
## -----------------------------------------------------------------------------
clave_registro <- function(datos) {
  campos <- c("fuente_anio", "fuente_archivo", "fuente_fila_excel")
  stopifnot(all(campos %in% names(datos)), !anyNA(datos[, campos]))
  paste(datos$fuente_anio, datos$fuente_archivo, datos$fuente_fila_excel,
        sep = "::")
}

## -----------------------------------------------------------------------------
## 8. Load helpers (99_utils.R)
## -----------------------------------------------------------------------------
if (!file.exists("99_utils.R")) {
  stop("99_utils.R is not found in the R/ folder.", call. = FALSE)
}
source("99_utils.R", encoding = "UTF-8")