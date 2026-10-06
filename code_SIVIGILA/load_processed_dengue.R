## -----------------------------------------------------------------------------
## load_processed_dengue.R
##
## INTERACTIVE UTILITY — NOT PART OF THE PIPELINE
##
## Defines cargar_dengue_procesado(), a helper to load into the console the
## processed annual dengue partitions without writing long paths.
##
## Typical usage:
##   source("00_config.R")                    # if not already loaded
##   source("load_processed_dengue.R")
##   dengue_2023_2025 <- cargar_dengue_procesado(2023:2025)
##   dengue_todos     <- cargar_dengue_procesado(SIVIGILA_ANIOS, combinar = FALSE)
##
## It is not sourced from 00_config.R to avoid a circular dependency,
## because this script depends on 00_config.R.
## -----------------------------------------------------------------------------

source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## LOAD FUNCTION FOR THE HISTORICAL DENGUE DATABASE (EVENT 210)
## Purpose: Safely load processed dengue annual partitions.
## Allows loading specific years or ranges, aligning columns automatically
## and preventing RAM overflows.
## -----------------------------------------------------------------------------

cargar_dengue_procesado <- function(
    anios = max(SIVIGILA_ANIOS),
    tipo = c("analitica", "canonica"),
    sin_duplicados = FALSE,  ## Adjusted to FALSE by default temporarily
    combinar = TRUE,
    verbose = TRUE) {
  
  tipo <- match.arg(tipo)
  anios <- sort(unique(as.integer(anios)))
  
  if (any(is.na(anios)) || any(!anios %in% SIVIGILA_ANIOS)) {
    stop("Years must be between 2007 and 2025.", call. = FALSE)
  }
  
  if (isTRUE(sin_duplicados)) {
    subfolder <- if (tipo == "analitica") {
      "analitica_sin_duplicados_por_anio"
    } else {
      "canonica_sin_duplicados_por_anio"
    }
    prefix <- if (tipo == "analitica") {
      "dengue_analitica_sin_duplicados_"
    } else {
      "dengue_canonica_sin_duplicados_"
    }
  } else {
    subfolder <- if (tipo == "analitica") "analitica_por_anio" else "canonica_por_anio"
    prefix <- if (tipo == "analitica") "dengue_analitica_" else "dengue_canonica_"
  }
  
  paths <- file.path(
    SIVIGILA_PATH, "resultados", "dengue", "datos_procesados", subfolder,
    paste0(prefix, anios, ".rds")
  )
  
  if (any(!file.exists(paths))) {
    stop("Missing partitions: ", paste(basename(paths[!file.exists(paths)]), collapse = ", "))
  }
  
  if (length(anios) > 6L && isTRUE(combinar)) {
    warning(
      "Will combine ", length(anios), " years in memory. ",
      "Use combinar = FALSE if you only need to iterate over them.", call. = FALSE
    )
  }
  
  tables <- lapply(seq_along(paths), function(i) {
    if (isTRUE(verbose)) message("Loading DENGUE ", anios[[i]], "...")
    readRDS(paths[[i]])
  })
  names(tables) <- as.character(anios)
  
  if (!isTRUE(combinar)) return(tables)
  
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(x) {
    missing_cols <- setdiff(columns, names(x))
    for (variable in missing_cols) x[[variable]] <- NA
    x[, columns, drop = FALSE]
  })
  
  output <- do.call(rbind, tables)
  rownames(output) <- NULL
  output
}

## -----------------------------------------------------------------------------
## Usage examples:
## dengue_2025 <- cargar_dengue_procesado(2025)
## dengue_2020_2025 <- cargar_dengue_procesado(2020:2025)
## dengue_por_anio <- cargar_dengue_procesado(SIVIGILA_ANIOS, combinar = FALSE)
## -----------------------------------------------------------------------------
