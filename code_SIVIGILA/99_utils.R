# 99_utils.R
# Canonical helpers for the dengue pipeline. No project dependencies.
# Automatically loaded from 00_config.R at the end.

## -----------------------------------------------------------------------------
## 1. Missing values and sentinels detection
## -----------------------------------------------------------------------------
# Textual sentinels that SIVIGILA uses as "no data" and are not real NAs.
.SENTINELAS_VACIO <- c(
  "sin dato", "sin datos", "sin informacion", "sin información",
  "no aplica", "no informa", "no reportado", "no reportada",
  "desconocido", "desconocida", "se desconoce", "sin definir",
  "na", "n/a", "null", "nan", "-", "--", "---", "*", "**"
)

es_vacio <- function(x, sentinelas = .SENTINELAS_VACIO) {
  if (inherits(x, "Date") || is.numeric(x) || is.integer(x) || is.logical(x)) {
    return(is.na(x))
  }
  text <- trimws(as.character(x))
  text_norm <- tolower(iconv(text, from = "", to = "ASCII//TRANSLIT", sub = ""))
  is_na <- is.na(text) | !nzchar(text)
  is_sentinel <- !is_na & text_norm %in% sentinelas
  # Also: strings starting with "*" are usually SIVIGILA comments
  is_comment <- !is_na & grepl("^\\*", text)
  is_na | is_sentinel | is_comment
}

## -----------------------------------------------------------------------------
## 2. Robust date parsing (unifies 02, 03, 07, 10)
## -----------------------------------------------------------------------------
# Accepts:
#   - ISO Dates:              2007-03-20
#   - Latin American dates:   30/03/2025  /  3/06/2025
#   - Excel serials:          44197
#   - Dates with attached time: 30/03/2025 12:00:00 a. m.
# The 2025 files deliver the date with the time marker "a. m.", and without it,
# the full string doesn't match any lubridate order, which is why the time 
# suffix is discarded before parsing.
parsear_fecha <- function(x) {
  text <- trimws(as.character(x))
  text[es_vacio(text)] <- NA_character_
  
  # Remove the time part when it is attached to the date.
  text <- sub("\\s+\\d{1,2}:\\d{2}(:\\d{2})?.*$", "", text)
  text <- trimws(text)
  text[!nzchar(text)] <- NA_character_
  
  output <- as.Date(rep(NA_character_, length(text)))
  
  # Excel serials (e.g. 44197)
  number <- suppressWarnings(as.numeric(gsub(",", ".", text, fixed = TRUE)))
  is_serial <- !is.na(number) & number >= 1 & number <= 100000
  output[is_serial] <- as.Date(number[is_serial], origin = "1899-12-30")
  
  # ISO and textual variants via lubridate if available
  if (requireNamespace("lubridate", quietly = TRUE)) {
    remaining <- which(is.na(output) & !is.na(text))
    if (length(remaining)) {
      parsed <- suppressWarnings(lubridate::parse_date_time(
        text[remaining],
        orders = c("Ymd", "Ymd HMS", "Ymd HM", "dmY", "mdY",
                   "Y/m/d", "d/m/Y", "d-m-Y", "Y-m-d"),
        quiet = TRUE
      ))
      parsed <- as.Date(parsed)
      output[remaining[!is.na(parsed)]] <- parsed[!is.na(parsed)]
    }
  }
  output
}

## -----------------------------------------------------------------------------
## 3. Code normalization (DIVIPOLA, providers, etc.)
## -----------------------------------------------------------------------------
normalizar_codigo <- function(x, ancho) {
  text <- trimws(as.character(x))
  text[is.na(text) | !nzchar(text)] <- NA_character_
  indices <- which(!is.na(text) & grepl("^[0-9]+(?:\\.0+)?$", text))
  digits <- sub("\\.0+$", "", text[indices])
  digits <- sub("^0+(?=[0-9])", "", digits, perl = TRUE)
  zeros <- vapply(pmax(0L, ancho - nchar(digits)),
                  function(n) paste0(rep("0", n), collapse = ""), character(1))
  text[indices] <- paste0(zeros, digits)
  text
}

# Unifies notification DIVIPOLA where COD_MUN_N can come as:
#   - 5 digits (Full DIVIPOLA, 2007–2010 cases)
#   - 3 digits (municipality only; joined with COD_DPTO_N)
#   - 4 digits (rare; treated as full)
#   - "000" (sentinel: unknown municipality → NA)
unir_divipola_notificacion <- function(cod_dpto, cod_mun) {
  dpto <- normalizar_codigo(cod_dpto, 2L)
  original <- trimws(as.character(cod_mun))
  original[es_vacio(original)] <- NA_character_
  digits <- sub("\\.0+$", "", original)
  
  is_numeric <- !is.na(digits) & grepl("^[0-9]+$", digits)
  is_pure_zero <- is_numeric & grepl("^0+$", digits)
  
  is_full <- is_numeric & nchar(digits) %in% 4:5 & !is_pure_zero
  output <- rep(NA_character_, length(original))
  output[is_full] <- sprintf("%05d", as.integer(digits[is_full]))
  
  is_municipality <- is_numeric & nchar(digits) <= 3L & !is_pure_zero & !is.na(dpto)
  output[is_municipality] <- paste0(
    dpto[is_municipality], sprintf("%03d", as.integer(digits[is_municipality]))
  )
  output
}

# Old variant (03) that always assumes dpto=2 and mun=3
unir_divipola_simple <- function(dpto, mun) {
  d <- normalizar_codigo(dpto, 2L)
  m <- normalizar_codigo(mun, 3L)
  ifelse(is.na(d) | is.na(m), NA_character_, paste0(d, m))
}

## -----------------------------------------------------------------------------
## 4. Code labeling (numeric map → text)
## -----------------------------------------------------------------------------
etiquetar_codigo <- function(x, mapa) {
  output <- unname(mapa[as.character(x)])
  unknown <- is.na(output) & !is.na(x) & !es_vacio(x)
  output[unknown] <- paste0("Code ", x[unknown])
  output
}

## -----------------------------------------------------------------------------
## 5. Deduplication: digital footprint and group marking
## -----------------------------------------------------------------------------
crear_clave <- function(datos, columnas) {
  columns <- intersect(columnas, names(datos))
  values <- lapply(datos[, columns, drop = FALSE], function(x) {
    text <- trimws(as.character(x))
    text[es_vacio(text)] <- "<MISSING>"
    text
  })
  do.call(paste, c(values, sep = "\r"))
}

marcar_grupos_repetidos <- function(clave) {
  mark <- duplicated(clave) | duplicated(clave, fromLast = TRUE)
  group <- rep(NA_integer_, length(clave))
  if (any(mark)) group[mark] <- match(clave[mark], unique(clave[mark]))
  list(marca = mark, grupo = group)
}

## -----------------------------------------------------------------------------
## 6. SIVIGILA column names cleaning
## -----------------------------------------------------------------------------
limpiar_nombres_sivigila <- function(nombres) {
  x <- tolower(trimws(nombres))
  x <- chartr("áéíóúñÁÉÍÓÚÑ", "aeiounAEIOUN", x)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

## -----------------------------------------------------------------------------
## 7. Age in years (unifies 07, 10 and 03)
## -----------------------------------------------------------------------------
calcular_edad_anios <- function(edad, uni_med) {
  edad <- suppressWarnings(as.numeric(edad))
  uni_med <- suppressWarnings(as.integer(uni_med))
  out <- rep(NA_real_, length(edad))
  out[uni_med == 1L] <- edad[uni_med == 1L]
  out[uni_med == 2L] <- edad[uni_med == 2L] / 12
  out[uni_med == 3L] <- edad[uni_med == 3L] / 365.25
  out[uni_med == 4L] <- edad[uni_med == 4L] / (365.25 * 24)
  out[uni_med == 5L] <- edad[uni_med == 5L] / (365.25 * 24 * 60)
  out
}

grupo_edad_tesis <- function(edad_anios) {
  cut(
    edad_anios,
    breaks = c(-Inf, 1, 5, 15, 25, 45, 65, Inf),
    labels = c("<1", "1-4", "5-14", "15-24", "25-44", "45-64", "65+"),
    right = FALSE
  )
}

## -----------------------------------------------------------------------------
## 8. Safe figure overwriting
## -----------------------------------------------------------------------------
guardar_manifesto <- function(manifesto, ruta) {
  if (length(manifesto)) {
    write.csv(do.call(rbind, manifesto), ruta, row.names = FALSE,
              fileEncoding = "UTF-8")
  }
}