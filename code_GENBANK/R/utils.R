# Shared helpers used across the pipeline.

#' Standard NULL coalescing.
#' Returns `default` when `x` is NULL or has length zero.
`%||%` <- function(x, default) {
  if (is.null(x) || length(x) == 0L) default else x
}

#' String default: returns `default` for NULL, NA, or blank strings.
str_or <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0L) return(default)
  v <- trimws(as.character(x[[1L]]))
  if (is.na(v) || !nzchar(v)) default else v
}

#' Create a directory (recursive) and return its normalized path.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

#' Path relative to `root`.
rel_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else basename(path)
}

#' Sanitize a value for inclusion in a FASTA header or filename.
clean_header_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return("NA")
  v <- trimws(as.character(x[[1L]]))
  if (is.na(v) || !nzchar(v)) return("NA")
  gsub("[\r\n\t| ]+", "_", v)
}

#' Redact secrets from strings before logging them.
redact_secrets <- function(text, api_key = "", email = "") {
  if (nzchar(api_key)) text <- gsub(api_key, "[API_KEY]", text, fixed = TRUE)
  if (nzchar(email))   text <- gsub(email,   "[EMAIL]",   text, fixed = TRUE)
  text
}

#' Write an Excel workbook if openxlsx is available.
write_excel <- function(path, sheets) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(FALSE)
  wb <- openxlsx::createWorkbook()
  for (sheet in names(sheets)) {
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, sheets[[sheet]])
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  TRUE
}

#' Timestamped logger. Returns a list of closures (info/warn/error/path).
new_logger <- function(path = NULL) {
  state <- new.env(parent = emptyenv())
  state$path <- path
  
  emit <- function(level, fmt, ...) {
    line <- sprintf("%s [%s] %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    level,
                    sprintf(fmt, ...))
    base::message(line)
    if (!is.null(state$path)) {
      cat(line, "\n", sep = "", file = state$path, append = TRUE)
    }
    invisible(line)
  }
  list(
    info  = function(fmt, ...) emit("INFO",  fmt, ...),
    warn  = function(fmt, ...) emit("WARN",  fmt, ...),
    error = function(fmt, ...) emit("ERROR", fmt, ...),
    path  = function() state$path
  )
}

#' Source every R/*.R helper other than utils.R itself.
#' Call this once at the top of each script, after sourcing utils.R.
load_project_helpers <- function(project_root) {
  r_dir <- file.path(project_root, "R")
  if (!dir.exists(r_dir)) stop("Missing R/ directory under ", project_root)
  utils_file <- normalizePath(file.path(r_dir, "utils.R"),
                              winslash = "/", mustWork = FALSE)
  files <- sort(list.files(r_dir, "\\.R$", full.names = TRUE))
  files <- files[normalizePath(files, winslash = "/", mustWork = FALSE) != utils_file]
  for (f in files) source(f, local = FALSE)
  invisible(project_root)
}

#' Determine the project root from the running script or cwd.
#' Precedence: $DENGUE_SURAMERICA_ROOT, --file= argument, cwd with DESCRIPTION.
detect_project_root <- function() {
  env_root <- Sys.getenv("DENGUE_SURAMERICA_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  file_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_flag) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_flag[[1L]]),
                                 winslash = "/", mustWork = FALSE)
    return(dirname(dirname(script_path)))
  }
  probe <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in seq_len(10L)) {
    if (file.exists(file.path(probe, "DESCRIPTION"))) return(probe)
    parent <- dirname(probe)
    if (identical(parent, probe)) break
    probe <- parent
  }
  stop("Cannot locate project root. Set DENGUE_SURAMERICA_ROOT.")
}