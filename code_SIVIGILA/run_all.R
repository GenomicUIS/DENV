# run_all.R
# Sequential execution of the complete dengue pipeline.
# Usage:
#   source("run_all.R")                 # defines correr_pipeline()
#   correr_pipeline()                   # executes everything
#   correr_pipeline(hasta = "10")       # executes up to script 10 inclusive
#   correr_pipeline(desde = "13")       # executes from script 13 onwards

# Directory of this script (robust to previous setwd)
.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
if (!file.exists(file.path(.script_dir, "00_config.R"))) {
  .script_dir <- getwd()
}
setwd(.script_dir)

## -----------------------------------------------------------------------------
## Canonical execution order
## -----------------------------------------------------------------------------
.SCRIPTS_ORDEN <- c(
  "00_config.R",
  "01_load_databases.R",
  "02_initial_diagnosis.R",
  # --- Dengue mortality ---
  "03_prepare_dengue_mortality.R",
  "04_deduplicate_dengue_mortality.R",
  "05_generate_dengue_mortality_graphics.R",
  # --- Severe Dengue ---
  "06_prepare_severe_dengue.R",
  "07_deduplicate_severe_dengue.R",
  "08_generate_severe_dengue_graphics.R",
  # --- Dengue ---
  "09_prepare_dengue.R",
  "10_document_dengue.R",
  "11_deduplicate_dengue.R",
  "12_generate_dengue_graphics.R",
  # --- Integration ---
  "13_integrate_dengue_events.R",
  "14_integrated_care_figure.R"
)

## -----------------------------------------------------------------------------
## Main function
## -----------------------------------------------------------------------------
correr_pipeline <- function(desde = NULL, hasta = NULL, detener_en_error = TRUE) {
  scripts <- .SCRIPTS_ORDEN
  
  if (!is.null(desde)) {
    idx <- which(startsWith(scripts, desde))
    if (length(idx) == 0) stop("No script found starting with: ", desde)
    scripts <- scripts[idx[1]:length(scripts)]
  }
  if (!is.null(hasta)) {
    idx <- which(startsWith(scripts, hasta))
    if (length(idx) == 0) stop("No script found starting with: ", hasta)
    scripts <- scripts[1:idx[1]]
  }
  
  log_dir <- file.path(getwd(), "..", "salidas", "log_run_all")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  
  registro <- data.frame(
    script  = character(), inicio = character(), fin = character(),
    estado  = character(), mensaje = character(),
    stringsAsFactors = FALSE
  )
  
  loguear <- function(...) {
    line <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(...))
    cat(line, "\n")
    cat(line, "\n", file = log_file, append = TRUE)
  }
  
  loguear("=== START RUN_ALL ===")
  loguear("Directory: ", getwd())
  loguear("Total scripts: ", length(scripts))
  
  # Single shared environment for the whole run. Scripts are still isolated
  # from the user's global workspace, but they can see each other's objects
  # (SIVIGILA_PATH, helper functions, inventories, ...).
  pipeline_env <- new.env(parent = globalenv())
  
  for (s in scripts) {
    loguear("\n", strrep("=", 70))
    loguear(">>> ", s)
    loguear(strrep("=", 70))
    
    start_time <- Sys.time()
    status     <- "OK"
    error_msg  <- ""
    
    tryCatch(
      {
        source(s, encoding = "UTF-8", echo = FALSE, local = pipeline_env)
      },
      error = function(e) {
        status    <<- "ERROR"
        error_msg <<- conditionMessage(e)
        loguear("!!! ERROR in ", s, ": ", error_msg)
        if (detener_en_error) stop(e)
      },
      warning = function(w) {
        loguear("WARNING: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    loguear(sprintf("<<< %s | %s | %.1f s", s, status, duration))
    
    registro <- rbind(registro, data.frame(
      script = s, inicio = format(start_time), fin = format(end_time),
      estado = status, mensaje = error_msg, stringsAsFactors = FALSE
    ))
  }
  
  write.csv(registro, file.path(log_dir, "last_run.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  
  loguear("\n=== END RUN_ALL ===")
  loguear("Scripts OK: ", sum(registro$estado == "OK"), " / ", nrow(registro))
  invisible(registro)
}

# Automatic execution is disabled to prevent accidental runs.
# Uncomment the next line if you want `source("run_all.R")` to run everything.
# if (interactive() || !exists(".RUN_ALL_NO_AUTO", inherits = FALSE)) {
#   correr_pipeline()
# }