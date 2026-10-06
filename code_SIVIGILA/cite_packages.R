# Export the citation metadata from the environment used for the analysis.
# This is an optional documentation utility, not a pipeline stage.
# Run it in the SAME R installation and library configuration as the pipeline.
# No packages are installed or updated by this script.

paquetes <- c(
  "dplyr", "stringr", "fs", "purrr", "readxl", "janitor", "tidyr",
  "ggplot2", "scales", "patchwork", "sf", "lwgeom", "lubridate",
  "readr", "tibble"
)

destino <- file.path(getwd(), "referencias_software")
dir.create(destino, recursive = TRUE, showWarnings = FALSE)
instalados <- vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)

versiones <- data.frame(
  paquete = c("R", paquetes),
  version = c(
    as.character(getRversion()),
    vapply(paquetes, function(p) {
      if (instalados[[p]]) as.character(utils::packageVersion(p)) else NA_character_
    }, character(1))
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(versiones, file.path(destino, "VERSIONES_SOFTWARE.csv"),
                 row.names = FALSE, na = "", fileEncoding = "UTF-8")

referencias_txt <- character()
referencias_bib <- character()
for (p in c("R", paquetes)) {
  disponible <- identical(p, "R") || isTRUE(instalados[[p]])
  if (!disponible) {
    referencias_txt <- c(referencias_txt, paste0("## ", p),
                          "Not installed: no local citation could be retrieved.", "")
    next
  }
  cita <- if (identical(p, "R")) utils::citation() else utils::citation(p)
  referencias_txt <- c(referencias_txt, paste0("## ", p),
                        capture.output(print(cita)), "")
  referencias_bib <- c(referencias_bib, paste0("% ", p),
                        as.character(utils::toBibtex(cita)), "")
}
writeLines(enc2utf8(referencias_txt),
           file.path(destino, "CITAS_RECOMENDADAS.txt"), useBytes = TRUE)
writeLines(enc2utf8(referencias_bib),
           file.path(destino, "REFERENCIAS_SOFTWARE.bib"), useBytes = TRUE)
writeLines(capture.output(utils::sessionInfo()),
           file.path(destino, "SESSION_INFO.txt"))
if (isTRUE(instalados[["sf"]])) {
  writeLines(capture.output(sf::sf_extSoftVersion()),
             file.path(destino, "VERSIONES_BIBLIOTECAS_ESPACIALES.txt"))
}
message("Citation metadata exported to: ", normalizePath(destino))
if (any(!instalados)) {
  message("Missing packages: ", paste(paquetes[!instalados], collapse = ", "))
}
