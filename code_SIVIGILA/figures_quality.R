# Spanish numerical conventions, common to all three collections.
# Output functions translated to English equivalents to match plot text.
scales::number_options(decimal.mark = ",", big.mark = ".")
porcentaje_es <- function(x, ...) scales::percent(x, ..., decimal.mark = ",", big.mark = ".", suffix = " %")
porcentaje_formato_es <- function(...) scales::label_percent(..., decimal.mark = ",", big.mark = ".", suffix = " %")
escala_miles_es <- function() setNames(c(1, 1000, 1000000), c("", " thousand", " million"))

guardar_evidencia_figura <- function(grafico, ruta, ancho, alto) {
  # The PDF version preserves text and vector strokes; TIFF continues at 600 dpi.
  ggplot2::ggsave(sub("\\.tiff$", ".pdf", ruta), grafico,
                  device = grDevices::cairo_pdf, width = ancho, height = alto, units = "in",
                  bg = "white", limitsize = FALSE)
  
  event_name <- basename(dirname(dirname(dirname(ruta))))
  dest_dir <- file.path(SIVIGILA_PATH, "salidas", "auditoria_integral", "objetos_figuras", event_name)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(grafico, file.path(dest_dir, sub("\\.tiff$", ".rds", basename(ruta))))
}