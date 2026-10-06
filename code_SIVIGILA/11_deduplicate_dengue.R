source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## CONSERVATIVE DEDUPLICATION OF DENGUE (EVENT 210)
## Purpose: Eliminate exact analytical redundancy while keeping the most 
## recent version of each duplicated record.
## -----------------------------------------------------------------------------

results_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue")
canonical_dir <- file.path(results_dir, "datos_procesados", "canonica_por_anio")
analytical_dir <- file.path(results_dir, "datos_procesados", "analitica_por_anio")
canonical_sd_dir <- file.path(results_dir, "datos_procesados", "canonica_sin_duplicados_por_anio")
analytical_sd_dir <- file.path(results_dir, "datos_procesados", "analitica_sin_duplicados_por_anio")
audit_dir <- file.path(results_dir, "reportes", "deduplicacion")

invisible(lapply(
  c(canonical_sd_dir, analytical_sd_dir, audit_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

years <- SIVIGILA_ANIOS
canonical_paths <- file.path(canonical_dir, sprintf("dengue_canonica_%d.rds", years))
analytical_paths <- file.path(analytical_dir, sprintf("dengue_analitica_%d.rds", years))

if (any(!file.exists(canonical_paths)) || any(!file.exists(analytical_paths))) {
  stop("Missing annual canonical or analytical partitions for dengue.")
}

## 1. Extract records marked as duplicates directly from the RDS files
cat("Extracting duplicates metadata from annual partitions...\n")
exact_matches_list <- lapply(canonical_paths, function(ruta) {
  d <- readRDS(ruta)
  d[d$flag_duplicado_contenido_analitico & !is.na(d$flag_duplicado_contenido_analitico), ]
})
exact_matches <- do.call(rbind, exact_matches_list)
rownames(exact_matches) <- NULL

exact_matches$fec_aju_orden <- as.Date(exact_matches$fec_aju)
exact_matches$fec_arc_xl_orden <- as.Date(exact_matches$fec_arc_xl)

## 2. Group duplicates and select the most recent representative
groups <- split(seq_len(nrow(exact_matches)), exact_matches$grupo_duplicado_contenido_analitico)

select_representative <- function(indices) {
  ordering <- order(
    is.na(exact_matches$fec_aju_orden[indices]),
    -as.numeric(exact_matches$fec_aju_orden[indices]),
    is.na(exact_matches$fec_arc_xl_orden[indices]),
    -as.numeric(exact_matches$fec_arc_xl_orden[indices]),
    as.character(exact_matches$consecutive[indices]),
    na.last = TRUE
  )
  indices[ordering[[1L]]]
}

representative_indices <- vapply(groups, select_representative, integer(1))
excluded_indices <- unlist(Map(
  function(indices, representative) setdiff(indices, representative),
  groups, representative_indices
), use.names = FALSE)

representative_keys <- clave_registro(exact_matches)[representative_indices]
excluded_keys <- clave_registro(exact_matches)[excluded_indices]

## 3. Generate detailed audit table
audit <- do.call(rbind, lapply(names(groups), function(group_id) {
  indices <- groups[[group_id]]
  representative <- select_representative(indices)
  excluded <- setdiff(indices, representative)
  data.frame(
    grupo_duplicado = as.integer(group_id),
    tamano_grupo = length(indices),
    consecutive_conservado = as.character(exact_matches$consecutive[representative]),
    archivo_conservado = exact_matches$fuente_archivo[representative],
    anio_conservado = exact_matches$fuente_anio[representative],
    fila_excel_conservada = exact_matches$fuente_fila_excel[representative],
    fec_aju_conservado = as.character(exact_matches$fec_aju_orden[representative]),
    fec_arc_xl_conservado = as.character(exact_matches$fec_arc_xl_orden[representative]),
    consecutive_excluido = as.character(exact_matches$consecutive[excluded]),
    archivo_excluido = exact_matches$fuente_archivo[excluded],
    anio_excluido = exact_matches$fuente_anio[excluded],
    fila_excel_excluida = exact_matches$fuente_fila_excel[excluded],
    fec_aju_excluido = as.character(exact_matches$fec_aju_orden[excluded]),
    fec_arc_xl_excluido = as.character(exact_matches$fec_arc_xl_orden[excluded]),
    criterio = "Matching analytical content; the most recent FEC_AJU is kept, then the most recent FEC_ARC_XL.",
    stringsAsFactors = FALSE
  )
}))
rownames(audit) <- NULL

## 4. Process and save each year without redundant records
process_partition <- function(input_path, output_path) {
  datos <- readRDS(input_path)
  keys <- clave_registro(datos)
  exclude_mask <- keys %in% excluded_keys
  output <- datos[!exclude_mask, , drop = FALSE]
  output$representante_duplicado_conservado <- clave_registro(output) %in% representative_keys
  
  saveRDS(output, output_path, compress = "xz")
  
  c(
    original = nrow(datos),
    excluidos = sum(exclude_mask),
    deduplicada = nrow(output),
    representantes = sum(output$representante_duplicado_conservado),
    claves_registro_duplicadas = anyDuplicated(clave_registro(output))
  )
}

control <- vector("list", length(years))
for (i in seq_along(years)) {
  year <- years[[i]]
  output_can <- file.path(canonical_sd_dir, sprintf("dengue_canonica_sin_duplicados_%d.rds", year))
  output_ana <- file.path(analytical_sd_dir, sprintf("dengue_analitica_sin_duplicados_%d.rds", year))
  
  res_can <- process_partition(canonical_paths[[i]], output_can)
  res_ana <- process_partition(analytical_paths[[i]], output_ana)
  
  if (!identical(unname(res_can[1:4]), unname(res_ana[1:4]))) {
    stop("The canonical and analytical partitions differ in ", year, ".")
  }
  
  control[[i]] <- data.frame(
    ano = year,
    registros_originales = unname(res_can[["original"]]),
    redundantes_excluidos = unname(res_can[["excluidos"]]),
    registros_deduplicados = unname(res_can[["deduplicada"]]),
    representantes_conservados = unname(res_can[["representantes"]]),
    claves_registro_duplicadas = unname(res_can[["claves_registro_duplicadas"]]),
    stringsAsFactors = FALSE
  )
  message("Dengue ", year, ": ", res_can[["deduplicada"]], " records kept.")
}
control <- do.call(rbind, control)

## 5. Final reports
summary_data <- data.frame(
  indicador = c(
    INDICADORES_DEDUP[["base_canonica"]],
    INDICADORES_DEDUP[["grupos_redundantes"]],
    INDICADORES_DEDUP[["pertenecen_grupos"]],
    INDICADORES_DEDUP[["excluidos"]],
    INDICADORES_DEDUP[["representantes"]],
    INDICADORES_DEDUP[["base_deduplicada"]]
  ),
  valor = c(
    sum(control$registros_originales), length(groups), nrow(exact_matches),
    length(excluded_keys), length(representative_keys),
    sum(control$registros_deduplicados)
  ),
  stringsAsFactors = FALSE
)

write.csv(audit, file.path(audit_dir, "registros_excluidos_duplicados.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(summary_data, file.path(audit_dir, "resumen_deduplicacion.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(control, file.path(audit_dir, "deduplicacion_por_anio.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## Validations
stopifnot(
  nrow(exact_matches) == length(excluded_keys) + length(representative_keys),
  length(groups) == length(representative_keys),
  length(excluded_keys) == nrow(audit),
  !anyDuplicated(excluded_keys),
  !any(excluded_keys %in% representative_keys),
  sum(control$registros_originales) == filas_originales_esperadas("dengue"),
  sum(control$redundantes_excluidos) == length(excluded_keys),
  all(control$claves_registro_duplicadas == 0L)
)

cat("DENGUE deduplication completed successfully.\n")
