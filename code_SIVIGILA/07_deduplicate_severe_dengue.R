source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## DEDUPLICATION OF THE HISTORICAL SEVERE DENGUE DATABASE
## Purpose: Identify and exclude redundant records (identical analytical 
## content), keeping the most recent version according to the adjustment date.
## -----------------------------------------------------------------------------

## 1. Define input and output paths for clean data and auditing.
input_path <- file.path(
  SIVIGILA_PATH, "resultados", "dengue_grave", "datos_procesados",
  "dengue_grave_canonica.rds"
)
data_dir <- file.path(SIVIGILA_PATH, "resultados", "dengue_grave", "datos_procesados")
audit_dir <- file.path(
  SIVIGILA_PATH, "resultados", "dengue_grave", "reportes", "deduplicacion"
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

## Verify that the canonical base exists before continuing.
if (!file.exists(input_path)) {
  stop("The canonical base does not exist. Execute the severe dengue preparation script first.")
}

## 2. Load the data and verify that they contain duplicate marks.
datos <- readRDS(input_path)
stopifnot(
  "flag_duplicado_contenido_analitico" %in% names(datos),
  "grupo_duplicado_contenido_analitico" %in% names(datos)
)

## 3. Group the records marked as duplicates by analytical content.
candidate_indices <- which(datos$flag_duplicado_contenido_analitico)
groups <- split(
  candidate_indices,
  datos$grupo_duplicado_contenido_analitico[candidate_indices]
)

## Function to select which record stays from each duplicate group.
## Criteria: The one with the most recent adjustment date (fec_aju) is preferred,
## then the most recent extraction (fec_arc_xl) and finally the consecutive.
select_representative <- function(indices) {
  order_idx <- order(
    is.na(datos$fec_aju[indices]),
    -as.numeric(datos$fec_aju[indices]),
    is.na(datos$fec_arc_xl[indices]),
    -as.numeric(datos$fec_arc_xl[indices]),
    as.character(datos$consecutive[indices]),
    na.last = TRUE
  )
  indices[order_idx[[1L]]]
}

## 4. Identify which ones are kept and which ones are excluded.
indices_to_keep <- vapply(groups, select_representative, integer(1))
indices_to_exclude <- setdiff(candidate_indices, indices_to_keep)

## 5. Generate the audit report detailing what was excluded and why.
audit <- do.call(rbind, lapply(names(groups), function(group_id) {
  indices <- groups[[group_id]]
  kept <- select_representative(indices)
  excluded <- setdiff(indices, kept)
  data.frame(
    grupo_duplicado = as.integer(group_id),
    consecutive_conservado = as.character(datos$consecutive[kept]),
    archivo_conservado = as.character(datos$fuente_archivo[kept]),
    anio_conservado = datos$fuente_anio[kept],
    fila_excel_conservada = datos$fuente_fila_excel[kept],
    fec_aju_conservado = as.character(datos$fec_aju[kept]),
    fec_arc_xl_conservado = as.character(datos$fec_arc_xl[kept]),
    consecutive_excluido = as.character(datos$consecutive[excluded]),
    archivo_excluido = as.character(datos$fuente_archivo[excluded]),
    anio_excluido = datos$fuente_anio[excluded],
    fila_excel_excluida = datos$fuente_fila_excel[excluded],
    fec_aju_excluido = as.character(datos$fec_aju[excluded]),
    fec_arc_xl_excluido = as.character(datos$fec_arc_xl[excluded]),
    criterio = paste(
      "Matching analytical content; the record with the most recent FEC_AJU",
      "and then the most recent FEC_ARC_XL is kept."
    ),
    stringsAsFactors = FALSE
  )
}))
rownames(audit) <- NULL

## 6. Create the clean database without redundant records.
data_without_duplicates <- datos[!seq_len(nrow(datos)) %in% indices_to_exclude, , drop = FALSE]
data_without_duplicates$representante_duplicado_conservado <-
  clave_registro(data_without_duplicates) %in%
  clave_registro(datos)[indices_to_keep]

## Save the canonical base without duplicates (in RDS and CSV).
saveRDS(
  data_without_duplicates,
  file.path(data_dir, "dengue_grave_canonica_sin_duplicados.rds"),
  compress = "xz"
)
write.csv(
  data_without_duplicates,
  file.path(data_dir, "dengue_grave_canonica_sin_duplicados.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

## 7. Create the analytical sub-table without duplicates.
analytical_vars <- unique(c(
  "fuente_archivo", "fuente_anio", "fuente_fila_excel", "consecutive",
  "representante_duplicado_conservado", "ano", "semana", "fec_not", "ini_sin",
  "fec_con", "fec_hos", "fec_def", "fec_aju", "edad", "uni_med", "edad_anios",
  "grupo_edad_tesis", "sexo", "sexo_etiqueta", "area", "area_etiqueta",
  "ocupacion", "tip_ss", "regimen_salud_etiqueta", "cod_ase", "per_etn",
  "pertenencia_etnica_etiqueta", "gru_pob", "nom_grupo", "estrato",
  grep("^gp_", names(data_without_duplicates), value = TRUE), "sem_ges",
  "cod_dpto_o", "cod_mun_o", "divipola_ocurrencia", "departamento_ocurrencia",
  "municipio_ocurrencia", "cod_dpto_r", "cod_mun_r", "divipola_residencia",
  "departamento_residencia", "municipio_residencia", "cod_dpto_n", "cod_mun_n",
  "divipola_notificacion", "departamento_notificacion", "municipio_notificacion",
  "cod_pre", "cod_sub", "llave_upgd", "nom_upgd", "tip_cas",
  "tipo_caso_inicial_etiqueta", "pac_hos", "hospitalizado_etiqueta", "con_fin",
  "condicion_final_etiqueta", "ajuste", "ajuste_etiqueta", "estado_final_de_caso",
  "nom_est_f_caso", "cer_def", "cbmte",
  grep("^dias_", names(data_without_duplicates), value = TRUE),
  grep("^flag_", names(data_without_duplicates), value = TRUE),
  grep("^grupo_duplicado|^grupo_posible", names(data_without_duplicates), value = TRUE)
))
analytical_vars <- intersect(analytical_vars, names(data_without_duplicates))
analytical_data <- data_without_duplicates[, analytical_vars, drop = FALSE]

## Save the analytical base without duplicates.
saveRDS(
  analytical_data,
  file.path(data_dir, "dengue_grave_analitica_sin_duplicados.rds"),
  compress = "xz"
)
write.csv(
  analytical_data,
  file.path(data_dir, "dengue_grave_analitica_sin_duplicados.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

## 8. Export the audit files and the general deduplication summary.
write.csv(
  audit,
  file.path(audit_dir, "registros_excluidos_duplicados.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

summary_data <- data.frame(
  indicador = c(
    INDICADORES_DEDUP[["base_canonica"]],
    INDICADORES_DEDUP[["grupos_redundantes"]],
    INDICADORES_DEDUP[["excluidos"]],
    INDICADORES_DEDUP[["base_deduplicada"]],
    INDICADORES_DEDUP[["representantes"]],
    INDICADORES_DEDUP[["grupos_clave_estricta"]]
  ),
  valor = c(
    nrow(datos), length(groups), length(indices_to_exclude),
    nrow(data_without_duplicates), length(indices_to_keep),
    length(unique(na.omit(data_without_duplicates$grupo_posible_duplicado_clave_estricta[
      data_without_duplicates$flag_posible_duplicado_clave_estricta &
        !data_without_duplicates$flag_duplicado_contenido_analitico
    ])))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_data, file.path(audit_dir, "resumen_deduplicacion.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

## 9. Final validations to ensure clean base integrity.
stopifnot(
  length(indices_to_keep) == length(groups),
  length(indices_to_exclude) == length(candidate_indices) - length(groups),
  nrow(data_without_duplicates) == nrow(datos) - length(indices_to_exclude),
  nrow(analytical_data) == nrow(datos) - length(indices_to_exclude),
  anyDuplicated(clave_registro(data_without_duplicates)) == 0L,
  sum(
    duplicated(data_without_duplicates$grupo_duplicado_contenido_analitico) &
      !is.na(data_without_duplicates$grupo_duplicado_contenido_analitico)
  ) == 0L
)

cat("Severe dengue deduplication completed.\n")
cat("Original records:", nrow(datos), "\n")
cat("Redundant records excluded:", length(indices_to_exclude), "\n")
cat("Base for analysis and figures:", nrow(data_without_duplicates), "\n")
