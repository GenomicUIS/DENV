# 04_deduplicate_dengue_mortality.R
# Elimination of duplicated records and creation of a strict sample for graphics

# ## 1. CONFIGURATION AND PACKAGES ##
source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fs)
})

# Input and output directories
input_path <- fs::path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "datos_procesados", "mortalidad_dengue_analitica.rds")
data_dir <- fs::path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "datos_procesados")
audit_dir <- fs::path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "reportes", "deduplicacion")

fs::dir_create(c(data_dir, audit_dir))

# ## 2. DATA LOADING ##
message("Loading the mortality analytical database...")
datos <- readRDS(input_path)

# ## 3. IDENTIFICATION AND DEDUPLICATION ##
message("Identifying clinically identical records...")

# A. Define technical columns that should NOT be compared
# If the file or Excel row changes, but the clinic is the same, it IS a duplicate
technical_columns <- c(
  "fuente_evento", "fuente_archivo", "fuente_anio", "fuente_fila_excel",
  "separacion", "particion", "consecutive", "consecutive_origen", 
  "fec_arc_xl", "fec_aju", "confirmados", "va_sispro", 
  "estado_final_de_caso", "nom_est_f_caso", "nombre_evento", "cod_eve"
)

# B. Group and mark duplicates
grouped_data <- datos %>%
  group_by(across(-any_of(technical_columns))) %>%
  mutate(
    es_duplicado = n() > 1,
    grupo_duplicado_id = cur_group_id() # R assigns a unique number to each identical group
  ) %>%
  ungroup()

# C. Separate clean records from duplicates
unique_data <- grouped_data %>% filter(!es_duplicado)
raw_duplicates <- grouped_data %>% filter(es_duplicado)

# D. Select the "Representative" (the most recent one)
message("Resolving duplicate conflicts...")
kept_representatives <- raw_duplicates %>%
  group_by(grupo_duplicado_id) %>%
  arrange(desc(fec_aju), desc(fec_arc_xl), consecutive) %>%
  slice(1) %>%
  ungroup()

# E. Final clean database
data_without_duplicates <- bind_rows(unique_data, kept_representatives)

# ## 4. AUDIT ##
excluded <- raw_duplicates %>%
  anti_join(kept_representatives, by = c("fuente_anio", "fuente_archivo", "fuente_fila_excel"))

audit <- excluded %>%
  select(
    grupo_duplicado = grupo_duplicado_id,
    consecutive_excluido = consecutive,
    archivo_excluido = fuente_archivo,
    anio_excluido = fuente_anio,
    fec_aju_excluido = fec_aju
  ) %>%
  left_join(
    kept_representatives %>%
      select(
        grupo_duplicado = grupo_duplicado_id,
        consecutive_conservado = consecutive,
        archivo_conservado = fuente_archivo,
        anio_conservado = fuente_anio,
        fec_aju_conservado = fec_aju
      ),
    by = "grupo_duplicado"
  ) %>%
  mutate(criterio = "Identical content; the record with the most recent FEC_AJU was kept.")

# ## 5. FINAL BASE FOR GRAPHICS ##
graphics_data <- data_without_duplicates %>%
  filter(con_fin == "2") %>%
  select(-es_duplicado, -grupo_duplicado_id) # We clean temporary variables

data_without_duplicates <- data_without_duplicates %>% select(-es_duplicado, -grupo_duplicado_id)

# ## 5.5. DEDUPLICATION SUMMARY ##
# Single source of truth for downstream readers (e.g. 05_generate_dengue_mortality_graphics.R).
# Written BEFORE section 6 so that auditors see it alongside the audit table.
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
    nrow(datos),
    length(unique(raw_duplicates$grupo_duplicado_id)),
    nrow(raw_duplicates),
    nrow(excluded),
    nrow(kept_representatives),
    nrow(data_without_duplicates)
  ),
  stringsAsFactors = FALSE
)

# ## 6. EXPORT ##
message("Saving results...")
saveRDS(data_without_duplicates, fs::path(data_dir, "mortalidad_dengue_canonica_sin_duplicados.rds"), compress = "xz")
write_csv(data_without_duplicates, fs::path(data_dir, "mortalidad_dengue_canonica_sin_duplicados.csv"))

saveRDS(graphics_data, fs::path(data_dir, "mortalidad_dengue_muestra_graficas.rds"), compress = "xz")
write_csv(audit, fs::path(audit_dir, "registros_excluidos_duplicados.csv"))
write_csv(summary_data, fs::path(audit_dir, "resumen_deduplicacion.csv"))

# ## 7. CONSOLE REPORT ##
cat("\n--- DEDUPLICATION SUMMARY ---\n")
cat("Original evaluated records:", nrow(datos), "\n")
cat("Redundant records excluded:", nrow(excluded), "\n")
cat("Final deduplicated base:", nrow(data_without_duplicates), "\n")
cat("Strict sample for graphics (CON_FIN = 2):", nrow(graphics_data), "\n")
