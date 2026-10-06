# 03_prepare_dengue_mortality.R
# Harmonization, cleaning and analytical profiling of Dengue Mortality

source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

source("01_load_databases.R", encoding = "UTF-8")

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(readr)
})

## 1. Working directories
base_dir     <- fs::path(SIVIGILA_PATH, "resultados", "mortalidad_dengue")
data_dir    <- fs::path(base_dir, "datos_procesados")
reports_dir <- fs::path(base_dir, "reportes")
fs::dir_create(c(data_dir, reports_dir))

## 2. Database import
# Auxiliary functions (es_vacio, parsear_fecha, normalizar_codigo,
# etiquetar_codigo, unir_divipola_notificacion, calcular_edad_anios,
# grupo_edad_tesis) live in 99_utils.R and are already loaded by 00_config.R.
message("Loading the historical Dengue Mortality database...")
text_data <- cargar_evento("mortalidad_dengue")

## 3. Clinical and epidemiological cleaning
message("Executing clinical and epidemiological cleaning...")

analytical_data <- text_data %>%
  mutate(
    across(any_of(c("ano", "semana", "edad", "uni_med", "sem_ges")), as.integer),
    across(any_of(c("fec_not", "fec_con", "ini_sin", "fec_hos",
                    "fec_def", "fecha_nto", "fec_aju",
                    "fec_arc_xl")), parsear_fecha)
  ) %>%
  mutate(
    divipola_ocurrencia   = unir_divipola_notificacion(cod_dpto_o, cod_mun_o),
    divipola_residencia   = unir_divipola_notificacion(cod_dpto_r, cod_mun_r),
    divipola_notificacion = unir_divipola_notificacion(cod_dpto_n, cod_mun_n),
    llave_upgd = if_else(
      is.na(cod_pre) | is.na(cod_sub),
      NA_character_,
      paste0(cod_pre, "-", cod_sub)
    )
  ) %>%
  mutate(
    edad_anios       = calcular_edad_anios(edad, uni_med),
    grupo_edad_tesis = grupo_edad_tesis(edad_anios),
    dias_inicio_consulta          = as.numeric(difftime(fec_con, ini_sin, units = "days")),
    dias_consulta_hospitalizacion = as.numeric(difftime(fec_hos, fec_con, units = "days")),
    dias_inicio_defuncion         = as.numeric(difftime(fec_def, ini_sin, units = "days")),
    dias_inicio_notificacion      = as.numeric(difftime(fec_not, ini_sin, units = "days"))
  ) %>%
  mutate(
    flag_condicion_final_no_muerto    = is.na(con_fin) | con_fin != "2",
    flag_fecha_defuncion_faltante     = is.na(fec_def),
    flag_cronologia_inconsistente     = (dias_inicio_consulta < 0) |
      (dias_consulta_hospitalizacion < 0) |
      (dias_inicio_defuncion < 0),
    flag_edad_fuera_rango             = is.na(edad_anios) | edad_anios < 0 | edad_anios > 130,
    # Flags derived from raw files review:
    # 1. CER_DEF comes with textual "SIN DATO", which is not a real NA.
    # 2. COD_MUN_R = "000" is a sentinel for unknown municipality.
    # 3. Municipality name can bring "* ... MUNICIPIO DESCONOCIDO".
    flag_cert_def_sentinel            = con_fin == "2" & !is.na(cer_def) & es_vacio(cer_def),
    flag_residencia_municipio_cero    = !is.na(cod_mun_r) & cod_mun_r == "000",
    flag_municipio_texto_desconocido  = !is.na(municipio_residencia) &
      grepl("DESCONOCIDO|SIN DATO|^\\*", municipio_residencia)
  )

## 4. Export
message("Saving clean databases and reports...")

saveRDS(text_data,      fs::path(data_dir, "mortalidad_dengue_original.rds"))
saveRDS(analytical_data, fs::path(data_dir, "mortalidad_dengue_analitica.rds"))

quality_summary <- analytical_data %>%
  summarise(
    total_registros               = n(),
    alertas_no_muertos            = sum(flag_condicion_final_no_muerto, na.rm = TRUE),
    alertas_fechas_def_faltantes  = sum(flag_fecha_defuncion_faltante, na.rm = TRUE),
    alertas_cronologias_negativas = sum(flag_cronologia_inconsistente, na.rm = TRUE),
    alertas_edades_incoherentes   = sum(flag_edad_fuera_rango, na.rm = TRUE),
    alertas_cer_def_sentinel      = sum(flag_cert_def_sentinel, na.rm = TRUE),
    alertas_municipio_cero        = sum(flag_residencia_municipio_cero, na.rm = TRUE),
    alertas_municipio_desconocido = sum(flag_municipio_texto_desconocido, na.rm = TRUE)
  )

write_csv(quality_summary, fs::path(reports_dir, "resumen_calidad.csv"))

message("\nMortality preparation successfully completed!")