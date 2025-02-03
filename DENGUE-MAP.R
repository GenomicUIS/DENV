# 1. Análisis descriptivo de casos

# 1.1. Tendencia temporal

# Cargar los datos con read.csv
data <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv", 
                 header = TRUE, stringsAsFactors = FALSE, sep = ",")

#Revisar si se cargaron correctamente los datos y verificar los datos:
head(data)
str(data)

# Reemplazar valores faltantes con NA (opcional)
data[data == ""] <- NA

# Ver resumen de valores faltantes
summary(data)

#llamar las librerias
library(ggplot2)
library(dplyr)

dengue_summary <- data %>%
  group_by(ANO) %>%
  summarize(casos = n())

ggplot(dengue_summary, aes(x = ANO, y = casos)) +
  geom_line(color = "#0072B2", linewidth = 1) +  # Línea con color y grosor ajustado
  geom_point(color = "#D55E00", size = 2) +  # Puntos con color y tamaño ajustado
  labs(
    title = "Annual cases of DENV in Colombia",
    x = "Year",
    y = "Number of cases"
  ) +
  scale_x_continuous(breaks = dengue_summary$ANO) +  # Mostrar todos los años en el eje X
  theme_minimal() +  # Tema minimalista
  theme(
    panel.grid = element_blank(),         # Eliminar cuadrícula
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centrar y ajustar el título
    axis.title.x = element_text(size = 14),  # Ajustar tamaño del título del eje X
    axis.title.y = element_text(size = 14),  # Ajustar tamaño del título del eje Y
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)  # Rotar etiquetas del eje X
  )


# 1.2. Distribución por sexo

# Filtrar los datos por sexo M y F
data_filtrada <- data %>%
  filter(SEXO %in% c("M", "F")) %>%
  group_by(ANO, SEXO) %>%
  summarize(casos = n(), .groups = "drop")

# Crear el gráfico mejorado
ggplot(data_filtrada, aes(x = as.factor(ANO), y = casos, fill = SEXO)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras agrupadas por sexo
  scale_fill_manual(values = c("M" = "#0072B2", "F" = "#D55E00"),  # Colores personalizados
                    labels = c("M" = "Male", "F" = "Female")) +
  labs(
    title = "Distribution of DENV cases in Colombia",
    x = "Year",
    y = "Number of cases",
    fill = "Gender"
  ) +
  scale_x_discrete(labels = unique(data_filtrada$ANO)) +  # Mostrar todos los años en el eje X
  theme_minimal() +  # Tema minimalista
  theme(
    panel.grid = element_blank(),         # Quitar líneas de cuadrícula
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centrar y ajustar el título
    axis.title.x = element_text(size = 14),  # Ajustar tamaño del título del eje X
    axis.title.y = element_text(size = 14),  # Ajustar tamaño del título del eje Y
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotar etiquetas del eje X
    legend.title = element_text(size = 12),  # Ajustar título de la leyenda
    legend.text = element_text(size = 10)    # Ajustar texto de la leyenda
  )

# 1.3. Casos según grupo etario

# Agrupación por año y grupo etario
data_etario <- data %>%
  mutate(grupo_etario = case_when(
    EDAD <= 5 ~ "0-5 years",
    EDAD <= 17 ~ "6-17 years",
    EDAD <= 60 ~ "18-60 years",
    TRUE ~ "60+ years"
  )) %>%
  group_by(ANO, grupo_etario) %>%
  summarize(casos = n(), .groups = "drop")


# Crear el gráfico mejorado
ggplot(data_etario, aes(x = as.factor(ANO), y = casos / 1, fill = grupo_etario)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Barras agrupadas con bordes
  scale_fill_manual(values = c(
    "0-5 years" = "#0072B2",
    "6-17 years" = "#D55E00",
    "18-60 years" = "#009E73",
    "60+ years" = "#CC79A7"
  )) +  # Colores personalizados
  labs(
    title = "Distribution of DENV cases in Colombia",
    x = "Year",
    y = "Number of cases",
    fill = "Age group"
  ) +
  scale_x_discrete(labels = unique(data_etario$ANO)) +  # Mostrar todos los años en el eje X
  theme_minimal() +  # Tema minimalista
  theme(
    panel.grid = element_blank(),         # Quitar líneas de cuadrícula
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Centrar y ajustar el título
    axis.title.x = element_text(size = 14),  # Ajustar tamaño del título del eje X
    axis.title.y = element_text(size = 14),  # Ajustar tamaño del título del eje Y
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotar etiquetas del eje X
    legend.title = element_text(size = 12),  # Ajustar título de la leyenda
    legend.text = element_text(size = 10)    # Ajustar texto de la leyenda
  )

# 2. Mapas de distribución geográfica

# 2.1. Mapa nacional por departamentos.

# Instalar y cargar los paquetes necesario

library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(tidyr)

# 1. Cargar los datos desde un archivo CSV
datos <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv")

# 2. Sumarizar los casos por departamento
casos_por_dpto <- datos %>%
  group_by(Departamento_ocurrencia) %>%
  summarize(casos = n(), .groups = "drop") %>%
  mutate(Departamento_ocurrencia = toupper(trimws(Departamento_ocurrencia)))

# Verificar nombres únicos en los datos
print("Departamentos en los datos de casos:")
print(unique(casos_por_dpto$Departamento_ocurrencia))

# 3. Cargar los datos geoespaciales de los departamentos de Colombia
colombia_sf <- st_read("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\SHAPE\\Colombia_departamentos_poblacion.geojson") %>%
  mutate(DPTO_CNMBR = toupper(trimws(DPTO_CNMBR)))

# Verificar nombres únicos en el archivo geoespacial
print("Departamentos en el archivo geoespacial:")
print(unique(colombia_sf$DPTO_CNMBR))

# 4. Unir los datos de casos con el archivo geoespacial
colombia_sf <- colombia_sf %>%
  left_join(casos_por_dpto, by = c("DPTO_CNMBR" = "Departamento_ocurrencia"))

# Verificar sistema de coordenadas
print("CRS del archivo geoespacial:")
print(st_crs(colombia_sf))

# 5. Crear el mapa principal con la cuadrícula y coordenadas
mapa_principal <- ggplot() +
  geom_sf(data = colombia_sf, aes(fill = casos), color = "black", size = 0.1) +
  scale_fill_gradientn(colors = c("#FCDACA", "#F49489", "#EA332F", "#E60D0C", "#7C000C"), na.value = "grey80") +
  labs(
    title = "DENV cases by department in Colombia (2007-2023)",
    fill = "Number of cases",
    caption = "Datum: WGS84. SIVIGILA - 2024."
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +  # Escala en km
  annotation_north_arrow(
    location = "tr", which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +
  coord_sf(crs = st_crs(4326), expand = TRUE) +  # Coordenadas WGS84 con expansión
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"), # Cuadrícula
    axis.text.x = element_text(size = 10, color = "black"),  # Coordenadas eje X
    axis.text.y = element_text(size = 10, color = "black"),  # Coordenadas eje Y
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Mostrar el mapa
print(mapa_principal)

# mapas individuales 2007 a 2023

# Instalar y cargar los paquetes necesarios
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(tidyr)

# 1. Cargar los datos desde un archivo CSV
datos <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv")

# 2. Convertir la columna de fechas y extraer el año
datos <- datos %>%
  mutate(
    FEC_NOT = as.Date(FEC_NOT, format = "%Y-%m-%d"),  # Asegurar formato de fecha
    Año = format(FEC_NOT, "%Y")  # Extraer el año de la fecha
  )

# 3. Cargar los datos geoespaciales de los departamentos de Colombia
colombia_sf <- st_read("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\SHAPE\\Colombia_departamentos_poblacion.geojson") %>%
  mutate(DPTO_CNMBR = toupper(trimws(DPTO_CNMBR)))

# 4. Crear una función para generar los mapas
generar_mapa <- function(año) {
  # Filtrar los datos por año y sumarizar por departamento
  casos_por_dpto <- datos %>%
    filter(Año == año) %>%
    group_by(Departamento_ocurrencia) %>%
    summarize(casos = n(), .groups = "drop") %>%
    mutate(Departamento_ocurrencia = toupper(trimws(Departamento_ocurrencia)))
  
  # Unir los datos geoespaciales con los casos
  mapa_sf <- colombia_sf %>%
    left_join(casos_por_dpto, by = c("DPTO_CNMBR" = "Departamento_ocurrencia"))
  
  # Crear el mapa
  mapa <- ggplot() +
    geom_sf(data = mapa_sf, aes(fill = casos), color = "black", size = 0.1) +
    scale_fill_gradientn(colors = c("#FCDACA", "#F49489", "#EA332F", "#E60D0C", "#7C000C"), na.value = "grey80") +
    labs(
      title = paste("DENV cases by department in Colombia (", año, ")", sep = ""),
      fill = "Number of cases",
      caption = "Datum: WGS84. SIVIGILA - 2024."
    ) +
    annotation_scale(location = "bl", width_hint = 0.5) +  # Escala en km
    annotation_north_arrow(
      location = "tr", which_north = "true",
      style = north_arrow_fancy_orienteering()
    ) +
    coord_sf(crs = st_crs(4326), expand = TRUE) +  # Coordenadas WGS84 con expansión
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Guardar el mapa como archivo PNG
  ggsave(
    filename = paste0("Mapa_Dengue_", año, ".TIFF"),
    plot = mapa,
    width = 10, height = 8, dpi = 300
  )
}

# 5. Generar los mapas para cada año entre 2007 y 2023
for (año in 2007:2023) {
  generar_mapa(as.character(año))
}

# Mensaje final
print("Mapas generados y guardados exitosamente en el directorio de trabajo.")