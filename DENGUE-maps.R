########################################################################################

# 1. Descriptive analysis of cases

# 1.1. Temporal trend

# Load the data with read.csv
data <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv", 
                 header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Check if the data was loaded correctly and verify the data:
head(data)
str(data)

# Replace missing values with NA (optional)
data[data == ""] <- NA

# View summary of missing values
summary(data)

# Load the libraries
library(ggplot2)
library(dplyr)

dengue_summary <- data %>%
  group_by(ANO) %>%
  summarize(casos = n())

ggplot(dengue_summary, aes(x = ANO, y = casos)) +
  geom_line(color = "#0072B2", linewidth = 1) +  # Line with adjusted color and thickness
  geom_point(color = "#D55E00", size = 2) +  # Points with adjusted color and size
  labs(
    title = "Annual cases of DENV in Colombia",
    x = "Year",
    y = "Number of cases"
  ) +
  scale_x_continuous(breaks = dengue_summary$ANO) +  # Show all years on the X axis
  theme_minimal() +  # Minimalist theme
  theme(
    panel.grid = element_blank(),         # Remove grid
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and adjust the title
    axis.title.x = element_text(size = 14),  # Adjust the size of the X axis title
    axis.title.y = element_text(size = 14),  # Adjust the size of the Y axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)  # Rotate X axis labels
  )


# 1.2. Distribution by gender

# Filter the data by gender M and F
data_filtered <- data %>%
  filter(SEXO %in% c("M", "F")) %>%
  group_by(ANO, SEXO) %>%
  summarize(casos = n(), .groups = "drop")

# Create the enhanced plot
ggplot(data_filtered, aes(x = as.factor(ANO), y = casos, fill = SEXO)) +
  geom_bar(stat = "identity", position = "dodge") +  # Bars grouped by gender
  scale_fill_manual(values = c("M" = "#0072B2", "F" = "#D55E00"),  # Custom colors
                    labels = c("M" = "Male", "F" = "Female")) +
  labs(
    title = "Distribution of DENV cases in Colombia",
    x = "Year",
    y = "Number of cases",
    fill = "Gender"
  ) +
  scale_x_discrete(labels = unique(data_filtered$ANO)) +  # Show all years on the X axis
  theme_minimal() +  # Minimalist theme
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and adjust the title
    axis.title.x = element_text(size = 14),  # Adjust the size of the X axis title
    axis.title.y = element_text(size = 14),  # Adjust the size of the Y axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate X axis labels
    legend.title = element_text(size = 12),  # Adjust the legend title
    legend.text = element_text(size = 10)    # Adjust the legend text
  )

# 1.3. Cases by age group

# Grouping by year and age group
data_age_group <- data %>%
  mutate(age_group = case_when(
    EDAD <= 5 ~ "0-5 years",
    EDAD <= 17 ~ "6-17 years",
    EDAD <= 60 ~ "18-60 years",
    TRUE ~ "60+ years"
  )) %>%
  group_by(ANO, age_group) %>%
  summarize(casos = n(), .groups = "drop")


# Create the enhanced plot
ggplot(data_age_group, aes(x = as.factor(ANO), y = casos / 1, fill = age_group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Bars grouped with borders
  scale_fill_manual(values = c(
    "0-5 years" = "#0072B2",
    "6-17 years" = "#D55E00",
    "18-60 years" = "#009E73",
    "60+ years" = "#CC79A7"
  )) +  # Custom colors
  labs(
    title = "Distribution of DENV cases in Colombia",
    x = "Year",
    y = "Number of cases",
    fill = "Age group"
  ) +
  scale_x_discrete(labels = unique(data_age_group$ANO)) +  # Show all years on the X axis
  theme_minimal() +  # Minimalist theme
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and adjust the title
    axis.title.x = element_text(size = 14),  # Adjust the size of the X axis title
    axis.title.y = element_text(size = 14),  # Adjust the size of the Y axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate X axis labels
    legend.title = element_text(size = 12),  # Adjust the legend title
    legend.text = element_text(size = 10)    # Adjust the legend text
  )

# 2. Geographic distribution maps

# 2.1. National map by departments.

# Install and load the necessary packages

library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(tidyr)

# 1. Load the data from a CSV file
data <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv")

# 2. Summarize the cases by department
cases_by_dept <- data %>%
  group_by(Departamento_ocurrencia) %>%
  summarize(casos = n(), .groups = "drop") %>%
  mutate(Departamento_ocurrencia = toupper(trimws(Departamento_ocurrencia)))

# Verify unique names in the data
print("Departments in the case data:")
print(unique(cases_by_dept$Departamento_ocurrencia))

# 3. Load the geospatial data of the departments of Colombia
colombia_sf <- st_read("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\SHAPE\\Colombia_departamentos_poblacion.geojson") %>%
  mutate(DPTO_CNMBR = toupper(trimws(DPTO_CNMBR)))

# Verify unique names in the geospatial file
print("Departments in the geospatial file:")
print(unique(colombia_sf$DPTO_CNMBR))

# 4. Join the case data with the geospatial file
colombia_sf <- colombia_sf %>%
  left_join(cases_by_dept, by = c("DPTO_CNMBR" = "Departamento_ocurrencia"))

# Verify coordinate system
print("CRS of the geospatial file:")
print(st_crs(colombia_sf))

# 5. Create the main map with grid and coordinates
main_map <- ggplot() +
  geom_sf(data = colombia_sf, aes(fill = casos), color = "black", size = 0.1) +
  scale_fill_gradientn(colors = c("#FCDACA", "#F49489", "#EA332F", "#E60D0C", "#7C000C"), na.value = "grey80") +
  labs(
    title = "DENV cases by department in Colombia (2007-2023)",
    fill = "Number of cases",
    caption = "Datum: WGS84. SIVIGILA - 2024."
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +  # Scale in km
  annotation_north_arrow(
    location = "tr", which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +
  coord_sf(crs = st_crs(4326), expand = TRUE) +  # WGS84 coordinates with expansion
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"), # Grid
    axis.text.x = element_text(size = 10, color = "black"),  # X axis coordinates
    axis.text.y = element_text(size = 10, color = "black"),  # Y axis coordinates
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Display the map
print(main_map)

# Individual maps from 2007 to 2023

# Install and load the necessary packages
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(tidyr)

# 1. Load the data from a CSV file
data <- read.csv("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\DENGUE - SIVIGILA\\DENGUE\\Datos_SIVIGILA_Combinados2.csv")

# 2. Convert the date column and extract the year
data <- data %>%
  mutate(
    FEC_NOT = as.Date(FEC_NOT, format = "%Y-%m-%d"),  # Ensure date format
    Year = format(FEC_NOT, "%Y")  # Extract the year from the date
  )

# 3. Load the geospatial data of the departments of Colombia
colombia_sf <- st_read("C:\\Users\\caden\\Documents\\Documents - PC UIS\\Convocatoria HUS\\SHAPE\\Colombia_departamentos_poblacion.geojson") %>%
  mutate(DPTO_CNMBR = toupper(trimws(DPTO_CNMBR)))

# 4. Create a function to generate the maps
generate_map <- function(year) {
  # Filter the data by year and summarize by department
  cases_by_dept <- data %>%
    filter(Year == year) %>%
    group_by(Departamento_ocurrencia) %>%
    summarize(casos = n(), .groups = "drop") %>%
    mutate(Departamento_ocurrencia = toupper(trimws(Departamento_ocurrencia)))
  
  # Join the geospatial data with the cases
  map_sf <- colombia_sf %>%
    left_join(cases_by_dept, by = c("DPTO_CNMBR" = "Departamento_ocurrencia"))
  
  # Create the map
  map <- ggplot() +
    geom_sf(data = map_sf, aes(fill = casos), color = "black", size = 0.1) +
    scale_fill_gradientn(colors = c("#FCDACA", "#F49489", "#EA332F", "#E60D0C", "#7C000C"), na.value = "grey80") +
    labs(
      title = paste("DENV cases by department in Colombia (", year, ")", sep = ""),
      fill = "Number of cases",
      caption = "Datum: WGS84. SIVIGILA - 2024."
    ) +
    annotation_scale(location = "bl", width_hint = 0.5) +  # Scale in km
    annotation_north_arrow(
      location = "tr", which_north = "true",
      style = north_arrow_fancy_orienteering()
    ) +
    coord_sf(crs = st_crs(4326), expand = TRUE) +  # WGS84 coordinates with expansion
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
  
  # Save the map as a PNG file
  ggsave(
    filename = paste0("Dengue_Map_", year, ".TIFF"),
    plot = map,
    width = 10, height = 8, dpi = 300
  )
}

# 5. Generate the maps for each year between 2007 and 2023
for (year in 2007:2023) {
  generate_map(as.character(year))
}

# Final message
print("Maps generated and successfully saved in the working directory")
