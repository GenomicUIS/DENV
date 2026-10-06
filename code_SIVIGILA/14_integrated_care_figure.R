source("00_config.R", encoding = "UTF-8")

if (!file.exists("00_config.R")) {
  stop("Execute from the project's R/ folder.", call. = FALSE)
}

## -----------------------------------------------------------------------------
## INTEGRATED CARE FIGURE (6-PLOT PANEL)
## Purpose: Calculate care and hospitalization time distributions 
## reading directly from the deduplicated analytical bases.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

out <- file.path(SIVIGILA_PATH, "resultados", "integracion_atencion")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

events <- c("dengue", "dengue_grave", "mortalidad_dengue")
labels <- c(dengue = "Dengue (210)", dengue_grave = "Severe Dengue (220)", mortalidad_dengue = "Mortality (580)")
palette <- c(dengue = "#0072B2", dengue_grave = "#D55E00", mortalidad_dengue = "#8E4585")
shapes <- c(dengue = 16, dengue_grave = 2, mortalidad_dengue = 0)
lines <- c(dengue = "solid", dengue_grave = "longdash", mortalidad_dengue = "dotted")

intervals <- list(
  consulta = c("fec_con", "ini_sin"),
  hospitalizacion = c("fec_hos", "ini_sin"),
  consulta_hospitalizacion = c("fec_hos", "fec_con"),
  notificacion = c("fec_not", "ini_sin"),
  defuncion = c("fec_def", "ini_sin")
)

## -----------------------------------------------------------------------------
## 1. Extraction and calculation of indicators
## -----------------------------------------------------------------------------
hos_list <- list()
times_list <- list()

process_data <- function(d, ev) {
  for (yr in sort(unique(d$ano))) {
    z <- d[d$ano == yr, , drop = FALSE]
    N <- nrow(z)
    if (N == 0) next
    
    si <- z$pac_hos %in% "1"
    
    # 1. Hospitalization
    hos_list[[paste(ev, yr, sep = "_")]] <<- data.frame(
      evento = ev, ano = yr, proporcion = sum(si) / N, stringsAsFactors = FALSE
    )
    
    # 2. Care times
    for (k in names(intervals)) {
      if (k == "defuncion" && ev != "mortalidad_dengue") next
      
      v <- intervals[[k]]
      # Ensure dates exist
      if (!v[1] %in% names(z) || !v[2] %in% names(z)) next
      
      x <- as.integer(z[[v[1]]] - z[[v[2]]])
      eligible <- if (k %in% c("hospitalizacion", "consulta_hospitalizacion")) si else rep(TRUE, N)
      
      valid <- eligible & !is.na(x) & x >= 0
      y <- x[valid]
      nv <- length(y)
      
      times_list[[paste(ev, yr, k, sep = "_")]] <<- data.frame(
        evento = ev, ano = yr, indicador = k, n = nv,
        mismo_dia = sum(y == 0),
        un_dia = sum(y == 1),
        dos_tres_dias = sum(y >= 2 & y <= 3),
        cuatro_mas_dias = sum(y >= 4),
        mediana = if (nv) median(y) else NA_real_,
        p25 = if (nv) unname(quantile(y, .25)) else NA_real_,
        p75 = if (nv) unname(quantile(y, .75)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}

cat("Processing databases for care figure...\n")

# A. Dengue (Iterative by years)
for (a in SIVIGILA_ANIOS) {
  den_path <- file.path(SIVIGILA_PATH, "resultados", "dengue", "datos_procesados", "analitica_sin_duplicados_por_anio", sprintf("dengue_analitica_sin_duplicados_%d.rds", a))
  if (file.exists(den_path)) process_data(readRDS(den_path), "dengue")
}

# B. Severe Dengue
severe_path <- file.path(SIVIGILA_PATH, "resultados", "dengue_grave", "datos_procesados", "dengue_grave_analitica_sin_duplicados.rds")
if (file.exists(severe_path)) process_data(readRDS(severe_path), "dengue_grave")

# C. Mortality (On-the-fly filtering)
mort_path <- file.path(SIVIGILA_PATH, "resultados", "mortalidad_dengue", "datos_procesados", "mortalidad_dengue_canonica_sin_duplicados.rds")
if (file.exists(mort_path)) {
  d_mort <- readRDS(mort_path)
  d_mort <- d_mort[!is.na(d_mort$con_fin) & d_mort$con_fin == "2", ]
  process_data(d_mort, "mortalidad_dengue")
  rm(d_mort)
}

hos <- bind_rows(hos_list)
times <- bind_rows(times_list)

## -----------------------------------------------------------------------------
## 2. Construction of plots
## -----------------------------------------------------------------------------
custom_colors <- function() scale_colour_manual(values = palette, breaks = events, labels = labels, name = NULL)

custom_theme <- function() {
  theme_minimal(base_size = 11, base_family = "sans") + 
    theme(
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 12, colour = "#17202A"),
      plot.subtitle = element_text(size = 9, colour = "#52616F", margin = margin(b = 8)),
      axis.title = element_text(size = 10), 
      axis.text = element_text(size = 9),
      plot.margin = margin(9, 15, 9, 9), 
      legend.position = "none"
    )
}

# Panel A: Hospitalization
p_h <- ggplot(hos, aes(ano, proporcion, colour = evento, shape = evento, linetype = evento)) +
  geom_line(linewidth = .75) + 
  geom_point(size = 2.5) + 
  custom_colors() + 
  scale_shape_manual(values = shapes) + 
  scale_linetype_manual(values = lines) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = percent_format(decimal.mark = ",")) +
  scale_x_continuous(breaks = seq(2007, 2025, 3)) +
  labs(title = "A) Registered hospitalization", subtitle = "2007-2025 · Hospitalized / annual total", x = "Year", y = "Percentage") + 
  custom_theme()

titles <- c(
  consulta = "B) Symptom onset → consultation",
  hospitalizacion = "C) Symptom onset → hospitalization",
  consulta_hospitalizacion = "D) Consultation → hospitalization",
  notificacion = "E) Symptom onset → notification",
  defuncion = "F) Symptom onset → death"
)

plot_panel <- function(k) {
  if (k == "consulta_hospitalizacion") {
    d <- times %>% filter(ano >= 2013, indicador == "consulta_hospitalizacion") %>%
      group_by(evento) %>% summarise(across(c(n, mismo_dia, un_dia, dos_tres_dias, cuatro_mas_dias), sum), .groups = "drop")
    
    cats <- c(mismo_dia = "Same day", un_dia = "1 day", dos_tres_dias = "2-3 days", cuatro_mas_dias = "4 or more days")
    t <- bind_rows(lapply(names(cats), function(c) data.frame(evento = d$evento, categoria = unname(cats[c]), n = d[[c]], total = d$n)))
    t$categoria <- factor(t$categoria, levels = unname(cats))
    t$evento <- factor(t$evento, levels = events)
    t$proporcion <- t$n / t$total
    
    return(ggplot(t, aes(categoria, proporcion, fill = evento)) +
             geom_col(position = position_dodge(.8), width = .72) +
             geom_text(aes(label = percent(proporcion, accuracy = 0.1, decimal.mark = ",")),
                       position = position_dodge(.8), vjust = -0.5, size = 3, colour = "#243443") +
             scale_fill_manual(values = palette) +
             scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, .25), labels = percent_format(decimal.mark = ",")) +
             labs(title = titles[[k]], subtitle = "2013-2025 · Distribution among hospitalized with valid interval", x = NULL, y = "Percentage") + 
             custom_theme())
  }
  
  t <- times %>% filter(indicador == k, ano >= 2013)
  sub <- if (k == "defuncion") "2013-2025 · Only dengue mortality records" else if (k %in% c("hospitalizacion", "consulta_hospitalizacion")) "2013-2025 · Only hospitalized with valid interval" else "2013-2025 · Records with valid interval"
  
  g <- ggplot(t, aes(ano, mediana, colour = evento, shape = evento, linetype = evento)) +
    geom_ribbon(aes(ymin = p25, ymax = p75, fill = evento), alpha = .10, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = .8) + geom_point(size = 2.5) + custom_colors() +
    scale_shape_manual(values = shapes) + scale_linetype_manual(values = lines) + scale_fill_manual(values = palette) +
    scale_x_continuous(breaks = seq(2013, 2025, 2)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .13)), labels = label_number(decimal.mark = ",")) +
    labs(title = titles[[k]], subtitle = sub, x = "Year", y = "Days") + 
    custom_theme()
  
  if (k == "notificacion") {
    g <- g + scale_y_continuous(trans = "log1p", breaks = c(0, 3, 7, 14, 30, 90, 180), limits = c(0, NA), expand = expansion(mult = c(0, .08))) +
      labs(y = "Days (logarithmic scale: log[1 + days])") +
      annotate("label", x = 2018.5, y = 80, label = "Dengue, 2013: median 101 days\nDate anomaly; does not imply care delay", hjust = .5, size = 3, colour = "#7C2D12", fill = "#FFF7ED", linewidth = .2)
  }
  return(g)
}

legend_plot <- ggplot(data.frame(evento = events, x = seq_along(events)), aes(x, 1, colour = evento)) +
  geom_segment(aes(x = x - .22, xend = x - .07, yend = 1, linetype = evento), linewidth = 1) +
  geom_point(aes(x = x - .145, shape = evento), size = 2.5) +
  geom_text(aes(label = unname(labels[evento])), hjust = 0, size = 3.8, colour = "#243443") +
  custom_colors() + scale_shape_manual(values = shapes) + scale_linetype_manual(values = lines) + 
  coord_cartesian(xlim = c(0.6, 3.95), clip = "off") +
  theme_void() + theme(legend.position = "none", plot.margin = margin(0, 20, 0, 20))

panels <- wrap_plots(c(list(p_h), lapply(names(titles), plot_panel)), ncol = 2)

fig <- (legend_plot / panels) + plot_layout(heights = c(.08, 3)) +
  plot_annotation(
    title = "Care, hospitalization, and evolution of dengue in Colombia",
    subtitle = "Comparison of SIVIGILA databases · Hospitalization, care times, and outcome",
    caption = "Source: SIVIGILA, clean bases 2007-2025. A: hospitalized/annual total. B-F: missing dates and negative intervals are excluded.\nB, C, E, and F: median and 25-75 percentile band. E uses logarithmic scale. Comparisons are descriptive.",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", colour = "#17202A"),
      plot.subtitle = element_text(size = 12, colour = "#52616F", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, hjust = 0, lineheight = 1.2, colour = "#52616F"),
      plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(16, 18, 14, 18)
    )
  )

ggsave(file.path(out, "Figura_2_atencion_tres_eventos.tiff"), fig, device = "tiff", width = 16, height = 12, dpi = 600, compression = "lzw", bg = "white")

cat("Integrated care figure successfully generated at:", out, "\n")
