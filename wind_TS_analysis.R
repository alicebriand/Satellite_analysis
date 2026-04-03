# TS analysis of wind series

# pathway : ~/Documents/Alice/Données_vents/

# library -----------------------------------------------------------------

library(tidyverse)
library(tidync)
library(ncdf4)    # For reading NetCDF files
library(lubridate) # For working with dates
library(reshape2) # For data reshaping
library(ggplot2)
library(climaemet)

# load --------------------------------------------------------------------

# "~/Documents/Alice/Données_vent/wind_daily_200801_202501.nc" retroune
#  cat("Longitude range:", range(wind_df$longitude, na.rm = TRUE), "\n")
# Longitude range: 3.5625 5.9375 
# > cat("Latitude range:", range(wind_df$latitude, na.rm = TRUE), "\n")
# Latitude range: 42.3125 43.9375 

# Donc la zone de Nice n'est pas couverte, on charge un autre fichier

filename <- "~/Documents/Alice/Données_vent/CMEMS-HR/cmems_obs-wind_glo_phy_my_l4_0.125deg_PT1H_1770803483672.nc"

nc_data <- nc_open(filename)
print(names(nc_data$dim))
# [1] "time"      "latitude"  "longitude"
print(names(nc_data$var))
# [1] "eastward_wind"  "northward_wind"

nc_close(nc_data)

basename(filename)


load_wind <- function(filename) {
  # Extraction de la date à partir du nom de fichier (optionnel, si nécessaire)
  file_caracter <- substr(basename(filename), start = 12, stop = 17)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  wind_df <- tidync(filename) %>%
    hyper_filter(longitude = dplyr::between(longitude, 6.8925000, 7.4200000),
                 latitude = dplyr::between(latitude, 43.2136389, 43.7300000)) %>%
    hyper_tibble() %>%
    dplyr::rename(u = eastward_wind, v = northward_wind, lon = longitude, lat = latitude) %>%
    mutate(date = as.Date(time)) %>%
    dplyr::select(date, lon, lat, u, v) %>%
    summarise(u = mean(u, na.rm = TRUE),
              v = mean(v, na.rm = TRUE),
              .by = "date")
}

wind_data <- load_wind(filename)

wind_data$date_int <- seq(1, nrow(wind_data))

wind_data <- wind_data |> 
  mutate(wind_speed = round(sqrt(u^2 + v^2), 2),
         wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))


# code Louis --------------------------------------------------------------

# Scale one value to another for tidier double-y-axis plots
sec_axis_adjustement_factors <- function(var_to_scale, var_ref) {
  
  index_to_keep <- which(is.finite(var_ref))
  var_ref <- var_ref[index_to_keep]
  
  index_to_keep <- which(is.finite(var_to_scale))
  var_to_scale <- var_to_scale[index_to_keep]
  
  max_var_to_scale <- max(var_to_scale, na.rm = T) 
  min_var_to_scale <- min(var_to_scale, na.rm = T) 
  max_var_ref <- max(var_ref, na.rm = T) 
  min_var_ref <- min(var_ref, na.rm = T) 
  
  diff_to_scale <- max_var_to_scale - min_var_to_scale
  diff_to_scale <- ifelse(diff_to_scale == 0, 1 , diff_to_scale)
  diff_ref <- max_var_ref - min_var_ref
  diff <- diff_ref / diff_to_scale
  
  adjust <- (max_var_ref - max_var_to_scale*diff) 
  
  return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (var_to_scale * diff) + adjust",
                    trans_axis_operation = "var_to_scale = {scaled_var - adjust} / diff)"))
}

scaling_factor <- sec_axis_adjustement_factors(var_to_scale = wind_data$wind_dir, 
                                               var_ref = wind)

your_data <- your_data %>% mutate(second_y-axis_value_scaled = second_y-axis_value  * scaling_factor$diff + scaling_factor$adjust)

# plotting ----------------------------------------------------------------

ggplot(data = wind_data) +
  geom_segment(aes(x= date_int, y = 0, xend = u, yend = v))

ggplot(wind_data, aes(x = date, y = wind_speed)) +
  geom_line(aes(color = wind_dir)) +
  labs(title = "Vitesse du vent entre 2021 et 2025",
       x = "Date",
       y = "Vitesse (m/s)") +
  theme_minimal()

ggplot(wind_data, aes(x = date, y = wind_dir)) +
  geom_line(color = "red") +
  labs(title = "Vitesse du vent entre 2021 et 2025",
       x = "Date",
       y = "Direction (degrés)") +
  theme_minimal()

ggplot(wind_mean, aes(x = date, y = wind_speed, color = wind_dir)) +
  geom_point() +
  scale_color_gradientn(colors = c("blue", "red", "yellow", "green"), name = "Direction (degrés)") +
  labs(title = "Vitesse et direction du vent (1994 - 2025)",
       x = "Date",
       y = "Vitesse (m/s)") +
  theme_minimal()

ggplot(wind_mean, aes(x = date, y = wind_speed, color = wind_dir)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_gradientn(
    colors = c("#2C3E7A", "#4A90D9", "#A8D8A8", "#F4D03F", "#E74C3C", "#2C3E7A"),
    values = scales::rescale(c(0, 90, 180, 270, 360)),
    limits = c(0, 360),
    breaks = c(0, 90, 180, 270, 360),
    labels = c("N (0°)", "E (90°)", "S (180°)", "O (270°)", "N (360°)"),
    name   = "Direction du vent",
    guide  = guide_colorbar(
      barwidth  = 12,
      barheight = 0.5,
      title.position = "top",
      title.hjust    = 0.5
    )
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Vitesse et direction du vent près de Nice (1994–2025)",
    x       = NULL,
    y       = "Vitesse du vent (m s⁻¹)",
    caption = "Source : CMEMS-HR — Global Ocean Hourly Reprocessed Sea Surface Wind and Stress from Scatterometer and Model"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y  = element_text(size = 11, margin = margin(r = 10)),
    axis.text     = element_text(size = 10, color = "grey30"),
    axis.ticks    = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9, color = "grey30")
  )



# 2007 - 2025 wind --------------------------------------------------------

load("~/Vent/data/wind_2007_2025.Rdata")

wind_2007_2025 <- wind_2007_2025 |> 
  mutate(date_time = as.POSIXct(time * 3600, tz = "UTC", origin = "2007-01-11 00:00:00"))

# on moyenne le eastward_wind et le northward_wind par le temps

wind_mean <- wind_2007_2025 |> 
  mutate(date = as.Date(date_time)) |> 
  summarise(
    u = as.numeric(mean(eastward_wind)),
    v = as.numeric(mean(northward_wind)),
         .by = "date"
         ) |> 
  mutate(wind_speed = round(sqrt(u^2 + v^2), 2),
         wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))


## plotting ----------------------------------------------------------------

speed <- wind_mean$wind_speed
direction <- wind_mean$wind_dir

ggwindrose(
  speed = speed,
  direction = direction,
  n_directions = 8,
  n_speeds = 5,
  speed_cuts = NA,
  col_pal = "GnBu",
  legend_title = "Wind speed (m/s)",
  calm_wind = 0,
  n_col = 1,
  facet = NULL,
  plot_title = "Direction et vitesse du vent entre 2007 et 2025 près de Nice",
  stack_reverse = TRUE) +
  labs(
    subtitle = "2007-2020",
    caption = "Source: CMEMS-HR"
  )

# model_wind <- lm(eastward_wind_mean ~ time, data = wind_mean)
model_wind <- lm(wind_speed ~ date, data = wind_mean)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(wind_mean, aes(x = date, y = wind_speed)) +
  geom_line(color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
        "text",
        x = max(wind_mean$date, na.rm = TRUE),
        y = max(wind_mean$wind_speed, na.rm = TRUE) * 0.9,
        label = paste0(
          "y = ", round(intercept_wind, 3), " + ", round(slope_wind, 7), " * x",
          "\n", "p = ", ifelse(p_value_wind < 0.001, "< 0.001", format(p_value_wind, digits = 3))
        ),
        hjust = 1,  # Alignement à droite
        vjust = 1,  # Alignement en haut
        size = 6
      ) +
  labs(title = "Vitesse du vent entre 2007 et 2025 près de Nice",
       x = "Date",
       y = "Vitesse (m/s)") +
  theme_minimal()

# wind_speed_dir <- wind_mean |>
#   mutate(
#     wind_speed = sqrt(eastward_wind_mean^2 + northward_wind_mean^2),
#     wind_dir   = (atan2(eastward_wind_mean, northward_wind_mean) * 180 / pi) %% 360
#   )
# 
# ggwindrose(
#   speed        = as.numeric(wind_speed_dir$wind_speed),
#   direction    = as.numeric(wind_speed_dir$wind_dir),
#   n_directions = 8,
#   n_speeds     = 5,
#   col_pal      = "GnBu",
#   legend_title = "Wind speed (m/s)",
#   calm_wind    = 0,
#   plot_title   = "Direction du vent"
# ) +
#   labs(
#     subtitle = "2007-2025",
#     caption  = "Source: CMEMS-MR"
#   )


# 1994 - 2025 wind --------------------------------------------------------

load("~/Vent/data/wind_1994_2007.Rdata")
load("~/Vent/data/wind_2008_2025.Rdata")

wind_1994_2007 <- wind_1994_2007 |> 
  mutate(date_time = as.POSIXct(time * 3600, tz = "UTC", origin = "1994-06-01 00:00:00"))

# pour 2008 - 2025, on créer une nouvelle colonne pour le time avec - 8520 pour 
# revenir à time = 1 au début de la série

wind_2008_2025 <- wind_2008_2025 |> 
  mutate(time_new = time - 8520)

wind_2008_2025 <- wind_2008_2025 |> 
  mutate(date_time = as.POSIXct(time_new * 3600, tz = "UTC", origin = "2008-01-01 00:00:00"))

# on moyenne le eastward_wind et le northward_wind par le temps pour la première série

wind_first <- wind_1994_2007 |> 
  mutate(date = as.Date(date_time)) |> 
  summarise(
    u = as.numeric(mean(eastward_wind)),
    v = as.numeric(mean(northward_wind)),
    .by = "date"
  ) |> 
  mutate(wind_speed = round(sqrt(u^2 + v^2), 2),
         wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))

# on moyenne le eastward_wind et le northward_wind par le temps pour la deuxième série

wind_second <- wind_2008_2025 |> 
  mutate(date = as.Date(date_time)) |> 
  summarise(
    u = as.numeric(mean(eastward_wind)),
    v = as.numeric(mean(northward_wind)),
    .by = "date"
  ) |> 
  mutate(wind_speed = round(sqrt(u^2 + v^2), 2),
         wind_dir = round((270-(atan2(v, u)*(180/pi)))%%360))

# on peut ensuite merge les deux df

wind_mean <- rbind(wind_first, wind_second)

## plotting ----------------------------------------------------------------

speed <- wind_mean$wind_speed
direction <- wind_mean$wind_dir

ggwindrose(
  speed = speed,
  direction = direction,
  n_directions = 8,
  n_speeds = 5,
  speed_cuts = NA,
  col_pal = "GnBu",
  legend_title = "Vitesse du vent (m/s)",
  calm_wind = 0,
  n_col = 1,
  facet = NULL,
  plot_title = "Direction et vitesse du vent entre 1994 et 2025 près de Nice",
  stack_reverse = TRUE) +
  labs(
    subtitle = "1994-2025",
    caption = "Source: CMEMS-HR"
  )

# model_wind <- lm(eastward_wind_mean ~ time, data = wind_mean)
model_wind <- lm(wind_speed ~ date, data = wind_mean)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(wind_mean, aes(x = date, y = wind_speed)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(wind_mean$date, na.rm = TRUE),
    y = max(wind_mean$wind_speed, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_wind, 3), " ", round(slope_wind, 7), " × x",
      "\np = ", ifelse(p_value_wind < 0.001, "< 0.001", format(p_value_wind, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "#2C3E7A",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Évolution de la vitesse du vent près de Nice (1994–2025)",
    x = NULL,
    y = "Vitesse du vent (m s⁻¹)",
    caption = "Source : CMEMS-HR — Global Ocean Hourly Reprocessed Sea Surface Wind and Stress from Scatterometer and Model"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y  = element_text(size = 11, margin = margin(r = 10)),
    axis.text     = element_text(size = 10, color = "grey30"),
    axis.ticks    = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70", linewidth = 0.5)
  )

# separate wind -----------------------------------------------------------

# West
West <- wind_mean |> 
  filter(wind_dir > 255, wind_dir < 285)

# Complete missing dates in the date range
West <- West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# East
East <- wind_mean |> 
  filter(wind_dir > 75, wind_dir < 105)

# Complete missing dates in the date range
East <- East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Nord Est
North_East <- wind_mean |> 
  filter(wind_dir > 30, wind_dir < 60)

# Complete missing dates in the date range
North_East <- North_East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Sud Ouest
South_West <- wind_mean |> 
  filter(wind_dir > 210, wind_dir < 240)

# Complete missing dates in the date range
South_West <- South_West %>%
  complete(date = seq(min(date), max(date), by = "day"))


## plotting ----------------------------------------------------------------

# model_wind <- lm(eastward_wind_mean ~ time, data = wind_mean)
model_wind <- lm(wind_speed ~ date, data = North_East)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(North_East, aes(x = date, y = wind_speed)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(North_East$date, na.rm = TRUE),
    y = max(North_East$wind_speed, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_wind, 3), " ", round(slope_wind, 7), " × x",
      "\np = ", ifelse(p_value_wind < 0.001, "< 0.001", format(p_value_wind, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    family = "serif",
    color = "#2C3E7A",
    fontface = "italic"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Évolution de la vitesse du vent du Nord Est près de Nice (1994–2025)",
    x = NULL,
    y = "Vitesse du vent (m s⁻¹)",
    caption = "Source : CMEMS-HR — Global Ocean Hourly Reprocessed Sea Surface Wind and Stress from Scatterometer and Model"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y  = element_text(size = 11, margin = margin(r = 10)),
    axis.text     = element_text(size = 10, color = "grey30"),
    axis.ticks    = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70", linewidth = 0.5)
  )

# Thomas wind data --------------------------------------------------------

Wind_T <- read.csv("~/Vent/data/Q_06_previous-1950-2024_RR-T-Vent.csv", 
                   header = TRUE, sep = ";")

Wind_T <- Wind_T |> 
  filter(NUM_POSTE == "6088001") |> 
  select("LAT", "LON", "NUM_POSTE", "FFM", "DXY", "HXI", "RR")

Wind_T <- Wind_T |> 
  mutate(date = seq(as.Date("1950-01-01"), as.Date("2024-12-31"), by = "day"))

Wind_T <- Wind_T |> 
  filter(date >= "1994-06-01")

## plotting ----------------------------------------------------------------

model_wind <- lm(FFM ~ date, data = Wind_T)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(Wind_T, aes(x = date, y = FFM)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(Wind_T$date, na.rm = TRUE),
    y = max(Wind_T$FFM, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_wind, 3), " ", round(slope_wind, 7), " × x",
      "\np = ", ifelse(p_value_wind < 0.001, "< 0.001", format(p_value_wind, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "#2C3E7A",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Évolution de la vitesse du vent près de Nice (1994–2025)",
    x = NULL,
    y = "Vitesse du vent (m s⁻¹)",
    caption = "Source : Archives Météo France"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y  = element_text(size = 11, margin = margin(r = 10)),
    axis.text     = element_text(size = 10, color = "grey30"),
    axis.ticks    = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70", linewidth = 0.5)
  )


# rose des vents

speed <- Wind_T$FFM
direction <- Wind_T$DXY

ggwindrose(
  speed = speed,
  direction = direction,
  n_directions = 8,
  n_speeds = 5,
  speed_cuts = NA,
  col_pal = "GnBu",
  legend_title = "Vitesse du vent (m/s)",
  calm_wind = 0,
  n_col = 1,
  facet = NULL,
  plot_title = "Direction et vitesse du vent entre 1994 et 2024 près de Nice",
  stack_reverse = TRUE) +
  labs(
    subtitle = "1994-2024",
    caption = "Source: Archives Météo France"
  )

