# météo_france_TS_analysis

# pathway : ~/Satellite_analysis/météo_france_TS_analysis/

# library -----------------------------------------------------------------

library(tidyverse)
library(tidync)
library(ncdf4)    # For reading NetCDF files
library(lubridate) # For working with dates
library(reshape2) # For data reshaping
library(ggplot2)
library(climaemet)
library(heatwaveR)
library(scales)
library(RColorBrewer)

# load --------------------------------------------------------------------

# "~/Documents/Alice/Données_vent/wind_daily_200801_202501.nc" retourne
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






# Météo France data --------------------------------------------------------

Wind_T <- read.csv("~/Vent/data/Q_06_previous-1950-2024_RR-T-Vent.csv", 
                   header = TRUE, sep = ";")

Wind_T <- Wind_T |> 
  filter(NUM_POSTE == "6088001") |> 
  select("LAT", "LON", "NUM_POSTE", "FFM", "DXY", "HXI", "RR", "TM")

Wind_T <- Wind_T |> 
  mutate(date = seq(as.Date("1950-01-01"), as.Date("2024-12-31"), by = "day"))

Wind_T <- Wind_T |> 
  filter(date >= "1991-01-01")

Wind_2015_2024 <- Wind_T |> 
  filter(date >= as.Date("2015-01-01"), date <= as.Date("2024-12-31"))
  
## separate wind -----------------------------------------------------------

# West
West <- Wind_T |> 
  filter(DXY >= 255, DXY <= 285)
West <- West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# East
East <- Wind_T |> 
  filter(DXY >= 75, DXY <= 105)
East <- East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# North
North <- Wind_T |> 
  filter(DXY >= 345 | DXY <= 15)
North <- North %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Nord Est
North_East <- Wind_T |> 
  filter(DXY >= 30, DXY <= 60)
North_East <- North_East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# North West
North_West <- Wind_T |> 
  filter(DXY >= 300, DXY <= 330)
North_West <- North_West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# South
South <- Wind_T |> 
  filter(DXY >= 165, DXY <= 195)
South <- South %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Sud Ouest
South_West <- Wind_T |> 
  filter(DXY >= 210, DXY <= 240)
South_West <- South_West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# South East
South_East <- Wind_T |> 
  filter(DXY >= 120, DXY <= 150)
South_East <- South_East %>%
  complete(date = seq(min(date), max(date), by = "day"))

## plotting ----------------------------------------------------------------

model_wind <- lm(FFM ~ date, data = Wind_2015_2024)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(Wind_2015_2024, aes(x = date, y = FFM)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(Wind_2015_2024$date, na.rm = TRUE),
    y = max(Wind_2015_2024$FFM, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_wind, 3), " + ", round(slope_wind, 7), " × x",
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
    title = "Évolution de la vitesse du vent près de Nice (2015–2024)",
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


# Wind vs plume datas -----------------------------------------------------

## wind speed vs plume area ------------------------------------------------

### plotting ---------------------------------

# In Gangloff 2017, they envestigated the relationship between the plume and the
# wind velocity and wind direction

adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$aire_panache_km2, Wind_T$FFM)

SEXTANT_1998_2025_spm_95$scaled_aire_panache_km2 <- SEXTANT_1998_2025_spm_95$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Wind_T, 
             aes(x = date, y = FFM, color = "Vitesse du vent"), size = 0.5) +
  geom_point(data = SEXTANT_1998_2025_spm_95, 
             aes(x = date, y = scaled_aire_panache_km2, color = "Aire des panaches"), size = 0.5) +
  scale_color_manual(values = c("Vitesse du vent" = "#4A90D9", "Aire des panaches" = "red3")) +
  scale_y_continuous(
    name = "Vitesse du vent (m/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution de la vitesse du vent et de l'aire des panaches selon le produit SEXTANT OC5",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "5 year",  
    date_labels = "%Y"       
  )

### wind vs plume area correlation ---------------------------------

Vent_SEXTANT_panache <- inner_join(Wind_T, SEXTANT_1998_2025_spm_95, by = "date")

cor.test(Vent_SEXTANT_panache$FFM, Vent_SEXTANT_panache$aire_panache_km2, method = "spearman")


## wind direction vs plume area --------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$aire_panache_km2, North_West$DXY)

SEXTANT_1998_2025_spm_95$scaled_aire_panache_km2 <- SEXTANT_1998_2025_spm_95$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = North_West, 
             aes(x = date, y = DXY, color = "Vent de Nord Ouest"), size = 0.5) +
  geom_point(data = SEXTANT_1998_2025_spm_95, 
             aes(x = date, y = scaled_aire_panache_km2, color = "Aire des panaches"), size = 0.5) +
  scale_color_manual(values = c("Vent de Nord Ouest" = "#4A90D9", "Aire des panaches" = "red3")) +
  scale_y_continuous(
    name = "Vent de Nord Ouest (en °)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Vent de Nord Ouest et aire des panaches selon le produit SEXTANT OC5",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


# wind climatology --------------------------------------------------------

## wind speed climatology --------------------------------------------------

Wind_1991_2020_TS <- Wind_T |> 
  filter(date >= as.Date("1991-01-01"), date <= as.Date("2020-12-31")) |> 
  mutate(year = year(date), 
         month = month(date), 
         doy = yday(date))

Wind_1991_2020_climatology <- Wind_1991_2020_TS %>% 
  dplyr::filter(date >= as.Date("1991-01-01"))

Wind_1991_2020_climatology_year <- Wind_1991_2020_climatology %>% 
  summarise(wind_year_clim = mean(FFM, na.rm = TRUE), .by = "year")

Wind_1991_2020_climatology_month <- Wind_1991_2020_climatology %>%
  group_by(month) %>%
  summarise(
    wind_month_clim = mean(FFM, na.rm = TRUE),
    wind_month_clim_std = sd(FFM, na.rm = TRUE)
  )

Wind_1991_2020_climatology_day <- Wind_1991_2020_climatology %>% 
  summarise(wind_doy_clim = mean(FFM, na.rm = TRUE), .by = "doy")

Wind_1991_2020_climatology_doy <- ts2clm(data = Wind_1991_2020_TS, x = date, 
                                      y = FFM, climatologyPeriod = c("1991-01-01", "2020-12-31"), 
                                      windowHalfWidth = 3, smoothPercentileWidth = 15 )

# Série récente 2015-2024 sur laquelle on calcule les anomalies
Wind_2015_2024_TS <- Wind_T |> 
  filter(date >= as.Date("2015-01-01"), date <= as.Date("2024-12-31")) |> 
  mutate(year  = year(date),
         month = month(date),
         doy   = yday(date))

# Anomalies mensuelles
Wind_monthly_anom <- Wind_2015_2024_TS |> 
  mutate(date = floor_date(date, "month")) |> 
  summarise(mean_FFM = mean(FFM, na.rm = TRUE), .by = c("date", "year", "month")) |> 
  left_join(Wind_1991_2020_climatology_month, by = "month") |> 
  mutate(
    wind_month_anomaly     = mean_FFM - wind_month_clim,
    wind_month_anomaly_std = wind_month_anomaly / wind_month_clim_std  # anomalie normalisée
  )
  # This rounds all dates to the first day of the month
  # That way we can calculate monthly averages, but still have the full
  # date values (e.g. 2023-11-14) that ggplot2 needs to plot the values correctly
  # mutate(date = floor_date(date, "month")) |> 
  # summarise(mean_FFM = mean(FFM, na.rm = TRUE), .by = c("date", "year", "month")) |> 
  # left_join(Wind_1991_2020_climatology_month, by = c("month")) |> 
  # mutate(wind_month_anomaly = mean_FFM - wind_month_clim)

# vent par mois entre 1991 et 2024
# Wind_month <- Wind_T |>
#   mutate(month = floor_date(date, "month")) |>
#   group_by(month) |>
#   summarise(
#     FFM_mean = mean(FFM, na.rm = TRUE),
#     # Moyenne circulaire pour la direction
#     DXY_mean_circ = atan2(
#       mean(sin(DXY * pi / 180), na.rm = TRUE),  # composante Sud-Nord
#       mean(cos(DXY * pi / 180), na.rm = TRUE)   # composante Ouest-Est
#     ) * 180 / pi,
#     n = n()
#   ) |>
#   # Ramener les valeurs négatives entre 0 et 360°
#   mutate(DXY_mean_circ = (DXY_mean_circ + 360) %% 360)
#  
# climatologie du vent 
# Wind_clim <- Wind_T |>
#   filter(date >= as.Date("1991-01-01") & date <= as.Date("2020-12-31")) |>
#   mutate(month = as.numeric(format(date, "%m"))) |>  # extraire le numéro du mois
#   group_by(month) |>
#   summarise(
#     FFM_mean = mean(FFM, na.rm = TRUE),
#     FFM_sd = sd(FFM, na.rm = TRUE),
#     DXY_mean_circ = atan2(
#       mean(sin(DXY * pi / 180), na.rm = TRUE),
#       mean(cos(DXY * pi / 180), na.rm = TRUE)
#     ) * 180 / pi,
#     n = n()
#   ) |>
#   mutate(
#     DXY_mean_circ = (DXY_mean_circ + 360) %% 360,
#     month_name = month.abb[month]  # ajouter le nom du mois en abbrégé
#   )


## wind direction climatology ----------------------------------------------

Wind_1991_2020_TS <- Wind_T |> 
  filter(date >= as.Date("1991-01-01"), date <= as.Date("2020-12-31")) |> 
  mutate(year = year(date), 
         month = month(date), 
         doy = yday(date))

Wind_1991_2020_climatology <- Wind_1991_2020_TS %>%
  dplyr::filter(date >= as.Date("1991-01-01")) %>%
  mutate(Direction =
           case_when(
             DXY >= 255 & DXY <= 285 ~ "Ouest",
             DXY >= 75 & DXY <= 105 ~ "Est",
             DXY >= 345 | DXY <= 15 ~ "Nord",
             DXY >= 30 & DXY <= 60 ~ "Nord - Est",
             DXY >= 300 & DXY <= 330 ~ "Nord - Ouest",
             DXY >= 165 & DXY <= 195 ~ "Sud",
             DXY >= 210 & DXY <= 240 ~ "Sud - Ouest",
             DXY >= 120 & DXY <= 150 ~ "Sud - Est",
             TRUE ~ NA_character_
           ))

# Calculer les proportions par mois
proportions_par_mois <- Wind_1991_2020_climatology %>%
  filter(!is.na(Direction)) %>%  # Exclure les NA
  group_by(month, Direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Proportion = n / sum(n))  # Proportion par mois


## plotting ----------------------------------------------------------------

model_wind_month <- lm(FFM_mean ~ month, data = Wind_month)
p_value_wind_month <- summary(model_wind_month)$coefficients[2, 4]  # p-value pour la pente
intercept_wind_month <- coef(model_wind_month)[1]
slope_wind_month <- coef(model_wind_month)[2]

ggplot(Wind_month, aes(x = month, y = FFM_mean)) +
  geom_line(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(Wind_month$month, na.rm = TRUE),
    y = max(Wind_month$FFM_mean, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_wind_month, 3), " ", round(slope_wind_month, 7), " × x",
      "\np = ", ifelse(p_value_wind_month < 0.001, "< 0.001", format(p_value_wind_month, digits = 3))
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
    title = "Évolution de la vitesse du vent mensuelle près de Nice (1991–2020)",
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


# plotting climatology -------------------------------------------------------------

# Vitesse du vent par mois
ggplot(Wind_clim, aes(x = month, y = FFM_mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = FFM_mean - FFM_sd, 
                    ymax = FFM_mean + FFM_sd), width = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Climatologie de la vitesse du vent (1991 - 2020)",
    x = "Mois", y = "Vitesse du vent (m/s)") +
  theme_bw()


# Direction du vent par mois
ggplot(Wind_clim, aes(x = month, y = DXY_mean_circ)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(limits = c(0, 360), 
                     breaks = c(0, 90, 180, 270, 360),
                     labels = c("N", "E", "S", "W", "N")) +
  labs(x = "Mois", y = "Direction du vent (°)") +
  theme_bw()

# anomalie de la vitesse du vent sur la période 2015 - 2026
model_anom <- lm(wind_month_anomaly ~ date, data = Wind_monthly_anom)
p_value_anom <- summary(model_anom)$coefficients[2, 4]
slope_anom <- coef(model_anom)[2]
intercept_anom <- coef(model_anom)[1]

ggplot(Wind_monthly_anom, aes(x = date, y = wind_month_anomaly, fill = wind_month_anomaly > 0)) +
  
  geom_col(alpha = 0.8, width = 25) +  # width en jours
  geom_smooth(
    aes(x = date, y = wind_month_anomaly),
    method    = "lm",
    se        = TRUE,
    color     = "grey20",
    fill      = "grey60",
    linewidth = 0.8,
    alpha     = 0.2,
    inherit.aes = FALSE  # important pour ne pas hériter du fill des barres
  ) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
  
  scale_fill_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "steelblue4"),
    labels = c("TRUE" = "Au-dessus de la normale", "FALSE" = "En dessous de la normale"),
    name   = NULL
  ) +
  annotate(
    "text",
    x     = min(Wind_monthly_anom$date, na.rm = TRUE),
    y     = max(Wind_monthly_anom$wind_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0("y = ", round(intercept_anom, 3), " + ", round(slope_anom, 7), " × x",
                   "\np = ", ifelse(p_value_anom < 0.001, "< 0.001",
                                    format(p_value_anom, scientific = TRUE, digits = 3))),
    hjust = 0, vjust = 1,
    size  = 4, color = "grey20", fontface = "italic", family = "serif"
  ) +
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  
  labs(
    title    = "Anomalie mensuelle de vitesse du vent (2015–2024)",
    subtitle = "Par rapport à la climatologie 1991–2020",
    x        = NULL,
    y        = expression("Anomalie de vent (m s"^-1*")"),
    caption  = "Source : Météo-France"
  ) +
  
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey96", linewidth = 0.2),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    legend.key       = element_blank()
  )


## climatologie de la direction du vent ------------------------------------

ggplot(proportions_par_mois, aes(x = month, y = Proportion, fill = Direction)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Proportions des directions de vent par mois (1991-2020)",
       x = "Mois",
       y = "Proportion",
       fill = "Direction") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

# pour un graph plus propre

proportions_par_mois$month <- factor(proportions_par_mois$month, levels = 1:12,
                                     labels = month.name)

couleurs <- brewer.pal(n = 8, name = "Set3")
names(couleurs) <- levels(proportions_par_mois$Direction)  # Associe les noms des catégories aux couleurs

p <- ggplot(proportions_par_mois, aes(x = month, y = Proportion, fill = Direction)) +
  # Barres empilées à 100% par mois
  geom_bar(stat = "identity", position = "fill", width = 0.8, color = "black", linewidth = 0.1) +
  # Échelle des y en pourcentages
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  # Palette de couleurs personnalisée
  scale_fill_manual(values = couleurs) +
  # Titres et labels
  labs(x = "Mois",
       y = "Proportion (%)",
       fill = "Direction du vent",
       title = "Proportions mensuelles des directions du vent (1991-2020)",
       caption = "Source : Données climatologiques (1991-2020)") +
  # Thème épuré et professionnel
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.box.background = element_rect(linewidth = 0.5, color = "black"),
    plot.caption = element_text(size = 9, color = "grey50"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Afficher le graphique
print(p)


## température --------------------------------------------------------------

model_temp <- lm(TM ~ date, data = Wind_T)
p_value_temp <- summary(model_temp)$coefficients[2, 4]  # p-value pour la pente
intercept_temp <- coef(model_temp)[1]
slope_temp <- coef(model_temp)[2]

ggplot(Wind_T, aes(x = date, y = TM)) +
  geom_line(color = "orangered", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "orangered", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(Wind_T$date, na.rm = TRUE),
    y = max(Wind_T$TM, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_temp, 3), " + ", round(slope_temp, 7), " × x",
      "\np = ", ifelse(p_value_temp < 0.001, "< 0.001", format(p_value_temp, digits = 3))
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
    title = "Évolution de la température journalière près de Nice (1991–2024)",
    x = NULL,
    y = "Température (°C)",
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


## climatologie de la température ------------------------------------------

Temp_1991_2020_TS <- Wind_T |> 
  filter(date >= as.Date("1991-01-01"), date <= as.Date("2020-12-31")) |> 
  mutate(year = year(date), 
         month = month(date), 
         doy = yday(date))

Temp_1991_2020_climatology <- Wind_1991_2020_TS %>%
  dplyr::filter(date >= as.Date("1991-01-01"))

Temp_1991_2020_climatology_year <- Temp_1991_2020_climatology %>% 
  summarise(temp_year_clim = mean(TM, na.rm = TRUE), .by = "year")

Temp_1991_2020_climatology_month <- Temp_1991_2020_climatology %>%
  group_by(month) %>%
  summarise(
    temp_month_clim = mean(TM, na.rm = TRUE),
    temp_month_clim_std = sd(TM, na.rm = TRUE)
  )

Temp_1991_2020_climatology_day <- Temp_1991_2020_climatology %>% 
  summarise(temp_doy_clim = mean(TM, na.rm = TRUE), .by = "doy")

Wind_1991_2020_climatology_doy <- ts2clm(data = Temp_1991_2020_TS, x = date, 
                                         y = TM, climatologyPeriod = c("1991-01-01", "2020-12-31"), 
                                         windowHalfWidth = 3, smoothPercentileWidth = 15 )

# Série récente 2015-2024 sur laquelle on calcule les anomalies
Temp_2015_2024_TS <- Wind_T |> 
  filter(date >= as.Date("2015-01-01"), date <= as.Date("2024-12-31")) |> 
  mutate(year  = year(date),
         month = month(date),
         doy   = yday(date))

# Anomalies mensuelles
Temp_monthly_anom <- Temp_2015_2024_TS |> 
  mutate(date = floor_date(date, "month")) |> 
  summarise(mean_TM = mean(TM, na.rm = TRUE), .by = c("date", "year", "month")) |> 
  left_join(Temp_1991_2020_climatology_month, by = "month") |> 
  mutate(
    temp_month_anomaly     = mean_TM - temp_month_clim,
    temp_month_anomaly_std = temp_month_anomaly / temp_month_clim_std  # anomalie normalisée
  )

## plotting temperature climatology ----------------------------------------

# plotting de la climatologie
ggplot(Temp_1991_2020_climatology_month, aes(x = month, y = temp_month_clim)) +
  geom_bar(stat = "identity", fill = "orangered") +
  geom_errorbar(aes(ymin = temp_month_clim - temp_month_clim_std, 
                    ymax = temp_month_clim + temp_month_clim_std), width = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Climatologie de la température moyenne (1991 - 2020)",
       x = "Mois", y = "Température (en °C)") +
  theme_bw()

# plotting de l'anomalie à la climatologie
model_anom_temp <- lm(temp_month_anomaly ~ date, data = Temp_monthly_anom)
p_value_anom_temp <- summary(model_anom_temp)$coefficients[2, 4]
slope_anom_temp <- coef(model_anom_temp)[2]
intercept_anom_temp <- coef(model_anom_temp)[1]

ggplot(Temp_monthly_anom, aes(x = date, y = temp_month_anomaly, fill = temp_month_anomaly > 0)) +
  
  geom_col(alpha = 0.8, width = 25) +  # width en jours
  geom_smooth(
    aes(x = date, y = temp_month_anomaly),
    method    = "lm",
    se        = TRUE,
    color     = "grey20",
    fill      = "grey60",
    linewidth = 0.8,
    alpha     = 0.2,
    inherit.aes = FALSE  # important pour ne pas hériter du fill des barres
  ) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
  
  scale_fill_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "steelblue4"),
    labels = c("TRUE" = "Au-dessus de la normale", "FALSE" = "En dessous de la normale"),
    name   = NULL
  ) +
  annotate(
    "text",
    x     = min(Temp_monthly_anom$date, na.rm = TRUE),
    y     = max(Temp_monthly_anom$temp_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0("y = ", round(intercept_anom_temp, 3), " + ", round(slope_anom_temp, 7), " × x",
                   "\np = ", ifelse(p_value_anom_temp < 0.001, "< 0.001",
                                    format(p_value_anom_temp, scientific = TRUE, digits = 3))),
    hjust = 0, vjust = 1,
    size  = 4, color = "grey20", fontface = "italic", family = "serif"
  ) +
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  
  labs(
    title    = "Anomalie mensuelle de la température (2015–2024)",
    subtitle = "Par rapport à la climatologie 1991–2020",
    x        = NULL,
    y        = expression("Température (en °C)"),
    caption  = "Source : Archives Météo-France"
  ) +
  
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey96", linewidth = 0.2),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    legend.key       = element_blank()
  )


# Gangloff et al. 2017 ----------------------------------------------------

# 1. Classifier les vents
Wind_T <- Wind_T |>
  mutate(wind_sector = case_when(
    DXY >= 295 | DXY <= 15  ~ "Offshore", # souffle vers le large
    DXY >= 80  & DXY <= 160 ~ "Onshore",  # souffle vers la côte
    TRUE ~ "Autre"
  ))

# 2. Garder seulement les deux secteurs principaux
data_filtered <- Wind_T |>
  filter(wind_sector != "Autre")

# quelle proportion de "offshore" et "onshore"
prop.table(table(all_data$wind_sector)) * 100

all_data <- data_filtered |>
  left_join(Y6442010_depuis_2000, by = "date") |>
  left_join(SEXTANT_1998_2025_spm_95, by = "date")

# 3. Visualiser aire du panache vs débit, coloré par secteur de vent
ggplot(all_data, aes(x = débit, y = aire_panache_km2, color = wind_sector)) +
  geom_point(alpha = 0.6) +
  scale_y_log10() +
  labs(x = "Débit du Var (m³/s)", 
       y = "Aire du panache (km²)",
       color = "Secteur de vent") +
  theme_bw()

# 4. Vitesse du vent vs aire du panache
ggplot(all_data, aes(x = FFM, y = aire_panache_km2, color = wind_sector)) +
  geom_point(alpha = 0.6) +
  labs(x = "Vitesse du vent (m/s)",
       y = "Aire du panache (km²)",
       color = "Secteur de vent") +
  theme_bw()


# 5. Comparer statistiquement l'influence des vents sur l'aire du panache

# Test de Shapiro
shapiro.test(all_data$aire_panache_km2[all_data$wind_sector == "Offshore"])
# p-value < 2.2e-16 : les données des aires de panache ne sont pas normalement 
# distribuées
shapiro.test(all_data$aire_panache_km2[all_data$wind_sector == "Onshore"])
# p-value < 2.2e-16 : les données  des aires de panache ne sont pas normalement 
# distribuées

# Données non normales --> test de Wilcoxon non paramétrique
wilcox.test(aire_panache_km2 ~ wind_sector, data = all_data)

ggplot(all_data, aes(x = wind_sector, y = aire_panache_km2, fill = wind_sector)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Secteur de vent",
       y = "Aire du panache (km²)") +
  theme_bw()
