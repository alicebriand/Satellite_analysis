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

ggplot(wind_data, aes(x = date, y = wind_speed, color = wind_dir)) +
  geom_point() +
  scale_color_gradientn(colors = c("blue", "red", "yellow", "green"), name = "Direction (degrés)") +
  labs(title = "Vitesse du vent entre 2021 et 2025",
       x = "Date",
       y = "Vitesse (m/s)") +
  theme_minimal()

# script de Robert

speed <- wind_data$wind_speed
direction <- wind_data$wind_dir

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
  plot_title = "direction du vent",
  stack_reverse = TRUE) +
    labs(
    subtitle = "2021-2020",
    caption = "Source: CMEMS-HR"
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


