# TS analysis SEXTANT

# 03/03/2026

# pathway : "~/Satellite_analysis/TS analysis SEXTANT.R


# This script will load Sextant satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the Var river


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2) # Detects cores automagically
library(heatwaveR)

# functions -----------------------------------------------------------------

## scaling function --------------------------------------------------------

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

# loading data ------------------------------------------------------------

load("data/SEXTANT/sextant_2015_2025_SPM.Rdata")
load("data/SEXTANT/all_spm_propre_sextant_2024.RData")
load("data/sextant_2015_2025_CHL.Rdata")
load("data/Hydro France/Y6442010_Hydro.Rdata")

# climatology -------------------------------------------------------------

sextant_2015_2025_spm_TS <- sextant_2015_2025_spm |> 
  filter(analysed_spim >= 0) |> # There appear to be some erroneuos negative values in the data
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = "date") |> 
  mutate(year = year(date), 
         month = month(date), 
         doy = yday(date))

sextant_2015_2025_SPM_climatology <- sextant_2015_2025_spm_TS %>% 
  dplyr::filter(date >= as.Date("2016-01-01"))

sextant_2015_2025_SPM_climatology_year <- sextant_2015_2025_SPM_climatology %>% 
  summarise(spm_year_clim = mean(mean_spm, na.rm = TRUE), .by = "year")

sextant_2015_2025_SPM_climatology_month <- sextant_2015_2025_SPM_climatology %>% 
  summarise(spm_month_clim = mean(mean_spm, na.rm = TRUE), .by = "month")

sextant_2015_2025_SPM_climatology_day <- sextant_2015_2025_SPM_climatology %>% 
  summarise(spm_doy_clim = mean(mean_spm, na.rm = TRUE), .by = "doy")

sextant_spm_climatology_doy <- ts2clm(data = sextant_2015_2025_spm_TS, x = date, 
                                      y = mean_spm, climatologyPeriod = c("2016-01-01", "2025-12-31"), 
                                      windowHalfWidth = 3, smoothPercentileWidth = 15 )

# faire une TS en faisant la climatologie mensuelle en spm moins la moyenne mensuelle en spm et voir comment ça rend par rapport à 0

# THis can cause problems when the date is not a full date, e.g. 2022-12-14
# sextant_2015_2025_SPM <- sextant_2015_2025_SPM %>% 
#   mutate(year = year(date), 
#          month = month(date), 
#          month_year = paste(month, year, sep = "-"))

# sextant_2015_2025_SPM <- sextant_2015_2025_SPM %>% 
#   mutate(month_year = as.Date(month_year, form"%m%Y"))

# mean_spm_by_year_month <- sextant_2015_2025_SPM %>%
#   group_by(year, month, year_month) %>%
#   summarise(mean_spm_year_month = mean(mean_spm, na.rm = TRUE))

# sextant_2015_2025_SPM <- sextant_2015_2025_SPM %>%
#   left_join(mean_spm_by_year_month, by = c("year", "month", "year_month"))

# i want to put a new column on the sextant_2015_2025_SPM data frame that is the 
# column spm mean of the sextant_2015_2025_SPM_climatology_month data frame group by the column year_motnh of the sextant_2015_2025_SPM data frame

# sextant_2015_2025_SPM <- sextant_2015_2025_SPM %>%
#   mutate(month = month(date)) %>%
#   left_join(, by = "month")
# 
# sextant_2015_2025_SPM <- sextant_2015_2025_SPM %>%
#   left_join(sextant_2015_2025_SPM_climatology_month, by = "month")

sextant_2015_2025_SPM_monthly_anom <- sextant_2015_2025_spm_TS |> 
  # This rounds all dates to the first day of the month
  # That way we can calculate monthly averages, but still have the full
  # date values (e.g. 2023-11-14) that ggplot2 needs to plot the values correctly
  mutate(date = floor_date(date, "month")) |> 
  summarise(mean_spm = mean(mean_spm, na.rm = TRUE), .by = c("date", "year", "month")) |> 
  left_join(sextant_2015_2025_SPM_climatology_month, by = c("month")) |> 
  mutate(spm_month_anomaly = mean_spm - spm_month_clim)

# plotting ----------------------------------------------------------------

## climatology -------------------------------------------------------------

# create a line plot of the annual climatology of spm
ggplot(sextant_2015_2025_SPM_climatology_year, aes(x = year, y = spm_year_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie annuelle de la concentration en matière particulaire en suspension entre 2015 et 2025 avec le produit Sextant",
       x = "Année",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the monthly climatology of spm
ggplot(sextant_2015_2025_SPM_climatology_month, aes(x = month, y = spm_month_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie mensuelle de la concentration en matière particulaire en suspension entre 2015 et 2025 avec le produit Sextant",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the daily climatology of spm
ggplot(sextant_2015_2025_SPM_climatology_day, aes(x = doy, y = spm_doy_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2015 et 2025 avec le produit Sextant",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the daily climatology of spm
ggplot(sextant_spm_climatology_doy, aes(x = doy, y = seas)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2015 et 2025 avec le produit Sextant et lissé sur une fenêtre de 7 jours",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a plot that shows the anomaly of spm of 10 years of data view by Sextant
ggplot(sextant_2015_2025_SPM_monthly_anom, aes(x = date, y = spm_month_anomaly)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Anomalie mensuelle de la concentration en matière particulaire en suspension entre 2015 et 2025 avec le produit Sextant",
     x = "Mois",
     y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## SPM ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec sextant
ggplot(data = sextant_2015_2025_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 2015 et 2025 avec le produit Sextant",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec sextant contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(sextant_2015_2025_SPM$mean_spm, Y6442010_Hydro_complete$débit)

sextant_2015_2025_SPM$scaled_mean_spm <- sextant_2015_2025_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = sextant_2015_2025_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension entre 2015 et 2025 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on fait un zoom sur 2024
# faire un data frame sur l'année 2024 pour le débit


Y6442010_Hydro_2024 <- Y6442010_Hydro_complete |> 
  dplyr::filter(Date >= as.Date("2024-01-01"), Date <= as.Date("2024-12-31"))

adjust_factors <- sec_axis_adjustement_factors(all_spm_propre_sextant_2024$mean_spm, Y6442010_Hydro_2024$débit)

all_spm_propre_sextant_2024$scaled_mean_spm <- all_spm_propre_sextant_2024$mean_spm * adjust_factors$diff + adjust_factors$adjust


ggplot() +
  geom_line(
    data = Y6442010_Hydro_2024,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = all_spm_propre_sextant_2024,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension en 2024 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


## CHL ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec sextant
ggplot(data = sextant_2015_2025_CHL, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "chartreuse3") +
  labs(title = "Evolution de la concentration en chlorophylle moyenne entre 2015 et 2025 avec le produit Sextant",
       x = "Date",
       y = "Concentration moyenne en chlorophylle (en µg/L)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec sextant contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(sextant_2015_2025_CHL$mean_chl, Y6442010_Hydro_complete$débit)

sextant_2015_2025_CHL$scaled_mean_chl <- sextant_2015_2025_CHL$mean_chl * adjust_factors$diff + adjust_factors$adjust


ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = sextant_2015_2025_CHL,
    aes(x = date, y = scaled_mean_chl, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "chartreuse3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Chlorophylle (en µg/L)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en chlorophylle entre 2015 et 2025 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )
