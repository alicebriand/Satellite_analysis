# TS analysis MODIS

# 03/03/2026

# pathway : "~/Documents/Alice/Code/MODIS_TS_analysis


# This script will load MODIS satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the 
# Var river at the Napoleon bridge


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(heatwaveR)
library(doParallel); registerDoParallel(cores = detectCores()-2)

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

load("data/MODIS/SPM/MODIS_2015_2024_SPM.Rdata")
load("data/MODIS/SPM/all_spm_propre_MODIS_2023.Rdata")
# load("data/MODIS/CHL/")
load("data/Hydro France/Y6442010_Hydro.Rdata")
load("data/Hydro France/Y6442010_2015_2024.Rdata")

# Y6442010_2015_2024 <- Y6442010_depuis_2000 |> 
#   filter(date >= as.Date("2015-01-01"), date <= as.Date("2024-12-28"))
# 
# save(Y6442010_2015_2024, file = "data/Hydro France/Y6442010_2015_2024.Rdata")

# climatology -------------------------------------------------------------

# MODIS_2015_2025_SPM_TS <- MODIS_2015_2024_SPM |> 
#   filter(analysed_spim >= 0) |> # There appear to be some erroneuos negative values in the data
#   summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = "date") |> 
#   mutate(year = year(date), 
#          month = month(date), 
#          doy = yday(date))

MODIS_2015_2024_SPM <- MODIS_2015_2024_SPM %>% 
  mutate(year = year(date), 
         month = month(date),
         doy = yday(date))

MODIS_2016_2024_SPM_climatology <- MODIS_2015_2024_SPM %>% 
  dplyr::filter(date >= as.Date("2016-01-01"))

MODIS_2016_2024_SPM_climatology_year <- MODIS_2016_2024_SPM_climatology %>% 
  summarise(spm_year_clim = mean(mean_spm, na.rm = TRUE), .by = "year")

MODIS_2016_2024_SPM_climatology_month <- MODIS_2016_2024_SPM_climatology %>% 
  summarise(spm_month_clim = mean(mean_spm, na.rm = TRUE), .by = "month")

MODIS_2016_2024_SPM_climatology_day <- MODIS_2016_2024_SPM_climatology %>% 
  summarise(spm_doy_clim = mean(mean_spm, na.rm = TRUE), .by = "doy")

MODIS_spm_climatology_doy <- ts2clm(data = MODIS_2016_2024_SPM_climatology, x = date, 
                                    y = mean_spm, climatologyPeriod = c("2016-01-01", "2024-12-24"), 
                                    windowHalfWidth = 3, smoothPercentileWidth = 15 )
MODIS_2016_2024_SPM_monthly_anom <- MODIS_2016_2024_SPM_climatology |> 
  # This rounds all dates to the first day of the month
  # That way we can calculate monthly averages, but still have the full
  # date values (e.g. 2024-11-14) that ggplot2 needs to plot the values correctly
  mutate(date = floor_date(date, "month")) |> 
  summarise(mean_spm = mean(mean_spm, na.rm = TRUE), .by = c("date", "year", "month")) |> 
  left_join(MODIS_2016_2024_SPM_climatology_month, by = c("month")) |> 
  mutate(spm_month_anomaly = mean_spm - spm_month_clim)


# plotting ----------------------------------------------------------------

## climatology -------------------------------------------------------------

# create a line plot of the annual climatology of spm
ggplot(MODIS_2016_2024_SPM_climatology_year, aes(x = year, y = spm_year_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie annuelle de la concentration en matière particulaire en suspension entre 2016 et 2024 avec le produit MODIS",
       x = "Année",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the monthly climatology of spm
ggplot(MODIS_2016_2024_SPM_climatology_month, aes(x = month, y = spm_month_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie mensuelle de la concentration en matière particulaire en suspension entre 2016 et 2024 avec le produit MODIS",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the daily climatology of spm
ggplot(MODIS_2016_2024_SPM_climatology_day, aes(x = doy, y = spm_doy_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2016 et 2024 avec le produit MODIS",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the daily climatology of spm
ggplot(MODIS_spm_climatology_doy, aes(x = doy, y = seas)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2016 et 2024 avec le produit MODIS et lissé sur une fenêtre de 7 jours",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a plot that shows the anomaly of spm of 10 years of data view by Sextant
ggplot(MODIS_2016_2024_SPM_monthly_anom, aes(x = date, y = spm_month_anomaly)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Anomalie mensuelle de la concentration en matière particulaire en suspension entre 2016 et 2024 avec le produit MODIS",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## SPM ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS

ggplot(data = MODIS_2015_2024_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_point(color = "red3", size = 0.5) +
  labs(title = "Évolution de la concentration en matière particulaire en suspension moyenne entre 2015 et 2023 avec le produit MODIS issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

ggplot(data = all_spm_propre_MODIS_2023, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne en 2023 avec le produit MODIS issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(MODIS_2015_2024_SPM$mean_spm, Y6442010_2015_2024$débit)

MODIS_2015_2024_SPM$scaled_mean_spm <- MODIS_2015_2024_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(
    data = Y6442010_2015_2024,
    aes(x = date, y = débit, color = "Débit")
  ) +
  geom_point(
    data = MODIS_2015_2024_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension entre 2015 et 2024 avec le produit MODIS",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# scatter plot ------------------------------------------------------------

ggplot(Var_MODIS_SPM, aes(x = débit, y = mean_spm)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", fill = "pink", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Relation entre débit et concentration en SPM en échelle log",
    x = "Débit (m³/s)",
    y = "Concentration moyenne en SPM (g/m³)"
  ) +
  theme_minimal()

# runoff vs SPM concentration correlation ---------------------------------

Var_MODIS_SPM <- inner_join(Y6442010_2015_2024, MODIS_2015_2024_SPM, by = "date")

cor.test(Var_MODIS_SPM$débit, Var_MODIS_SPM$mean_spm, method = "spearman")







## CHL ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS
ggplot(data = MODIS_2015_2023_CHL, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "chartreuse3") +
  labs(title = "Evolution de la concentration en chlorophylle moyenne entre 2015 et 2023 avec le produit Sextant",
       x = "Date",
       y = "Concentration moyenne en chlorophylle (en µg/L)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(MODIS_2015_2023_CHL$mean_chl, Y6442010_Hydro_complete$débit)

MODIS_2015_2023_CHL$scaled_mean_chl <- MODIS_2015_2023_CHL$mean_chl * adjust_factors$diff + adjust_factors$adjust


ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = MODIS_2015_2023_CHL,
    aes(x = date, y = scaled_mean_chl, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "chartreuse3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Chlorophylle (en µg/L)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en chlorophylle entre 2015 et 2023 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal()
