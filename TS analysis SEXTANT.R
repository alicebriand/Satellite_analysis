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
library(ggpmisc)
library(ggpubr)
library(patchwork)
library(seasonal)
library(sf)
library(zoo)

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

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# loading data ------------------------------------------------------------

load("data/SEXTANT/SPM/sextant_1998_2025_SPM.Rdata")
load("data/SEXTANT/SPM/all_spm_propre_sextant_2024.RData")
load("data/SEXTANT/CHL/sextant_1998_2025_CHL.Rdata")
load("data/Hydro France/Y6442010_depuis_2000.Rdata")
load("data/SEXTANT/SPM/sextant_2015_2025_SPM.Rdata")
load("data/SEXTANT/SPM/sextant_2001_2020_SPM.Rdata")
load("data/SEXTANT/CHL/SEXTANT_1998_2025_chl_pixels.RData")
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")

# climatology of MES -------------------------------------------------------------

# combien de valeurs négatives
sum(SEXTANT_1998_2025_spm_pixels$analysed_spim < 0, na.rm = TRUE)
# [1] 21080

# supprimer seulement les valeurs négatives
SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_pixels |>
  filter(analysed_spim >= 0 | is.na(analysed_spim))

SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# on choisit une période de longue (1998 - 2025)
SEXTANT_1998_2025 <- SEXTANT_1998_2025_spm_clean |> 
  filter(date >= as.Date("1998-01-01"), date <= as.Date("2025-12-31"))

SEXTANT_1998_2025_stat <- SEXTANT_1998_2025_spm_clean |>
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# climatologie annuelle
SEXTANT_1998_2025_spm_year <- SEXTANT_1998_2025_stat |> 
  group_by(year) |> 
  summarise(mean_spm_year_clim = mean(analysed_spim, na.rm = TRUE), 
            median_spm_year_clim = median(analysed_spim, na.rm = TRUE),
            sd_spm_year_clim = sd(analysed_spim, na.rm = TRUE))

# climatologie mensuelle
SEXTANT_1998_2025_spm_month <- SEXTANT_1998_2025_stat |> 
  group_by(month) |>
  summarise(mean_spm_month_clim = mean(analysed_spim, na.rm = TRUE), 
            median_spm_month_clim = median(analysed_spim, na.rm = TRUE),
            sd_spm_month_clim = sd(analysed_spim, na.rm = TRUE))

# climatologie journalière
SEXTANT_1998_2025_spm_doy <- SEXTANT_1998_2025_stat |>
  group_by(doy) |> 
  summarise(mean_spm_doy_clim = mean(analysed_spim, na.rm = TRUE), 
            median_spm_doy_clim = median(analysed_spim, na.rm = TRUE),
            sd_spm_doy_clim = sd(analysed_spim, na.rm = TRUE))

# créer une climatologie à partir d'une TS journalière
sextant_spm_climatology_doy <- ts2clm(data = SEXTANT_1998_2025_spm_mean, x = date, 
                                      y = mean_spm, climatologyPeriod = c("1998-01-01", "2017-12-31"), 
                                      windowHalfWidth = 3, smoothPercentileWidth = 15 )

# anomalie mensuelle
sextant_1998_2025_spm_monthly_anom <- SEXTANT_1998_2025_spm_clean |> 
  mutate(date = floor_date(date, "month")) |> 
  filter(date >= as.character.Date("2015-01-01"), date <= as.Date ("2025-12-31")) |>
  summarise(mean_spm_month = mean(analysed_spim, na.rm = TRUE), .by = c("date", "year", "month")) |>
  left_join(SEXTANT_1998_2025_spm_month, by = c("month")) |> 
  mutate(spm_month_anomaly = mean_spm_month - mean_spm_month_clim)


# climatology of plume area -------------------------------------------------------------

# on ajoute au df la date avec l'année, le mois et le jour de l'année
SEXTANT_1998_2025_spm_95 <- SEXTANT_1998_2025_spm_95 |> 
  # filter(analysed_spim >= 0) |> # There appear to be some erroneuos negative values in the data
  # summarise(mean_panache = mean(analysed_spim, na.rm = TRUE), .by = "date") |> 
  mutate(year = year(date),
         month = month(date),
         doy = yday(date))

# on crée la climatologie annuelle
SEXTANT_1998_2025_panache_year <- SEXTANT_1998_2025_spm_95 %>% 
  group_by(year) |> 
  summarise(mean_panache_year_clim = mean(aire_panache_km2, na.rm = TRUE),
            median_panache_year_clim = median(aire_panache_km2, na.rm = TRUE),
            sd_panache_year_clim = sd(aire_panache_km2, na.rm = TRUE))

# climatologie mensuelle
SEXTANT_1998_2025_panache_month <- SEXTANT_1998_2025_spm_95 %>%
  group_by(month) %>%
  summarise(mean_panache_month_clim = mean(aire_panache_km2, na.rm = TRUE), 
              median_panache_month_clim = median(aire_panache_km2, na.rm = TRUE),
              sd_panache_month_clim = sd(aire_panache_km2, na.rm = TRUE))

# climatologie journalière
SEXTANT_1998_2025_panache_doy <- SEXTANT_1998_2025_spm_95 |>
  group_by(doy) |> 
  summarise(mean_panache_doy_clim = mean(aire_panache_km2, na.rm = TRUE), 
            median_panache_doy_clim = median(aire_panache_km2, na.rm = TRUE),
            sd_panache_doy_clim = sd(aire_panache_km2, na.rm = TRUE))

# créer une climatologie à partir d'une TS journalière
# sextant_panache_climatology_doy <- ts2clm(data = SEXTANT_1998_2025_panache_TS, x = date,
#                                       y = aire_panache_km2, climatologyPeriod = c("2001-01-01", "2020-12-31"),
#                                       windowHalfWidth = 3, smoothPercentileWidth = 15 )


# anomalie mensuelle
sextant_1998_2025_panache_monthly_anom <- SEXTANT_1998_2025_spm_95 |> 
  mutate(date = floor_date(date, "month")) |> 
  filter(date >= as.character.Date("2015-01-01"), date <= as.Date ("2025-12-31")) |>
  summarise(mean_panache_month = mean(aire_panache_km2, na.rm = TRUE), .by = c("date", "year", "month")) |>
  left_join(SEXTANT_1998_2025_panache_month, by = c("month")) |> 
  mutate(panache_month_anomaly = mean_panache_month - mean_panache_month_clim)

# climatology of chl ------------------------------------------------------------

# on choisit une période de longue (1998 - 2025)
SEXTANT_1998_2025 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("1998-01-01"), date <= as.Date("2025-12-31"))

SEXTANT_1998_2025_stat <- SEXTANT_1998_2025_chl_clean |>
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# climatologie annuelle
SEXTANT_1998_2025_chl_year <- SEXTANT_1998_2025_stat |> 
  group_by(year) |> 
  summarise(mean_chl_year_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_year_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_year_clim = sd(analysed_chl_a, na.rm = TRUE))

# climatologie mensuelle
SEXTANT_1998_2025_chl_month <- SEXTANT_1998_2025_stat |> 
  group_by(month) |>
  summarise(mean_chl_month_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_month_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_month_clim = sd(analysed_chl_a, na.rm = TRUE))

# climatologie journanlière
SEXTANT_1998_2025_chl_doy <- SEXTANT_1998_2025_stat |>
  group_by(doy) |> 
  summarise(mean_chl_doy_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_doy_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_doy_clim = sd(analysed_chl_a, na.rm = TRUE))

# créer une climatologie à partir d'une TS journalière
sextant_chl_climatology_doy <- ts2clm(data = SEXTANT_1998_2025_chl_mean, x = date, 
                                      y = mean_chl, climatologyPeriod = c("1998-01-01", "2017-12-31"), 
                                      windowHalfWidth = 3, smoothPercentileWidth = 15 )

# anomalie mensuelle
sextant_1998_2025_chl_monthly_anom <- SEXTANT_1998_2025_chl_clean |> 
  mutate(date = floor_date(date, "month")) |> 
  filter(date >= as.character.Date("2015-01-01"), date <= as.Date ("2025-12-31")) |>
  summarise(mean_chl_month = mean(analysed_chl_a, na.rm = TRUE), .by = c("date", "year", "month")) |>
  left_join(SEXTANT_1998_2025_chl_month, by = c("month")) |> 
  mutate(chl_month_anomaly = mean_chl_month - mean_chl_month_clim)

# plotting ----------------------------------------------------------------

## climatology of MES -------------------------------------------------------------

# create a line plot of the annual climatology of spm
ggplot(SEXTANT_1998_2025_spm_year, aes(x = year, y = mean_spm_year_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie annuelle de la concentration en matière particulaire en suspension entre 2001 et 2020 avec le produit Sextant",
       x = "Année",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the monthly climatology of spm
ggplot(SEXTANT_1998_2025_spm_month, aes(x = month, y = mean_spm_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_spm_month_clim - sd_spm_month_clim,
      ymax = mean_spm_month_clim + sd_spm_month_clim
    ),
    fill = "steelblue", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "steelblue")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration en MES — Sextant OC5",
    x       = NULL,
    y       = "Concentration en MES (g/m³)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998-2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the daily climatology of spm
ggplot(SEXTANT_1998_2025_spm_doy, aes(x = doy, y = mean_spm_doy_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2001 et 2020 avec le produit Sextant",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

# create a line plot of the daily climatology of spm
ggplot(sextant_spm_climatology_doy, aes(x = doy, y = seas)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de la concentration en matière particulaire en suspension entre 2001 et 2020 avec le produit Sextant et lissé sur une fenêtre de 7 jours",
       x = "Mois",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal()

## climatology of plume -------------------------------------------------------------

# create a line plot of the annual climatology of spm
ggplot(SEXTANT_1998_2025_panache_year, aes(x = year, y = mean_panache_year_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie annuelle de l'extension des panaches turbides entre 2001 et 2020 avec le produit Sextant OC5",
       x = "Année",
       y = "Extension des panaches (en km²)") +
  theme_minimal()

# create a line plot of the monthly climatology of spm
ggplot(SEXTANT_1998_2025_panache_month, aes(x = month, y = mean_panache_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_panache_month_clim - sd_panache_month_clim,
      ymax = mean_panache_month_clim + sd_panache_month_clim
    ),
    fill = "steelblue", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "steelblue")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de l'extension des panaches turbides — Sextant OC5",
    x       = NULL,
    y       = "Extension des panaches (en km²)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998-2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the daily climatology of panache
ggplot(sextant_2001_2020_panache_climatology_day, aes(x = doy, y = panache_doy_clim)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de l'extension des panaches turbides entre 2001 et 2020 avec le produit Sextant OC5",
       x = "Mois",
       y = "Extension des panaches (en km²)") +
  theme_minimal()

# create a line plot of the daily climatology of panache
ggplot(sextant_panache_climatology_doy, aes(x = doy, y = seas)) +
  geom_line(color = "blue") +
  geom_point(color = "red3") +
  labs(title = "Climatologie journalière de l'extension des panaches turbides entre 2001 et 2020 avec le produit Sextant OC5 (lissage sur fenêtre de 7 jours)",
       x = "Mois",
       y = "Extension des panaches (en km²)") +
  theme_minimal()

## climatology chl -------------------------------------------------------------

# create a line plot of the yearly climatology of chl
ggplot(SEXTANT_1998_2025_chl_year, aes(x = year, y = mean_chl_year_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_year_clim - sd_chl_year_clim,
      ymax = mean_chl_year_clim + sd_chl_year_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie annuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie annuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie annuelle" = "chartreuse3")
    # ) +
    # scale_x_continuous(
    #   breaks = 1:12,
    #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
    #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie annuelle de la concentration en chlorophylle a — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the monthly climatology of chl
ggplot(SEXTANT_1998_2025_chl_month, aes(x = month, y = mean_chl_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_month_clim - sd_chl_month_clim,
      ymax = mean_chl_month_clim + sd_chl_month_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "chartreuse3")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration en chlorophylle a — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the daily climatology of chl
ggplot(SEXTANT_1998_2025_chl_doy, aes(x = doy, y = mean_chl_doy_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_doy_clim - sd_chl_doy_clim,
      ymax = mean_chl_doy_clim + sd_chl_doy_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie journalière"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie journalière"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie journalière" = "chartreuse3")
  ) +
  # scale_x_continuous(
  #   breaks = 1:12,
  #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
  #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  # ) +
  labs(
    title   = "Climatologie journalière de la concentration en chlorophylle a — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# anomalie mensuelle
model_sextant_chl_anom <- lm(chl_month_anomaly ~ date, data = sextant_1998_2025_chl_monthly_anom)
p_value_sextant_chl_anom <- summary(model_sextant_chl_anom)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_chl_anom <- coef(model_sextant_chl_anom)[1]
slope_sextant_chl_anom <- coef(model_sextant_chl_anom)[2]

# Créer le graphique
ggplot(sextant_1998_2025_chl_monthly_anom, aes(x = date, y = chl_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(chl_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(chl_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_chl_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_chl_monthly_anom$chl_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_chl_anom, 3), " + ", round(slope_sextant_chl_anom, 7), " × x",
      "\np = ", ifelse(p_value_sextant_chl_anom < 0.001, "< 0.001", format(p_value_sextant_chl_anom, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en chlorophylle a (1998–2025) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 1998–2025"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

## monthly anomaly of MES ---------------------------------------------------------

# Extraire le modèle linéaire
model_sextant_1998 <- lm(spm_month_anomaly ~ date, data = sextant_1998_2025_SPM_monthly_anom)
p_value_sextant_1998 <- summary(model_sextant_1998)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998 <- coef(model_sextant_1998)[1]
slope_sextant_1998 <- coef(model_sextant_1998)[2]

# Créer le graphique
ggplot(sextant_1998_2025_SPM_monthly_anom, aes(x = date, y = spm_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(spm_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(spm_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_SPM_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_SPM_monthly_anom$spm_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_1998, 3), " + ", round(slope_sextant_1998, 7), " × x",
      "\np ", ifelse(p_value_sextant_1998 < 0.001, "< 0.001", format(p_value_sextant_1998, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en matière particulaire en suspension (1998–2025) — Sextant OC5",
    x       = NULL,
    y       = "Concentration moyenne en matière particulaire en suspension (en g/m³)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 2001–2020"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

## monthly anomaly of turbid plumes ---------------------------------------------------------

# Extraire le modèle linéaire
model_sextant_1998 <- lm(panache_month_anomaly ~ date, data = sextant_1998_2025_SPM_monthly_anom)
p_value_sextant_1998 <- summary(model_sextant_1998)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998 <- coef(model_sextant_1998)[1]
slope_sextant_1998 <- coef(model_sextant_1998)[2]

# Créer le graphique
ggplot(sextant_1998_2025_SPM_monthly_anom, aes(x = date, y = spm_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(spm_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(spm_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_SPM_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_SPM_monthly_anom$spm_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_1998, 3), " + ", round(slope_sextant_1998, 7), " × x",
      "\np ", ifelse(p_value_sextant_1998 < 0.001, "< 0.001", format(p_value_sextant_1998, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en MES (1998–2025) — Sextant OC5",
    x       = NULL,
    y       = "Anomalie de concentration en MES (g/m³)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 2001–2020"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

## patchwork ----------------------------------------------------------------

### MES ---------------------------------------------------------------------

# --- Graphique 1 : climatologie mensuelle ---
p1 <- ggplot(SEXTANT_1998_2025_spm_month, aes(x = month, y = mean_spm_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_spm_month_clim - sd_spm_month_clim,
      ymax = mean_spm_month_clim + sd_spm_month_clim
    ),
    fill = "steelblue", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "steelblue")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration en matière particulaire en suspension — Sextant OC5",
    x       = NULL,
    y       = "Concentration moyenne en MES (en g/m³)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998-2025 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# Extraire le modèle linéaire
model_sextant_1998 <- lm(spm_month_anomaly ~ date, data = sextant_1998_2025_spm_monthly_anom)
p_value_sextant_1998 <- summary(model_sextant_1998)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998 <- coef(model_sextant_1998)[1]
slope_sextant_1998 <- coef(model_sextant_1998)[2]

# --- Graphique 2 : anomalie mensuelle ---
p2 <- ggplot(sextant_1998_2025_spm_monthly_anom, aes(x = date, y = spm_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(spm_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(spm_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_spm_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_spm_monthly_anom$spm_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_1998, 3), " + ", round(slope_sextant_1998, 7), " × x",
      "\np =", ifelse(p_value_sextant_1998 < 0.001, "< 0.001", format(p_value_sextant_1998, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 6,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en matière particulaire en suspension — Sextant OC5",
    x       = NULL,
    y       = "Concentration moyenne en MES (en g/m³)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 1998-2025"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

# --- Patchwork ---
p1 / p2 +
  plot_annotation(
    title   = "Concentration en Matière Particulaire en Suspension — Sextant OC5",
    caption = "Source : Sextant OC5",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )

### turbid plume ---------------------------------------------------------------------

# --- Graphique 1 : climatologie mensuelle ---
p1 <- ggplot(SEXTANT_1998_2025_panache_month, aes(x = month, y = mean_panache_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_panache_month_clim - sd_panache_month_clim,
      ymax = mean_panache_month_clim + sd_panache_month_clim
    ),
    fill = "steelblue", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(values = c("Climatologie mensuelle" = "steelblue")) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de l'extension des panaches turbides",
    x       = NULL,
    y       = "Extension des panaches (km²)",
    color   = NULL,
    caption = "Période de référence : 1998-2025 | Ruban : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 12, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 10, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 12, margin = margin(r = 10)),
    axis.text        = element_text(size = 11, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10)
  )

# --- Modèle linéaire ---
model_sextant_1998   <- lm(panache_month_anomaly ~ date, data = sextant_1998_2025_panache_monthly_anom)
p_value_sextant_1998  <- summary(model_sextant_1998)$coefficients[2, 4]
intercept_sextant_1998 <- coef(model_sextant_1998)[1]
slope_sextant_1998     <- coef(model_sextant_1998)[2]

# --- Graphique 2 : anomalie mensuelle ---
p2 <- ggplot(sextant_1998_2025_panache_monthly_anom, aes(x = date, y = panache_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(panache_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(panache_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_panache_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_panache_monthly_anom$panache_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_1998, 3), " + ", round(slope_sextant_1998, 7), " × x",
      "\np = ", ifelse(p_value_sextant_1998 < 0.001, "< 0.001", format(p_value_sextant_1998, digits = 3))
    ),
    hjust = 0, vjust = 1, size = 6, color = "grey20", fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de l'extension des panaches turbides",
    x       = NULL,
    y       = "Anomalie d'extension des panaches (km²)",
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 1998-2025"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 12, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 10, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 12, margin = margin(r = 10)),
    axis.text        = element_text(size = 11, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10)
  )

# --- Patchwork ---
(p1 / p2) +
  plot_annotation(
    title   = "Extension des panaches turbides — Sextant OC5",
    caption = "Source : Sextant OC5",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )

## SPM ----------------------------------------------------------------

# en échelle normale

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 1998 et 2025 avec sextant

model_sextant_1998 <- lm(mean_spm ~ date, data = sextant_1998_2025_SPM)
p_value_sextant_1998 <- summary(model_sextant_1998)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998 <- coef(model_sextant_1998)[1]
slope_sextant_1998 <- coef(model_sextant_1998)[2]

ggplot(data = sextant_1998_2025_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  geom_point(color = "red3", size = 1) +
  annotate(
    "text",
    x = max(sextant_1998_2025_SPM$date, na.rm = TRUE),
    y = max(sextant_1998_2025_SPM$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_sextant_1998, 3), " + ", round(slope_sextant_1998, 7), " * x",
      "\n", "p = ", ifelse(p_value_sextant_1998 < 0.001, "< 0.001", format(p_value_sextant_1998, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration en matière particulaire en suspension moyenne entre 1998 et 2025 avec le produit Sextant",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# en échelle log

# pour faire notre graph on doit transformer nos données "normales" en échelle log
# sauf qu'on a des valeurs NA et inférieures à 0

sextant_filtered <- sextant_1998_2025_SPM %>%
  filter(mean_spm > 0 & !is.na(mean_spm))

model_log_sextant_1998 <- lm(log10(mean_spm) ~ date, data = sextant_filtered)
p_value_log_sextant_1998 <- summary(model_log_sextant_1998)$coefficients[2, 4]

intercept_log_sextant_1998 <- coef(model_log_sextant_1998)[1]
slope_log_sextant_1998 <- coef(model_log_sextant_1998)[2]
p_value_log_sextant_1998 <- summary(model_log_sextant_1998)$coefficients[2, 4]

# Formater l'équation
equation_text_log <- paste0(
  "log10(y) = ", round(intercept_log_sextant_1998, 4),
  ifelse(sign(slope_log_sextant_1998) == 1, " + ", " - "),
  abs(round(slope_log_sextant_1998, 4)), " * x",
  "\n",  # Saut de ligne
  "p-value = ", format.pval(p_value_log_sextant_1998, digits = 3)
)




# ggplot(data = sextant_1998_2025_SPM, aes(x = date, y = mean_spm)) +
#   # geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
#   geom_line(color = "red3") +
#   labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 1998 et 2025 avec le produit Sextant",
#        x = "Date",
#        y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
#   theme_minimal() +
#   scale_x_date(
#     date_breaks = "1 year",
#     date_labels = "%Y"
#   ) +
#   scale_y_log10()

# Graphique
ggplot(data = sextant_filtered, aes(x = date, y = mean_spm)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  labs(
    title = "Evolution de la concentration en matière particulaire en suspension entre 1998 et 2025 (échelle log)",
    x = "Date",
    y = "Concentration moyenne (g/m³, échelle log)"
  ) +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_log10() +
  annotate(
    "text",
    x = as.Date("2010-01-01"),
    y = max(sextant_filtered$mean_spm, na.rm = TRUE) * 0.8,
    label = equation_text_log,
    hjust = 0,
    vjust = 1,
    size = 5,
    color = "black"
  )
# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec sextant contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(sextant_1998_2025_SPM$mean_spm, Y6442010_depuis_2000$débit)

sextant_1998_2025_SPM$scaled_mean_spm <- sextant_1998_2025_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

# en échelle normale
ggplot() +
  geom_line(
    data = Y6442010_depuis_2000,
    aes(x = date, y = débit, color = "Débit"), size = 0.5
  ) +
  geom_line(
    data = sextant_1998_2025_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM"), size = 0.5
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration moyenne en matière en suspension entre 1998 et 2025 avec le produit Sextant OC5",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# en échelle log

ggplot() +
  geom_line(
    data = Y6442010_depuis_2000,
    aes(x = date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = sextant_1998_2025_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_log10(
    name = "Débit (m³/s, log)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    sec.axis = sec_axis(
      ~ log10(. - adjust_factors$adjust) / log10(adjust_factors$diff),
      name = "Matière particulaire en suspension (g/m³, log)",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  ) +
  labs(
    title = "Débit du Var et concentration en matière en suspension (1998-2025, échelle log)",
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


# scatter plot ------------------------------------------------------------

ggplot(Var_SEXTANT_SPM, aes(x = débit, y = mean_spm)) +
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

Var_SEXTANT_SPM <- inner_join(Y6442010_depuis_2000, sextant_1998_2025_SPM, by = "date")

cor.test(Var_SEXTANT_SPM$débit, Var_SEXTANT_SPM$mean_spm, method = "spearman")

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


# décomposition X11 de l'aire des panaches ---------------------------------------

# pour ce faire on aggrège les données par mois et non par jour

# on ajoute au df la date avec l'année, le mois et le jour de l'année
SEXTANT_1998_2025_spm_95 <- SEXTANT_1998_2025_spm_95 |> 
  filter(date >= as.Date("2008-01-01"), date <= as.Date("2019-12-31")) |> 
  mutate(year = year(date),
         month = month(date),
         doy = yday(date))

# Agréger en mensuel
panache_mensuel <- SEXTANT_1998_2025_spm_95 |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(panache_mois = mean(aire_panache_km2, na.rm = TRUE))

# Créer la série temporelle sur la colonne débit uniquement
panache_ts <- ts(
  data      = panache_mensuel$panache_mois,  # ← juste la colonne
  start     = c(2008, 1),
  frequency = 12
)

# Appliquer X11
x11_result_panache <- seas(panache_ts, x11 = "")

# 3. Inspecter les résultats
summary(x11_result_panache)

# 4. Extraire les composantes
composantes_panache <- data.frame(
  date         = panache_mensuel$mois,
  observed     = as.numeric(original(x11_result_panache)),
  tendance     = as.numeric(trend(x11_result_panache)),
  saisonnalite = as.numeric(series(x11_result_panache, "d10")),  # facteurs saisonniers X11
  residus      = as.numeric(irregular(x11_result_panache))
)

# 5. Visualiser
composantes_long_panache <- composantes_panache |>
  pivot_longer(-date, names_to = "composante", values_to = "valeur") |>
  mutate(composante = factor(composante,
                             levels = c("observed", "tendance", "saisonnalite", "residus")))

# Graphique 1 — Signal brut + tendance
p1 <- ggplot(composantes_panache, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"), linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "steelblue", "Tendance" = "firebrick")) +
  labs(title = "a) Signal observé et tendance", x = NULL, y = "Aire des panaches (km²)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 2 — Saisonnalité
p2 <- ggplot(composantes_panache, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 0), ymax = 0), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(saisonnalite, 0)), fill = "chartreuse4", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière", x = NULL, y = "Aire des panaches (km²)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 3 — Résidus
p3 <- ggplot(composantes_panache, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)), fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)", x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# Assembler
(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 de l'aire des panaches — 2008–2019",
    subtitle = "Sextant OC5",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

# Moyenne de la saisonnalité par mois + min/max
saisonnalite_clim_panache <- composantes_panache |>
  mutate(month = month(date, label = TRUE, abbr = TRUE, locale = "fr_FR")) |>
  group_by(month) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE)
  )

# Plot
ggplot(saisonnalite_clim_panache, aes(x = month, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Saisonnalité X11 de l'aire des panaches",
    subtitle = "Moyenne mensuelle 2008–2019 (enveloppe = min/max)",
    x        = NULL,
    y        = "Facteur saisonnier"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# décomposition X11 de la concentration en chl ---------------------------------------

SEXTANT_1998_2025_chl_clean <- SEXTANT_1998_2025_chl_pixels |> 
  filter(analysed_chl_a >= 0) |> 
  filter(date >= as.Date("2008-01-01"), date <= as.Date("2019-12-31")) |> 
  mutate(year = year(date),
         month = month(date),
         doy = yday(date))

# Agréger en mensuel
chl_mensuel <- SEXTANT_1998_2025_chl_clean |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(CHL_mois = mean(analysed_chl_a, na.rm = TRUE))

# Créer la série temporelle sur la colonne débit uniquement
chl_ts <- ts(
  data      = chl_mensuel$CHL_mois,  # ← juste la colonne
  start     = c(2008, 1),
  frequency = 12
)

# Appliquer X11
x11_result_chl <- seas(chl_ts, x11 = "")

# 3. Inspecter les résultats
summary(x11_result_chl)

# 4. Extraire les composantes
composantes_chl <- data.frame(
  date         = chl_mensuel$mois,
  observed     = as.numeric(original(x11_result_chl)),
  tendance     = as.numeric(trend(x11_result_chl)),
  saisonnalite = as.numeric(series(x11_result_chl, "d10")),  # facteurs saisonniers X11
  residus      = as.numeric(irregular(x11_result_chl))
)

# 5. Visualiser
composantes_long_chl <- composantes_chl |>
  pivot_longer(-date, names_to = "composante", values_to = "valeur") |>
  mutate(composante = factor(composante,
                             levels = c("observed", "tendance", "saisonnalite", "residus")))

# Graphique 1 — Signal brut + tendance
p1 <- ggplot(composantes_chl, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"), linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "steelblue", "Tendance" = "firebrick")) +
  labs(title = "a) Signal observé et tendance", x = NULL, y = "Concentration en CHL (µg.L-1)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 2 — Saisonnalité
p2 <- ggplot(composantes_chl, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 0), ymax = 0), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(saisonnalite, 0)), fill = "chartreuse4", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière", x = NULL, y = "Concentration en CHL (µg.L-1)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 3 — Résidus
p3 <- ggplot(composantes_chl, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)), fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)", x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# Assembler
(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 de la concentration en chlorophylle a — 2008–2019",
    subtitle = "Sextant OC5",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

# Moyenne de la saisonnalité par mois + min/max
saisonnalite_clim_chl <- composantes_chl |>
  mutate(month = month(date, label = TRUE, abbr = TRUE, locale = "fr_FR")) |>
  group_by(month) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE)
  )

# Plot
ggplot(saisonnalite_clim_chl, aes(x = month, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "chartreuse3", alpha = 0.25) +
  geom_line(color = "chartreuse3", linewidth = 1) +
  geom_point(color = "chartreuse3", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Saisonnalité X11 de la concentration en chlorophylle a",
    subtitle = "Moyenne mensuelle 2008–2019 (enveloppe = min/max)",
    x        = NULL,
    y        = "Facteur saisonnier"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# décomposition X11 de la concentration en MES ---------------------------------------

SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_pixels |> 
  filter(analysed_spim >= 0) |> 
  filter(date >= as.Date("2008-01-01"), date <= as.Date("2019-12-31")) |> 
  mutate(year = year(date),
         month = month(date),
         doy = yday(date))

# Agréger en mensuel
spm_mensuel <- SEXTANT_1998_2025_spm_clean |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(spm_mois = mean(analysed_spim, na.rm = TRUE))

# Créer la série temporelle sur la colonne débit uniquement
spm_ts <- ts(
  data      = spm_mensuel$spm_mois,  # ← juste la colonne
  start     = c(2008, 1),
  frequency = 12
)

# Appliquer X11
x11_result_spm <- seas(spm_ts, x11 = "")

# 3. Inspecter les résultats
summary(x11_result_spm)

# 4. Extraire les composantes
composantes_spm <- data.frame(
  date         = spm_mensuel$mois,
  observed     = as.numeric(original(x11_result_spm)),
  tendance     = as.numeric(trend(x11_result_spm)),
  saisonnalite = as.numeric(series(x11_result_spm, "d10")),  # facteurs saisonniers X11
  residus      = as.numeric(irregular(x11_result_spm))
)

# 5. Visualiser
composantes_long_spm <- composantes_spm |>
  pivot_longer(-date, names_to = "composante", values_to = "valeur") |>
  mutate(composante = factor(composante,
                             levels = c("observed", "tendance", "saisonnalite", "residus")))

# Graphique 1 — Signal brut + tendance
p1 <- ggplot(composantes_spm, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"), linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "steelblue", "Tendance" = "firebrick")) +
  labs(title = "a) Signal observé et tendance", x = NULL, y = "Concentration en MES (g.m-3)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 2 — Saisonnalité
p2 <- ggplot(composantes_spm, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 0), ymax = 0), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(saisonnalite, 0)), fill = "red3", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière", x = NULL, y = "Concentration en MES (g.m-3)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 3 — Résidus
p3 <- ggplot(composantes_spm, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)), fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)", x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# Assembler
(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 de la concentration en MES — 2008–2019",
    subtitle = "Sextant OC5",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

# Moyenne de la saisonnalité par mois + min/max
saisonnalite_clim_spm <- composantes_spm |>
  mutate(month = month(date, label = TRUE, abbr = TRUE, locale = "fr_FR")) |>
  group_by(month) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE)
  )

# Plot
ggplot(saisonnalite_clim_spm, aes(x = month, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "red3", alpha = 0.25) +
  geom_line(color = "red3", linewidth = 1) +
  geom_point(color = "red3", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Saisonnalité X11 de la concentration en MES",
    subtitle = "Moyenne mensuelle 2008–2019 (enveloppe = min/max)",
    x        = NULL,
    y        = "Facteur saisonnier"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# décomposition X11 du débit ---------------------------------------

load("data/Hydro France/Y6442010_depuis_2000.Rdata")

# on choisit une période remplie donc entre 2008 et 2019

Y6442010_2008_2019 <- Y6442010_depuis_2000 |> 
  filter(date >= as.Date("2008-01-01"), date <= as.Date("2019-12-31"))

sum(is.na(Y6442010_2008_2019))

# Localiser et caractériser les trous
Y6442010_2008_2019 |>
  mutate(est_na = is.na(débit)) |>
  filter(est_na) |>
  mutate(
    groupe = cumsum(c(1, diff(as.numeric(date)) > 1))
  ) |>
  group_by(groupe) |>
  summarise(
    debut     = min(date),
    fin       = max(date),
    n_jours   = n()
  ) |>
  arrange(debut)

Y6442010_2008_2019 <- Y6442010_2008_2019 |>
  arrange(date) |>
  mutate(
    debit_interp = na.approx(débit, x = date, na.rm = FALSE)
  )

# Vérifier qu'il ne reste plus de NA
sum(is.na(Y6442010_2008_2019$debit_interp))

# Visualiser pour vérifier que l'interpolation est cohérente
Y6442010_2008_2019 |>
  mutate(est_interpole = is.na(débit) & !is.na(debit_interp)) |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = debit_interp), color = "steelblue", linewidth = 0.5) +
  geom_point(
    data = ~ filter(.x, est_interpole),
    aes(y = debit_interp),
    color = "red", size = 1.5
  ) +
  labs(
    title    = "Débit du Var interpolé — 2008–2019",
    subtitle = "Points rouges = valeurs interpolées",
    x        = NULL,
    y        = "Débit (m³/s)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# Agréger en mensuel
debit_mensuel <- Y6442010_2008_2019 |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(debit = mean(debit_interp, na.rm = TRUE))

# Vérifier qu'il n'y a plus de NA
sum(is.na(debit_mensuel$debit))

# Créer la série temporelle sur la colonne débit uniquement
debit_ts <- ts(
  data      = debit_mensuel$debit,  # ← juste la colonne
  start     = c(2008, 1),
  frequency = 12
)

# Appliquer X11
x11_result_debit <- seas(debit_ts, x11 = "")

# 3. Inspecter les résultats
summary(x11_result_debit)

composantes_debit <- data.frame(
  date         = debit_mensuel$mois,
  observed     = as.numeric(original(x11_result_debit)),
  tendance     = as.numeric(trend(x11_result_debit)),
  saisonnalite = as.numeric(series(x11_result_debit, "d10")),  # facteurs saisonniers X11
  residus      = as.numeric(irregular(x11_result_debit))
)

# 5. Visualiser
composantes_long_debit <- composantes_debit |>
  pivot_longer(-date, names_to = "composante", values_to = "valeur") |>
  mutate(composante = factor(composante,
                             levels = c("observed", "tendance", "saisonnalite", "residus")))

# Graphique 1 — Signal brut + tendance
p1 <- ggplot(composantes_debit, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"), linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "steelblue", "Tendance" = "firebrick")) +
  labs(title = "a) Signal observé et tendance", x = NULL, y = "Débit (m³/s)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 2 — Saisonnalité
p2 <- ggplot(composantes_debit, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 0), ymax = 0), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(saisonnalite, 0)), fill = "chartreuse4", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière", x = NULL, y = "Débit (m³/s)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 3 — Résidus
p3 <- ggplot(composantes_debit, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)), fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)", x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# Assembler
(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 du débit du Var — 2008–2019",
    subtitle = "Station Y6442010 — Agrégation mensuelle",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

# Moyenne de la saisonnalité par mois + min/max
saisonnalite_clim_débit <- composantes_debit |>
  mutate(month = month(date, label = TRUE, abbr = TRUE, locale = "fr_FR")) |>
  group_by(month) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE)
  )

# Plot
ggplot(saisonnalite_clim_débit, aes(x = month, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Saisonnalité X11 du débit du Var",
    subtitle = "Moyenne mensuelle 2008–2019 (enveloppe = min/max)",
    x        = NULL,
    y        = "Facteur saisonnier"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# décomposition X11 des précipitations ---------------------------------------

Wind_T <- read.csv("~/Vent/data/Q_06_previous-1950-2024_RR-T-Vent.csv", 
                   header = TRUE, sep = ";")

Wind_T <- Wind_T |> 
  filter(NUM_POSTE == "6088001") |> 
  select("LAT", "LON", "NUM_POSTE", "FFM", "DXY", "HXI", "RR", "TM")

Wind_T <- Wind_T |> 
  mutate(date = seq(as.Date("1950-01-01"), as.Date("2024-12-31"), by = "day"))

Wind_T <- Wind_T |> 
  mutate(
    annee = year(date),
    mois = month(date)
  )

Wind_T <- Wind_T |> 
  filter(date >= "2008-01-01", date <= "2019-12-31")

# Agréger en mensuel
pluie_mois <- Wind_T |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(rr_mois = mean(RR, na.rm = TRUE))

# Créer la série temporelle sur la colonne débit uniquement
pluie_ts <- ts(
  data      = pluie_mois$rr_mois,  # ← juste la colonne
  start     = c(2008, 1),
  frequency = 12
)

# Appliquer X11
x11_result_pluie <- seas(pluie_ts, x11 = "")

# 3. Inspecter les résultats
summary(x11_result_pluie)

# 4. Extraire les composantes
composantes_pluie <- data.frame(
  date         = pluie_mois$mois,
  observed     = as.numeric(original(x11_result_pluie)),
  tendance     = as.numeric(trend(x11_result_pluie)),
  saisonnalite = as.numeric(series(x11_result_pluie, "d10")),  # facteurs saisonniers X11
  residus      = as.numeric(irregular(x11_result_pluie))
)

# 5. Visualiser
composantes_long_pluie <- composantes_pluie |>
  pivot_longer(-date, names_to = "composante", values_to = "valeur") |>
  mutate(composante = factor(composante,
                             levels = c("observed", "tendance", "saisonnalite", "residus")))

# Graphique 1 — Signal brut + tendance
p1 <- ggplot(composantes_pluie, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"), linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "steelblue", "Tendance" = "firebrick")) +
  labs(title = "a) Signal observé et tendance", x = NULL, y = "Précipitations (mm)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 2 — Saisonnalité
p2 <- ggplot(composantes_pluie, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 0), ymax = 0), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(saisonnalite, 0)), fill = "red3", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière", x = NULL, y = "Précipitations (mm)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# Graphique 3 — Résidus
p3 <- ggplot(composantes_pluie, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1), fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)), fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)", x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# Assembler
(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 des précipitations — 2008–2019",
    subtitle = "Sextant OC5",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

# Moyenne de la saisonnalité par mois + min/max
saisonnalite_clim_pluie <- composantes_pluie |>
  mutate(month = month(date, label = TRUE, abbr = TRUE, locale = "fr_FR")) |>
  group_by(month) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE)
  )

# Plot
ggplot(saisonnalite_clim_pluie, aes(x = month, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "red3", alpha = 0.25) +
  geom_line(color = "red3", linewidth = 1) +
  geom_point(color = "red3", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Saisonnalité X11 des précipitations",
    subtitle = "Moyenne mensuelle 2008–2019 (enveloppe = min/max)",
    x        = NULL,
    y        = "Facteur saisonnier"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())


# superposition -----------------------------------------------------------

# ensuite on veut supperposer les quatres
ggplot() +
  # Aire des panaches
  geom_ribbon(data = saisonnalite_clim_panache,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "#00897B", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_panache,
            aes(x = month, y = mean_sais, color = "Aire des panaches", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_panache,
             aes(x = month, y = mean_sais, color = "Aire des panaches", group = 1),
             size = 2.5) +
  # Débit du Var
  geom_ribbon(data = saisonnalite_clim_débit,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "blue", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_débit,
            aes(x = month, y = mean_sais, color = "Débit du Var", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_débit,
             aes(x = month, y = mean_sais, color = "Débit du Var", group = 1),
             size = 2.5) +
  # Chlorophylle a
  geom_ribbon(data = saisonnalite_clim_chl,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "chartreuse3", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_chl,
            aes(x = month, y = mean_sais, color = "Concentration en chlorophylle a", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_chl,
             aes(x = month, y = mean_sais, color = "Concentration en chlorophylle a", group = 1),
             size = 2.5) +
  # ── AJOUT : Concentration en MES ──────────────────────────────────────────
  geom_ribbon(data = saisonnalite_clim_spm,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "red3", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_spm,
            aes(x = month, y = mean_sais, color = "Concentration en MES", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_spm,
             aes(x = month, y = mean_sais, color = "Concentration en MES", group = 1),
             size = 2.5) +
  # ──────────────────────────────────────────────────────────────────────────
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Aire des panaches"              = "#00897B",
      "Débit du Var"                   = "blue",
      "Concentration en chlorophylle a" = "chartreuse3",
      "Concentration en MES"           = "red3"   # ← ajout
    ),
    guide = guide_legend(override.aes = list(linewidth = 1.5, size = 3))
  ) +
  labs(
    title    = "Saisonnalité X11 — Débit du Var, aire des panaches turbides, MES et chlorophylle a",
    subtitle = "Moyenne mensuelle 2008–2019 · enveloppe = min/max interannuel",
    x        = NULL,
    y        = "Facteur saisonnier",
    color    = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey93"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(color = "grey50", size = 10, margin = margin(b = 10)),
    legend.position   = "top",
    legend.text       = element_text(size = 11),
    axis.text         = element_text(color = "grey30"),
    axis.text.x       = element_text(size = 11)
  )

ggplot() +
  # Aire des panaches
  geom_ribbon(data = saisonnalite_clim_panache,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "#00897B", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_panache,
            aes(x = month, y = mean_sais, color = "Aire des panaches", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_panache,
             aes(x = month, y = mean_sais, color = "Aire des panaches", group = 1),
             size = 2.5) +
  # Débit du Var
  geom_ribbon(data = saisonnalite_clim_débit,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "blue", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_débit,
            aes(x = month, y = mean_sais, color = "Débit du Var", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_débit,
             aes(x = month, y = mean_sais, color = "Débit du Var", group = 1),
             size = 2.5) +
  # Chlorophylle a
  geom_ribbon(data = saisonnalite_clim_chl,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "chartreuse3", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_chl,
            aes(x = month, y = mean_sais, color = "Concentration en chlorophylle a", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_chl,
             aes(x = month, y = mean_sais, color = "Concentration en chlorophylle a", group = 1),
             size = 2.5) +
  # Concentration en MES
  geom_ribbon(data = saisonnalite_clim_spm,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "red3", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_spm,
            aes(x = month, y = mean_sais, color = "Concentration en MES", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_spm,
             aes(x = month, y = mean_sais, color = "Concentration en MES", group = 1),
             size = 2.5) +
  # ── AJOUT : Précipitations ─────────────────────────────────────────────────
  geom_ribbon(data = saisonnalite_clim_pluie,
              aes(x = month, ymin = min_sais, ymax = max_sais, group = 1),
              fill = "steelblue", alpha = 0.15) +
  geom_line(data = saisonnalite_clim_pluie,
            aes(x = month, y = mean_sais, color = "Précipitations", group = 1),
            linewidth = 1.1) +
  geom_point(data = saisonnalite_clim_pluie,
             aes(x = month, y = mean_sais, color = "Précipitations", group = 1),
             size = 2.5) +
  # ──────────────────────────────────────────────────────────────────────────
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Aire des panaches"               = "#00897B",
      "Débit du Var"                    = "blue",
      "Concentration en chlorophylle a" = "chartreuse3",
      "Concentration en MES"            = "red3",
      "Précipitations"                  = "steelblue"  # ← ajout
    ),
    guide = guide_legend(override.aes = list(linewidth = 1.5, size = 3))
  ) +
  labs(
    title    = "Saisonnalité X11 — Débit du Var, panaches turbides, MES, chlorophylle a et précipitations",
    subtitle = "Moyenne mensuelle 2008–2019 · enveloppe = min/max interannuel",
    x        = NULL,
    y        = "Facteur saisonnier",
    color    = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey93"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(color = "grey50", size = 10, margin = margin(b = 10)),
    legend.position   = "top",
    legend.text       = element_text(size = 11),
    axis.text         = element_text(color = "grey30"),
    axis.text.x       = element_text(size = 11)
  )


# spatial climatology -----------------------------------------------------

## MES and panache extension -----------------------------------------------

SEXTANT_1998_2025_spm_pixels <- SEXTANT_1998_2025_spm_pixels |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  ) |> 
  filter(date >= as.Date("1998-01-01"), date <= as.Date("2025-12-31"))

# ── Climatologie spatiale mensuelle (moyenne par pixel et par mois) ──
clim_spatiale_spm_month_sextant <- SEXTANT_1998_2025_spm_pixels |>
  filter(analysed_spim >= 0) |>
  mutate(month = month(date)) |>
  group_by(lon, lat, month) |>
  summarise(
    mean_spm  = mean(analysed_spim, na.rm = TRUE),
    median_spm = median(analysed_spim, na.rm = TRUE),
    sd_spm    = sd(analysed_spim, na.rm = TRUE),
    .groups   = "drop"
  )

# ── Climatologie spatiale annuelle (moyenne par pixel et par année) ──
clim_spatiale_spm_year <- SEXTANT_2008_2019_spm_pixels |>
  filter(analysed_spim >= 0) |>
  mutate(year = year(date)) |>
  group_by(lon, lat, year) |>
  summarise(
    mean_spm = mean(analysed_spim, na.rm = TRUE),
    .groups  = "drop"
  )

# ── Climatologie spatiale globale (moyenne par pixel sur toute la période) ──
clim_spatiale_spm_total <- SEXTANT_2008_2019_spm_pixels |>
  filter(analysed_spim >= 0) |>
  group_by(lon, lat) |>
  summarise(
    mean_spm   = mean(analysed_spim, na.rm = TRUE),
    median_spm = median(analysed_spim, na.rm = TRUE),
    sd_spm     = sd(analysed_spim, na.rm = TRUE),
    .groups    = "drop"
  )

# plot mensuel de la concentration en MES
ggplot(clim_spatiale_spm_month_sextant, aes(x = lon, y = lat, fill = mean_spm)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_spm_month_sextant$lon),
    ylim = range(clim_spatiale_spm_month_sextant$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("MES (g. m"^{-3}*")"),
    option   = "turbo",
    na.value = "white",
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10")
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de la concentration en matières en suspension",
    subtitle = "Période de référence : 1998-2025 · Produit Sextant OC5",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,
    barheight      = 0.8,
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )



# plot mensuel de l'erreur standard à la concentration en MES
ggplot(clim_spatiale_spm_month_sextant, aes(x = lon, y = lat, fill = sd_spm)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_spm_month_sextant$lon),
    ylim = range(clim_spatiale_spm_month_sextant$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("MES (g. m"^{-3}*")"),
    option   = "turbo",
    na.value = "white",
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10")
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de l'erreur standard à la concentration en MES",
    subtitle = "Période de référence : 1998-2025 · Produit Sextant OC5",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,
    barheight      = 0.8,
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )

## chl -----------------------------------------------

# combien de valeurs négatives
sum(SEXTANT_1998_2025_chl_pixels$analysed_chl_a < 0, na.rm = TRUE)
# [1] 19808

# supprimer seulement les valeurs négatives
SEXTANT_1998_2025_chl_clean <- SEXTANT_1998_2025_chl_pixels |>
  filter(analysed_chl_a >= 0, analysed_chl_a <= 20 | is.na(analysed_chl_a))

SEXTANT_1998_2025_chl_clean <- SEXTANT_1998_2025_chl_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  ) |> 
  filter(date >= as.Date("1998-01-01"), date <= as.Date("2025-12-31"))

# ── Climatologie spatiale mensuelle (moyenne par pixel et par mois) ──
clim_spatiale_chl_month <- SEXTANT_1998_2025_chl_clean |>
  mutate(month = month(date)) |>
  group_by(lon, lat, month) |>
  summarise(
    mean_chl  = mean(analysed_chl_a, na.rm = TRUE),
    median_chl = median(analysed_chl_a, na.rm = TRUE),
    sd_chl    = sd(analysed_chl_a, na.rm = TRUE),
    .groups   = "drop"
  )

# plot mensuel
ggplot(clim_spatiale_chl_month, aes(x = lon, y = lat, fill = mean_chl)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_chl_month$lon),
    ylim = range(clim_spatiale_chl_month$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("Chl a (µg. L"^{-1}*")"),
    option   = "plasma",
    na.value = "white",
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10")
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de la concentration en chlorophylle a",
    subtitle = "Période de référence : 1998–2025 · Produit Sextant OC5",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,
    barheight      = 0.8,
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )

# sd chl
ggplot(clim_spatiale_chl_month, aes(x = lon, y = lat, fill = sd_chl)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_chl_month$lon),
    ylim = range(clim_spatiale_chl_month$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("Chl a (µg. L"^{-1}*")"),
    option   = "plasma",
    na.value = "white",
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10")
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de l'erreur standard à la concentration en chlorophylle a",
    subtitle = "Période de référence : 1998–2025 · Produit Sextant OC5",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,
    barheight      = 0.8,
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )



