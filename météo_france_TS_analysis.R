# météo_france_TS_analysis

# pathway : ~/Satellite_analysis/météo_france_TS_analysis/

# données météo france relevées in situ à l'aéroport de Nice

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
library(ggpubr)
library(circular)
library(Kendall)
library(openair)

# load --------------------------------------------------------------------

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

# Météo France data --------------------------------------------------------

Wind_T <- read.csv("~/Vent/data/Q_06_previous-1950-2024_RR-T-Vent.csv", 
                   header = TRUE, sep = ";")

Wind_T <- Wind_T |> 
  filter(NUM_POSTE == "6088001") |> 
  select("LAT", "LON", "NUM_POSTE", "AAAAMMJJ", "FFM", "DXY", "HXI", "RR", "TM") |> 
  mutate(date = as.Date(as.character(AAAAMMJJ), format = "%Y%m%d"))

Wind_T <- Wind_T |> 
  mutate(
    annee = year(date),
    mois = month(date)
  )

Wind_T <- Wind_T |> 
  mutate(date = seq(as.Date("1950-01-01"), as.Date("2024-12-31"), by = "day"))

Wind_T <- Wind_T |> 
  filter(date >= "1998-01-01", date <= "2024-12-31")

# other data --------------------------------------------------------------

load("data/Hydro France/All_debit.Rdata")

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

# ya t'il eu une baisse dans la vitesse des vents d'est

model_wind <- lm(FFM ~ date, data = East)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(East, aes(x = date, y = FFM)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(East$date, na.rm = TRUE),
    y = max(East$FFM, na.rm = TRUE) * 0.95,
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
    title = "Évolution de la vitesse du vent d'Est près de Nice (1998-2025)",
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

# y a t'il eu une baisse dans la vitesse des vents de Nord Ouest

model_wind <- lm(FFM ~ date, data = North_West)
p_value_wind <- summary(model_wind)$coefficients[2, 4]  # p-value pour la pente
intercept_wind <- coef(model_wind)[1]
slope_wind <- coef(model_wind)[2]

ggplot(North_West, aes(x = date, y = FFM)) +
  geom_point(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(North_West$date, na.rm = TRUE),
    y = max(North_West$FFM, na.rm = TRUE) * 0.95,
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
    title = "Évolution de la vitesse du vent de Nord Ouest près de Nice (2008-2020)",
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

# ya t'il plus/moins de vent de Nord Ouest qu'avant
# ya t'il plus/moins de vent d'e Nord Ouest'Est qu'avant


# 1. Assigner les secteurs ------------------------------------------------
Wind_T <- Wind_T |>
  mutate(
    # 8 secteurs de 45°
    secteur = case_when(
      DXY >= 337.5 | DXY < 22.5   ~ "N",
      DXY >= 22.5  & DXY < 67.5   ~ "NE",
      DXY >= 67.5  & DXY < 112.5  ~ "E",
      DXY >= 112.5 & DXY < 157.5  ~ "SE",
      DXY >= 157.5 & DXY < 202.5  ~ "S",
      DXY >= 202.5 & DXY < 247.5  ~ "SO",
      DXY >= 247.5 & DXY < 292.5  ~ "O",
      DXY >= 292.5 & DXY < 337.5  ~ "NO"
    ),
    # Classification offshore/onshore comme Gangloff et al.
    # Offshore = Mistral (N) + Tramontane (NO-N) = 295-15°
    # Onshore  = vents de mer (E à SE)           = 80-160°
    regime = case_when(
      (DXY >= 295 | DXY <= 15)    ~ "Offshore (Mistral/Tramontane)",
      (DXY >= 80  & DXY <= 160)   ~ "Onshore",
      TRUE                         ~ "Autre"
    ),
    annee    = year(date),
    mois     = month(date),
    periode  = ifelse(annee <= 2013, "2008–2013", "2014–2019")
  )

# 1. Vitesse moyenne annuelle par secteur
vitesse_secteur <- Wind_T |>
  group_by(annee, secteur) |>
  summarise(vitesse_moy = mean(FFM, na.rm = TRUE), .groups = "drop")  # adapte FXY

# 2. Mann-Kendall + Theil-Sen par secteur
resultats_vitesse <- vitesse_secteur |>
  group_by(secteur) |>
  arrange(annee) |>
  summarise(
    mk_pval     = mk.test(vitesse_moy)$p.value,
    mk_tau      = mk.test(vitesse_moy)$statistic,
    pente_an    = sens.slope(vitesse_moy)$estimates,  # m/s par an
    .groups = "drop"
  ) |>
  mutate(
    significatif = ifelse(mk_pval < 0.05, "*", "ns"),
    mk_pval_fmt  = ifelse(mk_pval < 2.2e-16, "< 2.2×10⁻¹⁶", round(mk_pval, 4))
  )

print(resultats_vitesse)

# 2. Proportion annuelle par secteur --------------------------------------
freq_secteur <- Wind_T |>
  group_by(annee, secteur) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(annee) |>
  mutate(freq = n / sum(n) * 100) |>
  ungroup()

# 3. Proportion annuelle offshore/onshore ---------------------------------
prop_regime <- Wind_T |>
  group_by(annee, regime) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(annee) |>
  mutate(prop = n / sum(n))

# 4. Test de tendance de Mann-Kendall sur chaque régime ------------------

resultats_secteur <- freq_secteur |>
  group_by(secteur) |>
  arrange(annee) |>
  summarise(
    mk_pval     = mk.test(freq)$p.value,
    mk_tau      = mk.test(freq)$statistic,
    pente_an    = sens.slope(freq)$estimates,  # % par an
    .groups = "drop"
  ) |>
  mutate(
    significatif = ifelse(mk_pval < 0.05, "*", "ns"),
    mk_pval_fmt  = ifelse(mk_pval < 2.2e-16, "< 2.2×10⁻¹⁶", round(mk_pval, 4))
  )

print(resultats_secteur)

# 5. Test du chi² : comparaison des proportions entre deux périodes ------
table_contingence <- Wind_T |>
  count(periode, secteur) |>
  pivot_wider(names_from = secteur, values_from = n, values_fill = 0) |>
  column_to_rownames("periode")

chisq.test(table_contingence)

# 6. Visualisation — proportion par régime au cours du temps -------------
ggplot(prop_regime |> filter(regime != "Autre"),
       aes(x = annee, y = prop, color = regime)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c(
    "Offshore (Mistral/Tramontane)" = "#2166ac",
    "Onshore"                       = "#d6604d"
  )) +
  labs(
    title    = "Évolution de la proportion des régimes de vent — Nice 2008–2019",
    subtitle = "Test de Mann-Kendall appliqué sur chaque régime",
    x        = NULL,
    y        = "Proportion annuelle",
    color    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey50", size = 10),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

# 7. Rose des vents comparative — début vs fin de série ------------------

windRose(
  Wind_T,
  ws   = "FFM",
  wd   = "DXY",
  type = "periode",
  cols = "YlOrRd",
  main = "Rose des vents — Nice 2008–2019"
)



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

p1 <- ggwindrose(
  speed         = speed,
  direction     = direction,
  n_directions  = 8,
  n_speeds      = 5,
  speed_cuts    = NA,
  col_pal       = "GnBu",
  legend_title  = "Vent (m/s)",
  calm_wind     = 0,
  n_col         = 1,
  facet         = NULL,
  plot_title    = "Direction et vitesse du vent près de Nice",
  stack_reverse = TRUE) +
  labs(subtitle = "1998–2025") +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 14, hjust = 0.5, color = "grey50"),
    plot.caption    = element_text(size = 10, color = "grey50", hjust = 0),
    axis.text       = element_text(size = 13),
    legend.position = "bottom",              # ← ici
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 13)
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
  group_by(doy) |> 
  summarise(wind_doy_clim = mean(FFM, na.rm = TRUE),
            wind_doy_clim_std = sd(FFM, na.rm = TRUE))

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

# Anomalies journalières
Wind_daily_anom <- Wind_2015_2024_TS |> 
  summarise(mean_FFM = mean(FFM, na.rm = TRUE), .by = c("date", "year", "month", "doy")) |> 
  left_join(Wind_1991_2020_climatology_day, by = "doy") |> 
  mutate(
    wind_daily_anomaly     = mean_FFM - wind_doy_clim,
    wind_daily_anomaly_std = wind_daily_anomaly / wind_doy_clim_std  # anomalie normalisée
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

Wind_2000_2024_TS <- Wind_T |> 
  filter(date >= as.Date("2000-01-01"), date <= as.Date("2025-12-31")) |> 
  mutate(year = year(date), 
         month = month(date), 
         doy = yday(date))

Wind_2000_2024_climatology <- Wind_2000_2024_TS %>%
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

# # Calculer les proportions par mois
# proportions_par_mois <- Wind_2000_2024_climatology %>%
#   filter(!is.na(Direction)) %>%  # Exclure les NA
#   group_by(month, Direction) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   mutate(Proportion = n / sum(n))  # Proportion par mois

proportions_par_mois <- Wind_2000_2024_climatology %>%
  filter(!is.na(Direction)) %>%
  group_by(month, Direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(month) %>%                        # ← regrouper par mois
  mutate(Proportion = n / sum(n)) %>%        # ← proportion dans chaque mois
  ungroup()

## plotting ----------------------------------------------------------------

# climatologie mensuelle
model_wind_month <- lm(wind_month_clim ~ month, data = Wind_1991_2020_climatology_month)
p_value_wind_month <- summary(model_wind_month)$coefficients[2, 4]  # p-value pour la pente
intercept_wind_month <- coef(model_wind_month)[1]
slope_wind_month <- coef(model_wind_month)[2]

ggplot(Wind_1991_2020_climatology_month, aes(x = month, y = wind_month_clim)) +
  geom_line(color = "#4A90D9", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, 
              color = "#2C3E7A", fill = "#4A90D9", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(Wind_1991_2020_climatology_month$month, na.rm = TRUE),
    y = max(Wind_1991_2020_climatology_month$wind_month_clim, na.rm = TRUE) * 0.95,
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

# mettre la colonne month en mois
Wind_1991_2020_climatology_month$month <- factor(Wind_1991_2020_climatology_month$month, levels = 1:12, labels = month.abb)

# ggplot de la climatologie de la la vitesse du vent par mois
ggplot(Wind_1991_2020_climatology_month, aes(x = month, y = wind_month_clim)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = wind_month_clim - wind_month_clim_std, 
                    ymax = wind_month_clim + wind_month_clim_std), width = 0.3) +
  scale_x_discrete() +   # ← discrete car month est un factor
  labs(
    title = "Climatologie de la vitesse du vent (1991–2020)",
    x     = "Mois", 
    y     = "Vitesse du vent (m/s)"
  ) +
  theme_bw()

# Boxplot de la vitesse du vent par mois
ggplot(Wind_1991_2020_climatology_month, aes(x = month, y = wind_month_clim, fill = month)) +
  geom_boxplot() +
  labs(
    title = "Climatologie de la vitesse du vent (1991 - 2020)",
    x = "Mois",
    y = "Vitesse du vent (m/s)"
  ) +
  scale_x_discrete(labels = month.abb) +  # Affiche les abréviations des mois
  theme_bw() +
  theme(legend.position = "none")  # Masque la légende si elle n'est pas nécessaire

# ggplot de la climatologie de la de la vitesse du vent par jour
# Jours correspondant au milieu de chaque mois (année non bissextile)
mois_breaks <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)

ggplot(Wind_1991_2020_climatology_day, aes(x = doy, y = wind_doy_clim)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = wind_doy_clim - wind_doy_clim_std, 
                    ymax = wind_doy_clim + wind_doy_clim_std), 
                width = 0.3) +
  scale_x_continuous(
    breaks = mois_breaks,
    labels = month.abb
  ) +
  labs(
    title = "Climatologie de la vitesse du vent moyenne journalière (1991–2020)",
    x     = NULL,
    y     = "Vitesse du vent (m/s)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13)  # ← augmentez cette valeur
  )

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

# anomalie à la climatologie mensuelle de la vitesse du vent sur la période 2015 - 2024
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
    size  = 8, color = "grey20", fontface = "italic", family = "serif"
  ) +
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  
  labs(
    title    = "Anomalie mensuelle à la climatologie de la vitesse du vent (2015–2024)",
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


# anomalie à la climatologie mensuelle de la vitesse du vent sur la période 2015 - 2024
model_anom <- lm(wind_daily_anomaly ~ date, data = Wind_daily_anom)
p_value_anom <- summary(model_anom)$coefficients[2, 4]
slope_anom <- coef(model_anom)[2]
intercept_anom <- coef(model_anom)[1]

ggplot(Wind_daily_anom, aes(x = date, y = wind_daily_anomaly, fill = wind_daily_anomaly > 0)) +
  
  geom_col(alpha = 0.8, width = 3) +  # width en jours
  geom_smooth(
    aes(x = date, y = wind_daily_anomaly),
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
    x     = min(Wind_daily_anom$date, na.rm = TRUE),
    y     = max(Wind_daily_anom$wind_daily_anomaly, na.rm = TRUE) * 0.95,
    label = paste0("y = ", round(intercept_anom, 3), " + ", round(slope_anom, 7), " × x",
                   "\np = ", ifelse(p_value_anom < 0.001, "< 0.001",
                                    format(p_value_anom, scientific = TRUE, digits = 3))),
    hjust = 0, vjust = 1,
    size  = 8, color = "grey20", fontface = "italic", family = "serif"
  ) +
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  
  labs(
    title    = "Anomalie journalière à la climatologie de vitesse du vent (2015–2024)",
    subtitle = "Par rapport à la climatologie journalière 1991–2020",
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
  labs(title = "Proportions des directions de vent par mois (1998-2024)",
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

p2 <- ggplot(proportions_par_mois, aes(x = month, y = Proportion, fill = Direction)) +
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
       title = "Proportions mensuelles des directions du vent",
       subtitle = "1998-2025") +
  # Thème épuré et professionnel
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 13),
    legend.box.background = element_rect(linewidth = 0.5, color = "black"),
    plot.caption = element_text(size = 9, color = "grey50"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Afficher le graphique
print(p2)

# patchwork
(p1 | p2) + 
  plot_layout(widths = c(1, 1.9)) &  # p2 prend plus de largeur
plot_annotation(
    # title      = "Rose des vents et climatologie saisonnière",
    caption    = "Source : Archives Météo France",
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    # legend.position = "right",
    theme      = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0),
      plot.margin = margin(5, 5, 5, 0),
      theme(legend.position = "bottom"))
  )

p1 <- p1 + theme(legend.position = "bottom")
p2 <- p2 + theme(legend.position = "bottom")

(p1 | p2) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    caption    = "Source : Archives Météo France",
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    theme      = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 14, color = "grey50", hjust = 0)
    )
  )

-## température --------------------------------------------------------------

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
  group_by(doy) %>%
  summarise(temp_doy_clim = mean(TM, na.rm = TRUE),
            temp_doy_clim_std = sd(TM, na.rm = TRUE))

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

# Anomalies journalières
Temp_daily_anom <- Temp_2015_2024_TS |> 
  summarise(mean_TM = mean(TM, na.rm = TRUE), .by = c("date", "year", "month", "doy")) |> 
  left_join(Temp_1991_2020_climatology_day, by = "doy") |> 
  mutate(
    temp_daily_anomaly     = mean_TM - temp_doy_clim,
    temp_daily_anomaly_std = temp_daily_anomaly / temp_doy_clim_std  # anomalie normalisée
  )

## plotting temperature climatology ----------------------------------------

# plotting de la climatologie mensuelle
ggplot(Temp_1991_2020_climatology_month, aes(x = month, y = temp_month_clim)) +
  geom_bar(stat = "identity", fill = "orangered") +
  geom_errorbar(aes(ymin = temp_month_clim - temp_month_clim_std, 
                    ymax = temp_month_clim + temp_month_clim_std), width = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Climatologie de la température moyenne (1991 - 2020)",
       x = "Mois", y = "Température (en °C)") +
  theme_bw()

# plotting de la climatologie journalière
# Jours correspondant au milieu de chaque mois (année non bissextile)
mois_breaks <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)

ggplot(Temp_1991_2020_climatology_day, aes(x = doy, y = temp_doy_clim)) +
  geom_bar(stat = "identity", fill = "orangered") +
  geom_errorbar(aes(ymin = temp_doy_clim - temp_doy_clim_std, 
                    ymax = temp_doy_clim + temp_doy_clim_std), 
                width = 0.3) +
  scale_x_continuous(
    breaks = mois_breaks,
    labels = month.abb
  ) +
  labs(
    title = "Climatologie de la température moyenne journalière (1991–2020)",
    x     = NULL,
    y     = "Température (°C)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13)  # ← augmentez cette valeur
  )

# plotting de l'anomalie à la climatologie mensuelle
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

# plotting de l'anomalie à la climatologie journalière
model_anom_temp <- lm(temp_daily_anomaly ~ date, data = Temp_daily_anom)
p_value_anom_temp <- summary(model_anom_temp)$coefficients[2, 4]
slope_anom_temp <- coef(model_anom_temp)[2]
intercept_anom_temp <- coef(model_anom_temp)[1]

ggplot(Temp_daily_anom, aes(x = date, y = temp_daily_anomaly, fill = temp_daily_anomaly > 0)) +
  
  geom_col(alpha = 0.8, width = 3) +  # width en jours
  geom_smooth(
    aes(x = date, y = temp_daily_anomaly),
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
    x     = min(Temp_daily_anom$date, na.rm = TRUE),
    y     = max(Temp_daily_anom$temp_daily_anomaly, na.rm = TRUE) * 0.95,
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
    title    = "Anomalie journalière de la température (2015–2024)",
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
prop.table(table(data_filtered$wind_sector)) * 100

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



# précipitations ----------------------------------------------------------

Wind_T <- Wind_T |> 
  mutate(
    annee = year(date),
    mois = month(date)
  )

# Option 2 — Totaux mensuels (plus de points, tendance saisonnière possible)
pluie_mensuelle <- Wind_T |>
  group_by(annee, mois) |>
  summarise(total = sum(RR, na.rm = TRUE), .groups = "drop") |>
  arrange(annee, mois)

# Conversion en série temporelle mensuelle
pluie_ts <- ts(pluie_mensuelle$total, 
               start = c(min(pluie_mensuelle$annee), 1), 
               frequency = 12)

smk.test(pluie_ts)

# Sauvegarder le résultat MK
mk_pluie <- mk.test(pluie_mensuelle$total)
mk_pval  <- mk_pluie$p.value
mk_label <- ifelse(mk_pval < 2.2e-16,
                   "Tendance pluie : p < 2.2×10⁻¹⁶",
                   paste0("Tendance pluie : p = ", round(mk_pval, 4),
                          ifelse(mk_pval < 0.05, " *", " (ns)")))

# les précipitations ont elles baisser ?
# oui entre 2009 et 2019 mais non significativement
# non entre 2000 et 2024, non significativement
# oui entre 1997 et 2024 mais non significativement

# Préparer les données journalières
merged_lag <- merge(
  All_debit,
  Wind_T |> select(date, RR),
  by = "date",
  all = FALSE
)

# Calculer la corrélation pour chaque lag (0 à 7 jours)
resultats_lag <- tibble(
  lag       = 0:7,
  rho       = NA_real_,
  p_value   = NA_real_
)

for (i in 0:7) {
  # Décaler les précipitations de i jours en avance sur le débit
  merged_lag_i <- merged_lag |>
    mutate(RR_lag = lag(RR, n = i))  # RR d'il y a i jours
  
  test <- cor.test(merged_lag_i$debit_cumule, 
                   merged_lag_i$RR_lag, 
                   method = "spearman", 
                   use = "complete.obs")
  
  resultats_lag$rho[i + 1]     <- test$estimate
  resultats_lag$p_value[i + 1] <- test$p.value
}

# Afficher les résultats
print(resultats_lag)

# Identifier le meilleur lag
meilleur_lag <- resultats_lag |> 
  filter(p_value < 0.05) |>   # seulement les significatifs
  slice_max(abs(rho), n = 1)

cat("Meilleur lag :", meilleur_lag$lag, "jours\n")
cat("Rho =", round(meilleur_lag$rho, 3), "\n")
cat("p =", round(meilleur_lag$p_value, 4), "\n")

# Visualisation
ggplot(resultats_lag, aes(x = lag, y = rho)) +
  geom_col(aes(fill = p_value < 0.05), width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "ns"),
                    name = "Significativité") +
  geom_text(aes(label = round(rho, 2)), vjust = -0.5, size = 4) +
  labs(title = "Corrélation de Spearman précipitations → débit selon le lag",
       x = "Lag (jours)", y = "Rho de Spearman") +
  theme_bw(base_size = 13)


# on veut mettre en lien avec le débit liquide du Var

load("data/Hydro France/Y6442010_depuis_2000.Rdata")
load("data/Hydro France/All_debit.Rdata")

Y6442010_2006_2024 <- Y6442010_depuis_2000 |>
  filter(date >= "2006-01-01", date <= "2024-12-31")

All_debit <- All_debit |> 
  filter(date >= "2014-01-01", date <= "2024-12-31")

Wind_T <- Wind_T |> 
  filter(date >= as.Date("2014-01-01"), date <= as.Date("2024-12-31"))

# mise à l'échelle
adjust_factors <- sec_axis_adjustement_factors(Wind_T$RR, All_debit$debit_cumule)
Wind_T$scaled_RR <- Wind_T$RR * adjust_factors$diff + adjust_factors$adjust

# Calcul de la corrélation entre débit et aire des panaches
merged_data <- merge(
  All_debit,
  Wind_T,
  by = "date",
  all = FALSE
) |> 
  mutate(RR_lag1 = lag(RR, n = 1))

correlation <- cor(merged_data$debit_cumule, merged_data$RR_lag1, 
                   method = "spearman", use = "complete.obs")
p_value <- cor.test(merged_data$debit_cumule, merged_data$RR_lag1, 
                    method = "spearman")$p.value
ggplot() +
  geom_line(
    data = All_debit,
    aes(x = date, y = debit_cumule, color = "Débit"),
    linewidth = 0.4
  ) +
  geom_line(
    data = Wind_T,
    aes(x = date, y = scaled_RR, color = "Précipitations"),
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("Précipitations" = "aquamarine", "Débit" = "darkolivegreen3"),
    name = "Légende"
  ) +
  scale_y_continuous(
    name = expression("Débit du Var (m"^{3}*".s"^{-1}*")"),
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Précipitations (mm)")
    )
  ) +
  labs(
    title = "Évolution des précipitations et du débit cumulé des fleuves niçois",
    caption = "Sources : Archives Météo France - Hydro Portail - MNCA",
    x = "Date",
    color = "Variable"
  ) +
  annotate(
    "text",
    x = min(c(All_debit$date, Wind_T$date), na.rm = TRUE),
    y = max(c(All_debit$debit_cumule, Wind_T$RR), na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    label = gsub(" \\(ns\\)", "", mk_label),  # supprime le (ns)
    size = 8,
    color = "aquamarine3",
    family = "serif",
    fontface = "italic"
  ) +
  annotate(
    "text",
    x = max(c(All_debit$date, Wind_T$date), na.rm = TRUE),
    y = max(c(All_debit$debit_cumule, Wind_T$RR), na.rm = TRUE),
    hjust = 1,
    vjust = 1,
    label = paste0(
      "r = ", round(correlation, 2),
      "\n(lag 1 jour)",
      "\np ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.caption  = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title    = element_text(face = "bold", family = "serif"),
    axis.text     = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title  = element_text(face = "bold"),
    plot.margin   = margin(1, 1.5, 1, 1, "cm")
  ) +
  scale_x_date(
    date_breaks = "2 year",
    date_labels = "%Y"
  )


# mettre en relation précipitation et débit -------------------------------

## estimate liquid flow rate trend -----------------------------------------

### Var ---------------------------------------------------------------------

load("data/Hydro France/Y6442010_depuis_2000.Rdata")

Y6442010_depuis_2006 <- Y6442010_depuis_2000 |> 
  filter(date >= "2006-01-01", date <= "2025-12-31")  # ← ajout

debit_clean <- Y6442010_depuis_2006 |>
  drop_na(débit)

mk_debit <- mk.test(debit_clean$débit)
print(mk_debit)

debit_clean <- debit_clean |>
  mutate(date_num = as.numeric(date - min(date)))

sen_debit       <- sens.slope(debit_clean$débit)
slope_debit_jan <- sen_debit$estimates * 365

intercept_debit <- median(
  debit_clean$débit - sen_debit$estimates * debit_clean$date_num,
  na.rm = TRUE
)

debit_clean <- debit_clean |>
  mutate(theilsen_fit_debit = intercept_debit + sen_debit$estimates * date_num)

cat("Tendance débit :", round(slope_debit_jan, 3), "m³/s/an\n")
cat("Mann-Kendall p =", round(mk_debit$p.value, 4), "\n")

debit_clean_var <- debit_clean
mk_var  <- mk_debit
sen_var <- sen_debit

### Paillon ---------------------------------------------------------------------

load("~/River_runoff_analysis/data/MNCA/Paillon_all_debit.Rdata")

Paillon_all_debit <- Paillon_all_debit |> 
  filter(date <= "2025-12-31")   # ← ajout

debit_clean <- Paillon_all_debit |>
  drop_na(ABA_debit_mean)

mk_debit <- mk.test(debit_clean$ABA_debit_mean)
print(mk_debit)

debit_clean <- debit_clean |>
  mutate(date_num = as.numeric(date - min(date)))

sen_debit       <- sens.slope(debit_clean$ABA_debit_mean)
slope_debit_jan <- sen_debit$estimates * 365

intercept_debit <- median(
  debit_clean$ABA_debit_mean - sen_debit$estimates * debit_clean$date_num,
  na.rm = TRUE
)

debit_clean <- debit_clean |>
  mutate(theilsen_fit_debit = intercept_debit + sen_debit$estimates * date_num)

cat("Tendance débit :", round(slope_debit_jan, 3), "m³/s/an\n")
cat("Mann-Kendall p =", round(mk_debit$p.value, 4), "\n")

debit_clean_paillon <- debit_clean
mk_paillon  <- mk_debit
sen_paillon <- sen_debit

### Magnan ---------------------------------------------------------------------

load("~/River_runoff_analysis/data/MNCA/Magnan_all_debit.Rdata")

Magnan_all_debit <- Magnan_all_debit |> 
  filter(date <= "2025-12-31")   # ← ajout

debit_clean <- Magnan_all_debit |>
  drop_na(AAM_debit_mean)

mk_debit <- mk.test(debit_clean$AAM_debit_mean)
print(mk_debit)

debit_clean <- debit_clean |>
  mutate(date_num = as.numeric(date - min(date)))

sen_debit       <- sens.slope(debit_clean$AAM_debit_mean)
slope_debit_jan <- sen_debit$estimates * 365

intercept_debit <- median(
  debit_clean$AAM_debit_mean - sen_debit$estimates * debit_clean$date_num,
  na.rm = TRUE
)

debit_clean <- debit_clean |>
  mutate(theilsen_fit_debit = intercept_debit + sen_debit$estimates * date_num)

cat("Tendance débit :", round(slope_debit_jan, 3), "m³/s/an\n")
cat("Mann-Kendall p =", round(mk_debit$p.value, 4), "\n")

debit_clean_magnan <- debit_clean
mk_magnan  <- mk_debit
sen_magnan <- sen_debit

### afficher les plots -------------------------------------------------------

plot_tendance <- function(data, date_col, debit_col,
                          mk_pval, slope_an, titre, ylim_max = NULL) {
  
  sig_label   <- ifelse(mk_pval < 2.2e-16,
                        "p < 2.2×10⁻¹⁶ *",
                        ifelse(mk_pval < 0.05,
                               paste0("p = ", round(mk_pval, 4), " *"),
                               paste0("p = ", round(mk_pval, 4), " (ns)")))
  slope_label <- paste0("Pente = ", round(slope_an, 3), " m³/s/an")
  
  y_max_visible <- if (!is.null(ylim_max)) ylim_max else max(data[[debit_col]], na.rm = TRUE)
  
  p <- ggplot(data, aes(x = .data[[date_col]])) +
    geom_line(aes(y = .data[[debit_col]]),
              color = "steelblue", alpha = 0.6, linewidth = 0.4) +
    annotate("text",
             x     = min(data[[date_col]]),
             y     = y_max_visible * 0.95,
             label = paste(sig_label, slope_label, sep = "\n"),
             hjust = 0, vjust = 1, size = 8,
             color = ifelse(mk_pval < 0.05, "firebrick", "gray40")) +
    labs(title = titre, x = NULL, y = "Débit (m³/s)") +
    theme_bw(base_size = 11)
  
  if (!is.null(ylim_max)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_max))
  }
  
  return(p)
}

# Graphiques
p1 <- plot_tendance(debit_clean_var,     "date", "débit",
                    mk_var$p.value,     sen_var$estimates     * 365, "Var")

p2 <- plot_tendance(debit_clean_paillon, "date", "ABA_debit_mean",
                    mk_paillon$p.value, sen_paillon$estimates * 365, "Paillon")

p3 <- plot_tendance(debit_clean_magnan,  "date", "AAM_debit_mean",
                    mk_magnan$p.value,  sen_magnan$estimates  * 365, "Magnan")

p1 / p2 / p3


# ── Panneau d : précipitations + corrélation avec débit du Var ──────────────

# Filtrer les précipitations sur la même période que le Var
Wind_T_filtered <- Wind_T |>
  filter(date >= "2006-01-01", date <= "2025-12-31") |>
  mutate(annee = year(date), mois = month(date))

# Totaux mensuels
pluie_mensuelle <- Wind_T_filtered |>
  group_by(annee, mois) |>
  summarise(total = sum(RR, na.rm = TRUE), .groups = "drop") |>
  arrange(annee, mois)

# Série temporelle mensuelle
pluie_ts <- ts(pluie_mensuelle$total,
               start     = c(min(pluie_mensuelle$annee), 1),
               frequency = 12)

# Test saisonnier de Mann-Kendall
smk_result <- smk.test(pluie_ts)
smk_pval   <- smk_result$p.value

# Label tendance
smk_label <- ifelse(smk_pval < 2.2e-16,
                    "p < 2.2×10⁻¹⁶ *",
                    ifelse(smk_pval < 0.05,
                           paste0("p = ", round(smk_pval, 4), " *"),
                           paste0("p = ", round(smk_pval, 4), " (ns)")))

# Corrélation Spearman débit Var ~ précipitations lag 1 jour
merged_pluie_debit <- debit_clean_var |>
  select(date, débit) |>
  inner_join(
    Wind_T_filtered |> select(date, RR),
    by = "date"
  ) |>
  mutate(RR_lag1 = lag(RR, n = 1)) |>
  drop_na(débit, RR_lag1)

cor_result <- cor.test(merged_pluie_debit$débit,
                       merged_pluie_debit$RR_lag1,
                       method = "spearman")

cor_label <- paste0(
  "r = ", round(cor_result$estimate, 2),
  "\np ", ifelse(cor_result$p.value < 0.001, "< 0.001",
                 paste0("= ", round(cor_result$p.value, 3)))
)

# Mise à l'échelle pour double axe
adjust_factors <- sec_axis_adjustement_factors(Wind_T_filtered$RR,
                                               debit_clean_var$débit)
Wind_T_filtered <- Wind_T_filtered |>
  mutate(RR_scaled = RR * adjust_factors$diff + adjust_factors$adjust)

# Graphique
p4 <- ggplot() +
  geom_line(data = debit_clean_var,
            aes(x = date, y = débit, color = "Débit Var"),
            linewidth = 0.4, alpha = 0.6) +
  geom_line(data = Wind_T_filtered,
            aes(x = date, y = RR_scaled, color = "Précipitations"),
            linewidth = 0.4, alpha = 0.6) +
  scale_color_manual(
    values = c("Débit Var" = "steelblue", "Précipitations" = "darkorange"),
    name   = NULL
  ) +
  scale_y_continuous(
    name     = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = "Précipitations (mm)"
    )
  ) +
  annotate("text",
           x = min(debit_clean_var$date),
           y = max(debit_clean_var$débit, na.rm = TRUE) * 0.95,
           hjust = 0, vjust = 1, size = 8, color = "darkorange",
           label = paste0("Tendance pluie : ", smk_label)) +
  annotate("text",
           x = min(debit_clean_var$date),
           y = max(debit_clean_var$débit, na.rm = TRUE) * 0.75,
           hjust = 0, vjust = 1, size = 8, color = "grey20",
           label = cor_label) +
  labs(title = "Précipitations et débit du Var (2006–2025)",
       x = NULL, y = "Débit (m³/s)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

# ── Assemblage final ────────────────────────────────────────────────────────

(p1 / p2 / p3 / p4) +
  plot_annotation(
    title   = "Évolution des débits fluviaux et des précipitations — 2006–2025",
    caption = "Sources : Hydro France, MNCA, Météo France",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )



# ── Assemblage final ────────────────────────────────────────────────────────

# Label tendance précipitations sans (ns)
smk_label <- ifelse(smk_pval < 2.2e-16,
                    "p < 2.2×10⁻¹⁶ *",
                    ifelse(smk_pval < 0.05,
                           paste0("p = ", round(smk_pval, 4), " *"),
                           paste0("p = ", round(smk_pval, 4))))  # ← (ns) supprimé

# Refaire p4 avec le label corrigé
p4 <- ggplot() +
  geom_line(data = debit_clean_var,
            aes(x = date, y = débit, color = "Débit Var"),
            linewidth = 0.4, alpha = 0.6) +
  geom_line(data = Wind_T_filtered,
            aes(x = date, y = RR_scaled, color = "Précipitations"),
            linewidth = 0.4, alpha = 0.6) +
  scale_color_manual(
    values = c("Débit Var" = "steelblue", "Précipitations" = "darkorange"),
    name   = NULL
  ) +
  scale_y_continuous(
    name     = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = "Précipitations (mm)"
    )
  ) +
  scale_x_date(limits = c(as.Date("2006-01-01"), as.Date("2025-12-31"))) +
  # Tendance pluie → à gauche
  annotate("text",
           x = as.Date("2006-01-01"),
           y = max(debit_clean_var$débit, na.rm = TRUE) * 0.95,
           hjust = 0, vjust = 1, size = 8, color = "darkorange",
           label = paste0("Tendance pluie : ", smk_label)) +
  # Corrélation → à droite
  annotate("text",
           x = as.Date("2025-12-31"),      # ← x à droite
           y = max(debit_clean_var$débit, na.rm = TRUE) * 0.95,
           hjust = 1, vjust = 1, size = 8, color = "grey20",  # ← hjust = 1
           label = cor_label) +
  labs(title = "d) Précipitations et débit du Var (2006–2025)",
       x = NULL, y = "Débit (m³/s)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

# Refaire p1, p2, p3 avec annotations a), b), c)
p1 <- plot_tendance(debit_clean_var,     "date", "débit",
                    mk_var$p.value,     sen_var$estimates     * 365, "a) Var",
                    date_min = "2006-01-01", date_max = "2025-12-31")

p2 <- plot_tendance(debit_clean_paillon, "date", "ABA_debit_mean",
                    mk_paillon$p.value, sen_paillon$estimates * 365, "b) Paillon",
                    date_min = "2013-01-01", date_max = "2025-12-31")

p3 <- plot_tendance(debit_clean_magnan,  "date", "AAM_debit_mean",
                    mk_magnan$p.value,  sen_magnan$estimates  * 365, "c) Magnan",
                    date_min = "2014-01-01", date_max = "2025-12-31")

# Et supprimer le (ns) dans la fonction plot_tendance
plot_tendance <- function(data, date_col, debit_col,
                          mk_pval, slope_an, titre,
                          date_min, date_max,
                          ylim_max = NULL) {
  
  sig_label   <- ifelse(mk_pval < 2.2e-16,
                        "p < 2.2×10⁻¹⁶ *",
                        ifelse(mk_pval < 0.05,
                               paste0("p = ", round(mk_pval, 4), " *"),
                               paste0("p = ", round(mk_pval, 4))))
  slope_label <- paste0("Pente = ", round(slope_an, 3), " m³/s/an")
  
  y_max_visible <- if (!is.null(ylim_max)) ylim_max else max(data[[debit_col]], na.rm = TRUE)
  
  p <- ggplot(data, aes(x = .data[[date_col]])) +
    geom_line(aes(y = .data[[debit_col]]),
              color = "steelblue", alpha = 0.6, linewidth = 0.4) +
    scale_x_date(limits = c(as.Date(date_min), as.Date(date_max))) +
    annotate("text",
             x     = as.Date(date_min),
             y     = y_max_visible * 0.95,
             label = paste(sig_label, slope_label, sep = "\n"),
             hjust = 0, vjust = 1, size = 8,
             color = ifelse(mk_pval < 0.05, "firebrick", "gray40")) +
    labs(title = titre, x = NULL, y = "Débit (m³/s)") +
    theme_bw(base_size = 11)
  
  if (!is.null(ylim_max)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_max))
  }
  
  return(p)
}

# Regénérer p1, p2, p3 avec la fonction corrigée
p1 <- plot_tendance(debit_clean_var,     "date", "débit",
                    mk_var$p.value,     sen_var$estimates     * 365, "a) Var",
                    date_min = "2006-01-01", date_max = "2025-12-31")

p2 <- plot_tendance(debit_clean_paillon, "date", "ABA_debit_mean",
                    mk_paillon$p.value, sen_paillon$estimates * 365, "b) Paillon",
                    date_min = "2013-01-01", date_max = "2025-12-31")

p3 <- plot_tendance(debit_clean_magnan,  "date", "AAM_debit_mean",
                    mk_magnan$p.value,  sen_magnan$estimates  * 365, "c) Magnan",
                    date_min = "2014-01-01", date_max = "2025-12-31")

# Assemblage
(p1 / p2 / p3 / p4) +
  plot_annotation(
    caption = "Sources : Hydro France, MNCA, Météo France",
    theme   = theme(
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )

cat("r de Spearman =", round(cor_result$estimate, 3), "\n")
print(cor_result$p.value)
