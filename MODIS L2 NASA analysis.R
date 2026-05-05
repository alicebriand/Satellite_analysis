# MODIS L2 NASA analysis

# 16/03/2026
# pathway : "~/Satellite_analysis/MODIS L2 NASA analysis.R

# This script will load MODIS L2 data from NASA and will try, thanks to the reflectance parameter,
# to retrieve SPM concentration of the study area with different formula. The 
# ultimate goal is to identify which algorithm is the best to retrieve SPM concentration. It
# will them perform TS analysis using the adapted algorithm

# Setup ------------------------------------------------------------------

# NB: When running a script with source(), it will performa ll of the actions therein
# This is not always ideal, especially if the script is downloading data or performing large analyses.
# In this case it is better to save the output of the first script, then load them in the second script
# source("~/Satellite_analysis/earth_data_access.R")

# The shared functions between scripts are however and excellent reason to use source()
# Therefore I have moved the shared functions to a script that is source()'d here
source("func.R")

load("~/River_runoff_analysis/data/Hydro France/Var_crues.Rdata")
Log_Var_Paillon_2016_2017 <- read_delim("~/Downloads/Var_Paillon_Mars2017/Var_Paillon_Mars2017/Log_Var_Paillon_2016_2017_propre.csv",
                                        delim = ",", locale = locale(decimal_mark = ","))
load("data/Hydro France/All_debit_2024.Rdata")
load("data/MODIS L2 NASA/study_area_df_2024")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(terra)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2)
library(heatwaveR)
library(ggpmisc)
library(sf)
library(ggpubr)
library(maptiles)
library(tidyterra)
library(patchwork)
library(giscoR) # Hi-res coastlines
library(ggspatial)  # pour la flèche nord et l'échelle

# function ---------------------------------------------------------------

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

# This function will convert reflectance into SPM concentration with the Var parameters :
Var_Morin <- function(study_area_df, sur_refl_b01_1 = "sur_refl_b01_1", A = 80, C = 0.1562, SPM_max = 500) {
  
  # Vérifier que la colonne de réflectance existe
  if (!(sur_refl_b01_1 %in% names(study_area_df))) {
    stop(paste("La colonne", sur_refl_b01_1, "n'existe pas dans le data frame."))
  }
  
  # Calculer SPM avec une approche conditionnelle
  study_area_df <- study_area_df %>%
    mutate(
      SPM = case_when(
        .data[[sur_refl_b01_1]] >= C ~ SPM_max,  # Si ρw >= C, SPM = SPM_max
        TRUE ~ (A * .data[[sur_refl_b01_1]]) / (1 - .data[[sur_refl_b01_1]] / C)  # Sinon, appliquer la formule
      )
    )
  
  return(study_area_df)
}

# This function will convert reflectance into SPM concentration with the Paillon parameters :
Paillon_Morin <- function(study_area_df, sur_refl_b01_1, A = 39, C = 0.2563) {
  # Vérifier que la colonne de réflectance existe
  if (!(sur_refl_b01_1 %in% names(study_area_df))) {
    stop(paste("La colonne", sur_refl_b01_1, "n'existe pas dans le data frame."))
  }
  
  # Calculer SPM avec la formule : SPM = A * ρw / (1 - ρw / C)
  study_area_df <- study_area_df %>%
    mutate(SPM = pmin((A * .data[[sur_refl_b01_1]]) / (1 - .data[[sur_refl_b01_1]] / C), 500)
    )
  
  return(study_area_df)
}

# NB: When writing a function it is a good idea (but not necessary) to name your arguments
# in a way that is not found in your code outside of the function
# E.g. Here I changed 'study_area_df' to 'df' and 'sur_refl_b01_1' to 'col_name'
MODIS_L2_SPM <- function(df, col_name, A = 80, C = 0.1562){
  
  # Vérifier que la colonne de réflectance existe
  if (!(col_name %in% names(df))) {
    stop(paste("La colonne", col_name, "n'existe pas dans le data frame."))
  }
  
  # Calculer SPM avec la formule : SPM = A * ρw / (1 - ρw / C)

  study_area_df <- study_area_df %>%
    mutate(SPM = pmin((A * .data[[sur_refl_b01_1]]) / (1 - .data[[sur_refl_b01_1]] / C), 500)
    )
  
  return(study_area_df)
}

# This is a formula seen in the paper of Doxaran et al., ... : 
Gironde_Doxaran <- function(study_area_df, sur_refl_b01_1, sur_refl_b02_1) {
  
  # Vérifier que les colonnes de réflectance existent
  if (!(sur_refl_b01_1 %in% names(study_area_df))) {
    stop(paste("La colonne", sur_refl_b01_1, "n'existe pas dans le data frame."))
  }
  if (!(sur_refl_b02_1 %in% names(study_area_df))) {
    stop(paste("La colonne", sur_refl_b02_1, "n'existe pas dans le data frame."))
  } 
  
  # Calculer SPM avec la formule : SPM = 12.996 * exp((R(B2)/R(B1)) / 0.189)
  study_area_df <- study_area_df %>%
    mutate(
      SPM = case_when(
        .data[[sur_refl_b01_1]] <= 0 ~ NA_real_,  # Évite la division par zéro ou des valeurs négatives
        TRUE ~ 12.996 * exp((.data[[sur_refl_b02_1]] / .data[[sur_refl_b01_1]]) / 0.189)
      )
    )

  df <- df |> 
    mutate(SPM = (A * .data[[col_name]]) / (1 - (.data[[col_name]] / C)))
  
  return(df)
}

# data analysis -----------------------------------------------------------

# we first start to exclude negative reflectance
study_area_df_2024 <- study_area_df_2024 |>  
  filter(sur_refl_b01_1 > 0, sur_refl_b02_1 > 0)

# It is a single equation, we can apply it directly to the data.frame with mutate()
# SPM = A * ρw / (1 - ρw / C); A = 80, C = 0.1562 # But where does this equation 
# and values come from? I do not find them in the literature?
study_area_df_2024 <- study_area_df_2024 |>  
  mutate(SPM_Morin_Var = (80 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.1562)),
         SPM_Morin_Paillon = (39 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.2563)))

# Équation Doxaran et al., 2009 (il faut les deux bandes réflectance)
study_area_df_2024 <- study_area_df_2024 |>
  mutate(SPM_Doxaran = 12.996 * exp((sur_refl_b02_1/sur_refl_b01_1)/0.189),
         SPM_Doxaran = ifelse(is.infinite(SPM_Doxaran), NA, SPM_Doxaran))

# A different algorithm based on Teng et al. 2025
# https://www.sciencedirect.com/science/article/pii/S003442572500149X
# SPM_org = a Rrs(lambda_RED)^b; a = 1992.2, b = 1.027
# NB: lambda_RED is taken here to be the MODIS band 1 waveband
study_area_df_2024 <- study_area_df_2024 |> 
mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
       Rrs_b02_01 = (sur_refl_b02_1/pi), # First convert Rhow_w to Rrs
       SPM_Teng_MO = 1992.2 * Rrs_b01_01^1.027,
       SPM_Teng_MM = 12662.7 * Rrs_b02_01^1.157,
       SPM_Teng_extrm_MM = 50556.7 * Rrs_b02_01^1.371)
# But these values are crazy high...

# So then this paper by Tsapanou et al. 2020
# http://www.teiath.gr/userfiles/pdrak/lab/coupling_remote_sensing_data.pdf
# Though this is for LandSat 8
# SPM = ((A * Rho_W)/(1-(Rhow_w/C)))+B
# A = 366,53 g m−3 , B = 0 g m−3 and C = 0.0324
study_area_df_2024 <- study_area_df_2024 |>
  mutate(SPM_Tsapanou = ((366.53 * sur_refl_b01_1)/(1-(sur_refl_b01_1/0.0324))))
# this produce a lot of negative values

# So we digress to the Nechad formula of 
# SPM = ((A * Rhow)/(1-(Rhow/C)))+B
# A = 289.29, C = 0.1686, B = 2.10 
study_area_df_2024 <- study_area_df_2024 |>
  mutate(SPM_Nechad = ((289.29 * sur_refl_b01_1)/(1-(sur_refl_b01_1/0.1686))) + 2.10)

# plotting -----------------------------------------------------------

# we have to choose some dates and look at what data look like : 
study_area_df_04_03_2024 <- study_area_df_2024 |> 
  filter(date == "2024-03-04")

# image avec des sortes de bandes
study_area_df_01_06_2024 <- study_area_df_2024 |> 
  filter(date == "2024-06-01")

# image claire
study_area_df_17_06_2024 <- study_area_df_2024 |> 
  filter(date == "2024-06-17")

study_area_df_25_11_2016 <- study_area_df_clean_2016 |> 
  filter(date == "2016-11-25")

# complexe pour 2017 car pas de crue
study_area_df_26_03_2017 <- study_area_df_clean_2017 |> 
  filter(date == "2017-03-26")

study_area_df_28_03_2017 <- study_area_df_clean_2017 |> 
  filter(date == "2017-03-28")

## Morin Var -------------------------------------------------------------------

# First, when we do the map there are some extreme values (until 8000), which is weird
# we want to erase them so maybe the map will be better
# study_area_df_04_03_2024_Morin <- study_area_df_04_03_2024 |> 
#   filter(SPM_Morin <= 100)

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# set the maximum values
max_spm <- max(study_area_df_17_06_2024$sur_refl_b01_1, na.rm = TRUE)

pl_map <- study_area_df_17_06_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = sur_refl_b01_1)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),  # ← écriture scientifique
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 17 juin 2024",
    subtitle = "Réflectance de MODIS à la bande 1 (645 nm)",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_17_06_2024$lon),
    ylim   = range(study_area_df_17_06_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/réflectance/sans filtrage des images contaminées/fig_MODIS_SPM_17_06_2024_reflectance_b1.png", pl_map, height = 9, width = 14)

## Morin Paillon -------------------------------------------------------------------

# First, when we do the map there are some extreme values (until 8000), which is weird
# we want to erase them so maybe the map will be better
# study_area_df_04_03_2024_Morin <- study_area_df_04_03_2024 |> 
#   filter(SPM_Morin <= 100)

# set the maximum values
max_spm <- max(study_area_df_04_03_2024$SPM_Morin_Paillon, na.rm = TRUE)

# Créer le graphique

pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Morin_Paillon)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),  # ← écriture scientifique
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Morin et al. (Paillon) appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Morin_Paillon.png", pl_map, height = 9, width = 14)

## Teng MO -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_MO, na.rm = TRUE)

pl_map1 <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_MO)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),  # ← écriture scientifique
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Teng et al. pour les eaux riches en matière organique appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Teng_MO.png", pl_map, height = 9, width = 14)

## Teng MM -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_MM, na.rm = TRUE)

# Créer le graphique
pl_map2 <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_MM)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Teng et al. pour les eaux riches en matière minérale appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Teng_MM.png", pl_map, height = 9, width = 14)

## Teng extrem MM -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_extrm_MM, na.rm = TRUE)

pl_map3 <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_extrm_MM)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Teng et al. pour les eaux très riches en matière minérale appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1 et 2/fig_MODIS_SPM_04_03_2024_Teng_extrem_MM.png", pl_map, height = 9, width = 14)

## patchwork Teng ----------------------------------------------------------

# Côte à côte
pl_map1 / pl_map2 | pl_map3 +
  plot_annotation(
    tag_levels = "A",
    title    = "Cartes de concentration en MES",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50")
    )
  )

## Tsapanou -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Tsapanou, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Tsapanou)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Tsapanou et al. appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )
# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Tsapanou_new.png", pl_map, height = 9, width = 14)

## Nechad -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Nechad, na.rm = TRUE)

pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Nechad)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Nechad appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Nechad.png", pl_map, height = 9, width = 14)

## Doxaran -------------------------------------------------------------------

# First, when we do the map there are some extreme values (until 8000), which is weird
# we want to erase them so maybe the map will be better
study_area_df_04_03_2024_1 <- study_area_df_04_03_2024 |>
  filter(SPM_Doxaran <= 2000)

# Utiliser le 99e percentile pour éviter que les valeurs extrêmes écrasent l'échelle
max_spm <- quantile(study_area_df_04_03_2024$SPM_Doxaran, 0.99, na.rm = TRUE)

# Optionnel : voir la distribution pour choisir le bon seuil
summary(study_area_df_04_03_2024$SPM_Doxaran)
hist(study_area_df_04_03_2024$SPM_Doxaran, breaks = 50)

pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Doxaran)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),
    limits = c(0, max_spm),
    oob    = scales::squish  # ← les valeurs > max_spm sont ramenées au max
    #   au lieu d'être mises en NA
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 04 mars 2024",
    subtitle = "Algorithme de Doxaran et al. appliqué aux données MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(study_area_df_04_03_2024$lon),
    ylim   = range(study_area_df_04_03_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/bande 1 et 2/fig_MODIS_SPM_04_03_2024_Doxaran.png", 
       pl_map, height = 9, width = 14)

# SPM prediction vs in situ data -------------------------------------------------------

## in situ data ------------------------------------------------------------

# we have to extract lat and lon from in situ data

# renommer le doc
Var_2016 <- Log_Var_Paillon_2016_2017

# garder seulement les paramètres d'intérêts
Var_2016 <- Var_2016 |> 
  select(StationID, Date, Latitude_dd, Longitude_dd, Location, `SPMmoyen(mg/l)`, 
         `SPMsd(mg/l)`)

# rename columns
Var_2016 <- Var_2016 |> 
  rename(date = Date,
         lat = Latitude_dd, 
         lon = Longitude_dd,
         SPM_moy = `SPMmoyen(mg/l)`, 
         SPM_sd = `SPMsd(mg/l)`)

# transformer la date en "vraie" date
Var_2016 <- Var_2016 %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"))

Var_2016 <- Var_2016 |>
  mutate(
    lat = as.numeric(gsub(",", ".", lat)),
    lon = as.numeric(gsub(",", ".", lon))
  )

Var_2016 <- Var_2016[-6,]

## MODIS data ------------------------------------------------------------

# we have to choose the speficical date
# study_area_df_25_11_2016 <- study_area_df_clean_2016 |>
#   filter(date == as.Date("2016-11-25"))

## merging data set --------------------------------------------------------

# we have to create a df with in situ data and SPM prediction
study_area_df_clean_2016_in_situ <- inner_join(study_area_df_25_11_2016, Var_2016, by = c("date", "lon", "lat"))

# If it worked we should have 7 rows, I hope it match

# it doesn't work

## spatial correspondence --------------------------------------------------------

# What if it doesn't match ? Then we would have to take the closer location to 
# our points, we should identify the 8 points around the point of interest (I did not
# do that finally)

# st_as_sf : convert foreign object to an sf object
# ce sont des df dont l'une des colonnes contient une géométrie

# Convertir les points in situ en objet spatial
pts_insitu <- st_as_sf(Var_2016, coords = c("lon", "lat"), crs = 4326) # crs = coordinate reference system

# Convertir les pixels MODIS en objet spatial (pour une date donnée)
pts_modis <- study_area_df_25_11_2016 |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Trouver le pixel MODIS le plus proche de chaque point in situ
idx_proche <- st_nearest_feature(pts_insitu, pts_modis)

# Extraire les coordonnées des pixels MODIS les plus proches
coords_modis <- st_coordinates(pts_modis[idx_proche, ])

# Extraire les valeurs MODIS correspondantes
Var_2016_modis <- Var_2016 |> 
  mutate(sur_refl_b01_modis = pts_modis$sur_refl_b01_1[idx_proche],
         sur_refl_b02_modis = pts_modis$sur_refl_b02_1[idx_proche])

Var_2016_modis <- Var_2016 |> 
  mutate(
    sur_refl_b01_modis = pts_modis$sur_refl_b01_1[idx_proche],
    sur_refl_b02_modis = pts_modis$sur_refl_b02_1[idx_proche],
    lon_modis = coords_modis[, "X"],   # <-- coordonnées du pixel MODIS associé
    lat_modis = coords_modis[, "Y"]    # <-- coordonnées du pixel MODIS associé
  )

# We have to convert reflectance in SPM concentration

## apply equations ---------------------------------------------------------

Var_2016_modis <- Var_2016_modis |>
  filter(sur_refl_b01_modis >= 0 ) |> 
  mutate(SPM_Morin_Var = (80 * sur_refl_b01_modis) / (1 - (sur_refl_b01_modis / 0.1562)),
         SPM_Morin_Paillon = (39 * sur_refl_b01_modis) / (1 - (sur_refl_b01_modis / 0.2563))) |> 
  mutate(SPM_Doxaran = 12.996 * exp((sur_refl_b02_modis/sur_refl_b01_modis)/0.189),
         SPM_Doxaran = ifelse(is.infinite(SPM_Doxaran), NA, SPM_Doxaran)) |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_modis/pi),
         Rrs_b02_01 = (sur_refl_b02_modis/pi), # First convert Rhow_w to Rrs
         SPM_Teng_MO = 1992.2 * Rrs_b01_01^1.027,
         SPM_Teng_MM = 12662.7 * Rrs_b02_01^1.157,
         SPM_Teng_extrm_MM = 50556.7 * Rrs_b02_01^1.371) |> 
  mutate(SPM_Tsapanou = ((289.29 * sur_refl_b01_modis)/(1-(sur_refl_b01_modis/0.1686)))+2.10) |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_modis/pi), # First convert Rhow_w to Rrs
         SPM_Nechad = ((200 * Rrs_b01_01)/(1-(Rrs_b01_01/0.17)))+0) |> 
  mutate(SPM_formule = 0 + 2.1 * sur_refl_b01_modis)

# We have to plot a scatter plot for each parameter

### ── Calcul des métriques ───────────────────────────────────────────────────────
metriques <- Var_2016_modis |>
  summarise(
    # Morin Var
    REQM_Morin_Var  = sqrt(mean((SPM_Morin_Var - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Morin_Var  = mean(abs((SPM_moy - SPM_Morin_Var) / SPM_moy) * 100, na.rm = TRUE),
    # Morin Paillon
    REQM_Morin_Paillon  = sqrt(mean((SPM_Morin_Paillon - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Morin_Paillon  = mean(abs((SPM_moy - SPM_Morin_Paillon) / SPM_moy) * 100, na.rm = TRUE),
    # Teng MO
    REQM_Teng_MO  = sqrt(mean((SPM_Teng_MO - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Teng_MO  = mean(abs((SPM_moy - SPM_Teng_MO) / SPM_moy) * 100, na.rm = TRUE),
    # Teng MM
    REQM_Teng_MM  = sqrt(mean((SPM_Teng_MM - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Teng_MM  = mean(abs((SPM_moy - SPM_Teng_MM) / SPM_moy) * 100, na.rm = TRUE),
    # Teng extrêmment riche en MM
    REQM_Teng_extrem_MM  = sqrt(mean((SPM_Teng_extrm_MM - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Teng_extrem_MM  = mean(abs((SPM_moy - SPM_Teng_extrm_MM) / SPM_moy) * 100, na.rm = TRUE),
    # Doxaran
    REQM_Doxaran  = sqrt(mean((SPM_Doxaran - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Doxaran  = mean(abs((SPM_moy - SPM_Doxaran) / SPM_moy) * 100, na.rm = TRUE),
    # Tsapanou
    REQM_Tsapanou  = sqrt(mean((SPM_Tsapanou - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Tsapanou  = mean(abs((SPM_moy - SPM_Tsapanou) / SPM_moy) * 100, na.rm = TRUE),
    # Nechad
    REQM_Nechad  = sqrt(mean((SPM_Nechad - SPM_moy)^2, na.rm = TRUE)),
    MAPE_Nechad  = mean(abs((SPM_moy - SPM_Nechad) / SPM_moy) * 100, na.rm = TRUE),
    # Formule Robert
    REQM_formule_Robert  = sqrt(mean((SPM_formule - SPM_moy)^2, na.rm = TRUE)),
    MAPE_formule_Robert  = mean(abs((SPM_moy - SPM_formule) / SPM_moy) * 100, na.rm = TRUE)
  )

## Morin Var -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Morin_Var,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Morin_Var,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Morin_Var, y = SPM_moy)) +
  
  theme_bw(base_size = 14) +
  
  geom_abline(intercept = 0, slope = 1, linewidth = 1,
              linetype = "dashed", color = "grey40") +
  
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue",
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  geom_point(size = 3, shape = 21,
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R²
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 5, color = "grey20") +
  
  # p-valeur
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3, size = 5, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  labs(
    title    = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Morin et al. pour le Var / mesures in situ du 25/11/2016",
    x        = expression(SPM[Morin] ~ (mg ~ L^{-1})),
    y        = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  theme(
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle    = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70")
  )

## Morin Paillon -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Morin_Paillon,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Morin_Paillon,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Morin_Paillon, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Morin et al. pour le Paillon / mesures in situ du 25/11/2016",
    x = expression(SPM[Morin] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Doxaran -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Doxaran,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Doxaran,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Doxaran, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison modèle développé par Doxaran et al., 2009 / mesures in situ du 25/11/2016",
    x = expression(SPM[Doxaran] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Teng MO -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Teng_MO,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Teng_MO,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Teng_MO, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Teng et al., 2025 (classification en eaux riches en MO) / mesures in situ du 25/11/2016",
    x = expression(SPM[Teng] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Teng MM -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Teng_MM,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Teng_MM,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Teng_MM, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Teng et al., 2025 (classification en eaux riches en MM) / mesures in situ du 25/11/2016",
    x = expression(SPM[Teng] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Teng extreme MM -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Teng_extrem_MM,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Teng_extrem_MM,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Teng_extrm_MM, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Teng et al., 2025 (classification en eaux extrêmement riches en MM) / mesures in situ du 25/11/2016",
    x = expression(SPM[Teng] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Tsapanou -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Tsapanou,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Tsapanou,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Tsapanou, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Tsapanou et al., 2020 / mesures in situ du 25/11/2016",
    x = expression(SPM[Tsapanou] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Nechad -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_Nechad,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_Nechad,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_Nechad, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de l'algorithme de Nechad et al., 2010 / mesures in situ du 25/11/2016",
    x = expression(SPM[Teng] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

## Formule -------------------------------------------------------------------

# Texte à afficher sur le graph
label_metriques <- paste0(
  "REQM = ",  round(metriques$REQM_formule_Robert,  2), " mg/L\n",
  # "Biais = ", round(metriques$Biais, 2), " mg/L\n",
  "MAPE = ",  round(metriques$MAPE_formule_Robert,  1), " %"
)

### ── Graphique ──────────────────────────────────────────────────────────────────

ggplot(Var_2016_modis, aes(x = SPM_formule, y = SPM_moy)) +
  
  # Fond et grille
  theme_bw(base_size = 14) +
  
  # Ligne 1:1
  geom_abline(intercept = 0, slope = 1, linewidth = 1, 
              linetype = "dashed", color = "grey40") +
  
  # Régression
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", 
              alpha = 0.15, colour = "steelblue", linewidth = 1.5) +
  
  # Points
  geom_point(size = 3, shape = 21, 
             fill = "steelblue", color = "white", stroke = 0.5, alpha = 0.85) +
  
  # R² seul (ligne du haut)
  stat_cor(aes(label = after_stat(rr.label)),
           label.x.npc = "left", label.y.npc = "top",
           size = 10, color = "grey20") +
  
  # p-valeur seule (ligne du bas, décalée vers le bas)
  stat_cor(aes(label = after_stat(p.label)),
           label.x.npc = "left", label.y.npc = "top",
           vjust = 3,        # <-- décale vers le bas
           size = 10, color = "grey20") +
  
  # RMSE, Biais, MAPE dans un encadré en bas à droite
  annotate("label",
           x = Inf, y = -Inf,
           hjust = 1.05, vjust = -0.3,
           label = label_metriques,
           size = 4.5, color = "grey20",
           fill = "white", label.size = 0.3,   # encadré avec bordure
           fontface = "plain", lineheight = 1.5) +
  
  # Labels
  labs(
    title = "Prédiction de la concentration en MES vs concentration in situ en MES",
    subtitle = "Comparaison de la formule ... / mesures in situ du 25/11/2016",
    x = expression(SPM[formule] ~ (mg ~ L^{-1})),
    y = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )

# graphs ------------------------------------------------------------------

# les résultats sont très nuls donc je me demande si c'est possible que la localisation des points n'est pas
# fonctionné --> c'est le capteur qui n'était pas bon (mais c'est bon j'ai remplacé par Terra)

pts <- st_as_sf(Var_2016_modis, coords = c("lon", "lat"), crs = 4326)
tuiles <- get_tiles(pts, provider = "Esri.WorldImagery", zoom = 10)

# Carte de vérification
ggplot() +
  geom_spatraster_rgb(data = tuiles) +
  # Pixels MODIS (grands carrés)
  geom_point(data = Var_2016_modis,
             aes(x = lon_modis, y = lat_modis, color = "Pixels MODIS"),
             size = 8, shape = 15, alpha = 0.5) +
  # Points in situ (petits ronds)
  geom_point(data = Var_2016_modis,
             aes(x = lon, y = lat, color = "Points in situ"),
             size = 3, shape = 16) +
  # Relier chaque point in situ à son pixel MODIS
  geom_segment(data = Var_2016_modis,
               aes(x = lon, y = lat, xend = lon_modis, yend = lat_modis),
               linetype = "dashed", color = "white", linewidth = 0.5) +
  scale_color_manual(values = c("Points in situ" = "yellow",
                                "Pixels MODIS"   = "red")) +
  labs(title = "Appariement points in situ / pixels MODIS (250 m)",
       x = "Longitude", y = "Latitude",
       color = "Type de mesure") +
  theme_minimal()


# graphique relations mathématiques (claude l'a fait)

# ── Plages de réflectance ──────────────────────────────────────────────────────
rhow_seq <- seq(0.001, 0.15, length.out = 500)   # ρw (645 nm)
Rrs_seq  <- rhow_seq / pi                          # Rrs (645 nm)
rhow_NIR <- seq(0.001, 0.10, length.out = 500)    # ρw NIR (approximation)
Rrs_NIR  <- rhow_NIR / pi                          # Rrs NIR

# ── GRAPHIQUE 1 : équations en ρw (Rhow) ──────────────────────────────────────
df_rhow <- tibble(
  rhow = rhow_seq,
  Morin_Var      = (80  * rhow_seq) / (1 - (rhow_seq / 0.1562)),
  Morin_Paillon  = (39  * rhow_seq) / (1 - (rhow_seq / 0.2563)),
  Tsapanou       = ((289.29 * rhow_seq) / (1 - (rhow_seq / 0.1686))) + 2.10,
  Formule_Robert = 0 + 2.1 * rhow_seq,
) |>
  # Supprimer les valeurs aberrantes (asymptotes)
  mutate(across(-rhow, ~ ifelse(. < 0 | . > 1000, NA, .))) |>
  pivot_longer(-rhow, names_to = "Algorithme", values_to = "SPM")

p1 <- ggplot(df_rhow, aes(x = rhow, y = SPM, color = Algorithme)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 300)) +
  theme_bw(base_size = 13) +
  labs(
    title    = expression(bold("SPM = f(ρ"[w]*"(620-670 nm))")),
    subtitle = "Morin, Tsapanou, Formule Robert",
    x        = expression(rho[w] ~ "(645 nm)"),
    y        = expression(SPM ~ (mg ~ L^{-1})),
    color    = "Algorithme"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ── GRAPHIQUE 2 : équations en Rrs ────────────────────────────────────────────
df_rrs <- tibble(
  Rrs_645 = Rrs_seq,
  Rrs_NIR = Rrs_NIR,
  rhow_645 = rhow_seq,
  rhow_NIR = rhow_NIR,
  Teng_MO          = 1992.2  * Rrs_seq^1.027,
  Teng_MM          = 12662.7 * Rrs_NIR^1.157,
  Teng_extrm_MM    = 50556.7 * Rrs_NIR^1.371,
  Nechad         = ((200    * Rrs_seq) / (1 - (Rrs_seq / 0.17)))   + 0
  ) |>
  mutate(across(c(Teng_MO, Teng_MM, Teng_extrm_MM, Nechad),
                ~ ifelse(. < 0 | . > 1000, NA, .)))

# Teng en format long (x = Rrs 645 pour MO, Rrs NIR pour MM)
df_teng <- bind_rows(
  tibble(x = Rrs_seq,  SPM = 1992.2  * Rrs_seq^1.027,  Algorithme = "Teng_MO"),
  tibble(x = Rrs_NIR,  SPM = 12662.7 * Rrs_NIR^1.157,  Algorithme = "Teng_MM"),
  tibble(x = Rrs_NIR,  SPM = 50556.7 * Rrs_NIR^1.371,  Algorithme = "Teng_extrm_MM")
) |> mutate(SPM = ifelse(SPM < 0 | SPM > 1000, NA, SPM))

p2 <- ggplot(df_teng, aes(x = x, y = SPM, color = Algorithme)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, 500)) +
  theme_bw(base_size = 13) +
  labs(
    title    = expression(bold("SPM = f(R"[rs]*")")),
    subtitle = "MO (x = Rrs 645 nm) · MM & Extrême MM (x = Rrs 858.5)",
    x        = expression(R[rs] ~ (sr^{-1})),
    y        = expression(SPM ~ (mg ~ L^{-1})),
    color    = "Algorithme"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ── GRAPHIQUE 3 : équations Doxaran ────────────────────────────────────────────

# Doxaran en fonction du ratio Rhow NIR/645
df_dox <- tibble(
  ratio = seq(0.01, 3, length.out = 500),
  SPM   = 12.996 * exp(ratio / 0.189),
  Algorithme = "Doxaran"
) |> mutate(SPM = ifelse(SPM > 1000, NA, SPM))

p3 <- ggplot(df_dox, aes(x = ratio, y = SPM, color = Algorithme)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 300)) +
  theme_bw(base_size = 13) +
  labs(
    title    = expression(bold("SPM = f(ρ"[w]*"(620-670 nm) et (841-876 nm))")),
    subtitle = "Doxaran et al., 2009",
    x        = expression(rho[w] ~ "(645 nm) et (858.5 nm)"),
    y        = expression(SPM ~ (mg ~ L^{-1})),
    color    = "Algorithme"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ── Affichage ──────────────────────────────────────────────────────────────────
print(p1)
print(p2)
print(p3)

# Côte à côte
p1 | p2 | p3 +
  plot_annotation(
    tag_levels = "A",
    title    = "Comparaison des algorithmes SPM",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50")
    )
  )

# TS analysis -------------------------------------------------------------

## 2024 --------------------------------------------------------------------

# on a un df avec des valeurs par jour en bcp de points,
# on veut le réduire et moyenné ces valeurs pas jours ainsi que calculer l'extension
# du panache

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------
f <- "~/Downloads/MODIS NASA/L2 2024 Aqua/MYD09GQ.A2024001.h18v04.061.2024003111813.hdf"

r <- rast(f, subds = "sur_refl_b01_1")

# Vérifier la projection
crs(r, describe = TRUE)   # confirme : Sinusoidal, unité = mètres

# La résolution est en MÈTRES
res_m <- res(r)[1]
cat("Résolution :", res_m, "m\n")          # → ~231.6 m ≈ 250 m (résolution nominale MYD09GQ)

# Aire pixel directement en m² puis km²
aire_pixel_m2  <- res_m^2
aire_pixel_km2 <- aire_pixel_m2 / 1e6
cat("Aire d'un pixel :", round(aire_pixel_km2, 6), "km²\n")   # → ~0.0536 km²

## calcul de l'aire --------------------------------------------------------

# Reprojection en WGS84
r_wgs84 <- project(r, "EPSG:4326")

# Vérifier la nouvelle résolution (maintenant en degrés)
res_deg <- res(r_wgs84)[1]
cat("Résolution :", res_deg, "°\n")   # → ~0.002° ≈ 250 m

# Convertir en dataframe
df <- as.data.frame(r_wgs84, xy = TRUE) |>
  rename(lon = x, lat = y, sur_refl_b01 = 3)

# Appliquer le facteur d'échelle MODIS (obligatoire !)
# Les valeurs brutes sont des entiers, il faut diviser par 10000
df <- df |>
  mutate(sur_refl_b01 = sur_refl_b01 / 10000)

# Calculer l'aire pixel à partir de la résolution reprojetée
lat_ref    <- 43
res_lon_km <- res(r_wgs84)[1] * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res(r_wgs84)[2] * 111
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 6), "km²\n")

# define 95ème percentile -------------------------------------------------

# Calculer le 95ème percentile
seuil_95_Teng_MO <- quantile(study_area_df_2024$SPM_Teng_MO, 0.95, na.rm = TRUE)
seuil_95_Teng_MM <- quantile(study_area_df_2024$SPM_Teng_MM, 0.95, na.rm = TRUE)
seuil_95_Teng_extrm_MM <- quantile(study_area_df_2024$SPM_Teng_extrm_MM, 0.95, na.rm = TRUE)
seuil_95_Doxaran <- quantile(study_area_df_2024$SPM_Doxaran, 0.95, na.rm = TRUE)
seuil_95_Tsapanou <- quantile(study_area_df_2024$SPM_Tsapanou, 0.95, na.rm = TRUE)
seuil_95_Nechad <- quantile(study_area_df_2024$SPM_Nechad, 0.95, na.rm = TRUE)
seuil_95_Morin_Var <- quantile(study_area_df_2024$SPM_Morin_Var, 0.95, na.rm = TRUE)
seuil_95_Morin_Paillon <- quantile(study_area_df_2024$SPM_Morin_Paillon, 0.95, na.rm = TRUE)

cat("Seuil 95ème percentile :", seuil_95_Teng_MO, "g/m³\n") # 43.21933 g/m³
cat("Seuil 95ème percentile :", seuil_95_Teng_MM, "g/m³\n") # 386.1905 g/m³
cat("Seuil 95ème percentile :", seuil_95_Teng_extrm_MM, "g/m³\n") # 808.5443 g/m³
cat("Seuil 95ème percentile :", seuil_95_Doxaran, "g/m³\n") # 207136.5 g/m³
cat("Seuil 95ème percentile :", seuil_95_Tsapanou, "g/m³\n") # 38.54169 g/m³
cat("Seuil 95ème percentile :", seuil_95_Nechad, "g/m³\n") # 5.552708 g/m³
cat("Seuil 95ème percentile :", seuil_95_Morin_Var, "g/m³\n") # 10.56593 g/m³
cat("Seuil 95ème percentile :", seuil_95_Morin_Paillon, "g/m³\n") # 4.028886 g/m³

# Stats du panache par jour
# study_area_df_2024_95 <- study_area_df_2024 |>
#   group_by(date) |>
#   summarise(
#     # Teng_MO
#     pixel_count_Teng_MO = sum(SPM_Teng_MO >= seuil_95_Teng_MO, na.rm = TRUE),
#     mean_spm_Teng_MO = mean(SPM_Teng_MO[SPM_Teng_MO >= seuil_95_Teng_MO], na.rm = TRUE),
#     median_spm_Teng_MO = median(SPM_Teng_MO[SPM_Teng_MO >= seuil_95_Teng_MO], na.rm = TRUE),
#     aire_panache_km2_Teng_MO = pixel_count_Teng_MO * aire_pixel_km2,
#     # Teng_MM
#     pixel_count_Teng_MM = sum(SPM_Teng_MM >= seuil_95_Teng_MM, na.rm = TRUE),
#     mean_spm_Teng_MM = mean(SPM_Teng_MM[SPM_Teng_MM >= seuil_95_Teng_MM], na.rm = TRUE),
#     median_spm_Teng_MM = median(SPM_Teng_MM[SPM_Teng_MM >= seuil_95_Teng_MM], na.rm = TRUE),
#     aire_panache_km2_Teng_MM = pixel_count_Teng_MM * aire_pixel_km2,
#     # Teng_extreme_MM
#     pixel_count_Teng_extrm_MM = sum(SPM_Teng_extrm_MM >= seuil_95_Teng_extrm_MM, na.rm = TRUE),
#     mean_spm_Teng_extrm_MM = mean(SPM_Teng_extrm_MM[SPM_Teng_extrm_MM >= seuil_95_Teng_extrm_MM], na.rm = TRUE),
#     median_spm_Teng_extrm_MM = median(SPM_Teng_extrm_MM[SPM_Teng_extrm_MM >= seuil_95_Teng_extrm_MM], na.rm = TRUE),
#     aire_panache_km2_Teng_extrm_MM = pixel_count_Teng_extrm_MM * aire_pixel_km2,
#     # Doxaran
#     pixel_count_Doxaran = sum(SPM_Doxaran >= seuil_95_Doxaran, na.rm = TRUE),
#     mean_spm_Doxaran = mean(SPM_Doxaran[SPM_Doxaran >= seuil_95_Doxaran], na.rm = TRUE),
#     median_spm_Doxaran = median(SPM_Doxaran[SPM_Doxaran >= seuil_95_Doxaran], na.rm = TRUE),
#     aire_panache_km2_Doxaran = pixel_count_Doxaran * aire_pixel_km2,
#     # Tsapanou
#     pixel_count_Tsapanou = sum(SPM_Tsapanou >= seuil_95_Tsapanou, na.rm = TRUE),
#     mean_spm_Tsapanou = mean(SPM_Tsapanou[SPM_Tsapanou >= seuil_95_Tsapanou], na.rm = TRUE),
#     median_spm_Tsapanou = median(SPM_Tsapanou[SPM_Tsapanou >= seuil_95_Tsapanou], na.rm = TRUE),
#     aire_panache_km2_Tsapanou = pixel_count_Tsapanou * aire_pixel_km2,
#     # Nechad
#     pixel_count_Nechad = sum(SPM_Nechad >= seuil_95_Nechad, na.rm = TRUE),
#     mean_spm_Nechad = mean(SPM_Nechad[SPM_Nechad >= seuil_95_Nechad], na.rm = TRUE),
#     median_spm_Nechad = median(SPM_Nechad[SPM_Nechad >= seuil_95_Nechad], na.rm = TRUE),
#     aire_panache_km2_Nechad = pixel_count_Nechad * aire_pixel_km2,
#     # Morin Var
#     pixel_count_Morin_Var = sum(SPM_Morin_Var >= seuil_95_Morin_Var, na.rm = TRUE),
#     mean_spm_Morin_Var = mean(SPM_Morin_Var[SPM_Morin_Var >= seuil_95_Morin_Var], na.rm = TRUE),
#     median_spm_Morin_Var = median(SPM_Morin_Var[SPM_Morin_Var >= seuil_95_Morin_Var], na.rm = TRUE),
#     aire_panache_km2_Morin_Var = pixel_count_Morin_Var * aire_pixel_km2,
#     # Morin Paillon
#     pixel_count_Morin_Paillon = sum(SPM_Morin_Paillon >= seuil_95_Morin_Paillon, na.rm = TRUE),
#     mean_spm_Morin_Paillon = mean(SPM_Morin_Paillon[SPM_Morin_Paillon >= seuil_95_Morin_Paillon], na.rm = TRUE),
#     median_spm_Morin_Paillon = median(SPM_Morin_Paillon[SPM_Morin_Paillon >= seuil_95_Morin_Paillon], na.rm = TRUE),
#     aire_panache_km2_Morin_Paillon = pixel_count_Morin_Paillon * aire_pixel_km2
#   )

pixels_panache_2024 <- study_area_df_2024 |>
  mutate(
    SPM_panache_Teng_MO       = ifelse(SPM_Teng_MO >= seuil_95_Teng_MO, SPM_Teng_MO, NA),
    SPM_panache_Teng_MM       = ifelse(SPM_Teng_MM >= seuil_95_Teng_MM, SPM_Teng_MM, NA),
    SPM_panache_Teng_extrm_MM = ifelse(SPM_Teng_extrm_MM >= seuil_95_Teng_extrm_MM, SPM_Teng_extrm_MM, NA),
    SPM_panache_Doxaran       = ifelse(SPM_Doxaran >= seuil_95_Doxaran, SPM_Doxaran, NA),
    SPM_panache_Tsapanou      = ifelse(SPM_Tsapanou >= seuil_95_Tsapanou, SPM_Tsapanou, NA),
    SPM_panache_Nechad        = ifelse(SPM_Nechad >= seuil_95_Nechad, SPM_Nechad, NA),
    SPM_panache_Morin_Var     = ifelse(SPM_Morin_Var >= seuil_95_Morin_Var, SPM_Morin_Var, NA),
    SPM_panache_Morin_Paillon = ifelse(SPM_Morin_Paillon >= seuil_95_Morin_Paillon, SPM_Morin_Paillon, NA)
  )

save(study_area_df_2024_95, file = "data/MODIS L2 NASA/study_area_df_2024_95.Rdata")

## graphiques --------------------------------------------------------------

# pour savoir si ça identifie correctement les panaches il faut que je plotte 
# le seuil au 95ème percentile pour une journée

# Filtrer un jour
# pixels_04_03_2024 <- pixels_panache_2024 |> filter(date == "2024-03-04")

pixels_01_04_2024 <- pixels_panache_2024 |> filter(date == "2024-04-01")

# Créer le graphique
pl_map <- pixels_01_04_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_panache_Nechad)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) + # ← ici
  scale_fill_viridis_c(
    option = "plasma",
    name = "SPM [mg/l]",
    limits = c(0, max_spm)  # Utilise max_spm calculé précédemment
  ) +
  guides(fill = guide_colorbar(
    barwidth = 20,
    barheight = 2,
    title.position = "top",
    title.hjust = 0.5
  )) +
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°N)",
    fill = "SPM [mg/l]"
  ) +
  coord_sf(                                   
    xlim = range(study_area_df_04_03_2024_95$lon), 
    ylim = range(study_area_df_04_03_2024_95$lat), 
    expand = FALSE
  ) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2 2024 Aqua/SPM/95ème percentile/fig_MODIS_SPM_01_04_2024_Nechad_95.png", pl_map, height = 9, width = 14)

### plotter tout en même temps ----------------------------------------------

#### mean_spm ----------------------------------------------------------------

ggplot() +
  geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Teng_MO, color = "Teng MO")) +
  geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Teng_MM, color = "Teng MM")) +
  geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Teng_extrm_MM, color = "Teng extrême MM")) +
  # geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Doxaran, color = "Doxaran")) +
  # geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Tsapanou, color = "Tsapanou")) +
  # geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Nechad, color = "Nechad")) +
  geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Morin_Var, color = "Morin Var")) +
  geom_point(data = study_area_df_2024_95, aes(x= date, y = mean_spm_Morin_Paillon, color = "Morin Paillon")) +
  labs(
    title = "Prédiction de la concentration en MES selon les différents algorithmes - 2024",
    subtitle = "Prédictions des MES - MODIS niveau 2",
    x = "date",
    y = expression(SPM["prédictions"] ~ (mg ~ L^{-1}))
  ) +
  
  # Thème
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey50"),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey70")
  )
  
  
#### Teng MO -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Teng_MO, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Teng_MO <- study_area_df_2024_95$aire_panache_km2_Teng_MO * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Teng_MO, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en g/m³)") +
      labs(
      title    = "Prédiction de l'extension des panaches turbides et du débit cumulé",
      subtitle = "Algorithme développé par Teng et al., 2025 pour les eaux riches en Matière Organique",
      x        = expression(SPM[Morin] ~ (mg ~ L^{-1})),
      y        = expression(SPM["in situ"] ~ (mg ~ L^{-1}))
    ) +
      theme(
        plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle    = element_text(size = 14, hjust = 0.5, color = "grey50"),
        axis.title       = element_text(face = "bold"),
        axis.text        = element_text(color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "grey70")
      ))
    

#### Teng MM -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Teng_MM, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Teng_MM <- study_area_df_2024_95$aire_panache_km2_Teng_MM * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Teng_MO, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, 
                        name = "Aire des panaches (en km²)")  # ← parenthèse fermée ici
  ) +                                                          # ← et ici
  labs(
    title    = "Prédiction de l'extension des panaches turbides et du débit cumulé",
    subtitle = "Algorithme développé par Teng et al., 2025 pour les eaux riches en Matière Organique",
    x        = NULL,
    color    = NULL
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(face = "bold", size = 13, margin = margin(b = 10)),
    plot.subtitle     = element_text(size = 11, color = "grey50", margin = margin(b = 10)),
    axis.title.y      = element_text(size = 13, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 13, margin = margin(l = 10)),
    axis.text         = element_text(size = 11, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 11)
  )

#### Teng extrêmement riche en MM -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Teng_extrm_MM, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Teng_extrm_MM <- study_area_df_2024_95$aire_panache_km2_Teng_extrm_MM * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Teng_extrm_MM, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Teng et al., 2025 pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

#### Doxaran -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Doxaran, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Doxaran <- study_area_df_2024_95$aire_panache_km2_Doxaran * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Doxaran, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Doxaran et al., 2009 pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

#### Tsapanou -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Tsapanou, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Tsapanou <- study_area_df_2024_95$aire_panache_km2_Tsapanou * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Tsapanou, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Tsapanou et al., 2020 pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

#### Nechad -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Nechad, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Nechad <- study_area_df_2024_95$aire_panache_km2_Nechad * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Nechad, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Nechad et al., ... pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

#### Morin Var -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Morin_Var, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Morin_Var <- study_area_df_2024_95$aire_panache_km2_Morin_Var * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Morin_Var, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Morin et al., ... pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

#### Morin Paillon -----------------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(study_area_df_2024_95$aire_panache_km2_Morin_Paillon, All_debit_2024$debit_cumule)

study_area_df_2024_95$scaled_aire_panache_km2_Morin_Paillon <- study_area_df_2024_95$aire_panache_km2_Morin_Paillon * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = All_debit_2024, 
            aes(x = date, y = debit_cumule, color = "Débit cumulé"), size = 0.3) +
  geom_line(data = study_area_df_2024_95, 
            aes(x = date, y = scaled_aire_panache_km2_Morin_Paillon, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit MMDIS L2",
       caption = "Algorithme développé par Morin et al., ... pour les eaux très riches en MM",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# Comparaison sextant vs MODIS data ---------------------------------------

## plume area --------------------------------------------------------------

SEXTANT_2024_spm_95 <- SEXTANT_1998_2025_spm_95 |> 
  filter(date >= as.Date("2024-01-01"), date <= as.Date("2024-12-31"))

ggplot() +
  geom_line(data = SEXTANT_2024_spm_95,
            aes(x = date, y = aire_panache_km2, color = "SEXTANT OC5"),
            linewidth = 0.5, alpha = 0.8) +
  geom_line(data = study_area_df_2024_95,
            aes(x = date, y = aire_panache_km2_Morin_Paillon, color = "MODIS — Morin Paillon"),
            linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("SEXTANT OC5" = "deepskyblue3", "MODIS — Morin Paillon" = "maroon1")
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b",
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title    = "Extension des panaches turbides en 2024",
    subtitle = "Comparaison SEXTANT OC5 vs MODIS L2 — Algorithme Morin et al.",
    x        = NULL,
    y        = expression("Aire des panaches (km"^2*")"),
    color    = NULL,
    caption  = "Source : SEXTANT OC5 | MODIS-Aqua | Seuil : 95ème percentile"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0,
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 12, color = "grey50", hjust = 0,
                                     margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 12, margin = margin(r = 10)),
    axis.text         = element_text(size = 11, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(1.5, "cm")  # allonge les traits dans la légende
  )


# reflectance analysis -------------------------------------------------------------------------

## 2024 --------------------------------------------------------------------

reflectance_2024 <- study_area_df_2024 |> 
  summarise(
    moy_reflec_b1 = mean(sur_refl_b01_1),
    moy_reflec_b2 = mean(sur_refl_b02_1),
    std_reflec_b1 = sd(sur_refl_b01_1),
    std_reflec_b2 = sd(sur_refl_b02_1),
    .by = date)

adjust_factors <- sec_axis_adjustement_factors(reflectance_2024$moy_reflec_b1, All_debit_2024$debit_cumule)

reflectance_2024$scaled_moy_reflec_b1 <- reflectance_2024$moy_reflec_b1 * adjust_factors$diff + adjust_factors$adjust

# Calcul de la corrélation entre débit et la réflectance
merged_data <- merge(
  All_debit_2024,
  reflectance_2024,
  by = "date",
  all = FALSE
)

correlation <- cor(merged_data$debit_cumule, merged_data$moy_reflec_b1, method = "spearman", use = "complete.obs")
p_value <- cor.test(merged_data$debit_cumule, merged_data$moy_reflec_b1, method = "spearman")$p.value


ggplot() +
  geom_line(data = reflectance_2024,
            aes(x = date, y = scaled_moy_reflec_b1, color = "Réflectance MODIS"),
            linewidth = 0.5, alpha = 0.8) +
  geom_line(data = All_debit_2024,
            aes(x = date, y = debit_cumule, color = "Débit cumulé"),
            linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("Réflectance MODIS" = "maroon1", "Débit cumulé" = "darkolivegreen3")
  ) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Réflectance (en %)")
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b",
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  # Annotation pour la corrélation (en haut à droite)
  annotate(
    "text",
    x = max(c(All_debit_2024$date, reflectance_2024$date), na.rm = TRUE),
    y = max(c(All_debit_2024$débit, reflectance_2024$moy_reflec_b1), na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "R = ", round(correlation, 2),
      "\n", "p ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  labs(
    title    = "Réflectance moyenne de la bande 1 MODIS contre débit liquide cumulé",
    subtitle = "Débit liquide et réflectance MODIS (645 nm)",
    x        = NULL,
    y        = "Réflectance (%)",
    color    = NULL,
    caption  = "Source : Hydro France | MODIS-Aqua"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0,
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 12, color = "grey50", hjust = 0,
                                     margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 12, margin = margin(r = 10)),
    axis.text         = element_text(size = 11, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(1.5, "cm")  # allonge les traits dans la légende
  )

unique(study_area_df_2024$date)



