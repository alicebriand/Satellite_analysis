# MODIS L2 NASA analysis

# 16/03/2026
# pathway : "~/Satellite_analysis/MODIS L2 NASA analysis.R

# This script will load MODIS L2 data from NASA and will try, thanks to the reflectance parameter,
# to retrieve SPM concentration of the study area with different formula. The 
# ultimate goal is to identify which algorithm is the best to retrieve SPM concentration

# Setup ------------------------------------------------------------------

# NB: When running a script with source(), it will performa ll of the actions therein
# This is not always ideal, especially if the script is downloading data or performing large analyses.
# In this case it is better to save the output of the first script, then load them in the second script
# source("~/Satellite_analysis/earth_data_access.R")

# The shared functions between scripts are however and excellent reason to use source()
# Therefore I have moved the shared functions to a script that is source()'d here
source("func.R")

load("~/River_runoff_analysis/data/Hydro France/Var_crues.Rdata")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(terra)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2)
library(heatwaveR)
library(ggpmisc)

# function ---------------------------------------------------------------

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

# It is a single equation, we can apply it directly to the data.frame with mutate()
# SPM = A * ρw / (1 - ρw / C); A = 80, C = 0.1562 # But where does this equation 
# and values come from? I do not find them in the literature?
study_area_df_clean <- study_area_df_clean |>  
  mutate(SPM_Morin = (80 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.1562)))

# A different algorithm based on Teng et al. 2025
# https://www.sciencedirect.com/science/article/pii/S003442572500149X
# SPM_org = a Rrs(lambda_RED)^b; a = 1992.2, b = 1.027
# NB: lambda_RED is taken here to be the MODIS band 1 waveband
study_area_df_clean <- study_area_df_clean |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
       SPM_Teng = 1992.2 * Rrs_b01_01^1.027)
# But these values are crazy high...

# So then this paper by Tsapanou et al. 2020
# http://www.teiath.gr/userfiles/pdrak/lab/coupling_remote_sensing_data.pdf
# Though this is for LandSat 8
# SPM = ((A * Rho_W)/(1-(Rhow_w/C)))+B
# A = 289.29 g m−3 , B = 2.10 g m−3 and C = 0.1686
study_area_df_clean <- study_area_df_clean |> 
  mutate(SPM_Tsapanou = ((289.29 * sur_refl_b01_1)/(1-(sur_refl_b01_1/0.1686)))+2.10)
# This produces too many negative values...

# So we digress to the Nechad formula of 
# SPM = ((A * Rrs)/(1-(Rrs/C)))+B
# A ≈ 200-230, C ≈ 0.15-⁣0.17 B ≈0# Just as a starting guess
study_area_df_clean <- study_area_df_clean |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
        SPM_Nechad = ((200 * Rrs_b01_01)/(1-(Rrs_b01_01/0.17)))+0)

# OR we can try
# SPM[mg l−1]= a + b * ρsurf,645
# ρsurf,645 = sur_refl_b01_1
# a=289.29,b=2.1 # For starting
study_area_df_clean <- study_area_df_clean |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(SPM_formule = 0 + 2.1 * sur_refl_b01_1)

# # We have to filter the data frame because there are some absurd values
# study_area_df <- study_area_df |> 
#   mutate(sur_refl_b01_1 = ifelse(sur_refl_b01_1 >= 0.5, 0.5, sur_refl_b01_1))

## plotting -----------------------------------------------------------

# we have to choose some dates and look at what data look like : 
study_area_df_04_03_2024 <- study_area_df_clean |> 
  filter(date == "2024-03-04")

## Morin -------------------------------------------------------------------

# First, when we do the map there are some extreme values (until 8000), which is weird
# we want to erase them so maybe the map will be better

study_area_df_04_03_2024_Morin <- study_area_df_04_03_2024 |> 
  filter(SPM_Morin <= 100)

max_spm <- max(study_area_df_04_03_2024_Morin$SPM_Morin, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_Morin %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Morin)) +
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
    xlim = range(study_area_df_04_03_2024$lon), 
    ylim = range(study_area_df_04_03_2024$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2024/SPM/fig_MODIS_SPM_04_03_2024_Morin1.png", pl_map, height = 9, width = 14)

## Teng -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_Morin %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng)) +
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
    xlim = range(study_area_df_04_03_2024$lon), 
    ylim = range(study_area_df_04_03_2024$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2024/SPM/fig_MODIS_SPM_04_03_2024_Teng.png", pl_map, height = 9, width = 14)

## Tsapanou -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Tsapanou, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_Morin %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Tsapanou)) +
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
    xlim = range(study_area_df_04_03_2024$lon), 
    ylim = range(study_area_df_04_03_2024$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2024/SPM/fig_MODIS_SPM_04_03_2024_Tsapanou.png", pl_map, height = 9, width = 14)

## Nechad -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Nechad, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_Morin %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Nechad)) +
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
    xlim = range(study_area_df_04_03_2024$lon), 
    ylim = range(study_area_df_04_03_2024$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2024/SPM/fig_MODIS_SPM_04_03_2024_Nechad.png", pl_map, height = 9, width = 14)

## Formule -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_formule, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_Morin %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_formule)) +
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
    xlim = range(study_area_df_04_03_2024$lon), 
    ylim = range(study_area_df_04_03_2024$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2024/SPM/fig_MODIS_SPM_04_03_2024_formule.png", pl_map, height = 9, width = 14)
