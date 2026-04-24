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
study_area_df_clean_2024 <- study_area_df_clean_2024 |>  
  mutate(SPM_Morin_Var = (80 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.1562)),
         SPM_Morin_Paillon = (39 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.2563)))

# Équation Doxaran et al., 2009 (il faut les deux bandes réflectance)
study_area_df_clean_2024 <- study_area_df_clean_2024 |>
  mutate(SPM_Doxaran = 12.996 * exp((sur_refl_b02_1/sur_refl_b01_1)/0.189),
         SPM_Doxaran = ifelse(is.infinite(SPM_Doxaran), NA, SPM_Doxaran))

# A different algorithm based on Teng et al. 2025
# https://www.sciencedirect.com/science/article/pii/S003442572500149X
# SPM_org = a Rrs(lambda_RED)^b; a = 1992.2, b = 1.027
# NB: lambda_RED is taken here to be the MODIS band 1 waveband
study_area_df_clean_2024 <- study_area_df_clean_2024 |> 
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
# A = 289.29 g m−3 , B = 2.10 g m−3 and C = 0.1686
study_area_df_clean_2024 <- study_area_df_clean_2024 |> 
  mutate(SPM_Tsapanou = ((289.29 * sur_refl_b01_1)/(1-(sur_refl_b01_1/0.1686)))+2.10)
# This produces too many negative values...

# So we digress to the Nechad formula of 
# SPM = ((A * Rrs)/(1-(Rrs/C)))+B
# A ≈ 200-230, C ≈ 0.15-⁣0.17 B ≈0# Just as a starting guess
study_area_df_clean_2024 <- study_area_df_clean_2024 |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
        SPM_Nechad = ((200 * Rrs_b01_01)/(1-(Rrs_b01_01/0.17)))+0)

# OR we can try : 
# SPM[mg l−1]= a + b * ρsurf,645
# ρsurf,645 = sur_refl_b01_1
# a=289.29, b=2.1 # For starting
study_area_df_clean_2024 <- study_area_df_clean_2024 |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(SPM_formule = 0 + 2.1 * sur_refl_b01_1)

# plotting -----------------------------------------------------------

# we have to choose some dates and look at what data look like : 
study_area_df_04_03_2024 <- study_area_df_clean_2024 |> 
  filter(date == "2024-03-04")

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

# set the maximum values
max_spm <- max(study_area_df_04_03_2024$SPM_Morin_Var, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Morin_Var)) +
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Morin_Var.png", pl_map, height = 9, width = 14)

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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Morin_Paillon.png", pl_map, height = 9, width = 14)

## Teng MO -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_MO, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_MO)) +
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Teng_MO.png", pl_map, height = 9, width = 14)

## Teng MM -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_MM, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_MM)) +
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Teng_MM.png", pl_map, height = 9, width = 14)

## Teng extrem MM -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Teng_extrm_MM, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Teng_extrm_MM)) +
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Teng_extrem_MM.png", pl_map, height = 9, width = 14)

## Tsapanou -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Tsapanou, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Tsapanou.png", pl_map, height = 9, width = 14)

## Nechad -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_Nechad, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_Nechad.png", pl_map, height = 9, width = 14)

## Formule -------------------------------------------------------------------

max_spm <- max(study_area_df_04_03_2024$SPM_formule, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024 %>%
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1/fig_MODIS_SPM_04_03_2024_formule.png", pl_map, height = 9, width = 14)

## Doxaran -------------------------------------------------------------------

# First, when we do the map there are some extreme values (until 8000), which is weird
# we want to erase them so maybe the map will be better
study_area_df_04_03_2024_1 <- study_area_df_04_03_2024 |>
  filter(SPM_Doxaran <= 2000)

# set the maximum values
max_spm <- max(study_area_df_04_03_2024_1$SPM_Doxaran, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df_04_03_2024_1 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM_Doxaran)) +
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
    xlim = range(study_area_df_04_03_2024_1$lon), 
    ylim = range(study_area_df_04_03_2024_1$lat), 
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
ggsave("~/Downloads/MODIS NASA/L2 2016 Terra/SPM/bande 1 et 2/fig_MODIS_SPM_04_03_2024_Doxaran_1.png", pl_map, height = 9, width = 14)

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
  
