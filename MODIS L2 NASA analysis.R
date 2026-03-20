# MODIS L2 NASA analysis

# 16/03/2026
# pathway : "~/Satellite_analysis/MODIS L2 NASA analysis.R


# This script will load MODIS L2 data from NASA (first on the 03/10/2020 to
# see if everything is fine) and will try, thanks to the reflectance parameter,
# to retrieve SPM concentration of the study area

# Setup ------------------------------------------------------------------

source("~/Satellite_analysis/earth_data_access.R")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2)
library(heatwaveR)
library(ggpmisc)

# function ---------------------------------------------------------------

# This function will convert reflectance into SPM concentration thanks to ...

# Var_Morin <- function(study_area_df, sur_refl_b01_1, A = 80, C = 0.1562) {
# 
#   # Vérifier que la colonne de réflectance existe
#   if (!(sur_refl_b01_1 %in% names(study_area_df))) {
#     stop(paste("La colonne", sur_refl_b01_1, "n'existe pas dans le data frame."))
#   }
# 
#   # Calculer SPM avec la formule : SPM = A * ρw / (1 - ρw / C)
#   study_area_df <- study_area_df %>%
#     mutate(SPM = pmin((A * .data[[sur_refl_b01_1]]) / (1 - .data[[sur_refl_b01_1]] / C), 50)
#     )
# 
#   return(study_area_df)
# }

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
  
  return(study_area_df)
}

# data analysis -----------------------------------------------------------

## apply function to data set -----------------------------------------------------------

# MODIS_L2_SPM <- function(study_area_df, sur_refl_b01_1 = "sur_refl_b01_1", A = 80, C = 0.1562) {
#   study_area_df <- study_area_df %>%
#     mutate(SPM = case_when(
#       .data[[sur_refl_b01_1]] >= C ~ NA_real_,  # Évite les valeurs infinies
#       TRUE ~ (A * .data[[sur_refl_b01_1]]) / pmax(1 - .data[[sur_refl_b01_1]] / C, 0.01)  # Évite la division par zéro
#     ))
#   return(study_area_df)
# }

# we have to filter the data frame because there are some absurd values
# study_area_df <- study_area_df %>%
#   mutate(sur_refl_b01_1 = ifelse(sur_refl_b01_1 >= 0.5, 0.5, sur_refl_b01_1))

study_area_df <- study_area_df |>
  filter(sur_refl_b01_1 <= 0.6)

C <- 0.1562  # Définis la valeur de C
study_area_df <- study_area_df %>%
  mutate(sur_refl_b01_1 = pmin(sur_refl_b01_1, C - 0.01))  # Limite la réflectance à C - 0.01


study_area_df <- Var_Morin(study_area_df, sur_refl_b01_1 = "sur_refl_b01_1")

# we have to filter the data frame because there are some absurd values
# study_area_df_filtered <- study_area_df %>%
#   filter(SPM >= 0,
#          SPM <= 10000)

## plotting -----------------------------------------------------------

max_spm <- max(study_area_df$SPM, na.rm = TRUE)

# Créer le graphique
pl_map <- study_area_df %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = SPM)) +
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
  coord_quickmap(xlim = range(study_area_df$lon, na.rm = TRUE),
                 ylim = range(study_area_df$lat, na.rm = TRUE)) +
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
ggsave("~/Downloads/MODIS NASA/L2/fig_MODIS_SPM.png", pl_map, height = 9, width = 14)

