# MODIS L2 NASA analysis

# 16/03/2026
# pathway : "~/Satellite_analysis/MODIS L2 NASA analysis.R


# This script will load MODIS L2 data from NASA (first on the 03/10/2020 to
# see if everything is fine) and will try, thanks to the reflectance parameter,
# to retrieve SPM concentration of the study area

# Setup ------------------------------------------------------------------

# NB: When running a script with source(), it will performa ll of the actions therein
# This is not always ideal, especially if the script is downloading data or performing large analyses.
# In this case it is better to save the output of the first script, then load them in the second script
# source("~/Satellite_analysis/earth_data_access.R")

# The shared functions between scripts are hwoever and excellent reason to use source()
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

# This function will convert reflectance into SPM concentration thanks to ...

<<<<<<< HEAD
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
=======
# NB: When writing a function it is a good idea (but not necessary) to name your arguments
# in a way that is not found in your code outside of the function
# E.g. Here I changed 'study_area_df' to 'df' and 'sur_refl_b01_1' to 'col_name'
MODIS_L2_SPM <- function(df, col_name, A = 80, C = 0.1562) {
>>>>>>> origin/main~
  
  # Vérifier que la colonne de réflectance existe
  if (!(col_name %in% names(df))) {
    stop(paste("La colonne", col_name, "n'existe pas dans le data frame."))
  }
  
  # Calculer SPM avec la formule : SPM = A * ρw / (1 - ρw / C)
<<<<<<< HEAD
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
=======
  df <- df |> 
    mutate(SPM = (A * .data[[col_name]]) / (1 - (.data[[col_name]] / C)))
>>>>>>> origin/main~
  
  return(df)
}

# data analysis -----------------------------------------------------------


## Load data ---------------------------------------------------------------

# Load the MODIS mask first
# Change the filename if this is not correct
# MODIS_mask <- rast("~/data/MODIS/study_area_MOD44W_2020-01-01.tif") # File path on Robert's computer...
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2/masks/study_area_MOD44W_2020-01-01.tif")

# Load one day of data
study_area_df <- load_MODIS_tif(file_name = "~/data/MODIS/study_area_MYD09GQ_2020-10-03.tif", MODIS_mask)


## apply function to data set -----------------------------------------------------------

# MODIS_L2_SPM <- function(study_area_df, sur_refl_b01_1 = "sur_refl_b01_1", A = 80, C = 0.1562) {
#   study_area_df <- study_area_df %>%
#     mutate(SPM = case_when(
#       .data[[sur_refl_b01_1]] >= C ~ NA_real_,  # Évite les valeurs infinies
#       TRUE ~ (A * .data[[sur_refl_b01_1]]) / pmax(1 - .data[[sur_refl_b01_1]] / C, 0.01)  # Évite la division par zéro
#     ))
#   return(study_area_df)

# Apply SDM equation
# study_area_df <- MODIS_L2_SPM(study_area_df, col_name = "sur_refl_b01_1") 

# Or, because it is a single equation, we can apply it directly to the data.frame with mutate()
# SPM = A * ρw / (1 - ρw / C); A = 80, C = 0.1562 # But where does this equation and values come from? I do not find them in the literature?
# study_area_df <- study_area_df |> 
#   mutate(SPM = (80 * sur_refl_b01_1) / (1 - (sur_refl_b01_1 / 0.1562)))

# A different algorithm based on Teng et al. 2025
# https://www.sciencedirect.com/science/article/pii/S003442572500149X
# SPM_org = a Rrs(lambda_RED)^b; a = 1992.2, b = 1.027
# NB: lambda_RED is taken here to be the MODIS band 1 waveband
study_area_df <- study_area_df |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
         SPM = 1992.2 * Rrs_b01_01^1.027)
# But these values are crazy high...

# So then this paper by Tsapanou et al. 2020
# http://www.teiath.gr/userfiles/pdrak/lab/coupling_remote_sensing_data.pdf
# Though this is for LandSat 8
# SPM = ((A * Rho_W)/(1-(Rhow_w/C)))+B
# A = 289.29 g m−3 , B = 2.10 g m−3 and C = 0.1686
study_area_df <- study_area_df |> 
  mutate(SPM = ((289.29 * sur_refl_b01_1)/(1-(sur_refl_b01_1/0.1686)))+2.10)
# This produces too many negative values...

# So we digress to the Nechad formula of 
# SPM = ((A * Rrs)/(1-(Rrs/C)))+B
# A ≈ 200-230, C ≈ 0.15-⁣0.17 B ≈0# Just as a starting guess
study_area_df <- study_area_df |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(Rrs_b01_01 = (sur_refl_b01_1/pi), # First convert Rhow_w to Rrs
         SPM = ((200 * Rrs_b01_01)/(1-(Rrs_b01_01/0.17)))+0)

# OR we can try
# SPM[mg l−1]= a + b * ρsurf,645
# ρsurf,645 = sur_refl_b01_1
# a=289.29,b=2.1 # For starting
study_area_df <- study_area_df |> 
  filter(sur_refl_b01_1 >= 0 ) |> 
  mutate(SPM = 0 + 2.1 * sur_refl_b01_1)

# We have to filter the data frame because there are some absurd values
# study_area_df <- study_area_df |> 
#   mutate(sur_refl_b01_1 = ifelse(sur_refl_b01_1 >= 0.5, 0.5, sur_refl_b01_1))

# we have to filter the data frame because there are some absurd values
<<<<<<< HEAD
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
=======
# study_area_df_filtered <- study_area_df 
study_area_df_filtered <- study_area_df |>
  filter(SPM >= 0,
         SPM <= 6000)
>>>>>>> origin/main~

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
pl_map

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2/fig_MODIS_SPM.png", pl_map, height = 9, width = 14)
