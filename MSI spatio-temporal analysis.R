# MSI spatio-temporal analysis

# 02/03/2026

# pathway : 

# This script will load MSI OCEANCOLOUR_MED_BGC_HR_L3_NRT_009_205 product and particulary
# the layer cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(stars)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# function ----------------------------------------------------------------

# Fonction pour extraire les données d'un objet stars
extract_stars_to_df <- function(stars_obj, year) {
  # Convertir l'objet stars en un objet sf
  sf_obj <- st_as_sf(stars_obj)
  
  # Convertir l'objet sf en un data frame
  df <- st_drop_geometry(sf_obj)
  
  # Ajouter les informations temporelles
  n_time <- dim(stars_obj)[3]
  doy <- 1:n_time  # Jour de l'année (1 à 365 ou 366)
  
  # Répéter chaque ligne pour chaque jour de l'année
  df <- df %>%
    slice(rep(1:n(), each = n_time)) %>%
    mutate(
      doy = rep(doy, each = n()),
      date = as.Date(paste(year, "-01-01", sep = "")) + (doy - 1),
      year = year,
      month = month(date),
      day = day(date)
    )
  
  return(df)
}

# loading files list ------------------------------------------------------

# loading files
load("~/Downloads/MSI/SPM/MSI_2020_SPM.RData")
load("~/Downloads/MSI/SPM/MSI_2021_SPM.RData")
load("~/Downloads/MSI/SPM/MSI_2022_SPM.RData")
load("~/Downloads/MSI/SPM/MSI_2023_SPM.RData")
load("~/Downloads/MSI/SPM/MSI_2024_SPM.RData")
load("~/Downloads/MSI/SPM/MSI_2025_SPM.RData")

# Créer une liste avec les objets chargés
# MSI_2020_2025 <- list(
#   MSI_2020 = download_2020_SPM,
#   MSI_2021 = download_2021_SPM,
#   MSI_2022 = download_2022_SPM,
#   MSI_2023 = download_2023_SPM,
#   MSI_2024 = download_2024_SPM,
#   MSI_2025 = download_2025_SPM
# )

load("~/Downloads/MSI/SPM/MSI_2020_2025.RData")

# Liste des années et des objets stars correspondants
years <- c(2020, 2021, 2022, 2023, 2024, 2025)
stars_list <- MSI_2020_2025  # Ta liste d'objets stars

# Extraire chaque année en data frame
MSI_SPM_2020_2025 <- bind_rows(
  mapply(extract_stars_to_df, stars_list, years, SIMPLIFY = FALSE)
)

# Vérifier le résultat
glimpse(download_2020_SPM)
summary(download_2020_SPM)







# temporal analysis -------------------------------------------------------

SPM_2025_2025 <- rbind(ma_liste_SPM, .id = "source") %>%
  mutate(date = as.Date(date))  # S'assurer que "date" est bien un objet Date

