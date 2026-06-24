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
library(terra)
library(sf)

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



# 04/03/2024 --------------------------------------------------------------

list.files("~/Downloads/S2A_MSIL2A_20240304T102921_N0510_R108_T32TLP_20240313T160152.SAFE", recursive = TRUE)

# Chemin de base
safe_dir <- "~/Downloads/S2A_MSIL2A_20240304T102921_N0510_R108_T32TLP_20240313T160152.SAFE"

# Chemins des bandes à 10m
b02_path <- file.path(safe_dir, "GRANULE/L2A_T32TLP_A045436_20240304T102921/IMG_DATA/R10m/T32TLP_20240304T102921_B02_10m.jp2")
b03_path <- file.path(safe_dir, "GRANULE/L2A_T32TLP_A045436_20240304T102921/IMG_DATA/R10m/T32TLP_20240304T102921_B03_10m.jp2")
b04_path <- file.path(safe_dir, "GRANULE/L2A_T32TLP_A045436_20240304T102921/IMG_DATA/R10m/T32TLP_20240304T102921_B04_10m.jp2")
b08_path <- file.path(safe_dir, "GRANULE/L2A_T32TLP_A045436_20240304T102921/IMG_DATA/R10m/T32TLP_20240304T102921_B08_10m.jp2")

# Charger les bandes
b02 <- rast(b02_path)
b03 <- rast(b03_path)
b04 <- rast(b04_path)
b08 <- rast(b08_path)

# Vérifier le CRS (UTM 32N)
crs(b02)

# Définir la zone d'étude et la reprojeter en UTM
zone_sf <- st_bbox(c(xmin = 7.110978, xmax = 7.360000,
                     ymin = 43.523182, ymax = 43.730000),
                   crs = 4326) |>
  st_as_sfc() |>
  st_transform(crs(b02))

# Recadrer les bandes
b02_crop <- crop(b02, zone_sf)
b03_crop <- crop(b03, zone_sf)
b04_crop <- crop(b04, zone_sf)
b08_crop <- crop(b08, zone_sf)

# Convertir en réflectance de surface (diviser par 10000 pour Sentinel-2 L2A)
b02_refl <- b02_crop / 10000
b03_refl <- b03_crop / 10000
b04_refl <- b04_crop / 10000
b08_refl <- b08_crop / 10000

# Reprojeter en WGS84 pour faciliter le plot
b02_refl <- project(b02_refl, "EPSG:4326")
b03_refl <- project(b03_refl, "EPSG:4326")
b04_refl <- project(b04_refl, "EPSG:4326")
b08_refl <- project(b08_refl, "EPSG:4326")

# Vérifier les valeurs
summary(b04_refl)

# Convertir en dataframe pour ggplot
df_b02 <- as.data.frame(b02_refl, xy = TRUE) |> rename(lon = x, lat = y, refl_b02 = 3)
df_b03 <- as.data.frame(b03_refl, xy = TRUE) |> rename(lon = x, lat = y, refl_b03 = 3)
df_b04 <- as.data.frame(b04_refl, xy = TRUE) |> rename(lon = x, lat = y, refl_b04 = 3)
df_b08 <- as.data.frame(b08_refl, xy = TRUE) |> rename(lon = x, lat = y, refl_b08 = 3)

# Assembler en un seul dataframe
df_MSI <- df_b02 |>
  left_join(df_b03, by = c("lon", "lat")) |>
  left_join(df_b04, by = c("lon", "lat")) |>
  left_join(df_b08, by = c("lon", "lat"))

head(df_MSI)

