# ODATIS_MR_expert_analysis

# pathway : ~/Satellite_analysis/ODATIS_MR_expert_analysis/

# This script will load Odatis MR expert satellite data
# Then 

# library -----------------------------------------------------------------

library(tidyverse)
library(tidync)
library(ncdf4)    # For reading NetCDF files
library(lubridate) # For working with dates
library(reshape2) # For data reshaping
library(ggplot2)
library(rnaturalearth)
library(ggspatial)

# load data ---------------------------------------------------------------

# raw data

# MERIS (2009)
load("data/ODATIS-MR_expert/MERIS_SPM_G_AC_sub.RData")
load("data/ODATIS-MR_expert/MERIS_SPM_G_PO_sub.RData")
load("data/ODATIS-MR_expert/MERIS_SPM_R_AC_sub.RData")
load("data/ODATIS-MR_expert/MERIS_SPM_R_PO_sub.RData")

# MODIS (2019)
load("data/ODATIS-MR_expert/MODIS_SPM_G_NS_sub.RData")
load("data/ODATIS-MR_expert/MODIS_SPM_G_PO_sub.RData")
load("data/ODATIS-MR_expert/MODIS_SPM_R_NS_sub.RData")
load("data/ODATIS-MR_expert/MODIS_SPM_R_PO_sub.RData")

# OLCI A
load("data/ODATIS-MR_expert/OLCIA_SPM_G_AC_sub.RData")
load("data/ODATIS-MR_expert/OLCIA_SPM_G_PO_sub.RData")
load("data/ODATIS-MR_expert/OLCIA_SPM_R_AC_sub.RData")
load("data/ODATIS-MR_expert/OLCIA_SPM_R_PO_sub.RData")

# OLCI B (2019)
load("data/ODATIS-MR_expert/OLCIB_SPM_G_AC_sub.RData")
load("data/ODATIS-MR_expert/OLCIB_SPM_G_PO_sub.RData")
load("data/ODATIS-MR_expert/OLCIB_SPM_R_AC_sub.RData")
load("data/ODATIS-MR_expert/OLCIB_SPM_R_PO_sub.RData")

# merge data

OLCI_SPM_G_AC_sub <- bind_rows(OLCIA_SPM_G_AC_sub, OLCIB_SPM_G_AC_sub)
OLCI_SPM_G_PO_sub <- bind_rows(OLCIA_SPM_G_PO_sub, OLCIB_SPM_G_PO_sub)
# OLCI_SPM_R_AC_sub <- bind_rows(OLCIA_SPM_R_AC_sub, OLCIB_SPM_R_AC_sub)
# OLCI_SPM_R_PO_sub <- bind_rows(OLCIA_SPM_R_PO_sub, OLCIB_SPM_R_PO_sub)

load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_PO_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_PO_sub_95.Rdata")

# 95 percentile data

# MERIS
load("data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_PO_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_PO_sub_95.Rdata")

# MODIS
load("data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_NS_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_PO_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_NS_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_PO_sub_95.Rdata")

# OLCI
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_PO_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_AC_sub_95.Rdata")
load("data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_PO_sub_95.Rdata")

# Hydrological data

load("data/Hydro France/Y6442010_depuis_2000.Rdata")
load("data/Hydro France/Y6442010_2009.Rdata")
load("data/Hydro France/All_debit.Rdata")
load("data/Hydro France/All_debit_2019.Rdata")

# function ----------------------------------------------------------------

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

# set a threshold ---------------------------------------------------------

## MERIS -------------------------------------------------------------------

### calcul de l'aire --------------------------------------------------------

# Trier et extraire les valeurs uniques
lons_uniques <- sort(unique(MERIS_SPM_G_AC_sub$lon))
lats_uniques <- sort(unique(MERIS_SPM_G_AC_sub$lat))

# Résolution en degrés
res_lon <- diff(lons_uniques)[1]  # écart entre 2 pixels voisins en longitude
res_lat <- diff(lats_uniques)[1]  # idem en latitude

cat("Résolution lon :", res_lon, "°\n")
cat("Résolution lat :", res_lat, "°\n")

# Conversion en km (pour ~43°N, zone Méditerranée/Atlantique Sud de France)
lat_ref <- 43 
res_lon_km <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res_lat * 111

cat("Résolution lon :", round(res_lon_km, 3), "km\n")
cat("Résolution lat :", round(res_lat_km, 3), "km\n")

# Aire d'un pixel
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 4), "km²\n")

### define 95ème percentile -------------------------------------------------

#### MERIS_SPM_G_AC_sub ------------------------------------------------------

seuil_95_MERIS_SPM_G_AC_sub <- quantile(MERIS_SPM_G_AC_sub$`SPM-G-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MERIS_SPM_G_AC_sub, "g/m³\n")

# seuil = 4.587191 g/m³

# Stats du panache par jour
MERIS_SPM_G_AC_sub_95 <- MERIS_SPM_G_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-AC_mean`>= seuil_95_MERIS_SPM_G_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_MERIS_SPM_G_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_MERIS_SPM_G_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MERIS_SPM_G_PO_sub ------------------------------------------------------

seuil_95_MERIS_SPM_G_PO_sub <- quantile(MERIS_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MERIS_SPM_G_PO_sub, "g/m³\n")

# seuil = 0.5030083 g/m³

# Stats du panache par jour
MERIS_SPM_G_PO_sub_95 <- MERIS_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_95_MERIS_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_MERIS_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_MERIS_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MERIS_SPM_R_AC_sub ------------------------------------------------------

seuil_95_MERIS_SPM_R_AC_sub <- quantile(MERIS_SPM_R_AC_sub$`SPM-R-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MERIS_SPM_R_AC_sub, "g/m³\n")

# seuil = 9.404465 g/m³

# Stats du panache par jour
MERIS_SPM_R_AC_sub_95 <- MERIS_SPM_R_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-AC_mean`>= seuil_95_MERIS_SPM_R_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_MERIS_SPM_R_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_MERIS_SPM_R_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MERIS_SPM_R_PO_sub ------------------------------------------------------

seuil_95_MERIS_SPM_R_PO_sub <- quantile(MERIS_SPM_R_PO_sub$`SPM-R-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MERIS_SPM_R_PO_sub, "g/m³\n")

# seuil = 1.570896 g/m³

# Stats du panache par jour
MERIS_SPM_R_PO_sub_95 <- MERIS_SPM_R_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-PO_mean`>= seuil_95_MERIS_SPM_R_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_MERIS_SPM_R_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_MERIS_SPM_R_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save
# save(MERIS_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_AC_sub_95.Rdata")
# save(MERIS_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_PO_sub_95.Rdata")
# save(MERIS_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_AC_sub_95.Rdata")
# save(MERIS_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_PO_sub_95.Rdata")


## MODIS -------------------------------------------------------------------

### calcul de l'aire --------------------------------------------------------

# Trier et extraire les valeurs uniques
lons_uniques <- sort(unique(MODIS_SPM_G_NS_sub$lon))
lats_uniques <- sort(unique(MODIS_SPM_G_NS_sub$lat))

# Résolution en degrés
res_lon <- diff(lons_uniques)[1]  # écart entre 2 pixels voisins en longitude
res_lat <- diff(lats_uniques)[1]  # idem en latitude

cat("Résolution lon :", res_lon, "°\n")
cat("Résolution lat :", res_lat, "°\n")

# Conversion en km (pour ~43°N, zone Méditerranée/Atlantique Sud de France)
lat_ref <- 43 
res_lon_km <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res_lat * 111

cat("Résolution lon :", round(res_lon_km, 3), "km\n")
cat("Résolution lat :", round(res_lat_km, 3), "km\n")

# Aire d'un pixel
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 4), "km²\n")

### define 95ème percentile -------------------------------------------------

#### MODIS_SPM_G_NS_sub ------------------------------------------------------

seuil_95_MODIS_SPM_G_NS_sub <- quantile(MODIS_SPM_G_NS_sub$`SPM-G-NS_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MODIS_SPM_G_NS_sub, "g/m³\n")

# seuil = 0.6975178 g/m³

# Stats du panache par jour
MODIS_SPM_G_NS_sub_95 <- MODIS_SPM_G_NS_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-NS_mean`>= seuil_95_MODIS_SPM_G_NS_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95_MODIS_SPM_G_NS_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95_MODIS_SPM_G_NS_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MODIS_SPM_G_PO_sub ------------------------------------------------------

seuil_95_MODIS_SPM_G_PO_sub <- quantile(MODIS_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MODIS_SPM_G_PO_sub, "g/m³\n")

# seuil = 1.057702 g/m³

# Stats du panache par jour
MODIS_SPM_G_PO_sub_95 <- MODIS_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_95_MODIS_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_MODIS_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_MODIS_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MODIS_SPM_R_NS_sub ------------------------------------------------------

seuil_95_MODIS_SPM_R_NS_sub <- quantile(MODIS_SPM_R_NS_sub$`SPM-R-NS_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MODIS_SPM_R_NS_sub, "g/m³\n")

# seuil = 2.169223 g/m³

# Stats du panache par jour
MODIS_SPM_R_NS_sub_95 <- MODIS_SPM_R_NS_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-NS_mean`>= seuil_95_MODIS_SPM_R_NS_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-NS_mean`[`SPM-R-NS_mean` >= seuil_95_MODIS_SPM_R_NS_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-NS_mean`[`SPM-R-NS_mean` >= seuil_95_MODIS_SPM_R_NS_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### MODIS_SPM_R_PO_sub ------------------------------------------------------

seuil_95_MODIS_SPM_R_PO_sub <- quantile(MODIS_SPM_R_PO_sub$`SPM-R-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_MODIS_SPM_R_PO_sub, "g/m³\n")

# seuil = 3.685347 g/m³

# Stats du panache par jour
MODIS_SPM_R_PO_sub_95 <- MODIS_SPM_R_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-PO_mean`>= seuil_95_MODIS_SPM_R_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_MODIS_SPM_R_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_MODIS_SPM_R_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save
# save(MODIS_SPM_G_NS_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_NS_sub_95.Rdata")
# save(MODIS_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_PO_sub_95.Rdata")
# save(MODIS_SPM_R_NS_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_NS_sub_95.Rdata")
# save(MODIS_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_PO_sub_95.Rdata")


# OLCI --------------------------------------------------------------------

### calcul de l'aire --------------------------------------------------------

# Trier et extraire les valeurs uniques
lons_uniques <- sort(unique(OLCI_SPM_G_AC_sub$lon))
lats_uniques <- sort(unique(OLCI_SPM_G_AC_sub$lat))

# Résolution en degrés
res_lon <- diff(lons_uniques)[1]  # écart entre 2 pixels voisins en longitude
res_lat <- diff(lats_uniques)[1]  # idem en latitude

cat("Résolution lon :", res_lon, "°\n")
cat("Résolution lat :", res_lat, "°\n")

# Conversion en km (pour ~43°N, zone Méditerranée/Atlantique Sud de France)
lat_ref <- 43 
res_lon_km <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res_lat * 111

cat("Résolution lon :", round(res_lon_km, 3), "km\n")
cat("Résolution lat :", round(res_lat_km, 3), "km\n")

# Aire d'un pixel
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 4), "km²\n")

### define 95ème percentile -------------------------------------------------

#### OLCI_SPM_G_AC_sub ------------------------------------------------------

seuil_95_OLCI_SPM_G_AC_sub <- quantile(OLCI_SPM_G_AC_sub$`SPM-G-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCI_SPM_G_AC_sub, "g/m³\n")

# seuil = 8.744929 g/m³

# Stats du panache par jour
OLCI_SPM_G_AC_sub_95 <- OLCI_SPM_G_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-AC_mean`>= seuil_95_OLCI_SPM_G_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCI_SPM_G_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCI_SPM_G_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_SPM_G_PO_sub ------------------------------------------------------

seuil_95_OLCI_SPM_G_PO_sub <- quantile(OLCI_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCI_SPM_G_PO_sub, "g/m³\n")

# seuil = 0.5430563 g/m³

# Stats du panache par jour
OLCI_SPM_G_PO_sub_95 <- OLCI_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_95_OLCI_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCI_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCI_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_SPM_R_AC_sub ------------------------------------------------------

seuil_95_OLCI_SPM_R_AC_sub <- quantile(OLCI_SPM_R_AC_sub$`SPM-R-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCI_SPM_R_AC_sub, "g/m³\n")

# seuil = 12.32171 g/m³

# Stats du panache par jour
OLCI_SPM_R_AC_sub_95 <- OLCI_SPM_R_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-AC_mean`>= seuil_95_OLCI_SPM_R_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCI_SPM_R_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCI_SPM_R_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_SPM_R_PO_sub ------------------------------------------------------

seuil_95_OLCI_SPM_R_PO_sub <- quantile(OLCI_SPM_R_PO_sub$`SPM-R-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCI_SPM_R_PO_sub, "g/m³\n")

# seuil = 1.741517 g/m³

# Stats du panache par jour
OLCI_SPM_R_PO_sub_95 <- OLCI_SPM_R_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-PO_mean`>= seuil_95_OLCI_SPM_R_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCI_SPM_R_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCI_SPM_R_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save
# save(OLCI_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_AC_sub_95.Rdata")
# save(OLCI_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCI_SPM_G_PO_sub_95.Rdata")
# save(OLCI_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_AC_sub_95.Rdata")
# save(OLCI_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCI_SPM_R_PO_sub_95.Rdata")

# plotting threshold 95 ------------------------------------------------------

# ── 1. Tableau récapitulatif de tous les seuils ────────────────────────────
seuils_df <- data.frame(
  capteur  = c(rep("MERIS", 4), rep("MODIS", 4), rep("OLCI", 4)),
  algo     = c("G","G","R","R",  "G","G","R","R",  "G","G","R","R"),
  correc   = c("AC","PO","AC","PO", "NS","PO","NS","PO", "AC","PO","AC","PO"),
  seuil    = c(
    seuil_95_MERIS_SPM_G_AC_sub,  # 4.587
    seuil_95_MERIS_SPM_G_PO_sub,  # 0.503
    seuil_95_MERIS_SPM_R_AC_sub,  # 9.404
    seuil_95_MERIS_SPM_R_PO_sub,  # 1.571
    seuil_95_MODIS_SPM_G_NS_sub,  # 0.698
    seuil_95_MODIS_SPM_G_PO_sub,  # 1.058
    seuil_95_MODIS_SPM_R_NS_sub,  # 2.169
    seuil_95_MODIS_SPM_R_PO_sub,  # 3.685
    seuil_95_OLCI_SPM_G_AC_sub,   # 8.744929
    seuil_95_OLCI_SPM_G_PO_sub,   # 0.5430563
    seuil_95_OLCI_SPM_R_AC_sub,   # 12.32171
    seuil_95_OLCI_SPM_R_PO_sub    # 1.741517
  )
) |>
  mutate(
    label    = paste0(algo, "-", correc),
    label_val = round(seuil, 2)
  )

# ── 2. Barplot groupé par capteur ─────────────────────────────────────────
couleurs_algo <- c(
  "G-AC" = "#1f77b4", "G-PO" = "#aec7e8",
  "G-NS" = "#6baed6",
  "R-AC" = "#d62728", "R-PO" = "#f7a8a8",
  "R-NS" = "#fc8d59"
)

ggplot(seuils_df, aes(x = capteur, y = seuil, fill = label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "white") +
  geom_text(
    aes(label = label_val),
    position = position_dodge(width = 0.75),
    vjust = -0.4, size = 3, fontface = "bold"
  ) +
  scale_fill_manual(values = couleurs_algo, name = "Algo - Correction") +
  labs(
    title    = "Seuils au 95ème percentile — comparaison capteurs & algorithmes",
    subtitle = "G = Général | R = Régional | AC = Acolite | PO = Polymer | NS = NirSwir",
    x        = NULL,
    y        = "Seuil SPM (g/m³)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey50"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# plotting  ---------------------------------------------------------------

## MERIS_SPM_G_AC_sub ------------------------------------------------------

# en échelle normale

model_MERIS_SPM_G_PO_sub_95 <- lm(mean_spm ~ date, data = MERIS_SPM_G_PO_sub_95)
p_value_MERIS_SPM_G_PO_sub_95 <- summary(model_MERIS_SPM_G_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MERIS_SPM_G_PO_sub_95 <- coef(model_MERIS_SPM_G_PO_sub_95)[1]
slope_MERIS_SPM_G_PO_sub_95 <- coef(model_MERIS_SPM_G_PO_sub_95)[2]

ggplot(data = MERIS_SPM_G_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(MERIS_SPM_G_PO_sub_95$date, na.rm = TRUE),
    y = max(MERIS_SPM_G_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_MERIS_SPM_G_PO_sub_95, 3), " ",
      round(slope_MERIS_SPM_G_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_MERIS_SPM_G_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_MERIS_SPM_G_PO_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — MERIS (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )


## MERIS_SPM_G_PO_sub ------------------------------------------------------
## MERIS_SPM_R_AC_sub ------------------------------------------------------
## MERIS_SPM_R_PO_sub ------------------------------------------------------

## MODIS_SPM_G_NS_sub ------------------------------------------------------

model_MODIS_SPM_G_NS_sub_95 <- lm(mean_spm ~ date, data = MODIS_SPM_G_NS_sub_95)
p_value_MODIS_SPM_G_NS_sub_95 <- summary(model_MODIS_SPM_G_NS_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_SPM_G_NS_sub_95 <- coef(model_MODIS_SPM_G_NS_sub_95)[1]
slope_MODIS_SPM_G_NS_sub_95 <- coef(model_MODIS_SPM_G_NS_sub_95)[2]

ggplot(data = MODIS_SPM_G_NS_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(MODIS_SPM_G_NS_sub_95$date, na.rm = TRUE),
    y = max(MODIS_SPM_G_NS_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_MODIS_SPM_G_NS_sub_95, 3), " ",
      round(slope_MODIS_SPM_G_NS_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_MODIS_SPM_G_NS_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_MODIS_SPM_G_NS_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — MODIS (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : NirSwir | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## MODIS_SPM_G_PO_sub ------------------------------------------------------

model_MODIS_SPM_G_PO_sub_95 <- lm(mean_spm ~ date, data = MODIS_SPM_G_PO_sub_95)
p_value_MODIS_SPM_G_PO_sub_95 <- summary(model_MODIS_SPM_G_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_SPM_G_PO_sub_95 <- coef(model_MODIS_SPM_G_PO_sub_95)[1]
slope_MODIS_SPM_G_PO_sub_95 <- coef(model_MODIS_SPM_G_PO_sub_95)[2]

ggplot(data = MODIS_SPM_G_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(MODIS_SPM_G_PO_sub_95$date, na.rm = TRUE),
    y = max(MODIS_SPM_G_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_MODIS_SPM_G_PO_sub_95, 3), " ",
      round(slope_MODIS_SPM_G_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_MODIS_SPM_G_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_MODIS_SPM_G_PO_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — MODIS (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## MODIS_SPM_R_NS_sub ------------------------------------------------------

model_MODIS_SPM_R_NS_sub_95 <- lm(mean_spm ~ date, data = MODIS_SPM_R_NS_sub_95)
p_value_MODIS_SPM_R_NS_sub_95 <- summary(model_MODIS_SPM_R_NS_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_SPM_R_NS_sub_95 <- coef(model_MODIS_SPM_R_NS_sub_95)[1]
slope_MODIS_SPM_R_NS_sub_95 <- coef(model_MODIS_SPM_R_NS_sub_95)[2]

ggplot(data = MODIS_SPM_R_NS_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(MODIS_SPM_R_NS_sub_95$date, na.rm = TRUE),
    y = max(MODIS_SPM_R_NS_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_MODIS_SPM_R_NS_sub_95, 3), " ",
      round(slope_MODIS_SPM_R_NS_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_MODIS_SPM_R_NS_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_MODIS_SPM_R_NS_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — MODIS (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : NirSwir | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## MODIS_SPM_R_PO_sub ------------------------------------------------------

model_MODIS_SPM_R_PO_sub_95 <- lm(mean_spm ~ date, data = MODIS_SPM_R_PO_sub_95)
p_value_MODIS_SPM_R_PO_sub_95 <- summary(model_MODIS_SPM_R_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_SPM_R_PO_sub_95 <- coef(model_MODIS_SPM_R_PO_sub_95)[1]
slope_MODIS_SPM_R_PO_sub_95 <- coef(model_MODIS_SPM_R_PO_sub_95)[2]

ggplot(data = MODIS_SPM_R_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(MODIS_SPM_R_PO_sub_95$date, na.rm = TRUE),
    y = max(MODIS_SPM_R_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_MODIS_SPM_R_PO_sub_95, 3), " + ",
      round(slope_MODIS_SPM_R_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_MODIS_SPM_R_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_MODIS_SPM_R_PO_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — MODIS (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## OLCI_SPM_G_AC_sub ------------------------------------------------------

model_OLCI_SPM_G_AC_sub_95 <- lm(mean_spm ~ date, data = OLCI_SPM_G_AC_sub_95)
p_value_OLCI_SPM_G_AC_sub_95 <- summary(model_OLCI_SPM_G_AC_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_SPM_G_AC_sub_95 <- coef(model_OLCI_SPM_G_AC_sub_95)[1]
slope_OLCI_SPM_G_AC_sub_95 <- coef(model_OLCI_SPM_G_AC_sub_95)[2]

ggplot(data = OLCI_SPM_G_AC_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCI_SPM_G_AC_sub_95$date, na.rm = TRUE),
    y = max(OLCI_SPM_G_AC_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCI_SPM_G_AC_sub_95, 3), " + ",
      round(slope_OLCI_SPM_G_AC_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCI_SPM_G_AC_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCI_SPM_G_AC_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## OLCI_SPM_G_PO_sub ------------------------------------------------------

model_OLCI_SPM_G_PO_sub_95 <- lm(mean_spm ~ date, data = OLCI_SPM_G_PO_sub_95)
p_value_OLCI_SPM_G_PO_sub_95 <- summary(model_OLCI_SPM_G_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_SPM_G_PO_sub_95 <- coef(model_OLCI_SPM_G_PO_sub_95)[1]
slope_OLCI_SPM_G_PO_sub_95 <- coef(model_OLCI_SPM_G_PO_sub_95)[2]

ggplot(data = OLCI_SPM_G_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCI_SPM_G_PO_sub_95$date, na.rm = TRUE),
    y = max(OLCI_SPM_G_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCI_SPM_G_PO_sub_95, 3), " + ",
      round(slope_OLCI_SPM_G_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCI_SPM_G_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCI_SPM_G_PO_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## OLCI_SPM_R_AC_sub ------------------------------------------------------

model_OLCI_SPM_R_AC_sub_95 <- lm(mean_spm ~ date, data = OLCI_SPM_R_AC_sub_95)
p_value_OLCI_SPM_R_AC_sub_95 <- summary(model_OLCI_SPM_R_AC_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_SPM_R_AC_sub_95 <- coef(model_OLCI_SPM_R_AC_sub_95)[1]
slope_OLCI_SPM_R_AC_sub_95 <- coef(model_OLCI_SPM_R_AC_sub_95)[2]

ggplot(data = OLCI_SPM_R_AC_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCI_SPM_R_AC_sub_95$date, na.rm = TRUE),
    y = max(OLCI_SPM_R_AC_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCI_SPM_R_AC_sub_95, 3), " + ",
      round(slope_OLCI_SPM_R_AC_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCI_SPM_R_AC_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCI_SPM_R_AC_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

## OLCI_SPM_R_PO_sub ------------------------------------------------------

model_OLCI_SPM_R_PO_sub_95 <- lm(mean_spm ~ date, data = OLCI_SPM_R_PO_sub_95)
p_value_OLCI_SPM_R_PO_sub_95 <- summary(model_OLCI_SPM_R_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_SPM_R_PO_sub_95 <- coef(model_OLCI_SPM_R_PO_sub_95)[1]
slope_OLCI_SPM_R_PO_sub_95 <- coef(model_OLCI_SPM_R_PO_sub_95)[2]

ggplot(data = OLCI_SPM_R_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCI_SPM_R_PO_sub_95$date, na.rm = TRUE),
    y = max(OLCI_SPM_R_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCI_SPM_R_PO_sub_95, 3), " + ",
      round(slope_OLCI_SPM_R_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCI_SPM_R_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCI_SPM_R_PO_sub_95, digits = 3))
    ),
    hjust = 1, vjust = 1,
    size = 8,
    color = "red3",
    family = "serif",
    fontface = "italic"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI (ODATIS-MR)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (mg m"^{-3}*")"),
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

# River runoff vs plume area ----------------------------------------------

# We have to create new df to load only one year of data (2019 and 2009)

# 2009 (only Var data)
Y6442010_2009 <- Y6442010_depuis_2000 |> 
  filter(date >= as.Date("2009-01-01"), date <= as.Date("2009-12-31"))

# save(Y6442010_2009, file = "data/Hydro France/Y6442010_2009.Rdata")

# 2019 (Var, Magnan and Paillon data)

All_debit_2019 <- All_debit |> 
  filter(date >= as.Date("2019-01-01"), date <= as.Date("2019-12-31"))

All_debit_2019 <- All_debit_2019[-1,] # on supprime la première ligne car doublon

# save(All_debit_2019, file = "data/Hydro France/All_debit_2019.Rdata")

## MERIS_SPM_G_AC_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MERIS_SPM_G_AC_sub_clean <- na.omit(MERIS_SPM_G_AC_sub_95)

# save(MERIS_SPM_G_AC_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MERIS_SPM_G_AC_sub_clean.Rdata")

MERIS_SPM_G_AC_sub_95_deb <- MERIS_SPM_G_AC_sub_clean |> 
  left_join(Y6442010_2009, by = "date")
  
adjust_factors <- sec_axis_adjustement_factors(MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2, MERIS_SPM_G_AC_sub_95_deb$débit)

MERIS_SPM_G_AC_sub_95_deb$scaled_aire_panache <- MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MERIS_SPM_G_AC_sub_95_deb$débit)
shapiro.test(MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MERIS_SPM_G_AC_sub_95_deb$débit, MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MERIS_SPM_G_AC_sub_95_deb$débit, MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MERIS_SPM_G_AC_sub_95_deb$débit, 
                       MERIS_SPM_G_AC_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MERIS_SPM_G_AC_sub_95_deb,
             aes(x = date, y = débit, color = "Débit"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MERIS_SPM_G_AC_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MERIS_SPM_G_AC_sub_95_deb$date, na.rm = TRUE),
    y = max(MERIS_SPM_G_AC_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var — MERIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10)
  )

## MERIS_SPM_G_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MERIS_SPM_G_PO_sub_clean <- na.omit(MERIS_SPM_G_PO_sub_95)

# save(MERIS_SPM_G_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MERIS_SPM_G_PO_sub_clean.Rdata")

MERIS_SPM_G_PO_sub_95_deb <- MERIS_SPM_G_PO_sub_clean |> 
  left_join(Y6442010_2009, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2, MERIS_SPM_G_PO_sub_95_deb$débit)

MERIS_SPM_G_PO_sub_95_deb$scaled_aire_panache <- MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MERIS_SPM_G_PO_sub_95_deb$débit)
shapiro.test(MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MERIS_SPM_G_PO_sub_95_deb$débit, MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MERIS_SPM_G_PO_sub_95_deb$débit, MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MERIS_SPM_G_PO_sub_95_deb$débit, 
                       MERIS_SPM_G_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MERIS_SPM_G_PO_sub_95_deb,
             aes(x = date, y = débit, color = "Débit"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MERIS_SPM_G_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MERIS_SPM_G_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(MERIS_SPM_G_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var — MERIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## MERIS_SPM_R_AC_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MERIS_SPM_R_AC_sub_clean <- na.omit(MERIS_SPM_R_AC_sub_95)

# save(MERIS_SPM_R_AC_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MERIS_SPM_R_AC_sub_clean.Rdata")

MERIS_SPM_R_AC_sub_95_deb <- MERIS_SPM_R_AC_sub_clean |> 
  left_join(Y6442010_2009, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2, MERIS_SPM_R_AC_sub_95_deb$débit)

MERIS_SPM_R_AC_sub_95_deb$scaled_aire_panache <- MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MERIS_SPM_R_AC_sub_95_deb$débit)
shapiro.test(MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MERIS_SPM_R_AC_sub_95_deb$débit, MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MERIS_SPM_R_AC_sub_95_deb$débit, MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MERIS_SPM_R_AC_sub_95_deb$débit, 
                       MERIS_SPM_R_AC_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MERIS_SPM_R_AC_sub_95_deb,
             aes(x = date, y = débit, color = "Débit"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MERIS_SPM_R_AC_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MERIS_SPM_R_AC_sub_95_deb$date, na.rm = TRUE),
    y = max(MERIS_SPM_R_AC_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var — MERIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## MERIS_SPM_R_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MERIS_SPM_R_PO_sub_clean <- na.omit(MERIS_SPM_R_PO_sub_95)

# save(MERIS_SPM_R_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MERIS_SPM_R_PO_sub_clean.Rdata")

MERIS_SPM_R_PO_sub_95_deb <- MERIS_SPM_R_PO_sub_clean |> 
  left_join(Y6442010_2009, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2, MERIS_SPM_R_PO_sub_95_deb$débit)

MERIS_SPM_R_PO_sub_95_deb$scaled_aire_panache <- MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MERIS_SPM_R_PO_sub_95_deb$débit)
shapiro.test(MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MERIS_SPM_R_PO_sub_95_deb$débit, MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MERIS_SPM_R_PO_sub_95_deb$débit, MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MERIS_SPM_R_PO_sub_95_deb$débit, 
                       MERIS_SPM_R_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MERIS_SPM_R_PO_sub_95_deb,
             aes(x = date, y = débit, color = "Débit"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MERIS_SPM_R_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MERIS_SPM_R_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(MERIS_SPM_R_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit" = "steelblue4", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var — MERIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )




## MODIS_SPM_G_NS_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MODIS_SPM_G_NS_sub_clean <- na.omit(MODIS_SPM_G_NS_sub_95)

# save(MODIS_SPM_G_NS_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MODIS_SPM_G_NS_sub_clean.Rdata")

MODIS_SPM_G_NS_sub_95_deb <- MODIS_SPM_G_NS_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2, MODIS_SPM_G_NS_sub_95_deb$debit_cumule)

MODIS_SPM_G_NS_sub_95_deb$scaled_aire_panache <- MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MODIS_SPM_G_NS_sub_95_deb$debit_cumule)
shapiro.test(MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MODIS_SPM_G_NS_sub_95_deb$debit_cumule, MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MODIS_SPM_G_NS_sub_95_deb$debit_cumule, MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MODIS_SPM_G_NS_sub_95_deb$debit_cumule, 
                       MODIS_SPM_G_NS_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MODIS_SPM_G_NS_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MODIS_SPM_G_NS_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MODIS_SPM_G_NS_sub_95_deb$date, na.rm = TRUE),
    y = max(MODIS_SPM_G_NS_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — MODIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : NirSwir | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## MODIS_SPM_G_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MODIS_SPM_G_PO_sub_clean <- na.omit(MODIS_SPM_G_PO_sub_95)

# save(MODIS_SPM_G_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MODIS_SPM_G_PO_sub_clean.Rdata")

MODIS_SPM_G_PO_sub_95_deb <- MODIS_SPM_G_PO_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2, MODIS_SPM_G_PO_sub_95_deb$debit_cumule)

MODIS_SPM_G_PO_sub_95_deb$scaled_aire_panache <- MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MODIS_SPM_G_PO_sub_95_deb$debit_cumule)
shapiro.test(MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MODIS_SPM_G_PO_sub_95_deb$debit_cumule, MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MODIS_SPM_G_PO_sub_95_deb$debit_cumule, MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MODIS_SPM_G_PO_sub_95_deb$debit_cumule, 
                       MODIS_SPM_G_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MODIS_SPM_G_PO_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MODIS_SPM_G_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MODIS_SPM_G_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(MODIS_SPM_G_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — MODIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : NirSwir | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## MODIS_SPM_R_NS_sub ------------------------------------------------------
# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MODIS_SPM_R_NS_sub_clean <- na.omit(MODIS_SPM_R_NS_sub_95)

# save(MODIS_SPM_R_NS_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MODIS_SPM_R_NS_sub_clean.Rdata")

MODIS_SPM_R_NS_sub_95_deb <- MODIS_SPM_R_NS_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2, MODIS_SPM_R_NS_sub_95_deb$debit_cumule)

MODIS_SPM_R_NS_sub_95_deb$scaled_aire_panache <- MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MODIS_SPM_R_NS_sub_95_deb$debit_cumule)
shapiro.test(MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MODIS_SPM_R_NS_sub_95_deb$debit_cumule, MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MODIS_SPM_R_NS_sub_95_deb$debit_cumule, MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MODIS_SPM_R_NS_sub_95_deb$debit_cumule, 
                       MODIS_SPM_R_NS_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MODIS_SPM_R_NS_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MODIS_SPM_R_NS_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MODIS_SPM_R_NS_sub_95_deb$date, na.rm = TRUE),
    y = max(MODIS_SPM_R_NS_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — MODIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : NirSwir | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## MODIS_SPM_R_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MODIS_SPM_R_PO_sub_clean <- na.omit(MODIS_SPM_R_PO_sub_95)

# save(MODIS_SPM_R_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MODIS_SPM_R_PO_sub_clean.Rdata")

MODIS_SPM_R_PO_sub_95_deb <- MODIS_SPM_R_PO_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2, MODIS_SPM_R_PO_sub_95_deb$debit_cumule)

MODIS_SPM_R_PO_sub_95_deb$scaled_aire_panache <- MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MODIS_SPM_R_PO_sub_95_deb$debit_cumule)
shapiro.test(MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MODIS_SPM_R_PO_sub_95_deb$debit_cumule, MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MODIS_SPM_R_PO_sub_95_deb$debit_cumule, MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MODIS_SPM_R_PO_sub_95_deb$debit_cumule, 
                       MODIS_SPM_R_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MODIS_SPM_R_PO_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MODIS_SPM_R_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MODIS_SPM_R_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(MODIS_SPM_R_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — MODIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )





## OLCI_SPM_G_AC_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
OLCI_SPM_G_AC_sub_clean <- na.omit(OLCI_SPM_G_AC_sub_95)

# save(OLCI_SPM_G_AC_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/OLCI_SPM_G_AC_sub_clean.Rdata")

OLCI_SPM_G_AC_sub_95_deb <- OLCI_SPM_G_AC_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2, OLCI_SPM_G_AC_sub_95_deb$debit_cumule)

OLCI_SPM_G_AC_sub_95_deb$scaled_aire_panache <- OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(OLCI_SPM_G_AC_sub_95_deb$debit_cumule)
shapiro.test(OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(OLCI_SPM_G_AC_sub_95_deb$debit_cumule, OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(OLCI_SPM_G_AC_sub_95_deb$debit_cumule, OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(OLCI_SPM_G_AC_sub_95_deb$debit_cumule, 
                       OLCI_SPM_G_AC_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = OLCI_SPM_G_AC_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.6) +
  geom_point(data = OLCI_SPM_G_AC_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.6) +
  annotate(
    "text",
    x = min(OLCI_SPM_G_AC_sub_95_deb$date, na.rm = TRUE),
    y = max(OLCI_SPM_G_AC_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — OLCI (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## OLCI_SPM_G_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
OLCI_SPM_G_PO_sub_clean <- na.omit(OLCI_SPM_G_PO_sub_95)

# save(OLCI_SPM_G_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/OLCI_SPM_G_PO_sub_clean.Rdata")

OLCI_SPM_G_PO_sub_95_deb <- OLCI_SPM_G_PO_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2, OLCI_SPM_G_PO_sub_95_deb$debit_cumule)

OLCI_SPM_G_PO_sub_95_deb$scaled_aire_panache <- OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(OLCI_SPM_G_PO_sub_95_deb$debit_cumule)
shapiro.test(OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(OLCI_SPM_G_PO_sub_95_deb$debit_cumule, OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(OLCI_SPM_G_PO_sub_95_deb$debit_cumule, OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(OLCI_SPM_G_PO_sub_95_deb$debit_cumule, 
                       OLCI_SPM_G_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = OLCI_SPM_G_PO_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.6) +
  geom_point(data = OLCI_SPM_G_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.6) +
  annotate(
    "text",
    x = min(OLCI_SPM_G_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(OLCI_SPM_G_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — OLCI (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## OLCI_SPM_R_AC_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
OLCI_SPM_R_AC_sub_clean <- na.omit(OLCI_SPM_R_AC_sub_95)

# save(OLCI_SPM_R_AC_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/OLCI_SPM_R_AC_sub_clean.Rdata")

OLCI_SPM_R_AC_sub_95_deb <- OLCI_SPM_R_AC_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2, OLCI_SPM_R_AC_sub_95_deb$debit_cumule)

OLCI_SPM_R_AC_sub_95_deb$scaled_aire_panache <- OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(OLCI_SPM_R_AC_sub_95_deb$debit_cumule)
shapiro.test(OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(OLCI_SPM_R_AC_sub_95_deb$debit_cumule, OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(OLCI_SPM_R_AC_sub_95_deb$debit_cumule, OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(OLCI_SPM_R_AC_sub_95_deb$debit_cumule, 
                       OLCI_SPM_R_AC_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = OLCI_SPM_R_AC_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.6) +
  geom_point(data = OLCI_SPM_R_AC_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.6) +
  annotate(
    "text",
    x = min(OLCI_SPM_R_AC_sub_95_deb$date, na.rm = TRUE),
    y = max(OLCI_SPM_R_AC_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — OLCI (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Acolite | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

## OLCI_SPM_R_PO_sub ------------------------------------------------------

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
OLCI_SPM_R_PO_sub_clean <- na.omit(OLCI_SPM_R_PO_sub_95)

# save(OLCI_SPM_R_PO_sub_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/OLCI_SPM_R_PO_sub_clean.Rdata")

OLCI_SPM_R_PO_sub_95_deb <- OLCI_SPM_R_PO_sub_clean |> 
  inner_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2, OLCI_SPM_R_PO_sub_95_deb$debit_cumule)

OLCI_SPM_R_PO_sub_95_deb$scaled_aire_panache <- OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(OLCI_SPM_R_PO_sub_95_deb$debit_cumule)
shapiro.test(OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(OLCI_SPM_R_PO_sub_95_deb$debit_cumule, OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(OLCI_SPM_R_PO_sub_95_deb$debit_cumule, OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(OLCI_SPM_R_PO_sub_95_deb$debit_cumule, 
                       OLCI_SPM_R_PO_sub_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = OLCI_SPM_R_PO_sub_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.6) +
  geom_point(data = OLCI_SPM_R_PO_sub_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.6) +
  annotate(
    "text",
    x = min(OLCI_SPM_R_PO_sub_95_deb$date, na.rm = TRUE),
    y = max(OLCI_SPM_R_PO_sub_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var, du Paillon et du Magnan — OLCI (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : R | Atmospheric correction : Polymer | Seuil au 95ème percentile"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# ODATIS MR Expert vs Sextant ----------------------------------------------

# we want to see if there is a good correlation between panache area of sextant
# and odatis mr expert data
# We took the best product of each product of expert data : 
# MERIS_SPM_G_PO
# MODIS_SPM_G_NS
# OLCI_SPM_G_PO


# load sextant data -------------------------------------------------------

load("data/SEXTANT/SPM/SEXTANT_2009_spm_95.Rdata")
load("data/SEXTANT/SPM/SEXTANT_2019_spm_95.Rdata")

## MERIS --------------------------------------------------------------------

ODATIS_EXPERT_SEXTANT_panache <- left_join(MERIS_SPM_G_PO_sub_clean, SEXTANT_2009_spm_95, by = "date")

cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, method = "spearman")

# 1. Stocker la corrélation
cor_result <- cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, 
                       ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, 
                       method = "spearman")
rho     <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 2. Graphique
ggplot() +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.x, color = "MERIS")) +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.y, color = "SEXTANT OC5")) +
  annotate(
    "text",
    x = min(ODATIS_EXPERT_SEXTANT_panache$date, na.rm = TRUE),
    y = max(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(values = c("MERIS" = "deeppink", "SEXTANT OC5" = "deepskyblue")) +
  scale_y_continuous(
    name = "Aire des panaches (km²)"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "Aire des panaches selon les produits SEXTANT OC5 et MERIS (2009)",
    x = NULL,
    color = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    axis.title.x       = element_text(size = 14, margin = margin(t = 10)),
    axis.text          = element_text(size = 10, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "bottom",
    legend.text        = element_text(size = 14),
    legend.key.size    = unit(1.2, "cm")
  )

## MODIS --------------------------------------------------------------------

ODATIS_EXPERT_SEXTANT_panache <- left_join(MODIS_SPM_G_NS_sub_clean, SEXTANT_2019_spm_95, by = "date")

cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, method = "spearman")

# 1. Stocker la corrélation
cor_result <- cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, 
                       ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, 
                       method = "spearman")
rho     <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 2. Graphique
ggplot() +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.x, color = "MODIS")) +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.y, color = "SEXTANT OC5")) +
  annotate(
    "text",
    x = min(ODATIS_EXPERT_SEXTANT_panache$date, na.rm = TRUE),
    y = max(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(values = c("MODIS" = "maroon", "SEXTANT OC5" = "deepskyblue")) +
  scale_y_continuous(
    name = "Aire des panaches (km²)"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title   = "Aire des panaches selon les produits SEXTANT OC5 et MODIS (2019)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : NirSwir | Seuil au 95ème percentile,
    SEXTANT OC5 | Seuil au 95ème percentile"
  ) +
  theme_minimal() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "bottom",             # ← bottom au lieu de top
    legend.text       = element_text(size = 14),  # ← plus gros
    legend.key.size   = unit(1.2, "cm")       # ← points de la légende plus gros
  )

## OLCI --------------------------------------------------------------------

ODATIS_EXPERT_SEXTANT_panache <- left_join(OLCI_SPM_G_PO_sub_clean, SEXTANT_2019_spm_95, by = "date")

cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, method = "spearman")

# 1. Stocker la corrélation
cor_result <- cor.test(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, 
                       ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.y, 
                       method = "spearman")
rho     <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 2. Graphique
ggplot() +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.x, color = "OLCI")) +
  geom_line(data = ODATIS_EXPERT_SEXTANT_panache, aes(x = date, y = aire_panache_km2.y, color = "SEXTANT OC5")) +
  annotate(
    "text",
    x = min(ODATIS_EXPERT_SEXTANT_panache$date, na.rm = TRUE),
    y = max(ODATIS_EXPERT_SEXTANT_panache$aire_panache_km2.x, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(values = c("OLCI" = "limegreen", "SEXTANT OC5" = "deepskyblue")) +
  scale_y_continuous(
    name = "Aire des panaches (km²)"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title   = "Aire des panaches selon les produits SEXTANT OC5 et OLCI (2019)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil au 95ème percentile,
    SEXTANT OC5 | Seuil au 95ème percentile"
  ) +
  theme_minimal() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),
    axis.text         = element_text(size = 10, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "bottom",             # ← bottom au lieu de top
    legend.text       = element_text(size = 14),  # ← plus gros
    legend.key.size   = unit(1.2, "cm")       # ← points de la légende plus gros
  )

# 98 percentile -----------------------------------------------------------

#### MERIS_SPM_G_PO_sub ------------------------------------------------------

seuil_98_MERIS_SPM_G_PO_sub <- quantile(MERIS_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.98, na.rm = TRUE)
cat("Seuil 98ème percentile :", seuil_98_MERIS_SPM_G_PO_sub, "g/m³\n")

# seuil = 0.745317 g/m³

# Stats du panache par jour
# MERIS_SPM_G_AC_sub_95 <- MERIS_SPM_G_AC_sub |> 
#   group_by(date) |> 
#   summarise(
#     pixel_count = sum(`SPM-G-AC_mean`>= seuil_95_MERIS_SPM_G_AC_sub, na.rm = TRUE),
#     mean_spm = mean(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_MERIS_SPM_G_AC_sub], na.rm = TRUE),
#     median_spm = median(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_MERIS_SPM_G_AC_sub], na.rm = TRUE),
#     aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
#   )

#### MODIS_SPM_G_NS_sub ------------------------------------------------------

seuil_98_MODIS_SPM_G_NS_sub <- quantile(MODIS_SPM_G_NS_sub$`SPM-G-NS_mean`, 0.98, na.rm = TRUE)
cat("Seuil 98ème percentile :", seuil_98_MODIS_SPM_G_NS_sub, "g/m³\n")

# seuil = 1.388686 g/m³

# Stats du panache par jour
# MODIS_SPM_G_NS_sub_98 <- MODIS_SPM_G_NS_sub |> 
#   group_by(date) |> 
#   summarise(
#     pixel_count = sum(`SPM-G-NS_mean`>= seuil_98_MODIS_SPM_G_NS_sub, na.rm = TRUE),
#     mean_spm = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_98_MODIS_SPM_G_NS_sub], na.rm = TRUE),
#     median_spm = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_98_MODIS_SPM_G_NS_sub], na.rm = TRUE),
#     aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
#   )

#### OLCI_SPM_G_PO_sub ------------------------------------------------------

seuil_98_OLCI_SPM_G_PO_sub <- quantile(OLCI_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.98, na.rm = TRUE)
cat("Seuil 98ème percentile :", seuil_98_OLCI_SPM_G_PO_sub, "g/m³\n")

# seuil = 1.027433 g/m³

# Stats du panache par jour
# OLCI_SPM_G_PO_sub_98 <- OLCI_SPM_G_PO_sub |> 
#   group_by(date) |> 
#   summarise(
#     pixel_count = sum(`SPM-G-PO_mean`>= seuil_98_OLCI_SPM_G_PO_sub, na.rm = TRUE),
#     mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_98_OLCI_SPM_G_PO_sub], na.rm = TRUE),
#     median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_98_OLCI_SPM_G_PO_sub], na.rm = TRUE),
#     aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
#   )

# plotting threshold 98 ------------------------------------------------------

# ── 1. Tableau récapitulatif de tous les seuils ────────────────────────────
seuils_df <- data.frame(
  capteur  = c(rep("MERIS", 1), rep("MODIS", 1), rep("OLCI", 1)),
  algo     = c("G", "G","G"),
  correc   = c("PO", "NS","PO"),
  seuil    = c(
    seuil_98_MERIS_SPM_G_PO_sub,  # 0.745317
    seuil_98_MODIS_SPM_G_NS_sub,  # 1.388686
    seuil_98_OLCI_SPM_G_PO_sub   # 1.027433
  )
) |>
  mutate(
    label    = paste0(algo, "-", correc),
    label_val = round(seuil, 2)
  )

# ── 2. Barplot groupé par capteur ─────────────────────────────────────────
couleurs_algo <- c(
  "G-PO" = "#aec7e8",
  "G-NS" = "#6baed6"
)

ggplot(seuils_df, aes(x = capteur, y = seuil, fill = label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "white") +
  geom_text(
    aes(label = label_val),
    position = position_dodge(width = 0.75),
    vjust = -0.4, size = 3, fontface = "bold"
  ) +
  scale_fill_manual(values = couleurs_algo, name = "Algo - Correction") +
  labs(
    title    = "Seuils au 98ème percentile — comparaison capteurs & algorithmes",
    subtitle = "G = Général | R = Régional | AC = Acolite | PO = Polymer | NS = NirSwir",
    x        = NULL,
    y        = "Seuil SPM (g/m³)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey50"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# Fixed threshold ---------------------------------------------------------

# to not being influenced by threshold we have to set a fix one
# how to choose it ?
# The SEXTANT product seems the best by the correlation between panache area 
# and the liquid flow
# so we can choose firstly the threshold of sextant OC5 product

## MERIS -------------------------------------------------------------------

### MERIS_SPM_G_PO_sub -------------------------------------------------------------------

seuil_sextant_MERIS_SPM_G_PO_sub <- 0.94

# Stats du panache par jour
MERIS_SPM_G_PO_sub_0.94 <- MERIS_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_sextant_MERIS_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_sextant_MERIS_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_sextant_MERIS_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MERIS_SPM_G_PO_sub_0.94_clean <- na.omit(MERIS_SPM_G_PO_sub_0.94)

# save(MERIS_SPM_G_PO_sub_0.94_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MERIS_SPM_G_PO_sub_0.94_clean.Rdata")

MERIS_SPM_G_PO_sub_0.94_deb <- MERIS_SPM_G_PO_sub_0.94_clean |> 
  left_join(Y6442010_2009, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2, MERIS_SPM_G_PO_sub_0.94_deb$débit)

MERIS_SPM_G_PO_sub_0.94_deb$scaled_aire_panache <- MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MERIS_SPM_G_PO_sub_0.94_deb$débit)
shapiro.test(MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MERIS_SPM_G_PO_sub_0.94_deb$débit, MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MERIS_SPM_G_PO_sub_0.94_deb$débit, MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MERIS_SPM_G_PO_sub_0.94_deb$débit, 
                       MERIS_SPM_G_PO_sub_0.94_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MERIS_SPM_G_PO_sub_0.94_deb,
             aes(x = date, y = débit, color = "Débit du Var"),
             size = 2, alpha = 0.6) +
  geom_point(data = MERIS_SPM_G_PO_sub_0.94_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.6) +
  annotate(
    "text",
    x = min(MERIS_SPM_G_PO_sub_0.94_deb$date, na.rm = TRUE),
    y = max(MERIS_SPM_G_PO_sub_0.94_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit du Var" = "blue", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit du Var" = "blue", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit liquide moyen journalier du Var — MERIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil à 0.94 g/m³"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 10, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

## MODIS -------------------------------------------------------------------

### MODIS_SPM_G_NS_sub ------------------------------------------------------

seuil_sextant_MODIS_SPM_G_NS_sub <- 0.94

# Stats du panache par jour
MODIS_SPM_G_NS_sub_0.94 <- MODIS_SPM_G_NS_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-NS_mean`>= seuil_sextant_MODIS_SPM_G_NS_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_sextant_MODIS_SPM_G_NS_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_sextant_MODIS_SPM_G_NS_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
MODIS_SPM_G_NS_sub_0.94_clean <- na.omit(MODIS_SPM_G_NS_sub_0.94)

# save(MODIS_SPM_G_NS_sub_0.94_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/MODIS_SPM_G_NS_sub_0.94_clean.Rdata")

MODIS_SPM_G_NS_sub_0.94_deb <- MODIS_SPM_G_NS_sub_0.94_clean |> 
  left_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2, MODIS_SPM_G_NS_sub_0.94_deb$debit_cumule)

MODIS_SPM_G_NS_sub_0.94_deb$scaled_aire_panache <- MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MODIS_SPM_G_NS_sub_0.94_deb$debit_cumule)
shapiro.test(MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MODIS_SPM_G_NS_sub_0.94_deb$debit_cumule, MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MODIS_SPM_G_NS_sub_0.94_deb$debit_cumule, MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MODIS_SPM_G_NS_sub_0.94_deb$debit_cumule, 
                       MODIS_SPM_G_NS_sub_0.94_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MODIS_SPM_G_NS_sub_0.94_deb,
             aes(x = date, y = débit, color = "Débit cumulé"),
             size = 2, alpha = 0.6) +
  geom_point(data = MODIS_SPM_G_NS_sub_0.94_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.6) +
  annotate(
    "text",
    x = min(MODIS_SPM_G_NS_sub_0.94_deb$date, na.rm = TRUE),
    y = max(MODIS_SPM_G_NS_sub_0.94_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé du Var, du Paillon et du Magnan — MODIS (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : NirSwir | Seuil à 0.94 g/m³"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 10, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

## OLCI -------------------------------------------------------------------

### OLCI_SPM_G_PO_sub ------------------------------------------------------

seuil_sextant_OLCI_SPM_G_PO_sub <- 0.94

# Stats du panache par jour
OLCI_SPM_G_PO_sub_0.94 <- OLCI_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_sextant_OLCI_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_sextant_OLCI_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_sextant_OLCI_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# To compare precisely river runoff and plume area we have to merge data set to
# only keep values in both data set

# first clean data
OLCI_SPM_G_PO_sub_0.94_clean <- na.omit(OLCI_SPM_G_PO_sub_0.94)

# save(OLCI_SPM_G_PO_sub_0.94_clean, file = "data/ODATIS-MR_expert/95 percentile/clean/OLCI_SPM_G_PO_sub_0.94_clean.Rdata")

OLCI_SPM_G_PO_sub_0.94_deb <- OLCI_SPM_G_PO_sub_0.94_clean |> 
  left_join(All_debit_2019, by = "date")

adjust_factors <- sec_axis_adjustement_factors(OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2, OLCI_SPM_G_PO_sub_0.94_deb$debit_cumule)

OLCI_SPM_G_PO_sub_0.94_deb$scaled_aire_panache <- OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(OLCI_SPM_G_PO_sub_0.94_deb$debit_cumule)
shapiro.test(OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(OLCI_SPM_G_PO_sub_0.94_deb$debit_cumule, OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(OLCI_SPM_G_PO_sub_0.94_deb$debit_cumule, OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(OLCI_SPM_G_PO_sub_0.94_deb$debit_cumule, 
                       OLCI_SPM_G_PO_sub_0.94_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = OLCI_SPM_G_PO_sub_0.94_deb,
             aes(x = date, y = débit, color = "Débit cumulé"),
             size = 2, alpha = 0.6) +
  geom_point(data = OLCI_SPM_G_PO_sub_0.94_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.6) +
  annotate(
    "text",
    x = min(OLCI_SPM_G_PO_sub_0.94_deb$date, na.rm = TRUE),
    y = max(OLCI_SPM_G_PO_sub_0.94_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé du Var, du Paillon et du Magnan — OLCI (ODATIS-MR)",
    x       = NULL,
    color   = NULL,
    caption = "Source : ODATIS — MR Expert Product | Algorithm : G | Atmospheric correction : Polymer | Seuil à 0.94 g/m³"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 10, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

# spatial analysis -------------------------------------------------------

# On identifie 3 dates : une date de crue forte, de crue moyenne / petite et 
# un jour sans crue
# Ensuite, à partir des produits bruts on sélectionne seulement ces dates
# Et pour chaque produit on plotte la carte spatiale

## identification des dates ------------------------------------------------

# MERIS (2009)

# forte crue : 09 , 08/02/2009
# moyenne / petite crue : 25/01/2009
# pas de crue : 16/02/2009

# MODIS (2019)

# forte crue : 23/12/2019
# moyenne / petite crue : 27/04/2019
# pas de crue : 20/07/2019

# OLCI (2019)

# forte crue : 23/12/2019
# moyenne / petite crue : 27/04/2019
# pas de crue : 20/07/2019

# filter raw data ---------------------------------------------------------

# MERIS_SPM_G_AC -------------------------------------------------------------------

MERIS_SPM_G_AC_sub_big_flood <- MERIS_SPM_G_AC_sub |> 
  filter(date >= as.Date("2009-02-09"), date <= as.Date("2009-02-09"))

MERIS_SPM_G_AC_sub_small_flood <- MERIS_SPM_G_AC_sub |> 
  filter(date >= as.Date("2009-01-25"), date <= as.Date("2009-01-25"))

MERIS_SPM_G_AC_sub_no_flood <- MERIS_SPM_G_AC_sub |> 
  filter(date >= as.Date("2009-02-16"), date <= as.Date("2009-02-16"))

# MERIS_SPM_G_PO -------------------------------------------------------------------

MERIS_SPM_G_PO_sub_big_flood <- MERIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2009-02-09"), date <= as.Date("2009-02-09"))

MERIS_SPM_G_PO_sub_small_flood <- MERIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2009-01-25"), date <= as.Date("2009-01-25"))

MERIS_SPM_G_PO_sub_no_flood <- MERIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2009-02-16"), date <= as.Date("2009-02-16"))

# MODIS_SPM_G_NS ----------------------------------------------------------

MODIS_SPM_G_NS_sub_big_flood <- MODIS_SPM_G_NS_sub |> 
  filter(date >= as.Date("2019-12-23"), date <= as.Date("2019-12-23"))

MODIS_SPM_G_NS_sub_small_flood <- MODIS_SPM_G_NS_sub |> 
  filter(date >= as.Date("2019-04-27"), date <= as.Date("2019-04-27"))

MODIS_SPM_G_NS_sub_no_flood <- MODIS_SPM_G_NS_sub |> 
  filter(date >= as.Date("2019-08-03"), date <= as.Date("2019-08-03"))

# MODIS_SPM_G_PO ----------------------------------------------------------

MODIS_SPM_G_PO_sub_big_flood <- MODIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-12-23"), date <= as.Date("2019-12-23"))

MODIS_SPM_G_PO_sub_small_flood <- MODIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-04-27"), date <= as.Date("2019-04-27"))

MODIS_SPM_G_PO_sub_no_flood <- MODIS_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-08-03"), date <= as.Date("2019-08-03"))

# OLCI_SPM_G_AC -----------------------------------------------------------

OLCI_SPM_G_AC_sub_big_flood <- OLCI_SPM_G_AC_sub |> 
  filter(date >= as.Date("2019-12-23"), date <= as.Date("2019-12-23"))

OLCI_SPM_G_AC_sub_small_flood <- OLCI_SPM_G_AC_sub |> 
  filter(date >= as.Date("2019-04-27"), date <= as.Date("2019-04-27"))

OLCI_SPM_G_AC_sub_no_flood <- OLCI_SPM_G_AC_sub |> 
  filter(date >= as.Date("2019-08-03"), date <= as.Date("2019-08-03"))

# OLCI_SPM_G_PO -----------------------------------------------------------

OLCI_SPM_G_PO_sub_big_flood <- OLCI_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-12-23"), date <= as.Date("2019-12-23"))

OLCI_SPM_G_PO_sub_small_flood <- OLCI_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-04-27"), date <= as.Date("2019-04-27"))

OLCI_SPM_G_PO_sub_no_flood <- OLCI_SPM_G_PO_sub |> 
  filter(date >= as.Date("2019-08-03"), date <= as.Date("2019-08-03"))

# graph -------------------------------------------------------------------

world_hr <- ne_countries(scale = "medium", returnclass = "sf")

# spatial plotting --------------------------------------------------------

## MERIS_SPM_G_AC -------------------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_AC_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 09 février 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_AC_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
    borders("world", colour = "grey30", linewidth = 0.4) +
    coord_fixed(
      ratio = 1.2,
      xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
      ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
    ) +
  labs(
    title    = "Concentration en MES — 25 janvier 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_AC_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 16 février 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

## MERIS_SPM_G_PO -------------------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_PO_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 09 février 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_PO_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 25 janvier 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MERIS_SPM_G_PO_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MERIS_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MERIS_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 16 février 2009",
    subtitle = "Produit MERIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

## MODIS_SPM_G_NS ----------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_NS_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-NS_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_NS_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_NS_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 23 décembre 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : NirSwir",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_NS_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-NS_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_NS_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_NS_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 27 avril 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : NirSwir",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_NS_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-NS_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_NS_sub_no_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_NS_sub_no_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 03 août 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : NirSwir",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

## MODIS_SPM_G_PO ----------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_PO_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_PO_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_PO_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 23 décembre 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_PO_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_PO_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_PO_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 27 avril 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = MODIS_SPM_G_PO_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(MODIS_SPM_G_PO_sub_no_flood$lon, na.rm = TRUE),
    ylim = range(MODIS_SPM_G_PO_sub_no_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 03 août 2019",
    subtitle = "Produit MODIS | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

## OLCI_SPM_G_AC -----------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_AC_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_AC_sub_big_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_AC_sub_big_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 23 décembre 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_AC_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_AC_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_AC_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 27 avril 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_AC_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-AC_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_AC_sub_no_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_AC_sub_no_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 03 août 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Acolite",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

## OLCI_SPM_G_PO -----------------------------------------------------------

### Big flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_PO_sub_big_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_PO_sub_big_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_PO_sub_big_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 23 décembre 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### Small flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_PO_sub_small_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_PO_sub_small_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_PO_sub_small_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 27 avril 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

### No flood ---------------------------------------------------------------

ggplot() +
  geom_tile(data = OLCI_SPM_G_PO_sub_no_flood, 
            aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  annotation_north_arrow(location = "tr") +
  scale_fill_viridis_c(
    name     = "MES (g/m³)",
    option   = "turbo",
    na.value = "transparent"
  ) +
  borders("world", colour = "grey30", linewidth = 0.4) +
  coord_fixed(
    ratio = 1.2,
    xlim = range(OLCI_SPM_G_PO_sub_no_flood$lon, na.rm = TRUE),
    ylim = range(OLCI_SPM_G_PO_sub_no_flood$lat, na.rm = TRUE)
  ) +
  labs(
    title    = "Concentration en MES — 03 août 2019",
    subtitle = "Produit OLCI | Algorithm : G | Correction atmosphérique : Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)",
    caption  = "Source : ODATIS — MR Expert Product"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9, color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey50", linewidth = 0.5),
    legend.ACsition  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

