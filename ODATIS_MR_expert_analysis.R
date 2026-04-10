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

# load data ---------------------------------------------------------------

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


## OLCI B -------------------------------------------------------------------

### calcul de l'aire --------------------------------------------------------

# Trier et extraire les valeurs uniques
lons_uniques <- sort(unique(OLCIB_SPM_G_AC_sub$lon))
lats_uniques <- sort(unique(OLCIB_SPM_G_AC_sub$lat))

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

#### OLCIB_SPM_G_AC_sub ------------------------------------------------------

seuil_95_OLCIB_SPM_G_AC_sub <- quantile(OLCIB_SPM_G_AC_sub$`SPM-G-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIB_SPM_G_AC_sub, "g/m³\n")

# seuil = 8.947746 g/m³

# Stats du panache par jour
OLCIB_SPM_G_AC_sub_95 <- OLCIB_SPM_G_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-AC_mean`>= seuil_95_OLCIB_SPM_G_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCIB_SPM_G_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCIB_SPM_G_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_B_SPM_G_PO_sub ------------------------------------------------------

seuil_95_OLCIB_SPM_G_PO_sub <- quantile(OLCIB_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIB_SPM_G_PO_sub, "g/m³\n")

# seuil = 0.5108674 g/m³

# Stats du panache par jour
OLCIB_SPM_G_PO_sub_95 <- OLCIB_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_95_OLCIB_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCIB_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCIB_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_B_SPM_R_AC_sub ------------------------------------------------------

seuil_95_OLCIB_SPM_R_AC_sub <- quantile(OLCIB_SPM_R_AC_sub$`SPM-R-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIB_SPM_R_AC_sub, "g/m³\n")

# seuil = 12.46739 g/m³

# Stats du panache par jour
OLCIB_SPM_R_AC_sub_95 <- OLCIB_SPM_R_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-AC_mean`>= seuil_95_OLCIB_SPM_R_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCIB_SPM_R_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCIB_SPM_R_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCI_B_SPM_R_PO_sub ------------------------------------------------------

seuil_95_OLCIB_SPM_R_PO_sub <- quantile(OLCIB_SPM_R_PO_sub$`SPM-R-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIB_SPM_R_PO_sub, "g/m³\n")

# seuil = 1.680244 g/m³

# Stats du panache par jour
OLCIB_SPM_R_PO_sub_95 <- OLCIB_SPM_R_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-PO_mean`>= seuil_95_OLCIB_SPM_R_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCIB_SPM_R_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCIB_SPM_R_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save
# save(OLCIB_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_G_AC_sub_95.Rdata")
# save(OLCIB_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_G_PO_sub_95.Rdata")
# save(OLCIB_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_R_AC_sub_95.Rdata")
# save(OLCIB_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_R_PO_sub_95.Rdata")

## OLCI A -------------------------------------------------------------------
### calcul de l'aire --------------------------------------------------------

# Trier et extraire les valeurs uniques
lons_uniques <- sort(unique(OLCIA_SPM_G_AC_sub$lon))
lats_uniques <- sort(unique(OLCIA_SPM_G_AC_sub$lat))

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

#### OLCIA_SPM_G_AC_sub ------------------------------------------------------

seuil_95_OLCIA_SPM_G_AC_sub <- quantile(OLCIA_SPM_G_AC_sub$`SPM-G-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIA_SPM_G_AC_sub, "g/m³\n")

# seuil = 8.549291 g/m³

# Stats du panache par jour
OLCIA_SPM_G_AC_sub_95 <- OLCIA_SPM_G_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-AC_mean`>= seuil_95_OLCIA_SPM_G_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCIA_SPM_G_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-AC_mean`[`SPM-G-AC_mean` >= seuil_95_OLCIA_SPM_G_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCIA_SPM_G_PO_sub ------------------------------------------------------

seuil_95_OLCIA_SPM_G_PO_sub <- quantile(OLCIA_SPM_G_PO_sub$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIA_SPM_G_PO_sub, "g/m³\n")

# seuil = 0.5696055 g/m³

# Stats du panache par jour
OLCIA_SPM_G_PO_sub_95 <- OLCIA_SPM_G_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean`>= seuil_95_OLCIA_SPM_G_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCIA_SPM_G_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95_OLCIA_SPM_G_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCIA_SPM_R_AC_sub ------------------------------------------------------

seuil_95_OLCIA_SPM_R_AC_sub <- quantile(OLCIA_SPM_R_AC_sub$`SPM-R-AC_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIA_SPM_R_AC_sub, "g/m³\n")

# seuil = 12.18402 g/m³

# Stats du panache par jour
OLCIA_SPM_R_AC_sub_95 <- OLCIA_SPM_R_AC_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-AC_mean`>= seuil_95_OLCIA_SPM_R_AC_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCIA_SPM_R_AC_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-AC_mean`[`SPM-R-AC_mean` >= seuil_95_OLCIA_SPM_R_AC_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

#### OLCIA_SPM_R_PO_sub ------------------------------------------------------

seuil_95_OLCIA_SPM_R_PO_sub <- quantile(OLCIA_SPM_R_PO_sub$`SPM-R-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95_OLCIA_SPM_R_PO_sub, "g/m³\n")

# seuil = 1.790228 g/m³

# Stats du panache par jour
OLCIA_SPM_R_PO_sub_95 <- OLCIA_SPM_R_PO_sub |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-R-PO_mean`>= seuil_95_OLCIA_SPM_R_PO_sub, na.rm = TRUE),
    mean_spm = mean(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCIA_SPM_R_PO_sub], na.rm = TRUE),
    median_spm = median(`SPM-R-PO_mean`[`SPM-R-PO_mean` >= seuil_95_OLCIA_SPM_R_PO_sub], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save
# save(OLCIA_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/OLCIA_SPM_G_AC_sub_95.Rdata")
# save(OLCIA_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/OLCIA_SPM_G_PO_sub_95.Rdata")
# save(OLCIA_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/OLCIA_SPM_R_AC_sub_95.Rdata")
# save(OLCIA_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/OLCIA_SPM_R_PO_sub_95.Rdata")

# plotting threshold ------------------------------------------------------

# ── 1. Tableau récapitulatif de tous les seuils ────────────────────────────
seuils_df <- data.frame(
  capteur  = c(rep("MERIS", 4), rep("MODIS", 4), rep("OLCI-B", 4), rep("OLCI-A", 4)),
  algo     = c("G","G","R","R",  "G","G","R","R",  "G","G","R","R", "G","G","R","R"),
  correc   = c("AC","PO","AC","PO", "NS","PO","NS","PO", "AC","PO","AC","PO", "AC", "PO", "AC", "PO"),
  seuil    = c(
    seuil_95_MERIS_SPM_G_AC_sub,  # 4.587
    seuil_95_MERIS_SPM_G_PO_sub,  # 0.503
    seuil_95_MERIS_SPM_R_AC_sub,  # 9.404
    seuil_95_MERIS_SPM_R_PO_sub,  # 1.571
    seuil_95_MODIS_SPM_G_NS_sub,  # 0.698
    seuil_95_MODIS_SPM_G_PO_sub,  # 1.058
    seuil_95_MODIS_SPM_R_NS_sub,  # 2.169
    seuil_95_MODIS_SPM_R_PO_sub,  # 3.685
    seuil_95_OLCIB_SPM_G_AC_sub,  # 8.948
    seuil_95_OLCIB_SPM_G_PO_sub,  # 0.511
    seuil_95_OLCIB_SPM_R_AC_sub,  # 12.467
    seuil_95_OLCIB_SPM_R_PO_sub,  # 1.680
    seuil_95_OLCIA_SPM_G_AC_sub,  # 8.549291
    seuil_95_OLCIA_SPM_G_PO_sub,  # 0.5696055
    seuil_95_OLCIA_SPM_R_AC_sub,  # 12.18402
    seuil_95_OLCIA_SPM_R_PO_sub  # 1.790228
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

## OLCIB_SPM_G_AC_sub ------------------------------------------------------

model_OLCIB_SPM_G_AC_sub_95 <- lm(mean_spm ~ date, data = OLCIB_SPM_G_AC_sub_95)
p_value_OLCIB_SPM_G_AC_sub_95 <- summary(model_OLCIB_SPM_G_AC_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCIB_SPM_G_AC_sub_95 <- coef(model_OLCIB_SPM_G_AC_sub_95)[1]
slope_OLCIB_SPM_G_AC_sub_95 <- coef(model_OLCIB_SPM_G_AC_sub_95)[2]

ggplot(data = OLCIB_SPM_G_AC_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCIB_SPM_G_AC_sub_95$date, na.rm = TRUE),
    y = max(OLCIB_SPM_G_AC_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCIB_SPM_G_AC_sub_95, 3), " + ",
      round(slope_OLCIB_SPM_G_AC_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCIB_SPM_G_AC_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCIB_SPM_G_AC_sub_95, digits = 3))
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
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI B (ODATIS-MR)",
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

## OLCI_B_SPM_G_PO_sub ------------------------------------------------------

model_OLCIB_SPM_G_PO_sub_95 <- lm(mean_spm ~ date, data = OLCIB_SPM_G_PO_sub_95)
p_value_OLCIB_SPM_G_PO_sub_95 <- summary(model_OLCIB_SPM_G_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCIB_SPM_G_PO_sub_95 <- coef(model_OLCIB_SPM_G_PO_sub_95)[1]
slope_OLCIB_SPM_G_PO_sub_95 <- coef(model_OLCIB_SPM_G_PO_sub_95)[2]

ggplot(data = OLCIB_SPM_G_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCIB_SPM_G_PO_sub_95$date, na.rm = TRUE),
    y = max(OLCIB_SPM_G_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCIB_SPM_G_PO_sub_95, 3), " ",
      round(slope_OLCIB_SPM_G_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCIB_SPM_G_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCIB_SPM_G_PO_sub_95, digits = 3))
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
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI B (ODATIS-MR)",
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

## OLCI_B_SPM_R_AC_sub ------------------------------------------------------

model_OLCIB_SPM_R_AC_sub_95 <- lm(mean_spm ~ date, data = OLCIB_SPM_R_AC_sub_95)
p_value_OLCIB_SPM_R_AC_sub_95 <- summary(model_OLCIB_SPM_R_AC_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCIB_SPM_R_AC_sub_95 <- coef(model_OLCIB_SPM_R_AC_sub_95)[1]
slope_OLCIB_SPM_R_AC_sub_95 <- coef(model_OLCIB_SPM_R_AC_sub_95)[2]

ggplot(data = OLCIB_SPM_R_AC_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCIB_SPM_R_AC_sub_95$date, na.rm = TRUE),
    y = max(OLCIB_SPM_R_AC_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCIB_SPM_R_AC_sub_95, 3), " + ",
      round(slope_OLCIB_SPM_R_AC_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCIB_SPM_R_AC_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCIB_SPM_R_AC_sub_95, digits = 3))
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
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI B (ODATIS-MR)",
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

## OLCI_B_SPM_R_PO_sub ------------------------------------------------------

model_OLCIB_SPM_R_PO_sub_95 <- lm(mean_spm ~ date, data = OLCIB_SPM_R_PO_sub_95)
p_value_OLCIB_SPM_R_PO_sub_95 <- summary(model_OLCIB_SPM_R_PO_sub_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCIB_SPM_R_PO_sub_95 <- coef(model_OLCIB_SPM_R_PO_sub_95)[1]
slope_OLCIB_SPM_R_PO_sub_95 <- coef(model_OLCIB_SPM_R_PO_sub_95)[2]

ggplot(data = OLCIB_SPM_R_PO_sub_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkslateblue", fill = "red3", alpha = 0.15,
              linewidth = 0.8) +
  annotate(
    "text",
    x = max(OLCIB_SPM_R_PO_sub_95$date, na.rm = TRUE),
    y = max(OLCIB_SPM_R_PO_sub_95$mean_spm, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_OLCIB_SPM_R_PO_sub_95, 3), " + ",
      round(slope_OLCIB_SPM_R_PO_sub_95, 7), " × x",
      "\np = ", ifelse(p_value_OLCIB_SPM_R_PO_sub_95 < 0.001, 
                       "< 0.001", 
                       format(p_value_OLCIB_SPM_R_PO_sub_95, digits = 3))
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
    title   = "Concentration moyenne en MES dans les panaches de la baie des Anges — OLCI B (ODATIS-MR)",
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




