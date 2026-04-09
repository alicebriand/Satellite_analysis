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
save(MERIS_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_AC_sub_95.Rdata")
save(MERIS_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_G_PO_sub_95.Rdata")
save(MERIS_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_AC_sub_95.Rdata")
save(MERIS_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MERIS_SPM_R_PO_sub_95.Rdata")


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
save(MODIS_SPM_G_NS_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_NS_sub_95.Rdata")
save(MODIS_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_G_PO_sub_95.Rdata")
save(MODIS_SPM_R_NS_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_NS_sub_95.Rdata")
save(MODIS_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/MODIS_SPM_R_PO_sub_95.Rdata")


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
save(OLCIB_SPM_G_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_G_AC_sub_95.Rdata")
save(OLCIB_SPM_G_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_G_PO_sub_95.Rdata")
save(OLCIB_SPM_R_AC_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_R_AC_sub_95.Rdata")
save(OLCIB_SPM_R_PO_sub_95, file = "data/ODATIS-MR_expert/95 percentile/OLCIB_SPM_R_PO_sub_95.Rdata")

# plotting threshold ------------------------------------------------------

# Tableau récapitulatif des seuils
seuils_df <- data.frame(
  dataset  = c("SPM-G-AC", "SPM-G-PO", "SPM-R-AC", "SPM-R-PO"),
  seuil    = c(seuil_95_OLCIB_SPM_G_AC_sub,
               seuil_95_OLCIB_SPM_G_PO_sub,
               seuil_95_OLCIB_SPM_R_AC_sub,
               seuil_95_OLCIB_SPM_R_PO_sub),
  algo     = c("G", "G", "R", "R"),
  correc   = c("AC", "PO", "AC", "PO")
)

ggplot(seuils_df, aes(x = dataset, y = seuil, fill = algo)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = round(seuil, 2)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("G" = "#1f77b4", "R" = "#d62728"),
                    labels = c("G" = "Gordon", "R" = "Roesler")) +
  labs(
    title    = "Seuils au 95ème percentile par dataset",
    subtitle = "AC = correction atmosphérique | PO = port-out",
    x        = NULL,
    y        = "Seuil SPM (g/m³)",
    fill     = "Algorithme"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))



