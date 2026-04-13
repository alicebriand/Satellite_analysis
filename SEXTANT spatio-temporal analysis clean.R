# SEXTANT spatio-temporal analysis new
# 27/02/2026

# new script for spatio-temporal analysis of SEXTANT data

# This script will load it in bite sized pieces
# Then perform an time series analysis


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(sf)
library(rnaturalearth)
library(ggpmisc)
library(doParallel); registerDoParallel(cores = 14)

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)

# Functions ---------------------------------------------------------------

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

# to load sextant data
load_SEXTANT_spm_pixels <- function(file_name, lon_range, lat_range){
  
  # Find the date
  sextant_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
  
  # The necessary code
  sextant_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = sextant_one_date) |> 
    dplyr::select(lon, lat, date, analysed_spim)  # tous les pixels, sans filtre
  
  # Exit
  return(sextant_one)
}

# Load data ---------------------------------------------------------------


load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

load("data/Hydro France/Y6442010_depuis_2000.Rdata")

## SPM ---------------------------------------------------------------------

# SEXTANT_1998_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_1999_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1999/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2000_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2000/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2001_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2001/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2002_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2002/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2003_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2003/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2004_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2004/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2005_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2005/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2006_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2006/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2007_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2007/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2008_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2008/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2009_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2009/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2010_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2010/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2011_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2012_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2011_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2012_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2013_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2013/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2014_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2014/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2015_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2015/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2016_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2017_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2018_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2019_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2020_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2021_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2022_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2023_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2024_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2025_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2025/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### to define threshold with percentile 95 --------------------------------------------------------

# Load and combine

# SEXTANT_1998_spm_pixels <- plyr::ldply(SEXTANT_1998_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_1999_spm_pixels <- plyr::ldply(SEXTANT_1999_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2000_spm_pixels <- plyr::ldply(SEXTANT_2000_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2001_spm_pixels <- plyr::ldply(SEXTANT_2001_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2002_spm_pixels <- plyr::ldply(SEXTANT_2002_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2003_spm_pixels <- plyr::ldply(SEXTANT_2003_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2004_spm_pixels <- plyr::ldply(SEXTANT_2004_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2005_spm_pixels <- plyr::ldply(SEXTANT_2005_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2006_spm_pixels <- plyr::ldply(SEXTANT_2006_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2007_spm_pixels <- plyr::ldply(SEXTANT_2007_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2008_spm_pixels <- plyr::ldply(SEXTANT_2008_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2009_spm_pixels <- plyr::ldply(SEXTANT_2009_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2010_spm_pixels <- plyr::ldply(SEXTANT_2010_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2011_spm_pixels <- plyr::ldply(SEXTANT_2011_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2012_spm_pixels <- plyr::ldply(SEXTANT_2012_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2013_spm_pixels <- plyr::ldply(SEXTANT_2013_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2014_spm_pixels <- plyr::ldply(SEXTANT_2014_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2015_spm_pixels <- plyr::ldply(SEXTANT_2015_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2016_spm_pixels <- plyr::ldply(SEXTANT_2016_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2017_spm_pixels <- plyr::ldply(SEXTANT_2017_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2018_spm_pixels <- plyr::ldply(SEXTANT_2018_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2019_spm_pixels <- plyr::ldply(SEXTANT_2019_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2020_spm_pixels <- plyr::ldply(SEXTANT_2020_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2021_spm_pixels <- plyr::ldply(SEXTANT_2021_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2022_spm_pixels <- plyr::ldply(SEXTANT_2022_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2023_spm_pixels <- plyr::ldply(SEXTANT_2023_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2024_spm_pixels <- plyr::ldply(SEXTANT_2024_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2025_spm_pixels <- plyr::ldply(SEXTANT_2025_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# 
# # Combine and save
# SEXTANT_1998_2025_spm_pixels <- rbind(SEXTANT_1998_spm_pixels, SEXTANT_1999_spm_pixels, SEXTANT_2000_spm_pixels, 
#                                       SEXTANT_2001_spm_pixels, SEXTANT_2002_spm_pixels, SEXTANT_2003_spm_pixels,
#                                       SEXTANT_2004_spm_pixels, SEXTANT_2005_spm_pixels, SEXTANT_2006_spm_pixels,
#                                       SEXTANT_2007_spm_pixels, SEXTANT_2008_spm_pixels, SEXTANT_2009_spm_pixels,
#                                       SEXTANT_2010_spm_pixels, SEXTANT_2011_spm_pixels, SEXTANT_2012_spm_pixels,
#                                       SEXTANT_2013_spm_pixels, SEXTANT_2014_spm_pixels, SEXTANT_2015_spm_pixels,
#                                       SEXTANT_2016_spm_pixels, SEXTANT_2017_spm_pixels, SEXTANT_2018_spm_pixels, 
#                                       SEXTANT_2019_spm_pixels, SEXTANT_2020_spm_pixels, SEXTANT_2021_spm_pixels, 
#                                       SEXTANT_2022_spm_pixels, SEXTANT_2023_spm_pixels,SEXTANT_2024_spm_pixels, 
#                                       SEXTANT_2025_spm_pixels)
# 
# save(SEXTANT_1998_2025_spm_pixels, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

# Lire les attributs du fichier pour trouver la résolution
tidync("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")[["attribute"]]

# Ou inspecter les coordonnées lon/lat directement
nc <- tidync("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")

# Vérifier d'abord le type des colonnes
test <- hyper_tibble(nc)
str(test)

coords <- hyper_tibble(nc) |> 
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat)) |> 
  summarise(
    res_lon = abs(mean(diff(sort(unique(lon))))),
    res_lat = abs(mean(diff(sort(unique(lat)))))
  )

print(coords)

tmp <- hyper_tibble(nc) |> 
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat))

res_lon <- diff(sort(unique(tmp$lon)))[1]  # prend juste le premier écart
res_lat <- diff(sort(unique(tmp$lat)))[1]

cat("Résolution lon :", res_lon, "°\n")
cat("Résolution lat :", res_lat, "°\n")

## calcul de l'aire --------------------------------------------------------

# Conversion en km (pour ~43°N, zone Méditerranée/Atlantique Sud de France)
lat_ref <- 43  

res_lon_km <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res_lat * 111

cat("Résolution lon :", round(res_lon_km, 3), "km\n")
cat("Résolution lat :", round(res_lat_km, 3), "km\n")

# Aire d'un pixel
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 4), "km²\n")

# define 95ème percentile -------------------------------------------------

# Calculer le 95ème percentile
seuil_95 <- quantile(SEXTANT_1998_2025_spm_pixels$analysed_spim, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.94 g/m³

# Stats du panache par jour
SEXTANT_1998_2025_spm_95 <- SEXTANT_1998_2025_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(analysed_spim >= seuil_95, na.rm = TRUE),
    mean_spm = mean(analysed_spim[analysed_spim >= seuil_95], na.rm = TRUE),
    median_spm = median(analysed_spim[analysed_spim >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(SEXTANT_1998_2025_spm_95, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# en échelle normale

# mean spm
model_sextant_1998_95 <- lm(mean_spm ~ date, data = SEXTANT_1998_2025_spm_95)
p_value_sextant_1998_95 <- summary(model_sextant_1998_95)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998_95 <- coef(model_sextant_1998_95)[1]
slope_sextant_1998_95 <- coef(model_sextant_1998_95)[2]

ggplot(data = SEXTANT_1998_2025_spm_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.5) +
  # geom_point(data = SEXTANT_1998_2025_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(SEXTANT_1998_2025_spm_95$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_spm_95$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_sextant_1998_95, 3), " ", round(slope_sextant_1998_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_sextant_1998_95 < 0.001, "< 0.001", format(p_value_sextant_1998_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration moyenne en MES dans les panaches de la baie des Anges vu par le produit SEXTANT",
       x = "Date",
       y = "Concentration moyenne en MES (en mg/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# median spm
model_sextant_1998_95 <- lm(median_spm ~ date, data = SEXTANT_1998_2025_spm_95)
p_value_sextant_1998_95 <- summary(model_sextant_1998_95)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998_95 <- coef(model_sextant_1998_95)[1]
slope_sextant_1998_95 <- coef(model_sextant_1998_95)[2]

ggplot(data = SEXTANT_1998_2025_spm_95, aes(x = date, y = median_spm)) +
  geom_point(color = "red3", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(SEXTANT_1998_2025_spm_95$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_spm_95$median_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_sextant_1998_95, 3), " ", round(slope_sextant_1998_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_sextant_1998_95 < 0.001, "< 0.001", format(p_value_sextant_1998_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches de la baie des Anges vu par le produit SEXTANT",
       x = "Date",
       y = "Concentration médiane en MES (en mg/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


# en échelle log

data_log_spm <- SEXTANT_1998_2025_spm_95 |> 
  filter(mean_spm > 0, median_spm > 0)

# mean spm
model_sextant_1998_95_log <- lm(log10(mean_spm) ~ date, data = data_log_spm)
p_value_sextant_1998_95_log <- summary(model_sextant_1998_95_log)$coefficients[2, 4]
intercept_sextant_1998_95_log <- coef(model_sextant_1998_95_log)[1]
slope_sextant_1998_95_log <- coef(model_sextant_1998_95_log)[2]

ggplot(data = data_log_spm, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE,
              formula = y ~ x,
              aes(y = 10^predict(model_sextant_1998_95_log)),  # droite dans l'espace log
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_sextant_1998_95_log, 3), " ",
      round(slope_sextant_1998_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_sextant_1998_95_log < 0.001, "< 0.001",
                           format(p_value_sextant_1998_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de la concentration moyenne en MES dans les panaches de la baie des Anges vu par le produit SEXTANT (échelle log)",
       x = "Date",
       y = "Concentration moyenne en MES (en mg/m³)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")


# comparison between liquid flow rate vs plume extension

adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$aire_panache_km2, Y6442010_depuis_2000$débit)

SEXTANT_1998_2025_spm_95$scaled_aire_panache <- SEXTANT_1998_2025_spm_95$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(data = Y6442010_depuis_2000, 
            aes(x = date, y = débit, color = "Débit"), size = 0.3) +
  geom_line(data = SEXTANT_1998_2025_spm_95, 
            aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"), size = 0.3) +
  scale_color_manual(values = c("Débit" = "blue", "Aire des panaches" = "darkcyan")) +
  scale_y_continuous(
    # limits = c(0, 250),   # ← min et max de l'axe Y
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Aire des panaches (en km²)")
  ) +
  labs(title = "Évolution des panaches et du débit du Var vu par le produit SEXTANT OC5",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# runoff vs SPM concentration correlation ---------------------------------

Var_SEXTANT_panache <- inner_join(Y6442010_depuis_2000, SEXTANT_1998_2025_spm_95, by = "date")

cor.test(Var_SEXTANT_panache$débit, Var_SEXTANT_panache$aire_panache_km2, method = "spearman")

