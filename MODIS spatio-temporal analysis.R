# TS analysis MODIS

# 03/03/2026

# pathway : "~/Satellite_analysis/"


# This script will load MODIS satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the 
# Var river at the Napoleon bridge


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(stars)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# problème cluster --------------------------------------------------------


# apparement problème de worker
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

clusterExport(cl, varlist = c("lon_range", "lat_range", "load_MODIS_spm_pixels"))

clusterEvalQ(cl, {
  library(tidyverse)
  library(tidync)
})

# Stopper l'ancien cluster s'il existe
if (exists("cl")) stopCluster(cl)

# Fonction qui gère son propre cluster
load_year <- function(year_dirs, lon_range, lat_range) {
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("lon_range", "lat_range", "load_MODIS_spm_pixels"),
                envir = environment())
  clusterEvalQ(cl, { library(tidyverse); library(tidync) })
  
  result <- plyr::ldply(year_dirs, load_MODIS_spm_pixels,
                        .parallel = TRUE,
                        lon_range = lon_range,
                        lat_range = lat_range)
  stopCluster(cl)
  return(result)
}

# setup -------------------------------------------------------------------

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)

# functions -----------------------------------------------------------------

## scaling function --------------------------------------------------------

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

## df treatment function --------------------------------------------------------

filename <- "~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2015_SPM/L3m_20150101__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc"

nc <- tidync(filename)
print(nc)

load_MODIS_spm_pixels <- function(file_name, lon_range, lat_range){
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  # The necessary code
  MODIS_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    tidync::hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = file_date) |> 
    dplyr::select(lon, lat, date, `SPM-G-NS_mean`)
  
  # Exit
  return(MODIS_one)
}

# load data ---------------------------------------------------------------
## Hydro France data ---------------------------------------------------------------------

load("data/Hydro France/Y6442010_depuis_2002.Rdata")

load("data/MODIS/SPM/MODIS_2002_2024_spm_95.Rdata")

load("data/MODIS/SPM/MODIS_2002_2024_spm_pixels.RData")

MODIS_03_10_2020 <- MODIS_2002_2024_spm_pixels |> 
  filter(date == "2020-10-03")

# loading data ------------------------------------------------------------
## SPM ---------------------------------------------------------------------

MODIS_2002_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2002_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2003_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2003_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2004_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2004_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2005_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2005_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2006_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2006_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2007_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2007_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2008_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2008_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2009_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2009_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2010_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2010_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2011_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2011_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2012_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2012_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2013_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2013_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2014_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2014_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2015_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2015_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2016_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2016_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2017_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2017_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2018_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2018_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2019_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2019_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2020_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2020_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2021_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2021_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2022_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2022_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2023_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2023_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2024_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2024_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### to define threshold with percentile 95 --------------------------------------------------------

# Load and combine
MODIS_2002_spm_pixels <- plyr::ldply(MODIS_2002_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2003_spm_pixels <- plyr::ldply(MODIS_2003_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2004_spm_pixels <- plyr::ldply(MODIS_2004_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2005_spm_pixels <- plyr::ldply(MODIS_2005_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2006_spm_pixels <- plyr::ldply(MODIS_2006_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2007_spm_pixels <- plyr::ldply(MODIS_2007_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2008_spm_pixels <- plyr::ldply(MODIS_2008_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2009_spm_pixels <- plyr::ldply(MODIS_2009_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2010_spm_pixels <- plyr::ldply(MODIS_2010_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2011_spm_pixels <- plyr::ldply(MODIS_2011_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2012_spm_pixels <- plyr::ldply(MODIS_2012_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2013_spm_pixels <- plyr::ldply(MODIS_2013_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2014_spm_pixels <- plyr::ldply(MODIS_2014_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2015_spm_pixels <- plyr::ldply(MODIS_2015_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2016_spm_pixels <- plyr::ldply(MODIS_2016_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2017_spm_pixels <- plyr::ldply(MODIS_2017_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2018_spm_pixels <- plyr::ldply(MODIS_2018_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2019_spm_pixels <- plyr::ldply(MODIS_2019_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2020_spm_pixels <- plyr::ldply(MODIS_2020_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2021_spm_pixels <- plyr::ldply(MODIS_2021_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2022_spm_pixels <- plyr::ldply(MODIS_2022_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2023_spm_pixels <- plyr::ldply(MODIS_2023_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2024_spm_pixels <- plyr::ldply(MODIS_2024_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# tout charger en un seul objet

# Créer une liste nommée de tous les répertoires
all_dirs <- list(
  "2002" = MODIS_2002_spm_dir,
  "2003" = MODIS_2003_spm_dir,
  "2004" = MODIS_2004_spm_dir,
  "2005" = MODIS_2005_spm_dir,
  "2006" = MODIS_2006_spm_dir,
  "2007" = MODIS_2007_spm_dir,
  "2008" = MODIS_2008_spm_dir,
  "2009" = MODIS_2009_spm_dir,
  "2010" = MODIS_2010_spm_dir,
  "2011" = MODIS_2011_spm_dir,
  "2012" = MODIS_2012_spm_dir,
  "2013" = MODIS_2013_spm_dir,
  "2014" = MODIS_2014_spm_dir,
  "2015" = MODIS_2015_spm_dir,
  "2016" = MODIS_2016_spm_dir,
  "2017" = MODIS_2017_spm_dir,
  "2018" = MODIS_2018_spm_dir,
  "2019" = MODIS_2019_spm_dir,
  "2020" = MODIS_2020_spm_dir,
  "2021" = MODIS_2021_spm_dir,
  "2022" = MODIS_2022_spm_dir,
  "2023" = MODIS_2023_spm_dir,
  "2024" = MODIS_2024_spm_dir
)

# Tout charger en un seul objet
MODIS_2002_2024_spm_pixels <- purrr::map_dfr(all_dirs, ~ plyr::ldply(.x, load_MODIS_spm_pixels,
                                                               .parallel = TRUE,
                                                               lon_range = lon_range,
                                                               lat_range = lat_range))

# Appliquer sur toutes les années
MODIS_2002_2024_spm_pixels <- purrr::map_dfr(all_dirs, load_year,
                                             lon_range = lon_range,
                                             lat_range = lat_range)

# Combine and save
MODIS_2002_2024_spm_pixels <- rbind(MODIS_2002_spm_pixels, MODIS_2003_spm_pixels, MODIS_2004_spm_pixels,
                                    MODIS_2005_spm_pixels, MODIS_2006_spm_pixels, MODIS_2007_spm_pixels,
                                    MODIS_2008_spm_pixels, MODIS_2009_spm_pixels, MODIS_2010_spm_pixels,
                                    MODIS_2011_spm_pixels, MODIS_2012_spm_pixels, MODIS_2013_spm_pixels,
                                    MODIS_2014_spm_pixels, MODIS_2015_spm_pixels, MODIS_2016_spm_pixels, 
                                    MODIS_2017_spm_pixels, MODIS_2018_spm_pixels, MODIS_2019_spm_pixels, 
                                    MODIS_2020_spm_pixels, MODIS_2021_spm_pixels, MODIS_2022_spm_pixels, 
                                    MODIS_2023_spm_pixels,MODIS_2024_spm_pixels)

save(MODIS_2002_2024_spm_pixels, file = "data/MODIS/SPM/MODIS_2002_2024_spm_pixels.RData")

load("data/MODIS/SPM/MODIS_2002_2024_spm_pixels.RData")

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

# Ou inspecter les coordonnées lon/lat directement
nc <- tidync("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2015_SPM/L3m_20150101__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc")

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

## define 95ème percentile -------------------------------------------------

# Calculer le 95ème percentile
seuil_95 <- quantile(MODIS_2002_2024_spm_pixels$`SPM-G-NS_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "mg/m³\n")

# Seuil 95ème percentile : 0.576409 mg/m³

# Stats du panache par jour
MODIS_2002_2024_spm_95 <- MODIS_2002_2024_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-NS_mean` >= seuil_95, na.rm = TRUE),
    mean_spm = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95], na.rm = TRUE),
    median_spm = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MODIS_2002_2024_spm_95, file = "data/MODIS/SPM/MODIS_2002_2024_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# mean / median SPM or panache area plot

# en échelle normale

# model_MODIS_2002_95 <- lm(aire_panache_km2 ~ date, data = MODIS_2002_2024_spm_95)
model_MODIS_2002_95 <- lm(median_spm ~ date, data = MODIS_2002_2024_spm_95)
p_value_MODIS_2002_95 <- summary(model_MODIS_2002_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_2002_95 <- coef(model_MODIS_2002_95)[1]
slope_MODIS_2002_95 <- coef(model_MODIS_2002_95)[2]

# ggplot(data = MODIS_2002_2024_spm_95, aes(x = date, y = aire_panache_km2)) +
#   geom_point(color = "darkcyan", size = 0.5) +
#   # geom_point(data = MODIS_2002_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   annotate(
#     "text",
#     x = max(MODIS_2002_2024_spm_95$date, na.rm = TRUE),
#     y = max(MODIS_2002_2024_spm_95$aire_panache_km2, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "y = ", round(intercept_MODIS_2002_95, 3), " + ", round(slope_MODIS_2002_95, 7), " * x",
#       "\n", "p = ", ifelse(p_value_MODIS_2002_95 < 0.001, "< 0.001", format(p_value_MODIS_2002_95, digits = 3))
#     ),
#     hjust = 1,  # Alignement à droite
#     vjust = 1,  # Alignement en haut
#     size = 6
#   ) +
#   labs(title = "Évolution de la taille des panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR)",
#        x = "Date",
#        y = "Aire des panaches (en km²)") +
#   theme_minimal() +
#   scale_x_date(
#     date_breaks = "1 year",  
#     date_labels = "%Y"       
#   )

ggplot(data = MODIS_2002_2024_spm_95, aes(x = date, y = median_spm)) +
  geom_point(color = "red3", size = 0.5) +
  # geom_point(data = MODIS_2002_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(MODIS_2002_2024_spm_95$date, na.rm = TRUE),
    y = max(MODIS_2002_2024_spm_95$median_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_MODIS_2002_95, 3), " ", round(slope_MODIS_2002_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_MODIS_2002_95 < 0.001, "< 0.001", format(p_value_MODIS_2002_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR)",
       x = "Date",
       y = "Concentration médiane en SPM (en mg/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


# en échelle log

data_log_spm <- MODIS_2002_2024_spm_95 |> 
  filter(mean_spm > 0)

# data_log_spm <- MODIS_2002_2024_spm_95 |> 
#   filter(aire_panache_km2 > 0)

# model_MODIS_2002_95_log <- lm(log10(mean_spm) ~ date, data = data_log_spm)
model_MODIS_2002_95_log <- lm(log10(median_spm) ~ date, data = data_log_spm)
p_value_MODIS_2002_95_log <- summary(model_MODIS_2002_95_log)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_2002_95_log <- coef(model_MODIS_2002_95_log)[1]
slope_MODIS_2002_95_log <- coef(model_MODIS_2002_95_log)[2]

ggplot(data = data_log_spm, aes(x = date, y = median_spm)) +  # utilise data_log_spm
  geom_point(color = "red3", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE,
              formula = y ~ x,
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$median_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_MODIS_2002_95_log, 3), " + ",
      round(slope_MODIS_2002_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_MODIS_2002_95_log < 0.001, "< 0.001",
                           format(p_value_MODIS_2002_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR) (échelle log)",
       x = "Date",
       y = "Concentration médiane en MES (en mg/m³)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

# ggplot(data = data_log_spm, aes(x = date, y = aire_panache_km2)) +  # utilise data_log_spm
#   geom_point(color = "darkcyan", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, 
#               formula = y ~ x,
#               color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
#   annotate(
#     "text",
#     x = max(data_log_spm$date, na.rm = TRUE),
#     y = max(data_log_spm$aire_panache_km2, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "log(y) = ", round(intercept_MODIS_2002_95_log, 3), " + ", 
#       round(slope_MODIS_2002_95_log, 7), " * x",
#       "\n", "p = ", ifelse(p_value_MODIS_2002_95_log < 0.001, "< 0.001", 
#                            format(p_value_MODIS_2002_95_log, digits = 3))
#     ),
#     hjust = 1, vjust = 1, size = 6
#   ) +
#   labs(title = "Évolution de la taille des panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR) (échelle log)",
#        x = "Date",
#        y = "Aire des panaches (en km²)") +
#   theme_minimal() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")


# comparison between liquid flow rate and panache extension / mean / median SPM

adjust_factors <- sec_axis_adjustement_factors(MODIS_2002_2024_spm_95$median_spm, Y6442010_depuis_2002$débit)

MODIS_2002_2024_spm_95$scaled_median_spm <- MODIS_2002_2024_spm_95$median_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Y6442010_depuis_2002, 
             aes(x = date, y = débit, color = "Débit"), size = 0.5) +
  geom_point(data = MODIS_2002_2024_spm_95, 
             aes(x = date, y = scaled_median_spm, color = "Concentration médiane en MES"), size = 0.5) +
  scale_color_manual(values = c("Débit" = "blue", "Concentration médiane en MES" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Concentration médiane en MES (en mg/m³)")
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches et du débit du Var vu par le produit MODIS (ODATIS-MR)",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# runoff vs SPM concentration correlation ---------------------------------

Var_MODIS <- Y6442010_depuis_2002 %>%
  select(date, débit) %>%
  inner_join(
    MODIS_2002_2024_spm_95 %>% select(date, mean_spm, aire_panache_km2, median_spm),
    by = "date"
  )

cor.test(Var_MODIS$débit, Var_MODIS$median_spm, method = "spearman")

# scatter plot

# Fusionner les données
Var_MODIS <- Y6442010_depuis_2002 %>% 
  select(date, débit) %>% 
  left_join(
    MODIS_2002_2024_spm_95 %>% select(date, aire_panache_km2, mean_spm),
    by = "date"
  )

ggplot(data = Var_MODIS_clean, aes(x = débit, y = aire_panache_km2)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 1) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x,
    parse = TRUE,
    colour = "red",
    size = 4,
    label.x = 0.05,  # position horizontale (0 = gauche, 1 = droite)
    label.y = 0.95   # position verticale (0 = bas, 1 = haut)
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis", name = "Nombre d'observations") +
  theme_bw() +
    labs(x = "Débit (m³/s)", y = "Aire du panache (en km²)", title = "Débit liquide du Var contre l'aire des panaches vue par MODIS (ODATIS-MR)") +
  theme_minimal()

ggplot(data = Var_MODIS_clean, aes(x = débit, y = mean_spm)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 1) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x,
    parse = TRUE,
    colour = "red",
    size = 4,
    label.x = 0.05,  # position horizontale (0 = gauche, 1 = droite)
    label.y = 0.95   # position verticale (0 = bas, 1 = haut)
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis", name = "Nombre d'observations") +
  theme_bw() +
  labs(x = "Débit (m³/s)", y = "Concentration moyenne en MES (en mg/m³)", title = "Débit liquide du Var contre la concentration en MES dans les panaches vue par MODIS (ODATIS-MR)") +
  theme_minimal()


# clean data --------------------------------------------------------------

# on trouve des R² très bas lorsqu'on compare les débits avec mean_spm et l'extension
# du panache, il faut donc enlever les artefacts et nettoyer nos données

seuil_debit_bas <- quantile(Var_MODIS$débit, 0.10, na.rm = TRUE)
seuil_debit_haut <- quantile(Var_MODIS$débit, 0.90, na.rm = TRUE)
seuil_spm_bas <- quantile(Var_MODIS$mean_spm, 0.10, na.rm = TRUE)
seuil_spm_haut <- quantile(Var_MODIS$mean_spm, 0.90, na.rm = TRUE)
seuil_panache_bas <- quantile(Var_MODIS$aire_panache_km2, 0.10, na.rm = TRUE)
seuil_panache_haut <- quantile(Var_MODIS$aire_panache_km2, 0.90, na.rm = TRUE)


Var_MODIS_clean <- Var_MODIS %>%
  filter(!is.na(débit), !is.na(mean_spm)) %>%
  filter(aire_panache_km2 > 0) %>%   # ← enlever les jours sans panache détecté
  filter(
    !(débit > seuil_debit_haut & aire_panache_km2 < seuil_panache_bas),
    !(débit < seuil_debit_bas & aire_panache_km2 > seuil_panache_haut)
  )

Var_MODIS_filtered <- Var_MODIS %>%
  filter(!is.na(débit), !is.na(mean_spm), aire_panache_km2 > 0)

seuil_debit_bas    <- quantile(Var_MODIS_filtered$débit, 0.10, na.rm = TRUE)
seuil_debit_haut   <- quantile(Var_MODIS_filtered$débit, 0.90, na.rm = TRUE)
seuil_panache_bas  <- quantile(Var_MODIS_filtered$aire_panache_km2, 0.10, na.rm = TRUE)
seuil_panache_haut <- quantile(Var_MODIS_filtered$aire_panache_km2, 0.90, na.rm = TRUE)

Var_MODIS_clean <- Var_MODIS_filtered %>%
  filter(
    !(débit > seuil_debit_haut & aire_panache_km2 < seuil_panache_bas),
    !(débit < seuil_debit_bas & aire_panache_km2 > seuil_panache_haut)
  )

# Vérifier
nrow(Var_MODIS_filtered)  # avant filtrage artefacts
nrow(Var_MODIS_clean)     # après








# Gangloff SPM ------------------------------------------------------------

## zones emboîtées ----------------------------------------------------------

# définir les coordonnées de l'embouchure du Var (comme SEXTANT OC5 chla)
lon_embouchure <- 7.199082
lat_embouchure <- 43.654709

# définir des zones emboîtées de tailles croissantes autour de l'embouchure
rayons_km <- c(5, 10, 20, 40, 70, 100)  # en km

# Convertir en degrés
rayons_deg <- rayons_km / 111

# Calculer le percentile 95 pour chaque zone emboîtée
seuils_zones <- lapply(rayons_deg, function(r) {
  
  pixels_zone <- MODIS_2002_2024_spm_pixels |>
    filter(
      lon >= lon_embouchure - r & lon <= lon_embouchure + r,
      lat >= lat_embouchure - r & lat <= lat_embouchure + r
    )
  
  aire_km2 <- nrow(distinct(pixels_zone, lon, lat)) * aire_pixel_km2
  seuil    <- quantile(pixels_zone$`SPM-G-NS_mean`, 0.95, na.rm = TRUE)
  
  data.frame(rayon_km = r * 111, aire_km2 = aire_km2, seuil_95 = seuil)
})

seuils_zones_df <- bind_rows(seuils_zones)
print(seuils_zones_df)

# Visualiser le plateau
ggplot(seuils_zones_df, aes(x = aire_km2, y = seuil_95)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line() +
  geom_hline(yintercept = seuil_95, linetype = "dashed", color = "red") +
  labs(
    title = "Détermination du seuil de détection du panache turbide",
    subtitle = "Percentile 95 par zone emboîtée autour de l'embouchure",
    x = "Aire de la zone (km²)",
    y = "Percentile 95 des MES (g/m³)"
  ) +
  theme_bw()

# Le seuil retenu est la valeur du plateau (zone > ~5000 km²)
seuil_retenu <- seuils_zones_df |>
  filter(aire_km2 > 1900) |>
  summarise(seuil = mean(seuil_95)) |>
  pull(seuil)

cat("Seuil retenu :", seuil_retenu, "g/m³\n", na.rm = TRUE)
# 0.5951114 g/m³

## ROPP --------------------------------------------------------------------

n_images_total <- n_distinct(MODIS_2002_2024_spm_pixels$date)
# nombre de jour avec des données : 5917

ROPP <- MODIS_2002_2024_spm_pixels |>
  group_by(lon, lat) |>
  summarise(
    freq_above_seuil = sum(`SPM-G-NS_mean` >= seuil_retenu, na.rm = TRUE) / n_images_total,
    .groups = "drop"
  ) |>
  filter(freq_above_seuil >= 0.05)

cat("Nombre de pixels dans la ROPP :", nrow(ROPP), "\n")
# 681

ggplot(ROPP) +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = freq_above_seuil)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Fréquence au-dessus du seuil",
    labels = scales::percent_format(accuracy = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Région d'occurrence des panaches turbides (ROPP)",
    # subtitle = "Pixels où la concentration en MES dépasse le seuil dans au moins 5% des images MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim        = range(ROPP$lon),
    ylim        = range(ROPP$lat),
    expand      = FALSE,
    default_crs = sf::st_crs(4326)
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

## stat panache ------------------------------------------------------------

## Métriques du panache — sans filtre de couverture -----------------------

MODIS_panache_metrics <- MODIS_2002_2024_spm_pixels |>
  semi_join(ROPP, by = c("lon", "lat")) |>
  group_by(date) |>
  summarise(
    pixel_count      = sum(`SPM-G-NS_mean` >= seuil_retenu, na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    mean_spm         = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    max_spm          = max(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_retenu],  na.rm = TRUE),
    median_spm       = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    lat_sud          = min(lat[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    lon_ouest        = min(lon[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    lon_est          = max(lon[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    centroid_lon     = mean(lon[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    centroid_lat     = mean(lat[`SPM-G-NS_mean` >= seuil_retenu], na.rm = TRUE),
    .groups = "drop"
  ) |>
  # Exclure les jours sans aucun pixel de panache détecté
  filter(pixel_count > 0)

cat("Jours avec panache détecté :", nrow(MODIS_panache_metrics), "\n")
# 3221

# plotting ----------------------------------------------------------------

# correlation river flow and plume extension ------------------------------

debit_lags <- Y6442010_depuis_2002 |>
  arrange(date) |>
  select(date, débit) |> 
  mutate(
    debit_j1     = lag(débit, 1),             # j-1 seulement
    debit_2j     = (lag(débit, 1) + lag(débit, 2)) / 2,       # moyenne j-1, j-2
    debit_3j     = (lag(débit, 1) + lag(débit, 2) + lag(débit, 3)) / 3
  )

MODIS_panache_metrics <- MODIS_panache_metrics |>
  inner_join(debit_lags, by = "date") |>
  filter(
    aire_panache_km2 > 0,
    débit > 0,
    is.finite(débit)
  ) |>
  mutate(
    panache_log  = log10(aire_panache_km2),
    debit_j0_log = log10(débit),
    debit_j1_log = log10(debit_j1),
    debit_2j_log = log10(debit_2j),
    debit_3j_log = log10(debit_3j)
  )

# Quel lag est le plus corrélé ?
MODIS_panache_metrics |>
  filter(aire_panache_km2 > 0) |>
  summarise(
    r_j0 = cor(panache_log, debit_j0_log, use = "complete.obs", method = "spearman"),
    r_j1 = cor(panache_log, debit_j1_log, use = "complete.obs", method = "spearman"),
    r_2j = cor(panache_log, debit_2j_log, use = "complete.obs", method = "spearman"),
    r_3j = cor(panache_log, debit_3j_log, use = "complete.obs", method = "spearman")
  )

# c'est le lag au jour deux qui semble avoir la corrélation la plus grande avec le débit

## modèle log - log --------------------------------------------------------

modele_log <- lm(panache_log ~ debit_j1_log, data = MODIS_panache_metrics)
r2       <- summary(modele_log)$r.squared
pente    <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

# Calculer la p-value
p_val <- summary(modele_log)$coefficients[2, 4]

label_eq <- paste0(
  "Aire = ", round(10^ordonnee, 3), " × Q^", round(pente, 2), "\n",
  "R² = ", round(r2, 2), "\n",
  "p = ", formatC(p_val, format = "e", digits = 2)
)

ggplot(MODIS_panache_metrics, aes(x = débit, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text",
           x        = 10^(log10(min(MODIS_panache_metrics$débit, na.rm = TRUE)) + 0.1),
           y        = 10^(log10(max(MODIS_panache_metrics$aire_panache_km2, na.rm = TRUE)) - 0.1),
           label    = label_eq,
           hjust    = 0,     # aligné à gauche
           vjust    = 1,     # aligné en haut
           size     = 8,
           color    = "black",
           fontface = "italic",
           family   = "serif") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = "Aire du panache (km²)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    plot.margin      = margin(1, 1.5, 1, 1, "cm")
  )

## Modèle semi-log : log10(aire) ~ débit --------------------------------------------------------

modele_semilog <- lm(panache_log ~ débit, data = SEXTANT_panache_metrics)
r2       <- summary(modele_semilog)$r.squared
print(r2)
pente    <- coef(modele_semilog)[2]
ordonnee <- coef(modele_semilog)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(SEXTANT_panache_metrics, aes(x = débit, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  # geom_smooth(method = "lm", formula = y ~ x,
  #             color = "black", se = FALSE, linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  annotate("text",
           x = max(SEXTANT_panache_metrics$débit, na.rm = TRUE) * 0.7,
           y = min(SEXTANT_panache_metrics$aire_panache_km2, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = "Aire du panache (km²)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title = element_text(face = "bold", family = "serif"),
    axis.text = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")  # Plus de marge à droite pour l'annotation
  )

# corrélation MES mean et débit -------------------------------------------

modele <- lm(débit ~ mean_spm, data = SEXTANT_panache_metrics)
r2       <- summary(modele)$r.squared
print(r2)
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(SEXTANT_panache_metrics, aes(x = débit, y = mean_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(SEXTANT_panache_metrics$débit, na.rm = TRUE) * 0.7,
           y = min(SEXTANT_panache_metrics$mean_spm, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, vjust = 2, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration en MES (g m"^{-3}*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title = element_text(face = "bold", family = "serif"),
    axis.text = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")  # Plus de marge à droite pour l'annotation
  )

# corrélation MES max et débit -------------------------------------------

modele <- lm(débit ~ max_spm, data = SEXTANT_panache_metrics)
r2       <- summary(modele)$r.squared
print(r2)
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(SEXTANT_panache_metrics, aes(x = débit, y = max_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(SEXTANT_panache_metrics$débit, na.rm = TRUE) * 0.7,
           y = min(SEXTANT_panache_metrics$mean_spm, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration maximale en MES (g m"^{-3}*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title = element_text(face = "bold", family = "serif"),
    axis.text = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")  # Plus de marge à droite pour l'annotation
  )



# cartographie ------------------------------------------------------------

max_spm <- max(MODIS_03_10_2020$`SPM-G-NS_mean`, na.rm = TRUE)

pl_map <- MODIS_03_10_2020 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = `SPM-G-NS_mean`)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),  # ← écriture scientifique
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension — 03 octobre 2020",
    subtitle = "Concentration en MES (g/m³)",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(MODIS_03_10_2020$lon),
    ylim   = range(MODIS_03_10_2020$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Satellite_analysis/Graphiques/MODIS/carto_03_10_2020.png", pl_map, height = 9, width = 14)



