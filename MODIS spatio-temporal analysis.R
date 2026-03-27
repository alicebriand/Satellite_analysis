# TS analysis MODIS

# 03/03/2026

# pathway : "~/Documents/Alice/Code/MODIS_TS_analysis


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
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MODIS_2002_2024_spm_95, file = "data/MODIS/SPM/MODIS_2002_2024_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# mean spm or panache area plot

# en échelle normale

# model_MODIS_2002_95 <- lm(aire_panache_km2 ~ date, data = MODIS_2002_2024_spm_95)
model_MODIS_2002_95 <- lm(mean_spm ~ date, data = MODIS_2002_2024_spm_95)
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

ggplot(data = MODIS_2002_2024_spm_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.5) +
  # geom_point(data = MODIS_2002_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(MODIS_2002_2024_spm_95$date, na.rm = TRUE),
    y = max(MODIS_2002_2024_spm_95$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_MODIS_2002_95, 3), " + ", round(slope_MODIS_2002_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_MODIS_2002_95 < 0.001, "< 0.001", format(p_value_MODIS_2002_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration moyenne en MES dans les panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR)",
       x = "Date",
       y = "Concentration moyenne en SPM (en mg/m³)") +
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
model_MODIS_2002_95_log <- lm(log10(aire_panache_km2) ~ date, data = data_log_spm)
p_value_MODIS_2002_95_log <- summary(model_MODIS_2002_95_log)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_2002_95_log <- coef(model_MODIS_2002_95_log)[1]
slope_MODIS_2002_95_log <- coef(model_MODIS_2002_95_log)[2]

# ggplot(data = data_log_spm, aes(x = date, y = mean_spm)) +  # utilise data_log_spm
#   geom_point(color = "red3", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, 
#               formula = y ~ x,
#               color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
#   annotate(
#     "text",
#     x = max(data_log_spm$date, na.rm = TRUE),
#     y = max(data_log_spm$mean_spm, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "log(y) = ", round(intercept_MODIS_2002_95_log, 3), " + ", 
#       round(slope_MODIS_2002_95_log, 7), " * x",
#       "\n", "p = ", ifelse(p_value_MODIS_2002_95_log < 0.001, "< 0.001", 
#                            format(p_value_MODIS_2002_95_log, digits = 3))
#     ),
#     hjust = 1, vjust = 1, size = 6
#   ) +
#   labs(title = "Évolution de la concentration en MES dans les panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR) (échelle log)",
#        x = "Date",
#        y = "Concentration moyenne en MES (en mg/m³)") +
#   theme_minimal() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")

ggplot(data = data_log_spm, aes(x = date, y = aire_panache_km2)) +  # utilise data_log_spm
  geom_point(color = "darkcyan", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, 
              formula = y ~ x,
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$aire_panache_km2, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_MODIS_2002_95_log, 3), " + ", 
      round(slope_MODIS_2002_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_MODIS_2002_95_log < 0.001, "< 0.001", 
                           format(p_value_MODIS_2002_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de la taille des panaches de la baie des Anges vu par le produit MODIS (ODATIS-MR) (échelle log)",
       x = "Date",
       y = "Aire des panaches (en km²)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")





Y6442010_depuis_2002 <- Y6442010_depuis_2000 |>
  filter(date >= as.Date("2002-07-04"))

save(Y6442010_depuis_2002, file = "data/Hydro France/Y6442010_depuis_2002.Rdata")

  