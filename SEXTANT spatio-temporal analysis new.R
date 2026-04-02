# OLCI spatio-temporal analysis new
# 27/02/2026

# new script for spatio-temporal analysis of OLCI A and B data

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


# Load and prep one day of data
# testers...
# file_name <- "~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
# lon_range <- c(1, 4)

# load_SEXTANT_spm <- function(file_name, lon_range, lat_range){
#   
#   # Find the date
#   sextant_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
#   
#   # The necessary code
#   sextant_one <- tidync(file_name) |> 
#     hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
#                  lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
#     hyper_tibble() |> 
#     mutate(lon = as.numeric(lon),
#            lat = as.numeric(lat),
#            date = sextant_one_date) |> 
#     dplyr::select(lon, lat, date, mask, analysed_spim) |> 
#   filter(analysed_spim >= 1.2) |> 
#   summarise(pixel_count = n(),
#             mean_spm = mean(analysed_spim, na.rm = TRUE), .by = "date")
#   
#   # Exit
#   return(sextant_one)
# }

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

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_spatial.Rdata")

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

load("data/Hydro France/Y6442010_depuis_2000.Rdata")

## SPM ---------------------------------------------------------------------

SEXTANT_1998_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_1999_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1999/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2000_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2000/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2001_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2001/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2002_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2002/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2003_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2003/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2004_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2004/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2005_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2005/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2006_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2006/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2007_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2007/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2008_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2008/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2009_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2009/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2010_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2010/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2011_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2012_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2011_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2012_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2013_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2013/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2014_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2014/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2015_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2015/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2016_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2017_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2018_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2019_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2020_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2021_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2022_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2023_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2024_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2025_dir <- dir("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2025/", pattern = ".nc", recursive = TRUE, full.names = TRUE)


### threshold of 1.2 --------------------------------------------------------

# Load and combine

SEXTANT_1998_spm <- plyr::ldply(SEXTANT_1998_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_1999_spm <- plyr::ldply(SEXTANT_1999_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2000_spm <- plyr::ldply(SEXTANT_2000_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2001_spm <- plyr::ldply(SEXTANT_2001_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2002_spm <- plyr::ldply(SEXTANT_2002_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2003_spm <- plyr::ldply(SEXTANT_2003_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2004_spm <- plyr::ldply(SEXTANT_2004_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2005_spm <- plyr::ldply(SEXTANT_2005_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2006_spm <- plyr::ldply(SEXTANT_2006_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2007_spm <- plyr::ldply(SEXTANT_2007_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2008_spm <- plyr::ldply(SEXTANT_2008_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2009_spm <- plyr::ldply(SEXTANT_2009_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2010_spm <- plyr::ldply(SEXTANT_2010_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2011_spm <- plyr::ldply(SEXTANT_2011_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2012_spm <- plyr::ldply(SEXTANT_2012_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2013_spm <- plyr::ldply(SEXTANT_2013_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2014_spm <- plyr::ldply(SEXTANT_2014_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2015_spm <- plyr::ldply(SEXTANT_2015_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2016_spm <- plyr::ldply(SEXTANT_2016_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2017_spm <- plyr::ldply(SEXTANT_2017_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2018_spm <- plyr::ldply(SEXTANT_2018_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2019_spm <- plyr::ldply(SEXTANT_2019_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2020_spm <- plyr::ldply(SEXTANT_2020_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2021_spm <- plyr::ldply(SEXTANT_2021_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2022_spm <- plyr::ldply(SEXTANT_2022_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2023_spm <- plyr::ldply(SEXTANT_2023_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2024_spm <- plyr::ldply(SEXTANT_2024_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2025_spm <- plyr::ldply(SEXTANT_2025_dir, load_SEXTANT_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)


# Combine and save
SEXTANT_1998_2025_spm_spatial <- rbind(SEXTANT_1998_spm, SEXTANT_1999_spm, SEXTANT_2000_spm, 
                                       SEXTANT_2001_spm, SEXTANT_2002_spm, SEXTANT_2003_spm,
                                       SEXTANT_2004_spm, SEXTANT_2005_spm, SEXTANT_2006_spm,
                                       SEXTANT_2007_spm, SEXTANT_2008_spm, SEXTANT_2009_spm,
                                       SEXTANT_2010_spm, SEXTANT_2011_spm, SEXTANT_2012_spm,
                                       SEXTANT_2013_spm, SEXTANT_2014_spm, SEXTANT_2015_spm,
                                       SEXTANT_2016_spm, SEXTANT_2017_spm, SEXTANT_2018_spm, 
                                       SEXTANT_2019_spm, SEXTANT_2020_spm, SEXTANT_2021_spm, 
                                       SEXTANT_2022_spm, SEXTANT_2023_spm,SEXTANT_2024_spm, 
                                       SEXTANT_2025_spm)

save(SEXTANT_1998_2025_spm_spatial, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_SPM_spatial.RData")

### to define threshold with percentile 95 --------------------------------------------------------

# Load and combine

SEXTANT_1998_spm_pixels <- plyr::ldply(SEXTANT_1998_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_1999_spm_pixels <- plyr::ldply(SEXTANT_1999_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2000_spm_pixels <- plyr::ldply(SEXTANT_2000_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2001_spm_pixels <- plyr::ldply(SEXTANT_2001_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2002_spm_pixels <- plyr::ldply(SEXTANT_2002_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2003_spm_pixels <- plyr::ldply(SEXTANT_2003_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2004_spm_pixels <- plyr::ldply(SEXTANT_2004_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2005_spm_pixels <- plyr::ldply(SEXTANT_2005_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2006_spm_pixels <- plyr::ldply(SEXTANT_2006_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2007_spm_pixels <- plyr::ldply(SEXTANT_2007_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2008_spm_pixels <- plyr::ldply(SEXTANT_2008_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2009_spm_pixels <- plyr::ldply(SEXTANT_2009_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2010_spm_pixels <- plyr::ldply(SEXTANT_2010_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2011_spm_pixels <- plyr::ldply(SEXTANT_2011_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2012_spm_pixels <- plyr::ldply(SEXTANT_2012_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2013_spm_pixels <- plyr::ldply(SEXTANT_2013_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2014_spm_pixels <- plyr::ldply(SEXTANT_2014_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2015_spm_pixels <- plyr::ldply(SEXTANT_2015_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2016_spm_pixels <- plyr::ldply(SEXTANT_2016_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2017_spm_pixels <- plyr::ldply(SEXTANT_2017_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2018_spm_pixels <- plyr::ldply(SEXTANT_2018_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2019_spm_pixels <- plyr::ldply(SEXTANT_2019_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2020_spm_pixels <- plyr::ldply(SEXTANT_2020_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2021_spm_pixels <- plyr::ldply(SEXTANT_2021_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2022_spm_pixels <- plyr::ldply(SEXTANT_2022_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2023_spm_pixels <- plyr::ldply(SEXTANT_2023_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2024_spm_pixels <- plyr::ldply(SEXTANT_2024_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2025_spm_pixels <- plyr::ldply(SEXTANT_2025_dir, load_SEXTANT_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)


# Combine and save
SEXTANT_1998_2025_spm_pixels <- rbind(SEXTANT_1998_spm_pixels, SEXTANT_1999_spm_pixels, SEXTANT_2000_spm_pixels, 
                                       SEXTANT_2001_spm_pixels, SEXTANT_2002_spm_pixels, SEXTANT_2003_spm_pixels,
                                       SEXTANT_2004_spm_pixels, SEXTANT_2005_spm_pixels, SEXTANT_2006_spm_pixels,
                                       SEXTANT_2007_spm_pixels, SEXTANT_2008_spm_pixels, SEXTANT_2009_spm_pixels,
                                       SEXTANT_2010_spm_pixels, SEXTANT_2011_spm_pixels, SEXTANT_2012_spm_pixels,
                                       SEXTANT_2013_spm_pixels, SEXTANT_2014_spm_pixels, SEXTANT_2015_spm_pixels,
                                       SEXTANT_2016_spm_pixels, SEXTANT_2017_spm_pixels, SEXTANT_2018_spm_pixels, 
                                       SEXTANT_2019_spm_pixels, SEXTANT_2020_spm_pixels, SEXTANT_2021_spm_pixels, 
                                       SEXTANT_2022_spm_pixels, SEXTANT_2023_spm_pixels,SEXTANT_2024_spm_pixels, 
                                       SEXTANT_2025_spm_pixels)

save(SEXTANT_1998_2025_spm_pixels, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

# Spatial analysis --------------------------------------------------------

# Look at the monthly average SPM 
sextant_1998_monthly <- all_spm_propre_sextant_1998 |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) |> 
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = c("lon", "lat", "year", "month"))

# Map of the 12 months
ggplot(data = sextant_1998_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  # annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  facet_wrap(~month)
# facet_grid(year~month)

# on veut faire une carte pour 1 jour en particulier
# Filtrer pour le jour voulu
df_jour <- SEXTANT_1998_2025_spm_pixels %>%
  filter(date == as.Date("2020-10-03"))

ggplot(df_jour, aes(x = lon, y = lat, fill = analysed_spim)) +
  geom_raster() +  # pour des données grillées (raster)
  scale_fill_viridis_c(
    name = "SPM (mg/L)",
    option = "turbo",   # turbo, magma, plasma selon tes préférences
    na.value = "transparent"
  ) +
  coord_fixed(ratio = 1.2) +  # ratio lat/lon adapté à la Méditerranée
  labs(
    title = "Concentration en SPM — 03 octobre 2020",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()


cotes <- ne_coastline(scale = "medium", returnclass = "sf")  

# ou ne_countries() pour les frontières terrestres

ggplot() +
  geom_raster(data = df_jour, aes(x = lon, y = lat, fill = analysed_spim)) +
  geom_sf(data = cotes, color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(
    name = "SPM (mg/L)",
    option = "turbo",
    na.value = "transparent"
  ) +
  coord_sf(
    xlim = c(min(df_jour$lon), max(df_jour$lon)),
    ylim = c(min(df_jour$lat), max(df_jour$lat))
  ) +
  labs(
    title = "Concentration en SPM — 03 octobre 2020",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

# maintenant, on regarde la carte après avoir appliqué les valeurs de filtre au
# 95ème percentile

df_jour <- SEXTANT_1998_2025_spm_pixels %>%
  filter(date == as.Date("2020-10-03")) %>%
  mutate(analysed_spim = ifelse(analysed_spim > 0.94, analysed_spim, NA))

# ensuite les valeurs au dessus de ce seuil sont moyennées

mean_spm_jour <- SEXTANT_1998_2025_spm_pixels %>%
  filter(date == as.Date("2020-10-03"), analysed_spim > 0.94) %>%
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE))



df_jour <- SEXTANT_1998_2025_spm_pixels %>%
  filter(date == as.Date("2020-10-03"), analysed_spim > 0.94)

# Valeur moyenne unique
mean_val <- mean(df_jour$analysed_spim, na.rm = TRUE)

# Créer un polygone convexe autour des pixels du panache
panache_sf <- df_jour %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_union() %>%
  st_convex_hull()  # enveloppe convexe du panache

# Plot
ggplot() +
  geom_raster(data = df_jour %>% mutate(mean_spm = mean_val),
              aes(x = lon, y = lat, fill = mean_spm)) +
  geom_sf(data = panache_sf, fill = NA, color = "red", linewidth = 0.5) +  # contour
  geom_sf(data = cotes, color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(name = "SPM moyen (mg/L)", option = "turbo") +
  coord_sf(
    xlim = c(min(df_jour$lon), max(df_jour$lon)),
    ylim = c(min(df_jour$lat), max(df_jour$lat))
  ) +
  labs(
    title = paste0("SPM moyen dans le panache — 03 octobre 2020"),
    subtitle = paste0("Moyenne = ", round(mean_val, 3), " mg/L | n = ", nrow(df_jour), " pixels"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

# Temporal analysis -------------------------------------------------------

# Time series analysis
sextant_1998_ts <- sextant_1998 |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) |> 
  summarise(mean_spm = mean(analysed_spim, na.rm = TRUE), .by = c("year", "month", "date"))

# Time series plot of the monthly average SPM
ggplot(data = sextant_1998_ts, aes(x = date, y = mean_spm)) +
  geom_line() +
  geom_point() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = "Date", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()

# Boxplot of the SPM per month
ggplot(data = sextant_1998_ts, aes(x = factor(month), y = mean_spm)) +
  geom_boxplot() +
  labs(x = "Month", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()

# Boxplots per year
ggplot(data = sextant_1998_ts, aes(x = factor(year), y = mean_spm)) +
  geom_boxplot() +
  labs(x = "Year", y = "Mean SPM", title = "Monthly Average SPM in 1998") +
  theme_minimal()


# Animation ---------------------------------------------------------------

# Once a brick of data is loaded, it is possible to walk through one day at a time

# Animate using gganimate to show the months as the time steps
p_map <- ggplot(data = sextant_1998_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  # facet_wrap(~month) +
  labs(title = "SEXTANT SPM data for 1998", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  scale_fill_viridis_c(option = "D")

# Add animation
animated_plot <- p_map +
  transition_states(
    month,
    transition_length = 1,
    state_length = 1
  ) +
  enter_fade() +
  exit_fade()

# Render the animation
animate(animated_plot, fps = 10, duration = 10, 
        renderer = gifski_renderer(file = "animations/sextant_1998.gif",
                                   width = 1200,    # Increase width in pixels
                                   height = 1000))   # Increase height in pixels



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

# Ajouter l'aire du panache dans ton df
SEXTANT_1998_2025_spm_pixels <- SEXTANT_1998_2025_spm_spatial |> 
  mutate(aire_panache_km2 = pixel_count * aire_pixel_km2)

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

# plot de la médiane

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



# comparison between liquid flow rate and panache extension

adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$median_spm, Y6442010_depuis_2000$débit)

SEXTANT_1998_2025_spm_95$scaled_median_spm <- SEXTANT_1998_2025_spm_95$median_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Y6442010_depuis_2000, 
             aes(x = date, y = débit, color = "Débit"), size = 0.5) +
  geom_point(data = SEXTANT_1998_2025_spm_95, 
             aes(x = date, y = scaled_median_spm, color = "MES"), size = 0.5) +
  scale_color_manual(values = c("Débit" = "blue", "MES" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Concentration médiane en MES (en mg/m3)")
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches et du débit du Var vu par le produit SEXTANT",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# runoff vs SPM concentration correlation ---------------------------------

Var_SEXTANT_panache <- inner_join(Y6442010_depuis_2000, SEXTANT_1998_2025_spm_95, by = "date")

cor.test(Var_SEXTANT_panache$débit, Var_SEXTANT_panache$median_spm, method = "spearman")


# échelle log

# Définir les limites log de chaque axe
log_debit_min <- min(log10(Y6442010_depuis_2000$débit), na.rm = TRUE)
log_debit_max <- max(log10(Y6442010_depuis_2000$débit), na.rm = TRUE)
log_mean_min  <- min(log10(SEXTANT_1998_2025_spm_95$mean_spm[SEXTANT_1998_2025_spm_95$mean_spm > 0]), na.rm = TRUE)
log_mean_max  <- max(log10(SEXTANT_1998_2025_spm_95$mean_spm), na.rm = TRUE)

# Fonction pour projeter l'aire sur l'échelle du débit (en log)
mean_to_debit_scale <- function(x) {
  (log10(x) - log_mean_min) / (log_mean_max - log_mean_min) *
    (log_debit_max - log_debit_min) + log_debit_min
}

SEXTANT_1998_2025_spm_95 <- SEXTANT_1998_2025_spm_95 %>%
  filter(mean_spm > 0) %>%
  mutate(mean_scaled = 10^mean_to_debit_scale(mean_spm))

ggplot() +
  geom_point(data = Y6442010_depuis_2000,
             aes(x = date, y = débit, color = "Débit"), size = 0.5) +
  geom_point(data = SEXTANT_1998_2025_spm_95,
             aes(x = date, y = mean_scaled, color = "Concentration en MES"), size = 0.5) +
  scale_color_manual(values = c("Débit" = "blue", "Concentration en MES" = "red3")) +
  scale_y_log10(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ 10^((log10(.) - log_debit_min) / (log_debit_max - log_debit_min) *
              (log_mean_max - log_mean_min) + log_mean_min),
      name = "Concentration en moyenne en MES (en mg/m3)"
    )
  ) +
  labs(
    title = "Évolution de la concentration moyenne en MES dans les panaches et du débit du Var vu par SEXTANT (échelle log)",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

# scatter plot

# Fusionner les données
Var_sextant <- Y6442010_depuis_2000 %>% 
  select(date, débit) %>% 
  inner_join(
    SEXTANT_1998_2025_spm_95 %>% select(date, aire_panache_km2, mean_spm),
    by = "date"
  )

ggplot(data = Var_sextant, aes(x = débit, y = aire_panache_km2)) +
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
  labs(x = "Débit (m³/s)", y = "Aire du panache (en km²)", title = "Débit liquide du Var contre l'aire des panaches vue par SEXTANT") +
  theme_minimal()





# différence pixels count SEXTANT vs OLCI ---------------------------------

# Fusionner par date
df_comparison <- inner_join(
  SEXTANT_1998_2025_spm_spatial |> dplyr::select(date, pixel_count, aire_panache_km2) |> 
    rename(pixel_count_sextant = pixel_count, aire_sextant = aire_panache_km2),
  OLCI_2016_2024_spm_spatial |> dplyr::select(date, pixel_count, aire_panache_km2) |> 
    rename(pixel_count_olci = pixel_count, aire_olci = aire_panache_km2),
  by = "date"  # seulement les dates communes aux deux produits
)

df_comparison |> 
  summarise(
    mean_sextant = mean(pixel_count_sextant, na.rm = TRUE),
    mean_olci = mean(pixel_count_olci, na.rm = TRUE),
    median_sextant = median(pixel_count_sextant, na.rm = TRUE),
    median_olci = median(pixel_count_olci, na.rm = TRUE),
    cor = cor(pixel_count_sextant, pixel_count_olci, use = "complete.obs")
  )

df_comparison |> 
  summarise(
    cor_pearson  = cor(pixel_count_sextant, pixel_count_olci, 
                       use = "complete.obs", method = "pearson"),
    cor_spearman = cor(pixel_count_sextant, pixel_count_olci, 
                       use = "complete.obs", method = "spearman")
  )

# Scatter plot de comparaison
ggplot(df_comparison, aes(x = pixel_count_sextant, y = pixel_count_olci)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", color = "red3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  stat_cor(method = "spearman") +
  labs(x = "Pixel count SEXTANT", y = "Pixel count OLCI",
       title = "Comparaison du nombre de pixels panache entre SEXTANT et OLCI") +
  theme_minimal()

# Séries temporelles superposées
ggplot(df_comparison, aes(x = date)) +
  geom_line(aes(y = pixel_count_sextant, color = "SEXTANT"), alpha = 0.7) +
  geom_line(aes(y = pixel_count_olci, color = "OLCI"), alpha = 0.7) +
  scale_color_manual(values = c("SEXTANT" = "steelblue", "OLCI" = "orange")) +
  labs(x = "Date", y = "Nombre de pixels panache", color = "Produit",
       title = "Comparaison temporelle SEXTANT vs OLCI") +
  theme_minimal()

# Vérifier la normalité visuellement
ggplot(df_comparison) +
  geom_histogram(aes(x = pixel_count_sextant), fill = "steelblue", alpha = 0.6, bins = 50) +
  geom_histogram(aes(x = pixel_count_olci), fill = "orange", alpha = 0.6, bins = 50) +
  theme_minimal()

# distribution asymétrique donc on peut utiliser spearman







# clean data --------------------------------------------------------------

# on trouve des R² très bas lorsqu'on compare les débits avec mean_spm et l'extension
# du panache, il faut donc enlever les artefacts et nettoyer nos données

seuil_debit_bas <- quantile(Var_sextant$débit, 0.10, na.rm = TRUE)
seuil_debit_haut <- quantile(Var_sextant$débit, 0.90, na.rm = TRUE)
seuil_spm_bas <- quantile(Var_sextant$mean_spm, 0.10, na.rm = TRUE)
seuil_spm_haut <- quantile(Var_sextant$mean_spm, 0.90, na.rm = TRUE)
seuil_panache_bas <- quantile(Var_sextant$aire_panache_km2, 0.10, na.rm = TRUE)
seuil_panache_haut <- quantile(Var_sextant$aire_panache_km2, 0.90, na.rm = TRUE)


# Var_sextant_clean <- Var_sextant %>%
#   filter(!is.na(débit), !is.na(aire_panache_km2)) %>%
#   filter(
#     # Enlever : débit fort MAIS spm faible (pas de panache malgré crue)
#     !(débit > seuil_debit_haut & aire_panache_km2 < seuil_panache_bas),
#     # Enlever : débit faible MAIS spm fort (panache sans crue)
#     !(débit < seuil_debit_bas & aire_panache_km2 > seuil_panache_haut)
#   )

Var_sextant_clean <- Var_sextant %>%
  filter(!is.na(débit), !is.na(mean_spm)) %>%
  filter(aire_panache_km2 > 0) %>%   # ← enlever les jours sans panache détecté
  filter(
    !(débit > seuil_debit_haut & aire_panache_km2 < seuil_panache_bas),
    !(débit < seuil_debit_bas & aire_panache_km2 > seuil_panache_haut)
  )

Var_sextant_filtered <- Var_sextant %>%
  filter(!is.na(débit), !is.na(mean_spm), aire_panache_km2 > 0)

seuil_debit_bas    <- quantile(Var_sextant_filtered$débit, 0.10, na.rm = TRUE)
seuil_debit_haut   <- quantile(Var_sextant_filtered$débit, 0.90, na.rm = TRUE)
seuil_panache_bas  <- quantile(Var_sextant_filtered$aire_panache_km2, 0.10, na.rm = TRUE)
seuil_panache_haut <- quantile(Var_sextant_filtered$aire_panache_km2, 0.90, na.rm = TRUE)

Var_sextant_clean <- Var_sextant_filtered %>%
  filter(
    !(débit > seuil_debit_haut & aire_panache_km2 < seuil_panache_bas),
    !(débit < seuil_debit_bas & aire_panache_km2 > seuil_panache_haut)
  )

# Vérifier
nrow(Var_sextant_filtered)  # avant filtrage artefacts
nrow(Var_sextant_clean)     # après

