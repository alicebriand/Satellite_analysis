# SEXTANT spatio-temporal analysis clean
# 27/02/2026

# new script for spatio-temporal analysis of SEXTANT data

# This script will load it in bite sized pieces
# Then perform an time series analysis


# Setup ------------------------------------------------------------------

# Install rnaturalearthhires as necessry
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(sf)
library(rnaturalearth)
library(giscoR) # Hi-res coastlines
library(ggpmisc)
library(doParallel); registerDoParallel(cores = parallel::detectCores()-2)
library(ggpubr)  # Pour stat_cor()
library(scales)


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

# to load sextant data
load_SEXTANT_chl_pixels <- function(file_name, lon_range, lat_range){
  
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
    dplyr::select(lon, lat, date, analysed_chl_a)  # tous les pixels, sans filtre
  
  # Exit
  return(sextant_one)
}

# Load data ---------------------------------------------------------------

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

load("data/Hydro France/Y6442010_depuis_2000.Rdata")

## direction ---------------------------------------------------------------

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

SEXTANT_1998_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/1998/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_1999_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/1999/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2000_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2000/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2001_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2001/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2002_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2002/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2003_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2003/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2004_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2004/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2005_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2005/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2006_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2006/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2007_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2007/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2008_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2008/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2009_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2009/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2010_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2010/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2011_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2012_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2011_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2012_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2013_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2013/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2014_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2014/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2015_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2015/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2016_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2017_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2018_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2019_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2020_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2021_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2022_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2023_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2024_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
SEXTANT_2025_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2025/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### to define threshold with percentile 95 --------------------------------------------------------

#### SPM ---------------------------------------------------------------------

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

#### CHL ---------------------------------------------------------------------

SEXTANT_1998_chl_pixels <- plyr::ldply(SEXTANT_1998_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_1999_chl_pixels <- plyr::ldply(SEXTANT_1999_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2000_chl_pixels <- plyr::ldply(SEXTANT_2000_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2001_chl_pixels <- plyr::ldply(SEXTANT_2001_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2002_chl_pixels <- plyr::ldply(SEXTANT_2002_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2003_chl_pixels <- plyr::ldply(SEXTANT_2003_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2004_chl_pixels <- plyr::ldply(SEXTANT_2004_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2005_chl_pixels <- plyr::ldply(SEXTANT_2005_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2006_chl_pixels <- plyr::ldply(SEXTANT_2006_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2007_chl_pixels <- plyr::ldply(SEXTANT_2007_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2008_chl_pixels <- plyr::ldply(SEXTANT_2008_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2009_chl_pixels <- plyr::ldply(SEXTANT_2009_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2010_chl_pixels <- plyr::ldply(SEXTANT_2010_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2011_chl_pixels <- plyr::ldply(SEXTANT_2011_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2012_chl_pixels <- plyr::ldply(SEXTANT_2012_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2013_chl_pixels <- plyr::ldply(SEXTANT_2013_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2014_chl_pixels <- plyr::ldply(SEXTANT_2014_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2015_chl_pixels <- plyr::ldply(SEXTANT_2015_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2016_chl_pixels <- plyr::ldply(SEXTANT_2016_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2017_chl_pixels <- plyr::ldply(SEXTANT_2017_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2018_chl_pixels <- plyr::ldply(SEXTANT_2018_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2019_chl_pixels <- plyr::ldply(SEXTANT_2019_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2020_chl_pixels <- plyr::ldply(SEXTANT_2020_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2021_chl_pixels <- plyr::ldply(SEXTANT_2021_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2022_chl_pixels <- plyr::ldply(SEXTANT_2022_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2023_chl_pixels <- plyr::ldply(SEXTANT_2023_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2024_chl_pixels <- plyr::ldply(SEXTANT_2024_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
SEXTANT_2025_chl_pixels <- plyr::ldply(SEXTANT_2025_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Combine and save
SEXTANT_1998_2025_chl_pixels <- rbind(SEXTANT_1998_chl_pixels, SEXTANT_1999_chl_pixels, SEXTANT_2000_chl_pixels,
                                      SEXTANT_2001_chl_pixels, SEXTANT_2002_chl_pixels, SEXTANT_2003_chl_pixels,
                                      SEXTANT_2004_chl_pixels, SEXTANT_2005_chl_pixels, SEXTANT_2006_chl_pixels,
                                      SEXTANT_2007_chl_pixels, SEXTANT_2008_chl_pixels, SEXTANT_2009_chl_pixels,
                                      SEXTANT_2010_chl_pixels, SEXTANT_2011_chl_pixels, SEXTANT_2012_chl_pixels,
                                      SEXTANT_2013_chl_pixels, SEXTANT_2014_chl_pixels, SEXTANT_2015_chl_pixels,
                                      SEXTANT_2016_chl_pixels, SEXTANT_2017_chl_pixels, SEXTANT_2018_chl_pixels,
                                      SEXTANT_2019_chl_pixels, SEXTANT_2020_chl_pixels, SEXTANT_2021_chl_pixels,
                                      SEXTANT_2022_chl_pixels, SEXTANT_2023_chl_pixels,SEXTANT_2024_chl_pixels,
                                      SEXTANT_2025_chl_pixels)

save(SEXTANT_1998_2025_chl_pixels, file = "data/SEXTANT/CHL/SEXTANT_1998_2025_chl_pixels.RData")

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

# SPM ---------------------------------------------------------------------

## define 95ème percentile -------------------------------------------------

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

## plotting ----------------------------------------------------------------

# en échelle normale

# Calcul du modèle linéaire
model_sextant_1998_95 <- lm(mean_spm ~ date, data = SEXTANT_1998_2025_spm_95)
p_value_sextant_1998_95 <- summary(model_sextant_1998_95)$coefficients[2, 4]
intercept_sextant_1998_95 <- coef(model_sextant_1998_95)[1]
slope_sextant_1998_95 <- coef(model_sextant_1998_95)[2]

# Création du graphique
ggplot(data = SEXTANT_1998_2025_spm_95, aes(x = date, y = mean_spm)) +
  # Points avec style épuré
  geom_point(
    size = 2,
    shape = 21,
    fill = "red3",
    color = "white",
    stroke = 0.5,
    # alpha = 0.85
  ) +
  # Ligne de régression avec intervalle de confiance
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "darkslateblue",
    fill = "darkslateblue",
    alpha = 0.15,
    linewidth = 1.5
  ) +
  # Annotation pour l'équation et la p-value (en haut à droite)
  annotate(
    "text",
    x = max(SEXTANT_1998_2025_spm_95$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_spm_95$mean_spm, na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "y = ", round(intercept_sextant_1998_95, 3), " + ", round(slope_sextant_1998_95, 7), " × x",
      "\n", "p = ", ifelse(p_value_sextant_1998_95 < 0.001, "< 0.001", format(p_value_sextant_1998_95, digits = 3))
    ),
    size = 8,
    color = "grey20",
    fontface = "italic",
    family = "serif"
  ) +
  # Titre et labels
  labs(
    title = "Évolution de la concentration moyenne en MES dans les panaches turbides de la baie des Anges",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    y = expression("Concentration en MES (g m"^{-3}*")")
  ) +
  # Thème sobre et élégant
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    plot.margin = margin(1, 1, 1, 1, "cm")  # Marges ajustées pour éviter le chevauchement
  ) +
  # Échelle des dates
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  )


# median spm
model_sextant_1998_95 <- lm(median_spm ~ date, data = SEXTANT_1998_2025_spm_95)
p_value_sextant_1998_95 <- summary(model_sextant_1998_95)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998_95 <- coef(model_sextant_1998_95)[1]
slope_sextant_1998_95 <- coef(model_sextant_1998_95)[2]

ggplot(data = SEXTANT_1998_2025_spm_95, aes(x = date, y = median_spm)) +
  # Points avec style épuré
  geom_point(
    size = 2,
    shape = 21,
    fill = "red3",
    color = "white",
    stroke = 0.5,
    # alpha = 0.85
  ) +
  # Ligne de régression avec intervalle de confiance
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "darkslateblue",
    fill = "darkslateblue",
    alpha = 0.15,
    linewidth = 1.5
  ) +
  # Annotation pour l'équation et la p-value (en haut à droite)
  annotate(
    "text",
    x = max(SEXTANT_1998_2025_spm_95$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_spm_95$mean_spm, na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "y = ", round(intercept_sextant_1998_95, 3), " + ", round(slope_sextant_1998_95, 7), " × x",
      "\n", "p = ", ifelse(p_value_sextant_1998_95 < 0.001, "< 0.001", format(p_value_sextant_1998_95, digits = 3))
    ),
    size = 8,
    color = "grey20",
    fontface = "italic",
    family = "serif"
  ) +
  # Titre et labels
  labs(
    title = "Évolution de la concentration médiane en MES dans les panaches turbides de la baie des Anges",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    y = expression("Concentration en MES (g m"^{-3}*")")
  ) +
  # Thème sobre et élégant
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    plot.margin = margin(1, 1, 1, 1, "cm")  # Marges ajustées pour éviter le chevauchement
  ) +
  # Échelle des dates
  scale_x_date(
    date_breaks = "5 years",
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


## Comparison --------------------------------------------------------------

### mean MES / liquid flow rate --------------------------------------------------

# mise à l'échelle
adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$mean_spm, Y6442010_depuis_2000$débit)
SEXTANT_1998_2025_spm_95$scaled_mean_spm <- SEXTANT_1998_2025_spm_95$mean_spm * adjust_factors$diff + adjust_factors$adjust

# Calcul de la corrélation entre débit et aire des panaches
merged_data <- merge(
  Y6442010_depuis_2000,
  SEXTANT_1998_2025_spm_95,
  by = "date",
  all = FALSE
)

correlation <- cor(merged_data$débit, merged_data$mean_spm, method = "spearman", use = "complete.obs")
p_value <- cor.test(merged_data$débit, merged_data$mean_spm, method = "spearman")$p.value

ggplot() +
  # Ligne pour le débit
  geom_line(
    data = Y6442010_depuis_2000,
    aes(x = date, y = débit, color = "Débit"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Ligne pour l'aire des panaches
  geom_line(
    data = SEXTANT_1998_2025_spm_95,
    aes(x = date, y = scaled_mean_spm, color = "Concentration en MES"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Couleurs personnalisées
  scale_color_manual(
    values = c("Concentration en MES" = "red3", "Débit" = "blue"),
    name = "Légende"
  ) +
  # Axes avec échelle secondaire
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Concentration en MES (g m"^{-3}*")")
    )
  ) +
  # Titre et labels
  labs(
    title = "Évolution de la concentration moyenne en Matière particulaire en suspension panaches turbides et du débit du Var",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    color = "Variable"
  ) +
  # Annotation pour la corrélation (en haut à droite)
  annotate(
    "text",
    x = max(c(Y6442010_depuis_2000$date, SEXTANT_1998_2025_spm_95$date), na.rm = TRUE),
    y = max(c(Y6442010_depuis_2000$débit, SEXTANT_1998_2025_spm_95$mean_spm), na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "R = ", round(correlation, 2),
      "\n", "p ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  # Thème sobre et élégant
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
  ) +
  # Échelle des dates
  scale_x_date(
    date_breaks = "5 year",
    date_labels = "%Y"
  )

### median MES / liquid flow rate --------------------------------------------------

# mise à l'échelle
adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$median_spm, Y6442010_depuis_2000$débit)
SEXTANT_1998_2025_spm_95$scaled_median_spm <- SEXTANT_1998_2025_spm_95$median_spm * adjust_factors$diff + adjust_factors$adjust

# Calcul de la corrélation entre débit et aire des panaches
merged_data <- merge(
  Y6442010_depuis_2000,
  SEXTANT_1998_2025_spm_95,
  by = "date",
  all = FALSE
)

correlation <- cor(merged_data$débit, merged_data$median_spm, method = "spearman", use = "complete.obs")
p_value <- cor.test(merged_data$débit, merged_data$median_spm, method = "spearman")$p.value

ggplot() +
  # Ligne pour le débit
  geom_line(
    data = Y6442010_depuis_2000,
    aes(x = date, y = débit, color = "Débit"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Ligne pour l'aire des panaches
  geom_line(
    data = SEXTANT_1998_2025_spm_95,
    aes(x = date, y = scaled_median_spm, color = "Concentration en MES"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Couleurs personnalisées
  scale_color_manual(
    values = c("Concentration en MES" = "red3", "Débit" = "blue"),
    name = "Légende"
  ) +
  # Axes avec échelle secondaire
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Concentration en MES (g m"^{-3}*")")
    )
  ) +
  # Titre et labels
  labs(
    title = "Évolution de la concentration médiane en Matière particulaire en suspension panaches turbides et du débit du Var",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    color = "Variable"
  ) +
  # Annotation pour la corrélation (en haut à droite)
  annotate(
    "text",
    x = max(c(Y6442010_depuis_2000$date, SEXTANT_1998_2025_spm_95$date), na.rm = TRUE),
    y = max(c(Y6442010_depuis_2000$débit, SEXTANT_1998_2025_spm_95$median_spm), na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "R = ", round(correlation, 2),
      "\n", "p ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  # Thème sobre et élégant
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
  ) +
  # Échelle des dates
  scale_x_date(
    date_breaks = "5 year",
    date_labels = "%Y"
  )

### Plume extension / liquid flow rate --------------------------------------------------

# mise à l'échelle
adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_spm_95$aire_panache_km2, Y6442010_depuis_2000$débit)
SEXTANT_1998_2025_spm_95$scaled_aire_panache <- SEXTANT_1998_2025_spm_95$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# Calcul de la corrélation entre débit et aire des panaches
merged_data <- merge(
  Y6442010_depuis_2000,
  SEXTANT_1998_2025_spm_95,
  by = "date",
  all = FALSE
)

correlation <- cor(merged_data$débit, merged_data$aire_panache_km2, method = "spearman", use = "complete.obs")
p_value <- cor.test(merged_data$débit, merged_data$aire_panache_km2, method = "spearman")$p.value

ggplot() +
  # Ligne pour le débit
  geom_line(
    data = Y6442010_depuis_2000,
    aes(x = date, y = débit, color = "Débit"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Ligne pour l'aire des panaches
  geom_line(
    data = SEXTANT_1998_2025_spm_95,
    aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
    size = 0.8,
    linewidth = 0.4
  ) +
  # Couleurs personnalisées
  scale_color_manual(
    values = c("Aire des panaches" = "darkcyan", "Débit" = "blue"),
    name = "Légende"
  ) +
  # Axes avec échelle secondaire
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = "Aire des panaches (km²)"
    )
  ) +
  # Titre et labels
  labs(
    title = "Évolution de l'extension des panaches turbides et du débit du Var",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    color = "Variable"
  ) +
  # Annotation pour la corrélation (en haut à droite)
  annotate(
    "text",
    x = max(c(Y6442010_depuis_2000$date, SEXTANT_1998_2025_spm_95$date), na.rm = TRUE),
    y = max(c(Y6442010_depuis_2000$débit, SEXTANT_1998_2025_spm_95$scaled_aire_panache), na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "R = ", round(correlation, 2),
      "\n", "p ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  # Thème sobre et élégant
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
  ) +
  # Échelle des dates
  scale_x_date(
    date_breaks = "5 year",
    date_labels = "%Y"
  )

# CHL ---------------------------------------------------------------------





# Example plot with hi-res coast ------------------------------------------

## This is just an example and should be delete after the hir-res coastline code has been integrated into your figure workflow

# Load one day of SEXTANT data as an example raster
SEXTANT_one <- load_SEXTANT_spm_pixels("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc", 
                                       lon_range = lon_range, lat_range = lat_range)

# Hi-res Mediterranean country and coastline shapes
coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# Plot with hi-res coast and SEXTANT raster
ggplot() +
  geom_raster(data = SEXTANT_one, aes(x = lon, y = lat, fill = analysed_spim)) +
  # geom_sf(data = coastline_giscoR, colour = "black", linewidth = 0.3) + # Coastlines
  geom_sf(data = countries_giscoR, colour = "black", linewidth = 0.3) + # Countries
  coord_sf(xlim = range(SEXTANT_one$lon), ylim = range(SEXTANT_one$lat), expand = FALSE) +
  scale_fill_viridis_c(option = "turbo", na.value = "transparent") +
  theme_minimal()

