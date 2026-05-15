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
library(ggspatial)


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
load("data/SEXTANT/CHL/SEXTANT_1998_2025_chl_pixels.RData")
load("data/SEXTANT/SPM/sextant_1998_2025_SPM.Rdata")

load("data/Hydro France/Y6442010_depuis_2000.Rdata")
load("~/River_runoff_analysis/data/Hydro France/Var_crues.Rdata")

## direction ---------------------------------------------------------------

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
# 
# SEXTANT_1998_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/1998/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_1999_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/1999/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2000_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2000/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2001_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2001/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2002_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2002/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2003_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2003/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2004_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2004/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2005_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2005/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2006_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2006/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2007_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2007/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2008_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2008/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2009_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2009/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2010_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2010/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2011_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2012_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2011_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2012_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2013_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2013/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2014_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2014/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2015_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2015/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2016_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2017_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2018_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2019_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2020_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2021_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2022_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2023_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2024_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# SEXTANT_2025_dir_chl <- dir("~/pCloudDrive/Stage/SEXTANT/CHLA/merged/Standard/DAILY/2025/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

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

# SEXTANT_1998_chl_pixels <- plyr::ldply(SEXTANT_1998_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_1999_chl_pixels <- plyr::ldply(SEXTANT_1999_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2000_chl_pixels <- plyr::ldply(SEXTANT_2000_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2001_chl_pixels <- plyr::ldply(SEXTANT_2001_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2002_chl_pixels <- plyr::ldply(SEXTANT_2002_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2003_chl_pixels <- plyr::ldply(SEXTANT_2003_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2004_chl_pixels <- plyr::ldply(SEXTANT_2004_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2005_chl_pixels <- plyr::ldply(SEXTANT_2005_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2006_chl_pixels <- plyr::ldply(SEXTANT_2006_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2007_chl_pixels <- plyr::ldply(SEXTANT_2007_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2008_chl_pixels <- plyr::ldply(SEXTANT_2008_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2009_chl_pixels <- plyr::ldply(SEXTANT_2009_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2010_chl_pixels <- plyr::ldply(SEXTANT_2010_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2011_chl_pixels <- plyr::ldply(SEXTANT_2011_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2012_chl_pixels <- plyr::ldply(SEXTANT_2012_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2013_chl_pixels <- plyr::ldply(SEXTANT_2013_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2014_chl_pixels <- plyr::ldply(SEXTANT_2014_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2015_chl_pixels <- plyr::ldply(SEXTANT_2015_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2016_chl_pixels <- plyr::ldply(SEXTANT_2016_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2017_chl_pixels <- plyr::ldply(SEXTANT_2017_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2018_chl_pixels <- plyr::ldply(SEXTANT_2018_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2019_chl_pixels <- plyr::ldply(SEXTANT_2019_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2020_chl_pixels <- plyr::ldply(SEXTANT_2020_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2021_chl_pixels <- plyr::ldply(SEXTANT_2021_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2022_chl_pixels <- plyr::ldply(SEXTANT_2022_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2023_chl_pixels <- plyr::ldply(SEXTANT_2023_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2024_chl_pixels <- plyr::ldply(SEXTANT_2024_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# SEXTANT_2025_chl_pixels <- plyr::ldply(SEXTANT_2025_dir_chl, load_SEXTANT_chl_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# 
# # Combine and save
# SEXTANT_1998_2025_chl_pixels <- rbind(SEXTANT_1998_chl_pixels, SEXTANT_1999_chl_pixels, SEXTANT_2000_chl_pixels,
#                                       SEXTANT_2001_chl_pixels, SEXTANT_2002_chl_pixels, SEXTANT_2003_chl_pixels,
#                                       SEXTANT_2004_chl_pixels, SEXTANT_2005_chl_pixels, SEXTANT_2006_chl_pixels,
#                                       SEXTANT_2007_chl_pixels, SEXTANT_2008_chl_pixels, SEXTANT_2009_chl_pixels,
#                                       SEXTANT_2010_chl_pixels, SEXTANT_2011_chl_pixels, SEXTANT_2012_chl_pixels,
#                                       SEXTANT_2013_chl_pixels, SEXTANT_2014_chl_pixels, SEXTANT_2015_chl_pixels,
#                                       SEXTANT_2016_chl_pixels, SEXTANT_2017_chl_pixels, SEXTANT_2018_chl_pixels,
#                                       SEXTANT_2019_chl_pixels, SEXTANT_2020_chl_pixels, SEXTANT_2021_chl_pixels,
#                                       SEXTANT_2022_chl_pixels, SEXTANT_2023_chl_pixels,SEXTANT_2024_chl_pixels,
#                                       SEXTANT_2025_chl_pixels)
# 
# save(SEXTANT_1998_2025_chl_pixels, file = "data/SEXTANT/CHL/SEXTANT_1998_2025_chl_pixels.RData")

## paramètres journaliers --------------------------------------------------

# combien de valeurs négatives
sum(SEXTANT_1998_2025_spm_pixels$analysed_spim < 0, na.rm = TRUE)
# [1] 21080

# supprimer seulement les valeurs négatives
SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_pixels |>
  filter(analysed_spim >= 0 | is.na(analysed_spim))

SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# paramètres journaliers
SEXTANT_1998_2025_spm_mean <- SEXTANT_1998_2025_spm_clean |>
  mutate(date = as.Date(date)) |>
  summarise(
    mean_spm   = mean(analysed_spim, na.rm = TRUE),
    median_spm = median(analysed_spim, na.rm = TRUE),
    sd_spm     = sd(analysed_spim, na.rm = TRUE),
    .by = date
  )

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

# save(SEXTANT_1998_2025_spm_95, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_spm_95.Rdata")

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
  
  pixels_zone <- SEXTANT_1998_2025_spm_clean |>
    filter(
      lon >= lon_embouchure - r & lon <= lon_embouchure + r,
      lat >= lat_embouchure - r & lat <= lat_embouchure + r
    )
  
  aire_km2 <- nrow(distinct(pixels_zone, lon, lat)) * aire_pixel_km2
  seuil    <- quantile(pixels_zone$analysed_spim, 0.95, na.rm = TRUE)
  
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
  filter(aire_km2 > 2459) |>
  summarise(seuil = mean(seuil_95)) |>
  pull(seuil)

cat("Seuil retenu :", seuil_retenu, "g/m³\n", na.rm = TRUE)
# 0.94 g/m³

## ROPP --------------------------------------------------------------------

# Pixels où le panache est présent dans au moins 5% des images
n_images_total <- n_distinct(SEXTANT_1998_2025_spm_clean$date)

ROPP <- SEXTANT_1998_2025_spm_clean |>
  group_by(lon, lat) |>
  summarise(
    freq_above_seuil = sum(analysed_spim >= seuil_retenu, na.rm = TRUE) / n_images_total,
    .groups = "drop"
  ) |>
  filter(freq_above_seuil >= 0.05)   # au moins 5% des images

cat("Nombre de clean dans la ROPP :", nrow(ROPP), "\n")

## filtrer les images avec trop de clean manquants ------------------------

# Garder seulement les images avec > 80% de clean valides sur la ROPP
clean_ROPP <- SEXTANT_1998_2025_spm_clean |>
  semi_join(ROPP, by = c("lon", "lat"))

images_valides <- clean_ROPP |>
  group_by(date) |>
  summarise(
    n_clean_valides = sum(!is.na(analysed_spim)),
    n_clean_total   = n(),
    pct_valide       = n_clean_valides / n_clean_total,
    .groups = "drop"
  ) |>
  filter(pct_valide >= 0.80)

cat("Images valides :", nrow(images_valides), "/", n_distinct(SEXTANT_1998_2025_spm_clean$date), "\n")


# Diagnostics à faire absolument
cat("Seuil retenu :", seuil_retenu, "\n")
cat("Nombre de pixels dans la ROPP :", nrow(ROPP), "\n")
cat("Images valides :", nrow(images_valides), "\n")
cat("Jours avec panache > 0 :", sum(SEXTANT_panache_metrics$aire_panache_km2 > 0), "\n")

# Distribution des aires de panache
summary(SEXTANT_panache_metrics$aire_panache_km2)
hist(SEXTANT_panache_metrics$aire_panache_km2, breaks = 50)

# Distribution du débit 3j associé
summary(SEXTANT_panache_metrics$debit_3j_mean)

## stat panache ------------------------------------------------------------

# Statistiques du panache par jour
SEXTANT_panache_metrics <- SEXTANT_1998_2025_spm_clean |>
  filter(date %in% images_valides$date) |>
  semi_join(ROPP, by = c("lon", "lat")) |>
  group_by(date) |>
  summarise(
    # Aire d'extension
    pixel_count = sum(analysed_spim >= seuil_retenu, na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    # Concentrations
    mean_spm         = mean(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    max_spm          = max(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    median_spm       = median(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    # Points extrêmes
    lat_sud          = min(lat[analysed_spim >= seuil_retenu], na.rm = TRUE),
    lon_ouest        = min(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    lon_est          = max(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    # Centroïde
    centroid_lon     = mean(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    centroid_lat     = mean(lat[analysed_spim >= seuil_retenu], na.rm = TRUE),
    .groups = "drop"
  )

# plotting ----------------------------------------------------------------

# correlation river flow and plume extension ------------------------------

debit_lags <- Y6442010_depuis_2000 |>
  arrange(date) |>
  select(date, débit) |> 
  mutate(
    debit_j1     = lag(débit, 1),             # j-1 seulement
    debit_2j     = (lag(débit, 1) + lag(débit, 2)) / 2,       # moyenne j-1, j-2
    debit_3j     = (lag(débit, 1) + lag(débit, 2) + lag(débit, 3)) / 3
  )

# Jointure et nettoyage
SEXTANT_panache_metrics <- SEXTANT_panache_metrics |>
  inner_join(debit_lags |> select(date, débit, debit_j1, debit_2j, debit_3j), by = "date") |>
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

# # Tout est déjà dans SEXTANT_panache_metrics, pas besoin de rejoindre !
SEXTANT_panache_metrics |>
  filter(aire_panache_km2 > 0) |>
  summarise(
    r_j0 = cor(log10(aire_panache_km2), debit_j0_log, use = "complete.obs", method = "spearman"),
    r_j1 = cor(log10(aire_panache_km2), debit_j1_log, use = "complete.obs", method = "spearman"),
    r_2j = cor(log10(aire_panache_km2), debit_2j_log, use = "complete.obs", method = "spearman"),
    r_3j = cor(log10(aire_panache_km2), debit_3j_log, use = "complete.obs", method = "spearman")
  )

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ debit_j0_log, data = SEXTANT_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1_log, data = SEXTANT_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j_log, data = SEXTANT_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j_log, data = SEXTANT_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ débit, data = SEXTANT_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1, data = SEXTANT_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j, data = SEXTANT_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j, data = SEXTANT_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

## modèle log - log --------------------------------------------------------

# Modèle linéaire en log-log → relation puissance
modele_log <- lm(panache_log ~ debit_j0_log, data = SEXTANT_panache_metrics)
r2     <- summary(modele_log)$r.squared
print(r2)
pente  <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

# Équation puissance : aire = 10^ordonnee * débit^pente
# à afficher dans le graphique
label_eq <- paste0(
  "Aire = ", round(10^ordonnee, 3), " × Q^", round(pente, 2),
  "\nR² = ", round(r2, 2)
)

# Graphique log-log avec droite de régression
ggplot(SEXTANT_panache_metrics, aes(x = débit, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,            # régression sur les axes log
              color = "black", se = FALSE, linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text", 
           x = min(SEXTANT_panache_metrics$debit_log, na.rm = TRUE) * 1.5,
           y = max(SEXTANT_panache_metrics$aire_panache_km2, na.rm = TRUE) * 0.7,
           label = label_eq, hjust = 1.1, vjust = 2, size = 8, color = "grey20",
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



# CHL ---------------------------------------------------------------------

## paramètres journaliers --------------------------------------------------

# combien de valeurs négatives
sum(SEXTANT_1998_2025_chl_pixels$analysed_chl_a < 0, na.rm = TRUE)
# [1] 19808

# supprimer seulement les valeurs négatives
SEXTANT_1998_2025_chl_clean <- SEXTANT_1998_2025_chl_pixels |>
  filter(analysed_chl_a >= 0, analysed_chl_a <= 20 | is.na(analysed_chl_a))

SEXTANT_1998_2025_chl_clean <- SEXTANT_1998_2025_chl_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# paramètres journaliers
SEXTANT_1998_2025_chl_mean <- SEXTANT_1998_2025_chl_clean |>
  mutate(date = as.Date(date)) |>
  summarise(
    mean_chl   = mean(analysed_chl_a, na.rm = TRUE),
    median_chl = median(analysed_chl_a, na.rm = TRUE),
    sd_chl     = sd(analysed_chl_a, na.rm = TRUE),
    .by = date
  )

## climatologie ------------------------------------------------------------

# on choisit une période de 20 ans (1998 - 2017)
SEXTANT_1998_2017 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("1998-01-01"), date <= as.Date("2017-12-31"))

SEXTANT_1998_2017_stat <- SEXTANT_1998_2017 |>
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# climatologie annuelle
SEXTANT_1998_2017_chl_year <- SEXTANT_1998_2017_stat |> 
  group_by(year) |> 
  summarise(mean_chl_year_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_year_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_year_clim = sd(analysed_chl_a, na.rm = TRUE))

# climatologie mensuelle
SEXTANT_1998_2017_chl_month <- SEXTANT_1998_2017_stat |> 
  group_by(month) |>
  summarise(mean_chl_month_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_month_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_month_clim = sd(analysed_chl_a, na.rm = TRUE))

# climatologie journanlière
SEXTANT_1998_2017_chl_doy <- SEXTANT_1998_2017_stat |>
  group_by(doy) |> 
  summarise(mean_chl_doy_clim = mean(analysed_chl_a, na.rm = TRUE), 
            median_chl_doy_clim = median(analysed_chl_a, na.rm = TRUE),
            sd_chl_doy_clim = sd(analysed_chl_a, na.rm = TRUE))

# créer une climatologie à partir d'une TS journalière
sextant_chl_climatology_doy <- ts2clm(data = SEXTANT_1998_2025_chl_mean, x = date, 
                                      y = mean_chl, climatologyPeriod = c("1998-01-01", "2017-12-31"), 
                                      windowHalfWidth = 3, smoothPercentileWidth = 15 )

# anomalie mensuelle
sextant_1998_2025_chl_monthly_anom <- SEXTANT_1998_2025_chl_clean |> 
  mutate(date = floor_date(date, "month")) |> 
  # filter(date >= as.character.Date("1998-01-01"), date <= as.Date ("2025-12-31")) |> 
  summarise(mean_chl_month = mean(analysed_chl_a, na.rm = TRUE), .by = c("date", "year", "month")) |>
  left_join(SEXTANT_1998_2017_chl_month, by = c("month")) |> 
  mutate(chl_month_anomaly = mean_chl_month - mean_chl_month_clim)

## plotting ----------------------------------------------------------------

### cartographie --------------------------------------------------

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# on choisit des jours précis (valeurs aberrantes, jours de crues, jours sans crues)

# valeur moyenne haute
SEXTANT_1998_2025_chl_2025_02_23 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("2025-02-23"), date <= as.Date("2025-02-23"))

# crues
SEXTANT_1998_2025_chl_2020_10_03 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("2020-10-03"), date <= as.Date("2020-10-03"))

# le 05-11-2011 très forte crue (débit de 910 m3/s du Var) + sur les images MODIS 
# nuages
SEXTANT_1998_2025_chl_2020_11_05 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("2020-11-05"), date <= as.Date("2020-11-05"))

# le 10-11-2011 on voit un panache sur les images MODIS
SEXTANT_1998_2025_chl_2020_11_10 <- SEXTANT_1998_2025_chl_clean |> 
  filter(date >= as.Date("2020-11-10"), date <= as.Date("2020-11-10"))

# plot
ggplot(data = SEXTANT_1998_2025_chl_2020_11_10) +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = analysed_chl_a)) +
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
    name   = expression("Concentration en chlorophylle a (µg.L"^{-1}*")")  # ← écriture scientifique
    # limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en chlorophylle a — 10 novembre 2011",
    subtitle = "Produit SEXTANT OC5",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(SEXTANT_1998_2025_chl_2020_11_10$lon),
    ylim   = range(SEXTANT_1998_2025_chl_2020_11_10$lat),
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

### paramètres journaliers --------------------------------------------------

# Calcul du modèle linéaire
model_sextant_chl_1998_2017 <- lm(median_chl ~ date, data = SEXTANT_1998_2025_chl_mean)
p_value_sextant_chl_1998_2017 <- summary(model_sextant_chl_1998_2017 )$coefficients[2, 4]
intercept_sextant_chl_1998_2017 <- coef(model_sextant_chl_1998_2017 )[1]
slope_sextant_chl_1998_2017 <- coef(model_sextant_chl_1998_2017 )[2]

# Création du graphique
ggplot(data = SEXTANT_1998_2025_chl_mean, aes(x = date, y = median_chl)) +
  # Points avec style épuré
  # geom_ribbon(
  #   aes(ymin = mean_chl - sd_chl, ymax = mean_chl + sd_chl),
  #   fill = "chartreuse3",
  #   alpha = 0.3
  # ) +
  geom_point(
    size = 2,
    shape = 21,
    fill = "chartreuse3",
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
    x = max(SEXTANT_1998_2025_chl_mean$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_chl_mean$median_chl, na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    label = paste0(
      "y = ", round(intercept_sextant_chl_1998_2017, 3), " ", round(slope_sextant_chl_1998_2017, 7), " × x",
      "\n", "p ", ifelse(p_value_sextant_chl_1998_2017 < 0.001, "< 0.001", format(p_value_sextant_chl_1998_2017, digits = 3))
    ),
    size = 8,
    color = "grey20",
    fontface = "italic",
    family = "serif"
  ) +
  # Titre et labels
  labs(
    title = "Évolution de la concentration médiane journalière en chlorophylle a dans la zone d'étude",
    subtitle = "Produit SEXTANT OC5 (1998-2025)",
    x = "Date",
    y = expression("Concentration en chlorophylle a (µg.L"^{-1}*")")
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
    date_labels = "%Y")
# +
  # # Zoom sur l'axe Y sans supprimer les données
  # coord_cartesian(ylim = c(-1, 2))

### climatologie --------------------------------------------------

# create a line plot of the yearly climatology of chl
ggplot(SEXTANT_1998_2017_chl_year, aes(x = year, y = mean_chl_year_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_year_clim - sd_chl_year_clim,
      ymax = mean_chl_year_clim + sd_chl_year_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie annuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie annuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie annuelle" = "chartreuse3")
  # ) +
  # scale_x_continuous(
  #   breaks = 1:12,
  #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
  #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie annuelle de la concentration en chlorophylle a (1998–2017) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2017 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the monthly climatology of chl
ggplot(SEXTANT_1998_2017_chl_month, aes(x = month, y = mean_chl_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_month_clim - sd_chl_month_clim,
      ymax = mean_chl_month_clim + sd_chl_month_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "chartreuse3")
    ) +
    scale_x_continuous(
      breaks = 1:12,
      labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
                 "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration en chlorophylle a (1998–2017) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2017 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the daily climatology of chl
ggplot(SEXTANT_1998_2017_chl_doy, aes(x = doy, y = mean_chl_doy_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_doy_clim - sd_chl_doy_clim,
      ymax = mean_chl_doy_clim + sd_chl_doy_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie journalière"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie journalière"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie journalière" = "chartreuse3")
  ) +
  # scale_x_continuous(
  #   breaks = 1:12,
  #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
  #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  # ) +
  labs(
    title   = "Climatologie journalière de la concentration en chlorophylle a (1998–2017) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2017 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# anomalie mensuelle
model_sextant_chl_anom <- lm(chl_month_anomaly ~ date, data = sextant_1998_2025_chl_monthly_anom)
p_value_sextant_chl_anom <- summary(model_sextant_chl_anom)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_chl_anom <- coef(model_sextant_chl_anom)[1]
slope_sextant_chl_anom <- coef(model_sextant_chl_anom)[2]

# Créer le graphique
ggplot(sextant_1998_2025_chl_monthly_anom, aes(x = date, y = chl_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(chl_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(chl_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_chl_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_chl_monthly_anom$chl_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_chl_anom, 3), " ", round(slope_sextant_chl_anom, 7), " × x",
      "\np = ", ifelse(p_value_sextant_chl_anom < 0.001, "< 0.001", format(p_value_sextant_chl_anom, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en chlorophylle a (1998–2025) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 1998–2017"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

### patchwork ----------------------------------------------------------------

# --- Graphique 1 : climatologie mensuelle ---
p1 <- ggplot(SEXTANT_1998_2017_chl_month, aes(x = month, y = mean_chl_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_chl_month_clim - sd_chl_month_clim,
      ymax = mean_chl_month_clim + sd_chl_month_clim
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "chartreuse3")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration en chlorophylle a (1998–2017) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Période de référence : 1998–2017 | Barres : ± 1 écart-type"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# Extraire le modèle linéaire
model_sextant_chl_anom <- lm(chl_month_anomaly ~ date, data = sextant_1998_2025_chl_monthly_anom)
p_value_sextant_chl_anom <- summary(model_sextant_chl_anom)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_chl_anom <- coef(model_sextant_chl_anom)[1]
slope_sextant_chl_anom <- coef(model_sextant_chl_anom)[2]

# --- Graphique 2 : anomalie mensuelle ---
p2 <- ggplot(sextant_1998_2025_chl_monthly_anom, aes(x = date, y = chl_month_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = pmin(chl_month_anomaly, 0), ymax = 0),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = pmax(chl_month_anomaly, 0)),
    fill = "tomato", alpha = 0.3
  ) +
  geom_line(aes(color = "Anomalie mensuelle"), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(
    aes(color = "Tendance linéaire", fill = "Tendance linéaire"),
    method = "lm", se = TRUE, alpha = 0.15, linewidth = 1
  ) +
  annotate(
    "text",
    x = min(sextant_1998_2025_chl_monthly_anom$date, na.rm = TRUE),
    y = max(sextant_1998_2025_chl_monthly_anom$chl_month_anomaly, na.rm = TRUE) * 0.95,
    label = paste0(
      "y = ", round(intercept_sextant_chl_anom, 3), " ", round(slope_sextant_chl_anom, 7), " × x",
      "\np = ", ifelse(p_value_sextant_chl_anom < 0.001, "< 0.001", format(p_value_sextant_chl_anom, digits = 3))
    ),
    hjust = 0, vjust = 1,
    size = 6,
    color = "grey20",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Anomalie mensuelle" = "grey30", "Tendance linéaire" = "firebrick")
  ) +
  scale_fill_manual(
    values = c("Tendance linéaire" = "firebrick"),
    guide  = "none"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title   = "Anomalie mensuelle de la concentration en chlorophylle a (1998–2025) — Sextant OC5",
    x       = NULL,
    y       = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    color   = NULL,
    caption = "Source : Sextant OC5 | Climatologie de référence : 1998–2017"
  ) +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption       = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y       = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 10)),
    axis.text          = element_text(size = 12, color = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.ticks         = element_line(color = "grey70"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", linewidth = 0.5),
    legend.position    = "top",
    legend.text        = element_text(size = 11)
  )

p1 / p2 +
  plot_annotation(
    title   = "Concentration en chlorophylle a — Sextant OC5",
    caption = "Source : Sextant OC5",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )

## correlation MES / Chl a -------------------------------------------------------------

# Fusionner tes deux jeux de données par date
data_merged <- inner_join(
  SEXTANT_1998_2025_chl_mean,
  sextant_1998_2025_SPM,
  by = "date"
)

# Corrélation globale
cor.test(data_merged$mean_chl, data_merged$mean_spm, 
         method = "spearman")  # Spearman car distributions asymétriques

# Visualisation
ggplot(data_merged, aes(x = mean_spm, y = mean_chl)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +  # loess pour capturer la non-linéarité
  labs(x = "MES (mg/L)", y = "Chlorophylle a (µg/L)")


## série temporelle --------------------------------------------------------

adjust_factors <- sec_axis_adjustement_factors(SEXTANT_1998_2025_chl_mean$mean_chl, sextant_1998_2025_SPM$mean_spm)

sextant_1998_2025_SPM$scaled_mean_spm <- sextant_1998_2025_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

correlation <- cor(data_merged$mean_chl, data_merged$mean_spm, method = "spearman", use = "complete.obs")
p_value <- cor.test(data_merged$mean_chl, data_merged$mean_spm, method = "spearman")$p.value

ggplot() +
  geom_ribbon(
    data = SEXTANT_1998_2025_chl_mean,
    aes(
      x    = date,
      ymin = (mean_chl - sd_chl),
      ymax = (mean_chl + sd_chl)
    ),
    fill = "chartreuse3", alpha = 0.2
  ) +
  geom_line(
    data = SEXTANT_1998_2025_chl_mean,
    aes(x = date, y = scaled_mean_chl, color = "CHL"), size = 0.5
  ) +
  # Ribbon SPM — x manquant + ymin/ymax pas scalés car SPM est sur axe principal
  geom_ribbon(
    data = sextant_1998_2025_SPM,
    aes(
      x    = date,
      ymin = (mean_spm - std_spm) * adjust_factors$diff + adjust_factors$adjust,,
      ymax = (mean_spm + std_spm) * adjust_factors$diff + adjust_factors$adjust,
    ),
    fill = "red3", alpha = 0.2
  ) +
  geom_line(
    data = sextant_1998_2025_SPM,
    aes(x = date, y = mean_spm, color = "SPM"), size = 0.5
  ) +
  scale_color_manual(values = c("CHL" = "chartreuse3", "SPM" = "red3")) +
  scale_y_continuous(
    name = expression("Concentration en MES (g.m"^{-3}*")"),
    limits = c(-5, 15),
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    )
  ) +
  labs(
    title = "Concentration en chlorophylle a et en matière en suspension (1998–2025) — Sextant OC5",
    x = "Date",
    caption = "Source : Sextant OC5"
  ) +
  annotate(
    "text",
    x = max(c(SEXTANT_1998_2025_chl_mean$date, sextant_1998_2025_SPM$date), na.rm = TRUE),
    y = max(c(SEXTANT_1998_2025_chl_mean$mean_chl, sextant_1998_2025_SPM$mean_spm), na.rm = TRUE),
    hjust = 1,  # Alignement à droite
    vjust = 0,  # Alignement en haut
    label = paste0(
      "R = ", round(correlation, 2),
      "\n", "p ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
    ),
    size = 8,
    color = "grey20",
    family = "serif",
    fontface = "italic"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x     = element_text(size = 13, margin = margin(t = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )


# Analyse par seuil de turbidité ------------------------------------------

data_merged <- data_merged |>
  mutate(
    Turbidité = case_when(
      mean_spm < quantile(mean_spm, 0.33, na.rm = TRUE) ~ "Faible",
      mean_spm < quantile(mean_spm, 0.66, na.rm = TRUE) ~ "Modérée",
      TRUE                                               ~ "Forte"
    ),
    # Définir l'ordre des niveaux
    Turbidité = factor(Turbidité, levels = c("Faible", "Modérée", "Forte"))
  )

# Boxplot comparatif
ggplot(data_merged, aes(x = Turbidité, y = mean_chl, fill = Turbidité)) +
  geom_boxplot() +
  labs(x = "Classe de turbidité", y = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"))


# Test statistique
kruskal.test(mean_chl ~ Turbidité, data = data_merged)  # non-paramétrique

# Test de Kruskal-Wallis
kw_test <- kruskal.test(mean_chl ~ Turbidité, data = data_merged)

# Extraire la statistique et la p-value
kw_stat <- kw_test$statistic
kw_pvalue <- kw_test$p.value


ggplot(data_merged, aes(x = Turbidité, y = mean_chl, fill = Turbidité)) +
  geom_boxplot() +
  labs(
    x = "Classe de turbidité",
    y = expression("Concentration en chlorophylle a (µg.L"^{-1}*")"),
    title = "Comparaison de la concentration moyenne en chlorophylle-a par classe de turbidité"
  ) +
  annotate(
    "text",
    x = 1.5,  # Position horizontale (ajustez selon le nombre de catégories)
    y = Inf,   # Position verticale (en haut du graphique)
    label = paste0(
      "Test de Kruskal-Wallis: χ² = ", round(kw_stat, 2),
      "\n", "p ", ifelse(kw_pvalue < 0.001, "< 0.001", format(kw_pvalue, digits = 3))
    ),
    vjust = 1.5,  # Ajuste la position verticale du texte
    hjust = 0.5,  # Centre le texte horizontalement
    size = 8,     # Taille du texte
    color = "black",
    family = "serif",
    fontface = "italic"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 11, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.title.x     = element_text(size = 13, margin = margin(t = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# analyse saisonnière ------------------------------------

data_merged <- data_merged |>
  mutate(
    month = month(date),
    saison = case_when(
      month %in% c(12, 1, 2)  ~ "Hiver",
      month %in% c(3, 4, 5)   ~ "Printemps",
      month %in% c(6, 7, 8)   ~ "Été",
      month %in% c(9, 10, 11) ~ "Automne"
    )
  )

# Corrélation par saison
data_merged |>
  group_by(saison) |>
  summarise(
    cor_spearman = cor(mean_chl, mean_spm, method = "spearman", use = "complete.obs"),
    p_value = cor.test(mean_chl, mean_spm, method = "spearman")$p.value
  )

# Visualisation
ggplot(data_merged, aes(x = mean_spm, y = mean_chl, color = saison)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~saison)


# lag correlation ---------------------------------------------------------

# Cross-corrélation pour identifier le délai optimal
ccf_result <- ccf(
  data_merged$mean_spm, 
  data_merged$mean_chl, 
  lag.max = 30,  # tester jusqu'à 30 jours de décalage
  na.action = na.pass,
  main = "Cross-corrélation MES – Chlorophylle a"
)

# Le lag avec la corrélation maximale indique le délai de réponse


# analyse spatiale --------------------------------------------------------

# comment la relation évolue avec la distance à l'embouchure

# Calculer la distance de chaque pixel à l'embouchure
# (coordonnées de l'embouchure à adapter)
lon_embouchure <- 7.199082
lat_embouchure <- 43.654709

SEXTANT_chl_spm_pixels <- SEXTANT_1998_2025_chl_clean |>
  mutate(date = as.Date(date)) |>
  inner_join(
    SEXTANT_1998_2025_spm_clean |> mutate(date = as.Date(date)),
    by = c("date", "lon", "lat")
  )

data_pixels <- SEXTANT_chl_spm_pixels |>
  mutate(
    distance_km = sqrt((lon - lon_embouchure)^2 + 
                         (lat - lat_embouchure)^2) * 111,
    zone = cut(distance_km, 
               breaks = c(0, 10, 20, 40, Inf),
               labels = c("0-10 km", "10-20 km", "20-40 km", ">40 km"))
  )

# Corrélation chl ~ MES par classe de distance
data_pixels |>
  group_by(zone) |>
  summarise(cor = cor(analysed_chl_a, analysed_spim, 
                      method = "spearman", use = "complete.obs"))

# Moyenne CHL par pixel (position fixe dans l'espace)
chl_mean_par_pixel <- data_pixels |>
  group_by(lon, lat, distance_km, zone) |>
  summarise(
    mean_chl = mean(analysed_chl_a, na.rm = TRUE),
    mean_spm = mean(analysed_spim,   na.rm = TRUE),   # ← adapter
    .groups = "drop"
  )

# Visualisation du gradient spatial
ggplot(chl_mean_par_pixel, aes(x = distance_km, y = mean_chl)) +
  geom_point(alpha = 0.3, color = "chartreuse3") +
  geom_smooth(method = "loess", color = "darkgreen") +
  labs(
    x = "Distance à l'embouchure (km)",
    y = expression("Chlorophylle a moyenne (µg.L"^{-1}*")")
  ) +
  theme_bw()

# test de la significativité de chaque corrélation
# data_pixels |>
#   group_by(zone) |>
#   summarise(
#     cor_spearman = cor(analysed_chl_a, analysed_spim,
#                        method = "spearman", use = "complete.obs"),
#     p_value = cor.test(analysed_chl_a, analysed_spim,
#                        method = "spearman")$p.value,
#     n = n()
#   )

data_pixels |>
  group_by(zone) |>
  summarise(
    cor_spearman = cor(analysed_chl_a, analysed_spim,
                       method = "spearman", use = "complete.obs"),
    # Calcul manuel de la p-value via approximation t (valide pour grands n)
    n = n(),
    t_stat  = cor_spearman * sqrt((n - 2) / (1 - cor_spearman^2)),
    p_value = 2 * pt(-abs(t_stat), df = n - 2),
    p_label = ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  ) |>
  select(zone, cor_spearman, n, p_label)

# Stocker les résultats
cor_par_zone <- data_pixels |>
  group_by(zone) |>
  summarise(
    cor_spearman = cor(analysed_chl_a, analysed_spim,
                       method = "spearman", use = "complete.obs"),
    n       = n(),
    t_stat  = cor_spearman * sqrt((n - 2) / (1 - cor_spearman^2)),
    p_value = 2 * pt(-abs(t_stat), df = n - 2),
    p_label = ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3)),
    .groups = "drop"
  )

ggplot(cor_par_zone, aes(x = zone, y = cor_spearman, fill = zone)) +
  geom_col() +
  geom_text(
    aes(label = paste0("r = ", round(cor_spearman, 3),
                       "\np ", p_label,
                       "\nn = ", formatC(n, format = "d", big.mark = " "))),
    vjust = -0.5, size = 4
  ) +
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits = c(0, 0.85)) +
  labs(
    title    = "Corrélation de Spearman CHL ~ MES par zone de distance",
    subtitle = "Décroissance de la corrélation avec l'éloignement de l'embouchure",
    x        = "Zone de distance",
    y        = "Corrélation de Spearman"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Visualiser le gradient de corrélation pixel par pixel

# Corrélation par pixel (chaque position lon/lat)
cor_par_pixel <- data_pixels |>
  group_by(lon, lat, distance_km) |>
  summarise(
    cor_spearman = cor(analysed_chl_a, analysed_spim,
                       method = "spearman", use = "complete.obs"),
    n = n(),
    .groups = "drop"
  ) |>
  filter(n >= 30)  # garder seulement les pixels avec assez d'observations

# Carte de la corrélation spatiale
ggplot(cor_par_pixel, aes(x = lon, y = lat, fill = cor_spearman)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("blue", "white", "red"),
    limits  = c(-1, 1),
    name    = "Corrélation\nSpearman"
  ) +
  labs(
    title = "Corrélation spatiale CHL ~ MES",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw() +
  coord_fixed()  # pour ne pas déformer la carte

# Visualiser la décroissance de la corrélation avec la distance
ggplot(cor_par_pixel, aes(x = distance_km, y = cor_spearman)) +
  geom_point(alpha = 0.3, color = "grey50") +
  geom_smooth(method = "loess", color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Distance à l'embouchure (km)",
    y = "Corrélation de Spearman (CHL ~ MES)",
    title = "Décroissance de la corrélation CHL–MES avec la distance"
  ) +
  theme_bw()

# analyse par saison
data_pixels |>
  mutate(
    month  = month(date),
    saison = case_when(
      month %in% c(12, 1, 2)  ~ "Hiver",
      month %in% c(3, 4, 5)   ~ "Printemps",
      month %in% c(6, 7, 8)   ~ "Été",
      month %in% c(9, 10, 11) ~ "Automne"
    )
  ) |>
  group_by(zone, saison) |>
  summarise(
    cor_spearman = cor(analysed_chl_a, analysed_spim,
                       method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) |>
  # Visualisation en heatmap
  ggplot(aes(x = saison, y = zone, fill = cor_spearman)) +
  geom_tile() +
  geom_text(aes(label = round(cor_spearman, 2)), size = 5) +
  scale_fill_gradientn(
    colours = c("white" , "blue", "red"),
    limits  = c(0, 1),
    name    = "Corrélation\nSpearman"
  ) +
  labs(
    title = "Corrélation CHL ~ MES par zone et par saison",
    x = "Saison", y = "Zone de distance"
  ) +
  theme_bw()

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

