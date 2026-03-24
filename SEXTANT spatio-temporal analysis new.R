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
library(doParallel); registerDoParallel(cores = 14)

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)


# Functions ---------------------------------------------------------------

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

# All of the SEXTANT 1998 files
# sextant_1998_dir <- dir("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1998", pattern = ".nc", recursive = TRUE, full.names = TRUE)
# sextant_1999_dir <- dir("~/pCloudDrive/data/SEXTANT/SPM/merged/Standard/DAILY/1999", pattern = ".nc", recursive = TRUE, full.names = TRUE)

# Load and combine
# system.time(
# sextant_1998 <- plyr::ldply(sextant_1998_dir, load_sextant, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# )
# sextant_1999 <- plyr::ldply(sextant_1999_dir, load_sextant, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Combine and save
# sextant_1998_1999 <- rbind(sextant_1998, sextant_1999)
# save(sextant_1998, file = "data/SEXTANT/sextant_1998.RData")

# Load the data
# load("data/SEXTANT/SPM/all_spm_propre_sextant_1998.RData")

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

load("data/SEXTANT/SPM/SEXTANT_1998_2025_SPM_spatial.RData")

### no threshold --------------------------------------------------------

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

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_spatial.Rdata")



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
SEXTANT_1998_2025_spm_spatial <- SEXTANT_1998_2025_spm_spatial |> 
  mutate(aire_panache_km2 = pixel_count * aire_pixel_km2)


# plotting ----------------------------------------------------------------

model_sextant_1998_spatial <- lm(mean_spm ~ date, data = SEXTANT_1998_2025_spm_spatial)
p_value_sextant_1998_spatial <- summary(model_sextant_1998_spatial)$coefficients[2, 4]  # p-value pour la pente
intercept_sextant_1998_spatial <- coef(model_sextant_1998_spatial)[1]
slope_sextant_1998_spatial <- coef(model_sextant_1998_spatial)[2]

ggplot(data = SEXTANT_1998_2025_spm_spatial, aes(x = date, y = mean_spm)) +
  geom_point(color = "deepskyblue", size = 0.5) +
  # geom_point(data = SEXTANT_1998_2025_spm_spatial, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(SEXTANT_1998_2025_spm_spatial$date, na.rm = TRUE),
    y = max(SEXTANT_1998_2025_spm_spatial$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_sextant_1998_spatial, 3), " + ", round(slope_sextant_1998_spatial, 7), " * x",
      "\n", "p = ", ifelse(p_value_sextant_1998_spatial < 0.001, "< 0.001", format(p_value_sextant_1998_spatial, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration moyenne en SPM dans les panaches de la baie des Anges en fonction du temps vu par le produit SEXTANT entre 1998 et 2025",
       x = "Date",
       y = "Aire du panache (en km²)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

save(SEXTANT_1998_2025_spm_spatial, file = "data/SEXTANT/SPM/SEXTANT_1998_2025_spm_spatial.Rdata")
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_spatial.Rdata")
