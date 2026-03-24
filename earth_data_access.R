# earth_data_access.R

# This script provides the functions, workflow, and examples 
# to download any data product from the NASA earth data server:
# https://www.earthdata.nasa.gov/

# The user guide for MODIS products may be found here :
# https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf

# NB: While this works well for MODIS and Landsat data, it cannot access Sentinel


# Libraries ---------------------------------------------------------------

# Luna package is used to access large spatial data products
# NB: It is not available on r-universe not CRAN
if (!"luna" %in% installed.packages()) {
  install.packages('luna', repos = 'https://rspatial.r-universe.dev')
}

# Check for missing libraries and install them if necessary
# We will use th plyr package, but we will not load it explicitly
if (!all(c("maps", "tidyverse", "ncdf4", "terra", "doParallel", "plyr") %in% installed.packages())) {
  install.packages(c("maps", "tidyverse", "ncdf4", "terra", "doParallel", "plyr"), repos = "https://cloud.r-project.org/")
}

library(tidyverse)
library(ncdf4)
library(terra)
library(luna) # Used to access NASA data
library(doParallel); registerDoParallel(cores = detectCores() - 2)

# Load shared functions
source("func.R")

# Setup -------------------------------------------------------------------

# Load username and password
# This is a two column csv file with column names: usrname, psswrd
# It contains one row of data, which is the username and password for an account on earth data
# An account can be created here: https://urs.earthdata.nasa.gov/
# Once you have your account details, create the .csv file shown here and store in a secure location
earth_up <- read_csv("~/pCloudDrive/Stage/password_earth_data_access.csv")


# MODIS data --------------------------------------------------------------

# Lists all products that are currently searchable on Earth Data and related servers
# Or see user guide: https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf
earth_data_catalogue <- luna::getProducts()

# Filter out just the MODIS products
MODIS_catalogue <- earth_data_catalogue[grepl("MOD|MYD", earth_data_catalogue$short_name),]

# MODIS directory that lists all products
# https://nrt3.modaps.eosdis.nasa.gov/archive/allData/61

# Uncomment and run any of these lines to open the product page in your web browser
# MODIS/Aqua Surface Reflectance Daily L2G Global 250m SIN Grid V061
# productInfo("MYD09GQ")

# Level 1
# NB: L1 products are often very difficult to work with in R
# productInfo("MYD01")
# productInfo("MYD021KM")
# productInfo("MYD02HKM")
# productInfo("MYD02QKM")

# Level 2
# productInfo("MYD09GHK") #L2G 500 m
# "MYD09"

# Level 3
# productInfo("MYD09Q1" # MODIS/Aqua Surface Reflectance 8-Day L3 Global 250m SIN Grid V006
# moproductInfo("MYD09A1") # MODIS/Aqua Surface Reflectance 8-Day L3 Global 500m SIN Grid V006

# MODIS/Terra Land Water Mask Derived from MODIS and SRTM L3 Global 250m SIN Grid V061
# https://lpdaac.usgs.gov/documents/1915/MOD44W_User_Guide_ATBD_V61.pdf
# productInfo("MOD44W")

# Sentinel-3 products
# NB: While this does list products, accessing them isn't working well
S3_catalogue <- earth_data_catalogue[grepl("OLCI|Sentinel|SENTINEL", earth_data_catalogue$short_name, fixed = FALSE),]


# Workflow ----------------------------------------------------------------

## 1) Setup ---------------------------------------------------------------

# Chose where you would like to save the files
dl_dir <- "~/Downloads/MODIS NASA/L2/"

# Chosen start and end dates for downloading
start_date <- "2020-10-03"; end_date <- "2020-10-03"

# Determine what you want your bounding box to be
# NB: The processing functions will fail if too much data are loaded at once
# Here is the area surrounding the Bay of Angels

# choose which one you want according if you want to study the Paillon or the
# Var area

# Paillon
# study_coords <- matrix(c(
#   7.2100000, 43.6000000,  # Bottom-left corner
#   7.3600000, 43.6000000, # Bottom-right corner
#   7.3600000, 43.7300000, # Top-right corner
#   7.2100000, 43.7300000,  # Top-left corner
#   7.2100000, 43.6000000   # Close the polygon (same as first point)
# ), ncol = 2, byrow = TRUE)
# Var
study_coords <- matrix(c(
  6.8925000, 43.2136389,  # Bottom-left corner
  7.4200000, 43.2136389, # Bottom-right corner
  7.4200000, 43.7300000, # Top-right corner
  6.8925000, 43.7300000,  # Top-left corner
  6.8925000, 43.2136389   # Close the polygon (same as first point)
), ncol = 2, byrow = TRUE)

# Turn it into the necessary SpatVector object type
study_bbox <- vect(study_coords, crs = "EPSG:4326", type = "polygons")

# Print the object to verify it worked - should be four points that make a box
plot(study_coords)
maps::map(add = TRUE)

# Chose the product ID you want to download
product_ID <- "MYD09GQ"

# Look at the server and version info for the product of choice
earth_data_catalogue[earth_data_catalogue$short_name == product_ID,]

# Then choose the server and version number
# NB: 'LPCLOUD' is the preferred server, but is not always available
# NB: One should generally choose the newest version, i.e. the biggest number
product_server <- "LPCLOUD" # NB: Change this if not shown in the output shown above
product_version <- "061" # NB: Change this if not shown in the output shown above


## 2) Download files -------------------------------------------------------

# First have a peak at what files exist
# NB: If this throws an error, it may be necessary to manually change the download server
# Remove the server and version arguments and run again
# It should show all of the possible servers and versions from which the desired product can be downloaded
luna::getNASA("MYD09GQ", start_date, end_date, aoi = study_bbox, download = FALSE)

# If that looks reasonable, download them
# NB: If this doesn't work, then the product ID, even if it is listed, may not be findable by the luna package
luna::getNASA(product = product_ID, start_date = start_date, end_date = end_date, aoi = study_bbox, 
              download = TRUE, overwrite = FALSE, server = product_server, version = product_version,
              path = dl_dir, username = earth_up$username, password = earth_up$password)

# To follow the rest of the examples below we also want to download the MODIS mask files
luna::getNASA(product = "MOD44W", start_date = start_date, end_date = end_date, 
              aoi = study_bbox, download = TRUE, overwrite = FALSE,
              path = dl_dir, username = earth_up$username, password = earth_up$password)


## 3) Process files --------------------------------------------------------

# Set file pathways
# NB: Change the directory to where you saved the files if it was changed
# NB: Change the pattern in rast_files to match the product ID you used if it is different
mask_files <- luna::modisDate(list.files(path = dl_dir, pattern = "MOD44W\\.", full.names = TRUE))
rast_files <- luna::modisDate(list.files(path = dl_dir, pattern = "MYD09GQ\\.", full.names = TRUE))

# Chose specific files
mask_files <- mask_files[1,] # Change accordingly
rast_files <- rast_files[1,] # Change accordingly

# Process all of the water mask files
# IF not, create it or change the directories below as desired
plyr::d_ply(.data = mask_files, .variables = c("date"), .fun = proc_MODIS_hdf, .parallel = FALSE,
            bbox = study_bbox, out_dir = dl_dir, layer_num = 2, land_mask = TRUE)

# Load the desired mask file
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2/study_area_MOD44W_2020-01-01.tif")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)
maps::map(add = TRUE)

# Prep one day of MODIS data
# NB: This requires that this folder exists: ~/data/MODIS
# IF not, create it or change the directories below to match 
plyr::d_ply(.data = rast_files, .variables = c("date"), .fun = proc_MODIS_hdf, .parallel = FALSE,
            bbox = study_bbox, out_dir = "~/Downloads/MODIS NASA/L2/", layer_num = 2, land_mask = FALSE)

# Load a file
MODIS_rast <- rast("~/Downloads/MODIS NASA/L2/study_area_MYD09GQ_2020-10-03.tif")

# Check that it looks correct
plot(MODIS_rast)
maps::map(add = TRUE)

# Project the 250 m mask to the same grid as the 500 m raster data
MODIS_mask_proj <- project(MODIS_mask, MODIS_rast)

# Check that it worked
plot(MODIS_mask_proj)
maps::map(add = TRUE)

# Mask the raster data
MODIS_water <- mask(MODIS_rast, MODIS_mask_proj)

# Check to see if it looks correct - should show white where land is
plot(MODIS_water)
maps::map(add = TRUE)


## 4) Load data ------------------------------------------------------------

# Load the MODIS mask first
# Change the filename if this is not correct
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2/masks/study_area_MOD44W_2020-01-01.tif")

# Filter out just the .tif files (i.e. not the HDF files)
tif_files <- list.files(path = dl_dir, pattern = "\\.tif$", full.names = TRUE)

# Ensure the correct product ID is being used
product_ID_files <- tif_files[grepl(product_ID, tif_files)]

# Select just a subset of files as desired
product_ID_files <- product_ID_files[1]

# Load all files
study_area_df <- map_dfr(product_ID_files, load_MODIS_tif, MODIS_mask)

# One can then apply whatever algorithms one wants to this dataframe


## 5) Plot data ------------------------------------------------------------

# Map
pl_map <- study_area_df |> 
  # Remove pixels that are too high (i.e. clouds)
  filter(sur_refl_b01_1 <= 3000) |> 
  # Select one date
  filter(date == "2020-10-03") |> 
  # Round all surface reflectance values greater than 0.1 down to 0.1 for better plotting
  # mutate(Rrs = case_when(Rrs > 0.1 ~ 0.1, 
  #                        Rrs < 0 ~ 0, TRUE ~ Rrs)) |> 
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = sur_refl_b01_1)) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
  # NB: Change fill label to correctly indicate which band width was used
  labs(x = "Longitude (°E)", y = "Latitude (°N)", fill = "Surface reflectance (620-670 nm) ") +
  coord_quickmap(xlim = range(study_area_df$lon), ylim = range(study_area_df$lat)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top", 
        legend.box = "vertical",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2/fig_MODIS.png", pl_map, height = 9, width = 14)
