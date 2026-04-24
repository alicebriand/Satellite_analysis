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
library(giscoR) # Hi-res coastlines
library(luna) # Used to access NASA data
library(doParallel); registerDoParallel(cores = detectCores() - 2)

# Setup -------------------------------------------------------------------

# Load username and password
# This is a two column csv file with column names: usrname, psswrd
# It contains one row of data, which is the username and password for an account on earth data
# An account can be created here: https://urs.earthdata.nasa.gov/
# Once you have your account details, create the .csv file shown here and store in a secure location
earth_up <- read_csv("~/pCloudDrive/Stage/password_earth_data_access.csv")

# Functions ---------------------------------------------------------------

# Process MODIS HDF files in a batch
# Doesn't work well with L1 products or very large study areas
proc_MODIS_hdf <- function(file_name_df, bbox, layer_num, out_dir, land_mask = FALSE){
  
  # Get date from filenames
  files_date <- luna::modisDate(file_name_df$filename)
  
  # Check that only one date of data is being processed
  if(length(unique(files_date$date)) > 1 ) stop("Please ensure only one date of data is being processed.")
  
  # get product ID from file_names
  file_product <- strsplit(file_name_df$filename, "\\/")
  file_product <- sapply(file_product, "[[", length(file_product[[1]]))
  file_product <- unique(sapply(strsplit(file_product, "\\."), "[[", 1))
  
  # Check that only one product type is being processed
  if(length(file_product) > 1 ) stop("Please ensure only one product ID is being processed.")
  
  # Check if file already exists and skip if so
  file_name_out <- file.path(out_dir, paste0("study_area_",file_product,"_",unique(files_date$date),".tif"))
  if(file.exists(file_name_out)){
    message(file_name_out, " already exists. Delete it if you want to reprocess the data. Otherwise, all good.")
    return(NULL)
  }
  
  # Load  and merge the files with desired layers etc.
  if(grepl("MYD01|MYD02", file_product)){
    
    stop("MODIS L1 data are not currently working with this processing workflow. Rather use a different software.")
    
  } else {
    
    # TODO: Get this to detect if it should use 'subds' or 'lyrs'
    data_layer <- rast(file_name_df$filename[1], lyrs = layer_num)
    # data_layer
    # plot(data_layer)
    
    # Project to EPSG:4326 and crop
    # data_rectify <- rectify(data_layers[[1]])
    data_proj <- project(data_layer, y = "EPSG:4326")
    # plot(data_proj)
    data_crop <- crop(data_proj, bbox)
    # plot(data_crop)
    
    if(length(file_name_df$filename) > 1){
      # Run and merge each individual file
      data_step <- 2
      while(data_step <= length(file_name_df$filename)){
        data_proj_i <- project(rast(file_name_df$filename[1], lyrs = layer_num), y = data_proj)
        data_crop_i <- crop(data_proj_i, bbox)
        data_crop <- terra::merge(data_crop, data_crop_i)
        data_step <- data_step+1
        # plot(data_base)
      }
      rm(data_base_i, data_proj_i, data_crop_i); gc()
    }
    data_base <- data_crop
  }
  
  # Remove unneeded mask layers if compiling the water/land mask
  if(land_mask){
    data_base <- terra::ifel(data_base %in% c(1, 2, 3, 4, 5), NA, data_base)
    # plot(data_base)
  }
  
  # Save and quit
  writeRaster(data_base, file_name_out, overwrite = TRUE)
}

# This function helps us to load the .tif files as data.frames more easily
load_MODIS_tif <- function(file_name, mask_rast){
  
  # Get date from file_name
  file_date <- as.Date(str_split(gsub("\\.tif", "", basename(file_name)), "_")[[1]][4])
  
  # Load all files
  study_area_rast <- terra::rast(file_name)
  
  # Project the mask to the same grid as raster data
  mask_proj <- project(mask_rast, study_area_rast)
  
  # Mask the raster data
  study_area_water <- mask(study_area_rast, mask_proj)
  
  # Convert to data.frame
  study_area_df <- as.data.frame(study_area_water, xy = TRUE, na.rm = TRUE) |> 
    dplyr::rename(lon = x, lat = y) |> 
    mutate(date = file_date, .before = lon)
  
  # Exit
  return(study_area_df)
}

# Fonction pour calculer le % de pixels valides (non nuages)
check_cloud_cover <- function(file_hdf, seuil_nuages = 0.5) {
  
  r <- rast(file_hdf, lyrs = 1)  # bande 1
  vals <- values(r)
  
  # Pixels valides = réflectance entre 0 et 0.6
  pct_valide <- sum(vals > 0 & vals <= 0.6, na.rm = TRUE) / sum(!is.na(vals))
  
  return(pct_valide >= seuil_nuages)  # TRUE = image utilisable
}

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
productInfo("MYD09GHK") #L2G 500 m
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
dl_dir <- "~/Downloads/MODIS NASA/L2/2024/02/11/"

# Chosen start and end dates for downloading
start_date <- "2024-02-11"; end_date <- "2024-02-11"

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
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2/2024/02/11/MOD44W.A2024001.h18v04.061.2025064072734.hdf")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)
maps::map(add = TRUE)

# Prep one day of MODIS data
# NB: This requires that this folder exists: ~/data/MODIS
# IF not, create it or change the directories below to match 
plyr::d_ply(.data = rast_files, .variables = c("date"), .fun = proc_MODIS_hdf, .parallel = FALSE,
            bbox = study_bbox, out_dir = dl_dir, layer_num = 2, land_mask = FALSE)

# Load a file
MODIS_rast_b1 <- rast("~/Downloads/MODIS NASA/L2/2024/02/11/MYD09GQ.A2024042.h18v04.061.2024044055403.hdf")

# Check that it looks correct
plot(MODIS_rast_b1)
maps::map(add = TRUE)

# Project the 250 m mask to the same grid as the 500 m raster data
MODIS_mask_proj <- project(MODIS_mask, MODIS_rast_b1)

# Check that it worked
plot(MODIS_mask_proj)
maps::map(add = TRUE)

# Mask the raster data
MODIS_water <- mask(MODIS_rast_b1, MODIS_mask_proj)

# Check to see if it looks correct - should show white where land is
plot(MODIS_water)
maps::map(add = TRUE)

## 4) Load data ------------------------------------------------------------

# Load the MODIS mask first
# Change the filename if this is not correct
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2/2024/02/11/study_area_MOD44W_2024-01-01.tif")

# Filter out just the .tif files (i.e. not the HDF files)
tif_files <- list.files(path = dl_dir, pattern = "\\.tif$", full.names = TRUE)

# Ensure the correct product ID is being used
product_ID_files <- tif_files[grepl(product_ID, tif_files)]

# Select just a subset of files as desired
product_ID_files <- product_ID_files[1]

# Load all files
study_area_df <- map_dfr(product_ID_files, load_MODIS_tif, MODIS_mask)

## 5) Plot data ------------------------------------------------------------

# Hi-res Mediterranean country and coastline shapes
# coastline_giscoR <- gisco_get_coastallines(resolution = "01")
# countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# Map
pl_map <- study_area_df |> 
  filter(sur_refl_b01_1 <= 3000) |> 
  filter(date == "2024-02-11") |> 
  ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = sur_refl_b01_1)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) + # ← ici
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", fill = "Surface reflectance (620-670 nm)") +
  coord_sf(                                                    # ← remplace coord_quickmap
    xlim = range(study_area_df$lon), 
    ylim = range(study_area_df$lat), 
    expand = FALSE
  ) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top", 
        legend.box = "vertical",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# Save as desired
ggsave("~/Downloads/MODIS NASA/L2/2024/02/11/fig_MODIS_2024-02-11.png", pl_map, height = 9, width = 14)





# one year by one year ----------------------------------------------------

# Workflow ----------------------------------------------------------------

## 2024 --------------------------------------------------------------------

## 1) Setup ---------------------------------------------------------------

# Chose where you would like to save the files
dl_dir <- "~/Downloads/MODIS NASA/L2 2024 Aqua/"

# Chosen start and end dates for downloading
start_date <- "2024-03-04"; end_date <- "2024-03-04"

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

# If that looks reasonable, download them
# NB: If this doesn't work, then the product ID, even if it is listed, may not be findable by the luna package
luna::getNASA(product = product_ID, start_date = start_date, end_date = end_date, aoi = study_bbox,
              download = TRUE, overwrite = FALSE, server = product_server, version = product_version,
              path = dl_dir, username = earth_up$username, password = earth_up$password)

# To follow the rest of the examples below we also want to download the MODIS mask files
# luna::getNASA(product = "MOD44W", start_date = start_date, end_date = end_date, 
#               aoi = study_bbox, download = TRUE, overwrite = FALSE,
#               path = dl_dir, username = earth_up$username, password = earth_up$password)

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
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2024/MOD44W.A2024001.h18v04.061.2025064072734.hdf")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)
maps::map(add = TRUE)

## MODIS data --------------------------------------------------------------

# Lister tous les fichiers HDF téléchargés
all_hdf <- list.files("~/Downloads/MODIS NASA/L2 2024/",
                      pattern = "MYD09GQ\\.",
                      full.names = TRUE, recursive = TRUE)

# Convertir la liste all_hdf en data.frame avec les dates (comme rast_files)
all_hdf_df <- luna::modisDate(all_hdf)

# Ajouter une colonne mois
all_hdf_df$mois <- format(all_hdf_df$date, "%m")

# Traiter les données mois par mois

# Pour la bande 1 = layer 2
for (m in unique(all_hdf_df$mois)) {
  
  cat("Traitement bande 1 - mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  # Traiter date par date avec gestion des erreurs
  for (d in unique(hdf_mois$date)) {
    
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    
    # tryCatch permet de continuer même si un fichier plante
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE,
                  bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2024/tif_bande_1/",
                  layer_num = 2,
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
      
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  
  gc()
  cat("Mois", m, "terminé\n\n")
}

# Pour la bande 2 = layer 3
for (m in unique(all_hdf_df$mois)) {
  cat("Traitement bande 2 - mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  for (d in unique(hdf_mois$date)) {
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE, bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2024/tif_bande_2/",
                  layer_num = 3,  # ← à vérifier, probablement 1 pour b02
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  gc()
}

# Load the MODIS mask first
# Change the filename if this is not correct
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2024/study_area_MOD44W_2024-01-01.tif")

# Lister tous les tif produits
tif_files_b1 <- list.files("~/Downloads/MODIS NASA/L2 2024/tif_bande_1/", 
                        pattern = "MYD09GQ.*\\.tif$", 
                        full.names = TRUE)

tif_files_b2 <- list.files("~/Downloads/MODIS NASA/L2 2024/tif_bande_2/", 
                           pattern = "MYD09GQ.*\\.tif$", 
                           full.names = TRUE)

# Ensure the correct product ID is being used
product_ID_files_b1 <- tif_files_b1[grepl(product_ID, tif_files_b1)]
product_ID_files_b2 <- tif_files_b2[grepl(product_ID, tif_files_b2)]

# Load all files
study_area_df_2024_b1 <- map_dfr(product_ID_files_b1, load_MODIS_tif, MODIS_mask)
study_area_df_2024_b2 <- map_dfr(product_ID_files_b2, load_MODIS_tif, MODIS_mask)

study_area_df_2024 <- left_join(study_area_df_2024_b1, study_area_df_2024_b2, by = c("date", "lon", "lat"))

save(study_area_df_2024, file = "data/MODIS L2 NASA/study_area_df_2024")

# on va enlever manuellement les jours qui sont contaminés par des nuages (ça va être long)

dates_a_exclure <- c("2024-01-01", "2024-01-02", "2024-01-03", "2024-01-05", "2024-01-07",
                     "2024-01-08", "2024-01-09", "2024-01-10", "2024-01-13", "2024-01-14", 
                     "2024-01-17", "2024-01-19", "2024-01-21", "2024-01-22", "2024-01-24",
                     "2024-01-26", "2024-01-27", "2024-01-28", "2024-01-29", "2024-01-30",
                     "2024-02-01", "2024-02-02", "2024-02-03", "2024-02-01", "2024-02-05",
                     "2024-02-06", "2024-02-09", "2024-02-10", "2024-02-12", "2024-02-15",
                     "2024-02-16", "2024-02-18", "2024-02-22", "2024-02-23", "2024-02-24",
                     "2024-02-25", "2024-02-27", "2024-02-28", "2024-03-03", "2024-03-08",
                     "2024-03-09", "2024-03-11", "2024-03-15", "2024-03-17", "2024-03-22", 
                     "2024-03-25", "2024-03-26", "2024-03-28", "2024-03-29", "2024-03-30", 
                     "2024-03-31", "2024-04-04", "2024-04-07", "2024-04-08", "2024-04-09",
                     "2024-04-10", "2024-04-15", "2024-04-18", "2024-04-22", "2024-04-26",
                     "2024-04-27", "2024-04-28", "2024-04-29", "2024-04-30", "2024-05-01",
                     "2024-05-01", "2024-05-02", "2024-05-03", "2024-05-05", "2024-05-06",
                     "2024-05-12", "2024-05-14", "2024-05-15", "2024-05-20", "2024-05-23", 
                     "2024-05-24", "2024-05-27", "2024-05-29", "2024-05-30", "2024-06-02",
                     "2024-06-05", "2024-06-06", "2024-06-08", "2024-06-09", "2024-06-19", 
                     "2024-06-20", "2024-06-22", "2024-06-23", "2024-06-24", "2024-06-26", 
                     "2024-06-27", "2024-07-02", "2024-07-03", "2024-07-06", "2024-07-12", 
                     "2024-07-19", "2024-07-21", "2024-08-05", "2024-08-15", "2024-08-19",
                     "2024-08-27", "2024-08-31", "2024-09-04", "2024-09-08", "2024-09-17", 
                     "2024-09-18", "2024-09-19", "2024-09-22", "2024-09-24", "2024-09-26",
                     "2024-10-02", "2024-10-04", "2024-10-06", "2024-10-07", "2024-10-09",
                     "2024-10-12", "2024-10-14", "2024-10-16", "2024-10-17", "2024-10-19",
                     "2024-10-22", "2024-10-24", "2024-10-26", "2024-10-27", "2024-10-29",
                     "2024-11-06", "2024-11-09", "2024-11-12", "2024-11-18", "2024-11-20",
                     "2024-11-21", "2024-11-24", "2024-11-26", "2024-11-28", "2024-12-03",
                     "2024-12-05", "2024-12-06", "2024-12-07", "2024-12-09", "2024-12-10", 
                     "2024-12-11", "2024-12-13", "2024-12-03", "2024-12-18", "2024-12-19")

study_area_df_clean_2024 <- study_area_df_2024 |> 
  filter(!date %in% as.Date(dates_a_exclure)) |> 
  filter(sur_refl_b01_1 >= 0, sur_refl_b02_1 >= 0) # on a des valeurs de réflectance négatives, on les supprime

### plotting --------------------------------------------------------------------

# Créer un dossier pour les figures
# dir.create("~/Downloads/MODIS NASA/L2 2024/figures/bande 2/", recursive = TRUE, showWarnings = FALSE)

# Hi-res Mediterranean country and coastline shapes
# coastline_giscoR <- gisco_get_coastallines(resolution = "01")
# countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# ensuite on peut plotter la bande 1 ou 2
for (d in unique(study_area_df_clean_2024$date)) {
  
  df_jour <- study_area_df_clean_2024 |> 
    filter(date == d, sur_refl_b02_1 <= 3000)
  
  # Sauter si pas assez de pixels (image trop nuageuse)
  # if (nrow(df_jour) < 100) {
  #   cat("✗ Trop nuageux :", as.character(d), "\n")
  #   next
  # }
  
  pl <- df_jour |> 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = sur_refl_b02_1)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c() +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)", 
         fill = "Surface reflectance (841-876 nm)",
         title = as.character(d)) +
    coord_sf(xlim = range(study_area_df_clean_2024$lon), ylim = range(study_area_df_clean_2024$lat), expand = FALSE) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
  
  ggsave(paste0("~/Downloads/MODIS NASA/L2 2024/figures/bande 2/fig_MODIS_", d, ".png"), 
         pl, height = 9, width = 14)
  
  cat("✓ Figure sauvegardée :", as.character(d), "\n")
}





## 2016 --------------------------------------------------------------------

# Chose where you would like to save the files
dl_dir <- "~/Downloads/MODIS NASA/L2 2016 Terra/"

# Chosen start and end dates for downloading
start_date <- "2016-01-01"; end_date <- "2016-12-31"

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
# Attention ici je dowload TERRA data
product_ID <- "MOD09GQ"

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
# luna::getNASA("MOD09GQ", start_date, end_date, aoi = study_bbox, download = FALSE)

# If that looks reasonable, download them
# NB: If this doesn't work, then the product ID, even if it is listed, may not be findable by the luna package
# luna::getNASA(product = product_ID, start_date = start_date, end_date = end_date, aoi = study_bbox,
#               download = TRUE, overwrite = FALSE, server = product_server, version = product_version,
#               path = dl_dir, username = earth_up$username, password = earth_up$password)

# To follow the rest of the examples below we also want to download the MODIS mask files
# luna::getNASA(product = "MOD44W", start_date = start_date, end_date = end_date,
#               aoi = study_bbox, download = TRUE, overwrite = FALSE,
#               path = dl_dir, username = earth_up$username, password = earth_up$password)

## 3) Process files --------------------------------------------------------

# Set file pathways
# NB: Change the directory to where you saved the files if it was changed
# NB: Change the pattern in rast_files to match the product ID you used if it is different
mask_files <- luna::modisDate(list.files(path = dl_dir, pattern = "MOD44W\\.", full.names = TRUE))
rast_files <- luna::modisDate(list.files(path = dl_dir, pattern = "MOD09GQ\\.", full.names = TRUE))

# Chose specific files
mask_files <- mask_files[1,] # Change accordingly
rast_files <- rast_files[1,] # Change accordingly

# Process all of the water mask files
# IF not, create it or change the directories below as desired
plyr::d_ply(.data = mask_files, .variables = c("date"), .fun = proc_MODIS_hdf, .parallel = FALSE,
            bbox = study_bbox, out_dir = dl_dir, layer_num = 2, land_mask = TRUE)

# Load the desired mask file
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2016 Terra/MOD44W.A2016001.h18v04.061.2024007133840.hdf")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)
maps::map(add = TRUE)

## MODIS data --------------------------------------------------------------

# Lister tous les fichiers HDF téléchargés
all_hdf <- list.files("~/Downloads/MODIS NASA/L2 2016 Terra/",
                      pattern = "MOD09GQ\\.",
                      full.names = TRUE, recursive = TRUE)

# Convertir la liste all_hdf en data.frame avec les dates (comme rast_files)
all_hdf_df <- luna::modisDate(all_hdf)

# Ajouter une colonne mois
all_hdf_df$mois <- format(all_hdf_df$date, "%m")

# Traiter les données mois par mois

# Pour la bande 1 = layer 2
for (m in unique(all_hdf_df$mois)) {
  
  cat("Traitement du mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  # Traiter date par date avec gestion des erreurs
  for (d in unique(hdf_mois$date)) {
    
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    
    # tryCatch permet de continuer même si un fichier plante
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE,
                  bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2016 Terra/tif_bande_1/",
                  layer_num = 2,
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
      
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  
  gc()
  cat("Mois", m, "terminé\n\n")
}

# Pour la bande 2 = layer 3
for (m in unique(all_hdf_df$mois)) {
  
  cat("Traitement du mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  # Traiter date par date avec gestion des erreurs
  for (d in unique(hdf_mois$date)) {
    
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    
    # tryCatch permet de continuer même si un fichier plante
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE,
                  bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2016 Terra/tif_bande_2/",
                  layer_num = 3,
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
      
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  
  gc()
  cat("Mois", m, "terminé\n\n")
}

# Load the MODIS mask first
# Change the filename if this is not correct
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2016 Terra/MOD44W.A2016001.h18v04.061.2024007133840.hdf")

# Lister tous les tif produits
tif_files_b1 <- list.files("~/Downloads/MODIS NASA/L2 2016 Terra/tif_bande_1/", 
                        pattern = "MOD09GQ.*\\.tif$", 
                        full.names = TRUE)

tif_files_b2 <- list.files("~/Downloads/MODIS NASA/L2 2016 Terra/tif_bande_2/", 
                           pattern = "MOD09GQ.*\\.tif$", 
                           full.names = TRUE)

# Ensure the correct product ID is being used
product_ID_files_b1 <- tif_files_b1[grepl(product_ID, tif_files_b1)]
product_ID_files_b2 <- tif_files_b2[grepl(product_ID, tif_files_b2)]

# Load all files
study_area_df_2016_b1 <- map_dfr(product_ID_files_b1, load_MODIS_tif, MODIS_mask)
study_area_df_2016_b2 <- map_dfr(product_ID_files_b2, load_MODIS_tif, MODIS_mask)

## terra --------------------------------------------------------------------

# atention à faire que pour MODIS Terra car les couches ont été mal chargées
r <- rast("~/Downloads/MODIS NASA/L2 2016 Terra/MOD09GQ.A2016002.h18v04.061.2021339185604.hdf")
names(r)   # noms des couches
nlyr(r)    # nombre de couches

study_area_df_2016_b1 <- study_area_df_2016_b1 |> 
  rename(sur_refl_b01_1 = lyr1) |>
  select(date, lon, lat, sur_refl_b01_1)

study_area_df_2016_b2 <- study_area_df_2016_b2 |> 
  rename(sur_refl_b02_1 = lyr1) |>
  select(date, lon, lat, sur_refl_b02_1)

## all good ----------------------------------------------------------------

# Vérifier
# unique(study_area_df_2016_b1$date)  # toutes les dates chargées
# nrow(study_area_df_2016_b1)         # nombre total de pixels
# 
# unique(study_area_df_2016_b2$date)  # toutes les dates chargées
# nrow(study_area_df_2016_b2)         # nombre total de pixels

study_area_df_2016 <- left_join(study_area_df_2016_b1, study_area_df_2016_b2, by = c("date", "lon", "lat"))

# save(study_area_df_2016, file ="data/MODIS L2 NASA/study_area_df_2016.Rdata")

# on va enlever manuellement les jours qui sont contaminés par des nuages (ça va être long)

# ATTENTION CE SONT LEZ DATES POUR TERRA (REFAIRE SI TU VEUX AQUA)
dates_a_exclure <- c("2016-01-01", "2016-01-02", "2016-01-03", "2016-01-04", "2016-01-05",
                     "2016-01-07", "2016-01-08", "2016-01-09", "2016-01-10", "2016-01-11", 
                     "2016-01-14", "2016-01-15", "2016-01-18", "2016-01-19", "2016-01-20",
                     "2016-01-21", "2016-01-22", "2016-01-23", "2016-01-24", "2016-01-25", 
                     "2016-01-26", "2016-01-27", "2016-01-28", "2016-01-29", "2016-01-31",
                     "2016-02-02", "2016-02-03", "2016-02-06", "2016-02-07", "2016-02-08", 
                     "2016-02-09", "2016-02-12", "2016-02-13", "2016-02-14", "2016-02-16",
                     "2016-02-17", "2016-02-18", "2016-02-19", "2016-02-20", "2016-02-21", 
                     "2016-02-22", "2016-02-23", "2016-02-24", "2016-02-25", "2016-02-26", 
                     "2016-02-27", "2016-02-28", "2016-02-29", "2016-03-04", "2016-03-05", 
                     "2016-03-08", "2016-03-10", "2016-03-11", "2016-03-13", "2016-03-14", 
                     "2016-03-16", "2016-03-17", "2016-03-20", "2016-03-23", "2016-03-25", 
                     "2016-03-27", "2016-03-29", "2016-03-30", "2016-03-31", "2016-04-01", 
                     "2016-04-02", "2016-04-03", "2016-04-04", "2016-04-05", "2016-04-06", 
                     "2016-04-09", "2016-04-15", "2016-04-17", "2016-04-20", "2016-04-21", 
                     "2016-04-23", "2016-04-24", "2016-04-28", "2016-04-30", "2016-05-01", 
                     "2016-05-02", "2016-05-07", "2016-05-08", "2016-05-09", "2016-05-10", 
                     "2016-05-11", "2016-05-13", "2016-05-18", "2016-05-22", "2016-05-25", 
                     "2016-05-26", "2016-05-27", "2016-05-28", "2016-05-29", "2016-06-01", 
                     "2016-06-02", "2016-06-03", "2016-06-04", "2016-06-11", "2016-06-12", 
                     "2016-06-14", "2016-06-15", "2016-06-16", "2016-06-19", "2016-06-21", 
                     "2016-06-30", "2016-07-08", "2016-07-12", "2016-07-14", "2016-07-23", 
                     "2016-07-24", "2016-07-30", "2016-08-09", "2016-08-11", "2016-08-16", 
                     "2016-08-17", "2016-08-18", "2016-08-22", "2016-08-30", "2016-09-03",
                     "2016-09-04", "2016-09-14", "2016-09-16", "2016-09-17", "2016-09-21", 
                     "2016-09-23", "2016-09-24", "2016-09-27", "2016-10-01", "2016-10-06", 
                     "2016-10-09", "2016-10-12", "2016-10-13", "2016-10-14", "2016-10-16", 
                     "2016-10-17", "2016-10-23", "2016-10-24", "2016-10-25", "2016-10-26", 
                     "2016-10-28", "2016-11-04", "2016-11-05", "2016-11-06", "2016-11-09", 
                     "2016-11-10", "2016-11-12", "2016-11-13", "2016-11-15", "2016-11-16", 
                     "2016-11-17", "2016-11-18", "2016-11-19", "2016-11-20", "2016-11-21", 
                     "2016-11-22", "2016-11-23", "2016-11-24", "2016-11-26", "2016-12-03", 
                     "2016-12-04", "2016-12-05", "2016-12-10", "2016-12-12", "2016-12-15", 
                     "2016-12-19", "2016-12-20", "2016-12-21", "2016-12-24", "2016-12-25", 
                     "2016-12-26", "2016-12-30")

study_area_df_clean_2016 <- study_area_df_2016 |> 
  filter(!date %in% as.Date(dates_a_exclure)) |> 
  filter(sur_refl_b01_1 >= 0, sur_refl_b02_1 >= 0) # on a des valeurs de réflectance négatives, on les supprime

### plotting --------------------------------------------------------------------

# Créer un dossier pour les figures
dir.create("~/Downloads/MODIS NASA/L2 2016 Terra/figures/bande 1/", recursive = TRUE, showWarnings = FALSE)
dir.create("~/Downloads/MODIS NASA/L2 2016 Terra/figures/bande 2/", recursive = TRUE, showWarnings = FALSE)

# Hi-res Mediterranean country and coastline shapes
# coastline_giscoR <- gisco_get_coastallines(resolution = "01")
# countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# ensuite on peut plotter la bande 1 ou 2
for (d in unique(study_area_df_clean_2016$date)) {
  
  df_jour <- study_area_df_clean_2016 |> 
    filter(date == d, sur_refl_b02_1 <= 3000)
  
  # Sauter si pas assez de pixels (image trop nuageuse)
  # if (nrow(df_jour) < 100) {
  #   cat("✗ Trop nuageux :", as.character(d), "\n")
  #   next
  # }
  
  pl <- df_jour |> 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = sur_refl_b02_1)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c() +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)", 
         fill = "Surface reflectance (841-876 nm)",
         title = as.character(d)) +
    coord_sf(xlim = range(study_area_df_clean_2016$lon), ylim = range(study_area_df_clean_2016$lat), expand = FALSE) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
  
  ggsave(paste0("~/Downloads/MODIS NASA/L2 2016 Terra/figures/bande 2/fig_MODIS_", d, ".png"), 
         pl, height = 9, width = 14)
  
  cat("✓ Figure sauvegardée :", as.character(d), "\n")
}

# ensuite on peut plotter la bande 1 ou 2
for (d in unique(study_area_df_clean_2016$date)) {
  
  df_jour <- study_area_df_clean_2016 |> 
    filter(date == d, sur_refl_b01_1 <= 3000)
  
  # Sauter si pas assez de pixels (image trop nuageuse)
  # if (nrow(df_jour) < 100) {
  #   cat("✗ Trop nuageux :", as.character(d), "\n")
  #   next
  # }
  
  pl <- df_jour |> 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = sur_refl_b01_1)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c() +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)", 
         fill = "Surface reflectance (620-670 nm)",
         title = as.character(d)) +
    coord_sf(xlim = range(study_area_df_clean_2016$lon), ylim = range(study_area_df_clean_2016$lat), expand = FALSE) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
  
  ggsave(paste0("~/Downloads/MODIS NASA/L2 2016 Terra/figures/bande 1/fig_MODIS_", d, ".png"), 
         pl, height = 9, width = 14)
  
  cat("✓ Figure sauvegardée :", as.character(d), "\n")
}


## 2017 --------------------------------------------------------------------

# Chose where you would like to save the files
dl_dir <- "~/Downloads/MODIS NASA/L2 2017 Terra/"

# Chosen start and end dates for downloading
start_date <- "2017-01-01"; end_date <- "2017-12-31"

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
# Attention ici je dowload TERRA data
product_ID <- "MOD09GQ"

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
luna::getNASA("MOD09GQ", start_date, end_date, aoi = study_bbox, download = FALSE)

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
rast_files <- luna::modisDate(list.files(path = dl_dir, pattern = "MOD09GQ\\.", full.names = TRUE))

# Chose specific files
mask_files <- mask_files[1,] # Change accordingly
rast_files <- rast_files[1,] # Change accordingly

# Process all of the water mask files
# IF not, create it or change the directories below as desired
plyr::d_ply(.data = mask_files, .variables = c("date"), .fun = proc_MODIS_hdf, .parallel = FALSE,
            bbox = study_bbox, out_dir = dl_dir, layer_num = 2, land_mask = TRUE)

# Load the desired mask file
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2017 Terra/MOD44W.A2017001.h18v04.061.2024007170505.hdf")

# Check that it looks correct - should show white where land would be
plot(MODIS_mask)
maps::map(add = TRUE)

## MODIS data --------------------------------------------------------------

# Lister tous les fichiers HDF téléchargés
all_hdf <- list.files("~/Downloads/MODIS NASA/L2 2017 Terra/",
                      pattern = "MOD09GQ\\.",
                      full.names = TRUE, recursive = TRUE)

# Convertir la liste all_hdf en data.frame avec les dates (comme rast_files)
all_hdf_df <- luna::modisDate(all_hdf)

# Ajouter une colonne mois
all_hdf_df$mois <- format(all_hdf_df$date, "%m")

# Traiter les données mois par mois

# Pour la bande 1 = layer 2
for (m in unique(all_hdf_df$mois)) {
  
  cat("Traitement du mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  # Traiter date par date avec gestion des erreurs
  for (d in unique(hdf_mois$date)) {
    
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    
    # tryCatch permet de continuer même si un fichier plante
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE,
                  bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2017 Terra/tif_bande_1/",
                  layer_num = 2,
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
      
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  
  gc()
  cat("Mois", m, "terminé\n\n")
}

# Pour la bande 2 = layer 3
for (m in unique(all_hdf_df$mois)) {
  
  cat("Traitement du mois :", m, "\n")
  hdf_mois <- all_hdf_df[all_hdf_df$mois == m, ]
  
  # Traiter date par date avec gestion des erreurs
  for (d in unique(hdf_mois$date)) {
    
    hdf_jour <- hdf_mois[hdf_mois$date == d, ]
    date_lisible <- as.Date(d, origin = "1970-01-01")
    
    # tryCatch permet de continuer même si un fichier plante
    tryCatch({
      plyr::d_ply(.data = hdf_jour, .variables = c("date"), .fun = proc_MODIS_hdf,
                  .parallel = FALSE,
                  bbox = study_bbox,
                  out_dir = "~/Downloads/MODIS NASA/L2 2017 Terra/tif_bande_2/",
                  layer_num = 3,
                  land_mask = FALSE)
      cat("  ✓", as.character(date_lisible), "\n")
      
    }, error = function(e) {
      cat("  ✗ ERREUR", as.character(date_lisible), ":", conditionMessage(e), "\n")
    })
  }
  
  gc()
  cat("Mois", m, "terminé\n\n")
}

# Load the MODIS mask first
# Change the filename if this is not correct
MODIS_mask <- rast("~/Downloads/MODIS NASA/L2 2017 Terra/MOD44W.A2017001.h18v04.061.2024007170505.hdf")

# Lister tous les tif produits
tif_files_b1 <- list.files("~/Downloads/MODIS NASA/L2 2017 Terra/tif_bande_1/", 
                           pattern = "MOD09GQ.*\\.tif$", 
                           full.names = TRUE)

tif_files_b2 <- list.files("~/Downloads/MODIS NASA/L2 2017 Terra/tif_bande_2/", 
                           pattern = "MOD09GQ.*\\.tif$", 
                           full.names = TRUE)

# Ensure the correct product ID is being used
product_ID_files_b1 <- tif_files_b1[grepl(product_ID, tif_files_b1)]
product_ID_files_b2 <- tif_files_b2[grepl(product_ID, tif_files_b2)]

# Load all files
study_area_df_2017_b1 <- map_dfr(product_ID_files_b1, load_MODIS_tif, MODIS_mask)
study_area_df_2017_b2 <- map_dfr(product_ID_files_b2, load_MODIS_tif, MODIS_mask)

## terra --------------------------------------------------------------------

# atention à faire que pour MODIS Terra car les couches ont été mal chargées
r <- rast("~/Downloads/MODIS NASA/L2 2017 Terra/MOD09GQ.A2017003.h18v04.061.2021363135727.hdf")
names(r)   # noms des couches
nlyr(r)    # nombre de couches

study_area_df_2017_b1 <- study_area_df_2017_b1 |> 
  rename(sur_refl_b01_1 = lyr1) |>
  select(date, lon, lat, sur_refl_b01_1)

study_area_df_2017_b2 <- study_area_df_2017_b2 |> 
  rename(sur_refl_b02_1 = lyr1) |>
  select(date, lon, lat, sur_refl_b02_1)

## all good ----------------------------------------------------------------

# Vérifier
# unique(study_area_df_2016_b1$date)  # toutes les dates chargées
# nrow(study_area_df_2016_b1)         # nombre total de pixels
# 
# unique(study_area_df_2016_b2$date)  # toutes les dates chargées
# nrow(study_area_df_2016_b2)         # nombre total de pixels

study_area_df_2017 <- left_join(study_area_df_2017_b1, study_area_df_2017_b2, by = c("date", "lon", "lat"))

# save(study_area_df_2017, file ="data/MODIS L2 NASA/study_area_df_2017.Rdata")

# on va enlever manuellement les jours qui sont contaminés par des nuages (ça va être long)

# ATTENTION CE SONT LEZ DATES POUR TERRA (REFAIRE SI TU VEUX AQUA)
dates_a_exclure <- c("2017-01-02", "2017-01-04", "2017-01-05", "2017-01-06",
                     "2017-01-08", "2017-01-10", "2017-01-11", "2017-01-12",
                     "2017-01-14", "2017-01-15", "2017-01-16", "2017-01-19",
                     "2017-01-21", "2017-01-22", "2017-01-23", "2017-01-25", 
                     "2017-01-26", "2017-01-27", "2017-01-28", "2017-01-30", 
                     "2017-01-31", "2017-02-01", "2017-02-02", "2017-02-03", 
                     "2017-02-04", "2017-02-05", "2017-02-07", "2017-02-08", 
                     "2017-02-09", "2017-02-10", "2017-02-11", "2017-02-12", 
                     "2017-02-13", "2017-02-14", "2017-02-17", "2017-02-22", 
                     "2017-02-24", "2017-02-28", "2017-03-01", "2017-03-04",
                     "2017-03-05", "2017-03-06", "2017-03-09", "2017-03-10", 
                     "2017-03-12", "2017-03-15", "2017-03-16", "2017-03-18", 
                     "2017-03-19", "2017-03-20", "2017-03-21", "2017-03-23", 
                     "2017-03-24", "2017-03-25", "2017-03-30", "2017-03-31", 
                     "2017-04-01", "2017-04-02", "2017-04-05", "2017-04-09", 
                     "2017-04-11", "2017-04-15", "2017-04-24", "2017-04-26", 
                     "2017-04-27", "2017-04-30", "2017-05-01", "2017-05-02", 
                     "2017-05-04", "2017-05-05", "2017-05-06", "2017-05-10", 
                     "2017-05-11", "2017-05-12", "2017-05-13", "2017-05-14", 
                     "2017-05-18", "2017-05-22", "2017-05-23", "2017-05-24", 
                     "2017-06-01", "2017-06-05", "2017-06-06", "2017-06-09",
                     "2017-06-14", "2017-06-15", "2017-06-23", "2017-06-24",
                     "2017-06-30", "2017-07-08", "2017-07-09", "2017-07-11",
                     "2017-07-20", "2017-07-21", "2017-07-24", "2017-07-27", 
                     "2017-08-06", "2017-08-10", "2017-08-13", "2017-08-16", 
                     "2017-08-17", "2017-08-19", "2017-08-24", "2017-08-28", 
                     "2017-08-29", "2017-08-30", "2017-08-31", "2017-09-02", 
                     "2017-09-06", "2017-09-09", "2017-09-10", "2017-09-14", 
                     "2017-09-15", "2017-09-17", "2017-09-18", "2017-09-19", 
                     "2017-09-22", "2017-09-23", "2017-09-24", "2017-09-25", 
                     "2017-09-26", "2017-09-29", "2017-09-30", "2017-10-01",
                     "2017-10-02", "2017-10-03", "2017-10-06", "2017-10-09", 
                     "2017-10-26", "2017-10-27", "2017-10-28", "2017-10-29", 
                     "2017-10-30", "2017-11-02", "2017-11-04", "2017-11-05", 
                     "2017-11-06", "2017-11-07", "2017-11-08", "2017-11-09", 
                     "2017-11-11", "2017-11-13", "2017-11-15", "2017-11-19",
                     "2017-11-22", "2017-11-23", "2017-11-24", "2017-11-25", 
                     "2017-12-02", "2017-12-06", "2017-12-07", "2017-12-08", 
                     "2017-12-10", "2017-12-11", "2017-12-12", "2017-12-14", 
                     "2017-12-15", "2017-12-18", "2017-12-20", "2017-12-24", 
                     "2017-12-25", "2017-12-26", "2017-12-27", "2017-12-29",
                     "2017-12-30", "2017-12-31")

study_area_df_clean_2017 <- study_area_df_2017 |> 
  filter(!date %in% as.Date(dates_a_exclure)) |> 
  filter(sur_refl_b01_1 >= 0, sur_refl_b02_1 >= 0) # on a des valeurs de réflectance négatives, on les supprime

### plotting --------------------------------------------------------------------

# Créer un dossier pour les figures
dir.create("~/Downloads/MODIS NASA/L2 2017 Terra/figures/bande 1/", recursive = TRUE, showWarnings = FALSE)
dir.create("~/Downloads/MODIS NASA/L2 2017 Terra/figures/bande 2/", recursive = TRUE, showWarnings = FALSE)

# Hi-res Mediterranean country and coastline shapes
# coastline_giscoR <- gisco_get_coastallines(resolution = "01")
# countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# ensuite on peut plotter la bande 2
for (d in unique(study_area_df_clean_2017$date)) {
  
  df_jour <- study_area_df_clean_2017 |> 
    filter(date == d, sur_refl_b02_1 <= 3000)
  
  # Sauter si pas assez de pixels (image trop nuageuse)
  # if (nrow(df_jour) < 100) {
  #   cat("✗ Trop nuageux :", as.character(d), "\n")
  #   next
  # }
  
  pl <- df_jour |> 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = sur_refl_b02_1)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c() +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)", 
         fill = "Surface reflectance (841-876 nm)",
         title = as.character(d)) +
    coord_sf(xlim = range(study_area_df_clean_2017$lon), ylim = range(study_area_df_clean_2017$lat), expand = FALSE) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
  
  ggsave(paste0("~/Downloads/MODIS NASA/L2 2017 Terra/figures/bande 2/fig_MODIS_", d, ".png"), 
         pl, height = 9, width = 14)
  
  cat("✓ Figure sauvegardée :", as.character(d), "\n")
}

# ensuite on peut plotter la bande 1
for (d in unique(study_area_df_clean_2017$date)) {
  
  df_jour <- study_area_df_clean_2017 |> 
    filter(date == d, sur_refl_b01_1 <= 3000)
  
  # Sauter si pas assez de pixels (image trop nuageuse)
  # if (nrow(df_jour) < 100) {
  #   cat("✗ Trop nuageux :", as.character(d), "\n")
  #   next
  # }
  
  pl <- df_jour |> 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = sur_refl_b01_1)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c() +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)", 
         fill = "Surface reflectance (620-670 nm)",
         title = as.character(d)) +
    coord_sf(xlim = range(study_area_df_clean_2017$lon), ylim = range(study_area_df_clean_2017$lat), expand = FALSE) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
  
  ggsave(paste0("~/Downloads/MODIS NASA/L2 2017 Terra/figures/bande 1/fig_MODIS_", d, ".png"), 
         pl, height = 9, width = 14)
  
  cat("✓ Figure sauvegardée :", as.character(d), "\n")
}












# create new repertories for each date ------------------------------------

# # Répertoire de base où créer les dossiers
repertoire_base <- "~/Downloads/MODIS NASA/L2 2024/"
# 
# # Extraire les dates uniques où débit > 200
# # dates_filtrees <-  %>%
# #   filter(debit > 200) %>%
# #   pull(date) %>%
# #   unique()
# 
# load("~/River_runoff_analysis/data/Hydro France/Var_crues.Rdata")
# 
# # on extrait la colonne date
# 
# Var_crues <- unique(as.character(Var_crues$date))
# 
# # Créer un dossier pour chaque date
# for (date in Var_crues) {
#   chemin_dossier <- file.path(repertoire_base, as.character(date))
#   if (!dir.exists(chemin_dossier)) {
#     dir.create(chemin_dossier, recursive = TRUE)
#   }
# }

# créer un répertoire pour toutes les dates de l'année 2024

# Définir l'année souhaitée
# annee <- 2024
# 
# dates <- seq(as.Date(paste0(annee, "-01-01")),
#              as.Date(paste0(annee, "-12-31")),
#              by = "day")
# 
# for (date in dates) {
#   date <- as.Date(date, origin = "1970-01-01")  # ← reconversion nécessaire
#   chemin <- file.path(
#     annee,
#     format(date, "%m"),
#     format(date, "%d")
#   )
#   dir.create(chemin, recursive = TRUE, showWarnings = FALSE)
# }
