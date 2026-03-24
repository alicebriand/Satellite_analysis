# func.R

# This script contains functions that are used across multiple different scripts
# It is meant to be loaded via source("func.R") by any script that needs these functions


# Loading functions -------------------------------------------------------

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


# Processing functions ----------------------------------------------------

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

