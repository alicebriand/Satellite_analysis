# MODIS L2 NASA analysis

# 16/03/2026
# pathway : "~/Satellite_analysis/MODIS L2 NASA analysis.R


# This script will load MODIS L2 data from NASA (first on the 03/10/2020 to
# see if everything is fine) and will try, thanks to the reflectance parameter,
# to retrieve SPM concentration of the study area

# Setup ------------------------------------------------------------------

# Bioconductor ------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install("rhdf5")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2)
library(heatwaveR)
library(ggpmisc)
library(rhdf5)

# load data ---------------------------------------------------------------
h5ls(fichier_hdf)

MODIS_03_10_2020.hdf <- "~/Downloads/MODIS L2 NASA/MYD09GQ.A2020277.h18v04.061.2020347193143.hdf"
datasets <- h5ls(MODIS_03_10_2020.hdf)

# function ----------------------------------------------------------------

filename <- "~/Downloads/MODIS L2 NASA/MYD09GQ.A2020277.h18v04.061.2020347193143.hdf"
basename(filename)

MODIS_L2_SPM <- function(filename){
  
}
