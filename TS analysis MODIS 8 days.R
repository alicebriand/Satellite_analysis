# TS analysis MODIS of 8 days products

# 03/03/2026

# pathway : "~/Downloads/satellite analysis/TS analysis MODIS 8 days


# This script will load Sextant satellite data on 8 days
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the Var river


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)
library(heatwaveR)

# load data ---------------------------------------------------------------

load("~/Downloads/MODIS/produit 8j/SPM/data/MODIS_2015_2024_SPM_8_j.Rdata")


# interpolation -----------------------------------------------------------



# plotting ----------------------------------------------------------------
## SPM ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS
ggplot(data = MODIS_2015_2024_SPM_8_j, aes(x = date, y = mean_spm)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 2015 et 2024 avec le produit MODIS issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )
