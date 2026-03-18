# TS analysis SEXTANT, MODIS, OLCI and more

# 03/03/2026

# pathway : "~/Downloads/satellite analysis/TS analysis SEXTANT, MODIS, OLCI and more


# This script will load Sextant, MODIS, LOCI and more satellite data
# Then compare the different product temporally

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# loading data ---------------------------------------------------------------

load("data/SEXTANT/sextant_2015_2016_SPM.RData")
# load("data/SEXTANT/sextant_2015_2025_CHL.Rdata")

load("data/MODIS/SPM/MODIS_2015_2023_SPM.Rdata")
# load("~/Downloads/MODIS/CHL/MODIS_2015_2023_CHL.Rdata")

load("data/OLCI/SPM/")
load("~/Downloads/OLCI/CHL/OLCI_A_B_2016_2024_CHL.Rdata")

# plotting ----------------------------------------------------------------
## SPM ----------------------------------------------------------------

ggplot() +
  geom_line(
    data = sextant_2015_2025_SPM,
    aes(x = date, y = mean_spm, color = "SEXTANT")
  ) +
  geom_line(
    data = MODIS_2015_2023_SPM,
    aes(x = date, y = mean_spm, color = "MODIS")
  ) +
  geom_line(
    data = OLCI_A_B_2016_2024_SPM,
    aes(x = date, y = mean_spm, color = "OLCI")
  ) +
  scale_color_manual(values = c("SEXTANT" = "deepskyblue", "MODIS" = "maroon1",  OLCI = "limegreen")) +
  # scale_y_continuous(
  #   name = "Débit (m³/s)",
  #   sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  # ) +
  labs(
    title = "Comparaison des concentrations en matière en suspension entre 2015 et 2025 avec les produit Sextant, MODIS et OLCI",
    x = "Date",
    y = " Concentration moyenne en matière particulaire en suspension"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# zoom sur une année, ex 2019

OLCI_A_B_2019_SPM <- OLCI_A_B_2016_2024_SPM %>%
  dplyr::filter(date >= as.Date("2019-01-01"), date <= as.Date("2019-12-31"))

MODIS_2019_SPM <- MODIS_2015_2023_SPM %>%
  dplyr::filter(date >= as.Date("2019-01-01"), date <= as.Date("2019-12-31"))

sextant_2019_SPM <- sextant_2015_2025_SPM %>%
  dplyr::filter(date >= as.Date("2019-01-01"), date <= as.Date("2019-12-31"))

ggplot() +
  geom_line(
    data = OLCI_A_B_2019_SPM,
    aes(x = date, y = mean_spm, color = "OLCI")
  ) +
  geom_line(
    data = MODIS_2019_SPM,
    aes(x = date, y = mean_spm, color = "MODIS")
  ) +
  geom_line(
    data = sextant_2019_SPM,
    aes(x = date, y = mean_spm, color = "SEXTANT")
  ) +
  scale_color_manual(values = c("SEXTANT" = "deepskyblue", "MODIS" = "maroon1",  OLCI = "limegreen")) +
  # scale_y_continuous(
  #   name = "Débit (m³/s)",
  #   sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  # ) +
  labs(
    title = "Comparaison des concentrations en matière en suspension en 2019 avec les produit Sextant, MODIS et OLCI",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

## CHL ----------------------------------------------------------------

ggplot() +
  geom_line(
    data = sextant_2015_2025_CHL,
    aes(x = Date, y = mean_chl, color = "SEXTANT")
  ) +
  geom_line(
    data = MODIS_2015_2023_CHL,
    aes(x = date, y = mean_chl, color = "MODIS")
  ) +
  geom_line(
    data = OLCI_A_B_2016_2024_CHL,
    aes(x = date, y = mean_chl, color = "OLCI")
  ) +
  scale_color_manual(values = c("SEXTANT" = "deepskyblue", "MODIS" = "maroon1", OLCI = "limegreen")) +
  # scale_y_continuous(
  #   name = "Débit (m³/s)",
  #   sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  # ) +
  labs(
    title = "Comparaison des concentrations en matière en suspension entre 2015 et 2025 avec les produit Sextant, MODIS et OLCI",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )