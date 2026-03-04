# TS analysis MODIS

# 03/03/2026

# pathway : "~/Documents/Alice/Code/MODIS_TS_analysis


# This script will load MODIS satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the 
# Var river at the Napoleon bridge


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# functions -----------------------------------------------------------------

## scaling function --------------------------------------------------------

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

# loading data ------------------------------------------------------------

load("~/Downloads/MODIS/SPM/MODIS_2015_2023_SPM.Rdata")
load("~/Downloads/MODIS/CHL/MODIS_2015_2023_CHL.Rdata")
load("~/Documents/Alice/Hydro France/Y6442010_Hydro.Rdata")

# plotting ----------------------------------------------------------------
## SPM ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS
ggplot(data = MODIS_2015_2023_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 2015 et 2023 avec le produit MODIS issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(MODIS_2015_2023_SPM$mean_spm, Y6442010_Hydro_complete$débit)

MODIS_2015_2023_SPM$scaled_mean_spm <- MODIS_2015_2023_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = MODIS_2015_2023_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension entre 2015 et 2023 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )




## CHL ----------------------------------------------------------------

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS
ggplot(data = MODIS_2015_2023_CHL, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "chartreuse3") +
  labs(title = "Evolution de la concentration en chlorophylle moyenne entre 2015 et 2023 avec le produit Sextant",
       x = "Date",
       y = "Concentration moyenne en chlorophylle (en µg/L)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec MODIS contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(MODIS_2015_2023_CHL$mean_chl, Y6442010_Hydro_complete$débit)

MODIS_2015_2023_CHL$scaled_mean_chl <- MODIS_2015_2023_CHL$mean_chl * adjust_factors$diff + adjust_factors$adjust


ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = MODIS_2015_2023_CHL,
    aes(x = date, y = scaled_mean_chl, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "chartreuse3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Chlorophylle (en µg/L)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en chlorophylle entre 2015 et 2023 avec le produit Sextant",
    x = "Date"
  ) +
  theme_minimal()
