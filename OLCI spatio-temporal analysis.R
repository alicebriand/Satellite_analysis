# OLCI spatio-temporal analysis
# 27/02/2026

# pathway : "~/Downloads/Satellite analysis/OLCI spatio-temporal analysis.Rdata


# This script will load OLCI satellite data
# Then perform a spatial analysis


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# Get satellite download function

source("~/Downloads/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)


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


## loading function -----------------------------------------------------------

file_name <- "~/Downloads/OLCI/CHL_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2016_CHL/2016/L3m_20160426__FRANCE_03_OLA_CHL-G-PO_DAY_00.nc"

basename(file_name)

# load_OLCI_spm <- function(file_name, lon_range, lat_range){ 
#   file_caracter <- substr(basename(file_name), start = 5, stop = 12)
#   file_date <- as.Date(file_caracter, format = "%Y%m%d")
#   
#   df_chl <- tidync(file_name) %>% 
#     hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
#                  lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
#     hyper_tibble() %>% 
#     mutate(lon = as.numeric(lon), 
#            lat = as.numeric(lat),
#            date = file_date) %>% 
#     summarise(
#       mean_spm = mean(`CHL-G-PO_mean`, na.rm = TRUE), 
#       min_spm = min(`CHL-G-PO_mean`, na.rm = TRUE),
#       max_spm = max(`CHL-G-PO_mean`, na.rm = TRUE),
#       std_spm = sd(`CHL-G-PO_mean`, na.rm = TRUE),
#       .by = "date"
#     )
# }

# load_OLCI_spm <- function(file_name, lon_range, lat_range){
#   
#   # Find the date
#   OLCI_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
#   
#   # The necessary code
#   OLCI_one <- tidync(file_name) |> 
#     hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
#                  lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
#     hyper_tibble() |> 
#     mutate(lon = as.numeric(lon),
#            lat = as.numeric(lat),
#            date = OLCI_one_date) |> 
#     dplyr::select(lon, lat, date, mask, `CHL-G-PO_mean`)
#   
#   # Exit
#   return(OLCI_one)
# }

load_OLCI_spm <- function(file_name, lon_range, lat_range) {
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  OLCI_one <- tidync(file_name) %>%
    hyper_filter(
      lon = lon >= lon_range[1] & lon <= lon_range[2],
      lat = lat >= lat_range[1] & lat <= lat_range[2]
    ) %>%
    hyper_tibble() %>%
    mutate(
      lon = as.numeric(lon),
      lat = as.numeric(lat),
      date = file_date
    ) %>%
    # Filtrer les pixels où le masque est valide (ex. mask == 1)
    filter(`SPM-G-PO_flags` == 1) %>%  # ou l2_flags == 0, selon le fichier
    select(lon, lat, date, `SPM-G-PO_mean`)
  return(OLCI_one)
}



load_OLCI_chl <- function(file_name, lon_range, lat_range) {
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  OLCI_one <- tidync(file_name) %>%
    hyper_filter(
      lon = lon >= lon_range[1] & lon <= lon_range[2],
      lat = lat >= lat_range[1] & lat <= lat_range[2]
    ) %>%
    hyper_tibble() %>%
    mutate(
      lon = as.numeric(lon),
      lat = as.numeric(lat),
      date = file_date
    ) %>%
    # Filtrer les pixels où le masque est valide (ex. mask == 1)
    filter(`CHL-G-PO_flags` == 1) %>%  # ou l2_flags == 0, selon le fichier
    select(lon, lat, date, `CHL-G-PO_mean`)
  return(OLCI_one)
}


# load data ---------------------------------------------------------------
### Hydro France data ---------------------------------------------------------------------

load("~/Documents/Alice/Hydro France/Y6442010_Hydro.Rdata")

### SPM ---------------------------------------------------------------------

# for temporal analysis we have already prepared our data frames
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2016.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2017.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2018.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2019.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2020.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2021.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2022.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2023.Rdata")
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/all_spm_propre_OLCI_2024.Rdata")
# 
# load("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_SPM_2016_2024.Rdata")

# however, for spatial analysis we need the position of our pixels, thus we have to
# load netcfd data

# All of the OLCI files SPM from 2015 to 2025
# you can add former data that you downloaded
OLCI_2016_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2016_SPM/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2017_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2017_SPM/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2018_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2018_SPM/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2019_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2019_SPM/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2020_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2020_SPM/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2021_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2021_SPM/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2022_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2022_SPM/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2023_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2023_SPM/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2024_dir <- dir("~/Downloads/OLCI/SPM_sans_vrai_bounding_box/OLCI_A_B_ODATIS_MR_2024_SPM/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

# Load and combine

OLCI_2016_spm <- plyr::ldply(OLCI_2016_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2017_spm <- plyr::ldply(OLCI_2017_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2018_spm <- plyr::ldply(OLCI_2018_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2019_spm <- plyr::ldply(OLCI_2019_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2020_spm <- plyr::ldply(OLCI_2020_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2021_spm <- plyr::ldply(OLCI_2021_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2022_spm <- plyr::ldply(OLCI_2022_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2023_spm <- plyr::ldply(OLCI_2023_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2024_spm <- plyr::ldply(OLCI_2024_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Combine and save
OLCI_2016_2024_spm <- rbind(OLCI_2015_spm, OLCI_2016_spm, OLCI_2017_spm,
                               OLCI_2015_spm, OLCI_2019_spm, OLCI_2020_spm,
                               OLCI_2021_spm, OLCI_2022_spm, OLCI_2023_spm,
                               OLCI_2024_spm)

save(OLCI_2016_2024_spm, file = "OLCI_2016_2024_SPM_spatial.RData")

load("OLCI_2016_2024_SPM.RData")


### CHL ---------------------------------------------------------------------

# for temporal analysis we have already prepared our data frames
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2016.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2017.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2018.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2019.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2020.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2021.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2022.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2023.Rdata")
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/all_chl_propre_OLCI_2024.Rdata")

# for temporal analysis
load("~/Downloads/OLCI/CHL_sans_vrai_bounding_box/OLCI_CHL_2016_2024.Rdata")


# Spatial analysis --------------------------------------------------------

# Look at the monthly average SPM for one year (ex : 2016)
OLCI_2016_spm_monthly <- OLCI_2016_spm|> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date))
summarise(mean_spm = mean(`SPM-G-PO_mean`, na.rm = TRUE), .by = c("lon", "lat", "year"))

# Map of the 12 months
ggplot(data = OLCI_2017_spm_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  # annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon, ylim = lat) +
  facet_wrap(~year)
# facet_grid(year~month)













# Temporal analysis -------------------------------------------------------

## cleaning data -----------------------------------------------------------

### CHL ---------------------------------------------------------------------


# OLCI A and B take photos for the same day sometimes and the dates is thus duplicate
# we want only one date so we take the mean of the 2 days

OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024 %>%
  group_by(date) %>%  # Grouper par la colonne "date"
  summarise(
    mean_spm = mean(mean_spm, na.rm = TRUE),
    min_spm = mean(min_spm, na.rm = TRUE),   # Moyenne des min (ou vous pouvez garder le min global)
    max_spm = mean(max_spm, na.rm = TRUE),   # Moyenne des max (ou vous pouvez garder le max global)
    std_spm = mean(std_spm, na.rm = TRUE)    # Moyenne des écarts-types
  )

# data still present some outliers that we have to remove (when we look at Nasa 
# World View we don't see any panaches)

OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024_new %>% 
  dplyr::filter(mean_spm <= 20
  )

# Complete missing dates in the date range
OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024_new %>%
  complete(date = seq(min(date), max(date), by = "day"))


ggplot(data = OLCI_CHL_2016_2024_new, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = min_spm, ymax = max_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_line() +
  labs(title = "Evolution de la concentration en CHL moyenne entre 2016 et 2024 avec le satellite OLCI",
       x = "Date",
       y = "concentration moyenne en CHL") +
  theme_minimal()

# we can plot this series against runoff the Var river to see if OLCI satellite 
# have a good estimate of spm 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(OLCI_CHL_2016_2024_new$mean_spm, Y6442010_Hydro_complete$débit)

OLCI_CHL_2016_2024_new$scaled_mean_spm <- OLCI_CHL_2016_2024_new$mean_spm * adjust_factors$diff + adjust_factors$adjust

# plotting
ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = OLCI_CHL_2016_2024_new,
    aes(x = date, y = scaled_mean_spm, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "red")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "CHL (mg/L)")
  ) +
  labs(
    title = "Débit et concentration en CHL entre 2016 et 2024 (OLCI)",
    x = "Date"
  ) +
  theme_minimal()







### CHL ---------------------------------------------------------------------


# OLCI A and B take photos for the same day sometimes and the dates is thus duplicate
# we want only one date so we take the mean of the 2 days

OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024 %>%
  group_by(date) %>%  # Grouper par la colonne "date"
  summarise(
    mean_chl = mean(mean_chl, na.rm = TRUE),
    min_chl = mean(min_chl, na.rm = TRUE),   # Moyenne des min (ou vous pouvez garder le min global)
    max_chl = mean(max_chl, na.rm = TRUE),   # Moyenne des max (ou vous pouvez garder le max global)
    std_chl = mean(std_chl, na.rm = TRUE)    # Moyenne des écarts-types
  )

# we can plot the TS 

ggplot(data = OLCI_CHL_2016_2024_new, mapping = aes(x = date, y = mean_chl)) +
  geom_line()  # ou geom_line(), selon ce que tu veux afficher

# we observe a strong seasonal pattern around the beginning of each year

# data still present some outliers that we have to remove (when we look at Nasa 
# World View we don't see any panaches)

# OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024_new %>% 
#   dplyr::filter(mean_chl <= 20
#   )

# Complete missing dates in the date range
OLCI_CHL_2016_2024_new <- OLCI_CHL_2016_2024_new %>%
  complete(date = seq(min(date), max(date), by = "day"))


ggplot(data = OLCI_CHL_2016_2024_new, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = min_chl, ymax = max_chl,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_line() +
  labs(title = "Evolution de la concentration en CHL moyenne entre 2016 et 2024 avec le satellite OLCI",
       x = "Date",
       y = "concentration moyenne en CHL") +
  theme_minimal()

# we can plot this series against runoff the Var river to see if OLCI satellite 
# have a good estimate of chl 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(OLCI_CHL_2016_2024_new$mean_chl, Y6442010_Hydro_complete$débit)

OLCI_CHL_2016_2024_new$scaled_mean_chl <- OLCI_CHL_2016_2024_new$mean_chl * adjust_factors$diff + adjust_factors$adjust

# plotting
ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = OLCI_CHL_2016_2024_new,
    aes(x = date, y = scaled_mean_chl, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "red")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "CHL (mg/L)")
  ) +
  labs(
    title = "Débit et concentration en CHL entre 2016 et 2024 (OLCI)",
    x = "Date"
  ) +
  theme_minimal()










# Animation ---------------------------------------------------------------

# Once a brick of data is loaded, it is possible to walk through one day at a time

# Animate using gganimate to show the months as the time steps
p_map <- ggplot(data = OLCI_2015_2025_spm_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  # facet_wrap(~month) +
  labs(title = "OLCI CHL data from 2015 to 2025", x = "Longitude", y = "Latitude") +
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
        renderer = gifski_renderer(file = "animations/OLCI_1998.gif",
                                   width = 1200,    # Increase width in pixels
                                   height = 1000))   # Increase height in pixels









# Spatial analysis --------------------------------------------------------

# Look at the monthly average CHL for one year (ex : 2016)
OLCI_2016_spm_monthly <- OLCI_2016_spm %>% 
  mutate(year = lubridate::year(date),
          month = lubridate::month(date))

# Map of the 12 months
ggplot(data = OLCI_2016_spm_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  # annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon, ylim = lat) +
  facet_wrap(~year)
# facet_grid(year~month)



# Look at the monthly average CHL for one year (ex : 2017)
OLCI_2017_spm_monthly <- all_spm_propre_OLCI_2017 |> 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date))
summarise(mean_spm = mean(`CHL-G-PO_mean`, na.rm = TRUE), .by = c("lon", "lat", "year"))

# Map of the 12 months
ggplot(data = OLCI_2017_spm_monthly, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_spm)) +
  # annotation_borders(fill = "black", colour = "lightgreen") +
  coord_quickmap(xlim = lon, ylim = lat) +
  facet_wrap(~year)
# facet_grid(year~month)
