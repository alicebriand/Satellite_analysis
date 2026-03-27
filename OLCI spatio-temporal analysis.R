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

# problème cluster --------------------------------------------------------

# Fonction qui gère son propre cluster
load_year <- function(year_dirs, lon_range, lat_range) {
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("lon_range", "lat_range", "load_OLCI_spm_pixels"),
                envir = environment())
  clusterEvalQ(cl, { library(tidyverse); library(tidync) })
  
  result <- plyr::ldply(year_dirs, load_OLCI_spm_pixels,
                        .parallel = TRUE,
                        lon_range = lon_range,
                        lat_range = lat_range)
  stopCluster(cl)
  return(result)
}

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

# filename <- "~/Downloads/OLCI/SPM/2016/OLCI_A_ODATIS_MR_2016_SPM/L3m_20160426__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc"

# load_OLCI_spm <- function(file_name, lon_range, lat_range) {
#   file_caracter <- substr(basename(file_name), start = 5, stop = 12)
#   file_date <- as.Date(file_caracter, format = "%Y%m%d")
#   
#   OLCI_one <- tidync(file_name) %>%
#     hyper_filter(
#       lon = lon >= lon_range[1] & lon <= lon_range[2],
#       lat = lat >= lat_range[1] & lat <= lat_range[2]
#     ) %>%
#     hyper_tibble() %>%
#     mutate(
#       lon = as.numeric(lon),
#       lat = as.numeric(lat),
#       date = file_date
#     ) %>%
#     dplyr::select(lon, lat, date,`SPM-G-PO_mean`) |> 
#     filter(`SPM-G-PO_mean` >= 1.2) |> 
#     summarise(pixel_count = n(),
#               mean_spm = mean(`SPM-G-PO_mean`, na.rm = TRUE), .by = "date")
#   
#   return(OLCI_one)
# }

load_OLCI_spm_pixels <- function(file_name, lon_range, lat_range){
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  # The necessary code
  OLCI_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = file_date) |> 
    dplyr::select(lon, lat, date, `SPM-G-PO_mean`)  # tous les pixels, sans filtre
  
  # Exit
  return(OLCI_one)
}

# load data ---------------------------------------------------------------
## Hydro France data ---------------------------------------------------------------------

load("data/Hydro France/Y6442010_2016_2024.Rdata")

load("data/OLCI/SPM/OLCI_2016_2024_spm_95.Rdata")

## SPM ---------------------------------------------------------------------
### threshold of 1.2 --------------------------------------------------------

OLCI_2016_dir <- dir("~/Downloads/OLCI/SPM/2016/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2017_dir <- dir("~/Downloads/OLCI/SPM/2017/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2018_dir <- dir("~/Downloads/OLCI/SPM/2018/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2019_dir <- dir("~/Downloads/OLCI/SPM/2019/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2020_dir <- dir("~/Downloads/OLCI/SPM/2020/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2021_dir <- dir("~/Downloads/OLCI/SPM/2021/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2022_dir <- dir("~/Downloads/OLCI/SPM/2022/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2023_dir <- dir("~/Downloads/OLCI/SPM/2023/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_2024_dir <- dir("~/Downloads/OLCI/SPM/2024/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

all_dirs <- list(
  "2016" = OLCI_2016_dir,
  "2017" = OLCI_2017_dir,
  "2018" = OLCI_2018_dir,
  "2019" = OLCI_2019_dir,
  "2020" = OLCI_2020_dir,
  "2021" = OLCI_2021_dir,
  "2022" = OLCI_2022_dir,
  "2023" = OLCI_2023_dir,
  "2024" = OLCI_2024_dir
)

### threshold of 1.2 --------------------------------------------------------

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
OLCI_2016_2024_spm_spatial <- rbind(OLCI_2016_spm, OLCI_2017_spm, OLCI_2018_spm,
                                    OLCI_2019_spm, OLCI_2020_spm, OLCI_2021_spm,
                                    OLCI_2022_spm, OLCI_2023_spm, OLCI_2024_spm)

save(OLCI_2016_2024_spm_spatial, file = "data/OLCI/SPM/OLCI_2016_2024_spm_spatial.Rdata")

load("data/OLCI/SPM/OLCI_2016_2024_spm_spatial.Rdata")

### define threshold with percentile 95 --------------------------------------------------------

OLCI_2016_spm_pixels <- plyr::ldply(OLCI_2016_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2017_spm_pixels <- plyr::ldply(OLCI_2017_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2018_spm_pixels <- plyr::ldply(OLCI_2018_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2019_spm_pixels <- plyr::ldply(OLCI_2019_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2020_spm_pixels <- plyr::ldply(OLCI_2020_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2021_spm_pixels <- plyr::ldply(OLCI_2021_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2022_spm_pixels <- plyr::ldply(OLCI_2022_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2023_spm_pixels <- plyr::ldply(OLCI_2023_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_2024_spm_pixels <- plyr::ldply(OLCI_2024_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Tout charger en un seul objet
# OLCI_2016_2024_spm_pixels <- purrr::map_dfr(all_dirs, load_year,
#                                             lon_range = lon_range,
#                                             lat_range = lat_range)

OLCI_2016_2024_spm_pixels <- rbind(OLCI_2016_spm_pixels, OLCI_2017_spm_pixels, OLCI_2018_spm_pixels,
                                   OLCI_2019_spm_pixels, OLCI_2020_spm_pixels, OLCI_2021_spm_pixels, 
                                   OLCI_2022_spm_pixels, OLCI_2023_spm_pixels, OLCI_2024_spm_pixels)

save(OLCI_2016_2024_spm_pixels, file = "data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

load("data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

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





# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

nc <- tidync("~/Downloads/OLCI/SPM/2016/OLCI_A_ODATIS_MR_2016_SPM/L3m_20160426__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc")

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

## define 95ème percentile -------------------------------------------------

# Calculer le 95ème percentile
seuil_95 <- quantile(OLCI_2016_2024_spm_pixels$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "mg/m³\n")

# Seuil 95ème percentile : 0.4883142 mg/m³

# Stats du panache par jour
OLCI_2016_2024_spm_95 <- OLCI_2016_2024_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean` >= seuil_95, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(OLCI_2016_2024_spm_95, file = "data/OLCI/SPM/OLCI_2016_2024_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# mean spm or panache area plot

# en échelle normale

# model_OLCI_2016_95 <- lm(aire_panache_km2 ~ date, data = OLCI_2016_2024_spm_95)
model_OLCI_2016_95 <- lm(mean_spm ~ date, data = OLCI_2016_2024_spm_95)
p_value_OLCI_2016_95 <- summary(model_OLCI_2016_95)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_2016_95 <- coef(model_OLCI_2016_95)[1]
slope_OLCI_2016_95 <- coef(model_OLCI_2016_95)[2]

# ggplot(data = OLCI_2016_2024_spm_95, aes(x = date, y = aire_panache_km2)) +
#   geom_point(color = "darkcyan", size = 0.5) +
#   # geom_point(data = OLCI_2016_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   annotate(
#     "text",
#     x = max(OLCI_2016_2024_spm_95$date, na.rm = TRUE),
#     y = max(OLCI_2016_2024_spm_95$aire_panache_km2, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "y = ", round(intercept_OLCI_2016_95, 3), " + ", round(slope_OLCI_2016_95, 7), " * x",
#       "\n", "p = ", ifelse(p_value_OLCI_2016_95 < 0.001, "< 0.001", format(p_value_OLCI_2016_95, digits = 3))
#     ),
#     hjust = 1,  # Alignement à droite
#     vjust = 1,  # Alignement en haut
#     size = 6
#   ) +
#   labs(title = "Évolution de l'aire des panaches de la baie des Anges vu par le produit OLCI (ODATIS-MR)",
#        x = "Date",
#        y = "Aire du panache (km²)") +
#   theme_minimal() +
#   scale_x_date(
#     date_breaks = "1 year",  
#     date_labels = "%Y"       
#   )

ggplot(data = OLCI_2016_2024_spm_95, aes(x = date, y = mean_spm)) +
  geom_point(color = "red3", size = 0.5) +
  # geom_point(data = OLCI_2016_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(OLCI_2016_2024_spm_95$date, na.rm = TRUE),
    y = max(OLCI_2016_2024_spm_95$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_OLCI_2016_95, 3), " + ", round(slope_OLCI_2016_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_OLCI_2016_95 < 0.001, "< 0.001", format(p_value_OLCI_2016_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration en MES dans les panaches de la baie des Anges vu par le produit OLCI (ODATIS-MR)",
       x = "Date",
       y = "Concentration moyenne en MES (en mg/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# en échelle log

data_log_spm <- OLCI_2016_2024_spm_95 |> 
  filter(aire_panache_km2> 0)

model_OLCI_2016_95_log <- lm(log10(aire_panache_km2) ~ date, data = data_log_spm)
p_value_OLCI_2016_95_log <- summary(model_OLCI_2016_95_log)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_2016_95_log <- coef(model_OLCI_2016_95_log)[1]
slope_OLCI_2016_95_log <- coef(model_OLCI_2016_95_log)[2]

ggplot(data = data_log_spm, aes(x = date, y = aire_panache_km2)) +  # utilise data_log_spm
  geom_point(color = "darkcyan", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, 
              formula = y ~ x,
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$aire_panache_km2, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_OLCI_2016_95_log, 3), " + ", 
      round(slope_OLCI_2016_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_OLCI_2016_95_log < 0.001, "< 0.001", 
                           format(p_value_OLCI_2016_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de l'aire des panaches de la baie des Anges vu par le produit OLCI (ODATIS-Mr)(échelle log)",
       x = "Date",
       y = "Aire du panache (km²)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")


# comparison between liquid flow rate and panache extension

adjust_factors <- sec_axis_adjustement_factors(OLCI_2016_2024_spm_95$aire_panache_km2, Y6442010_2016_2024$débit)

OLCI_2016_2024_spm_95$scaled_aire_panache_km2<- OLCI_2016_2024_spm_95$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Y6442010_2016_2024, 
             aes(x = date, y = débit, color = "Débit"), size = 0.5) +
  geom_point(data = OLCI_2016_2024_spm_95, 
             aes(x = date, y = scaled_aire_panache_km2, color = "Aire des panaches"), size = 0.5) +
  scale_color_manual(values = c("Débit" = "blue", "Aire des panaches" = "limegreen")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(title = "Évolution de la taille des panaches et du débit du Var vu par le produit OLCI ODATIS-MR",
  # labs(title = "Évolution de la concentration en MES dans les panaches et du débit du Var vu par le produit OLCI ODATIS-MR",
            
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# scatter plot

# Fusionner les données
Var_OLCI <- Y6442010_2016_2024 %>% 
  select(date, débit) %>% 
  left_join(
    OLCI_2016_2024_spm_95 %>% select(date, aire_panache_km2, mean_spm),
    by = "date"
  )

ggplot(data = Var_OLCI, aes(x = débit, y = mean_spm)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 1) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x,
    parse = TRUE,
    colour = "red",
    size = 6,
    label.x = 0.05,  # position horizontale (0 = gauche, 1 = droite)
    label.y = 0.90   # position verticale (0 = bas, 1 = haut)
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis", name = "Nombre d'observations") +
  theme_bw() +
  labs(x = "Débit (m³/s)", y = "Concentration moyenne en MES (en mg/m³)", 
       title = "Débit liquide du Var contre la concentration moyenne en MES dans les panaches vue par OLCI (ODATIS -MR)") +
  theme_minimal()


