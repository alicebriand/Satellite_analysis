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
library(sf)
library(rnaturalearth)
library(ggpmisc)
library(data.table)
library(giscoR) # Hi-res coastlines
library(doParallel); registerDoParallel(cores = 14)

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)

# problème cluster --------------------------------------------------------

# Fonction qui gère son propre cluster
# load_year <- function(year_dirs, lon_range, lat_range) {
#   cl <- makeCluster(detectCores() - 2)
#   registerDoParallel(cl)
#   
#   clusterExport(cl, varlist = c("lon_range", "lat_range", "load_OLCI_spm_pixels"),
#                 envir = environment())
#   clusterEvalQ(cl, { library(tidyverse); library(tidync) })
#   
#   result <- plyr::ldply(year_dirs, load_OLCI_spm_pixels,
#                         .parallel = TRUE,
#                         lon_range = lon_range,
#                         lat_range = lat_range)
#   stopCluster(cl)
#   return(result)
# }

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

load("data/Hydro France/All_debit.Rdata")

All_debit_2016_2024 <- All_debit |> 
  filter(date >= as.Date("2014-01-01"), date <= as.Date("2024-12-31"))


load("data/OLCI/SPM/OLCI_2016_2024_spm_95.Rdata")

load("data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

load("data/OLCI/SPM/all_spm_propre_OLCI_B_2020.Rdata")

load("data/OLCI/CHL/OLCI_A_2016_2024_CHL.Rdata")

## SPM ---------------------------------------------------------------------
### threshold of 1.2 --------------------------------------------------------

# on commence à 2018 car c'est là qu'on a les deux capteurs

OLCI_A_2016_dir <- dir("~/Downloads/OLCI/SPM/2016/OLCI_A_ODATIS_MR_2016_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2017_dir <- dir("~/Downloads/OLCI/SPM/2017/OLCI_A_ODATIS_MR_2017_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2018_dir <- dir("~/Downloads/OLCI/SPM/2018/OLCI_A_ODATIS_MR_2018_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2019_dir <- dir("~/Downloads/OLCI/SPM/2019/OLCI_A_ODATIS_MR_2019_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2020_dir <- dir("~/Downloads/OLCI/SPM/2020/OLCI_A_ODATIS_MR_2020_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2021_dir <- dir("~/Downloads/OLCI/SPM/2021/OLCI_A_ODATIS_MR_2021_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2022_dir <- dir("~/Downloads/OLCI/SPM/2022/OLCI_A_ODATIS_MR_2022_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2023_dir <- dir("~/Downloads/OLCI/SPM/2023/OLCI_A_ODATIS_MR_2023_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_A_2024_dir <- dir("~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

OLCI_B_2018_dir <- dir("~/Downloads/OLCI/SPM/2018/OLCI_B_ODATIS_MR_2018_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2019_dir <- dir("~/Downloads/OLCI/SPM/2019/OLCI_B_ODATIS_MR_2019_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2020_dir <- dir("~/Downloads/OLCI/SPM/2020/OLCI_B_ODATIS_MR_2020_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2021_dir <- dir("~/Downloads/OLCI/SPM/2021/OLCI_B_ODATIS_MR_2021_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2022_dir <- dir("~/Downloads/OLCI/SPM/2022/OLCI_B_ODATIS_MR_2022_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2023_dir <- dir("~/Downloads/OLCI/SPM/2023/OLCI_B_ODATIS_MR_2023_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
OLCI_B_2024_dir <- dir("~/Downloads/OLCI/SPM/2024/OLCI_B_ODATIS_MR_2024_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### threshold of 1.2 --------------------------------------------------------

# OLCI_2016_spm <- plyr::ldply(OLCI_2016_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2017_spm <- plyr::ldply(OLCI_2017_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2018_spm <- plyr::ldply(OLCI_2018_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2019_spm <- plyr::ldply(OLCI_2019_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2020_spm <- plyr::ldply(OLCI_2020_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2021_spm <- plyr::ldply(OLCI_2021_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2022_spm <- plyr::ldply(OLCI_2022_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2023_spm <- plyr::ldply(OLCI_2023_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# OLCI_2024_spm <- plyr::ldply(OLCI_2024_dir, load_OLCI_spm, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
# 
# # Combine and save
# OLCI_2016_2024_spm_spatial <- rbind(OLCI_2016_spm, OLCI_2017_spm, OLCI_2018_spm,
#                                     OLCI_2019_spm, OLCI_2020_spm, OLCI_2021_spm,
#                                     OLCI_2022_spm, OLCI_2023_spm, OLCI_2024_spm)
# 
# save(OLCI_2016_2024_spm_spatial, file = "data/OLCI/SPM/OLCI_2016_2024_spm_spatial.Rdata")
# 
# load("data/OLCI/SPM/OLCI_2016_2024_spm_spatial.Rdata")

### define threshold with percentile 95 --------------------------------------------------------

OLCI_A_2016_spm_pixels <- plyr::ldply(OLCI_A_2016_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2017_spm_pixels <- plyr::ldply(OLCI_A_2017_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2018_spm_pixels <- plyr::ldply(OLCI_A_2018_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2019_spm_pixels <- plyr::ldply(OLCI_A_2019_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2020_spm_pixels <- plyr::ldply(OLCI_A_2020_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2021_spm_pixels <- plyr::ldply(OLCI_A_2021_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2022_spm_pixels <- plyr::ldply(OLCI_A_2022_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2023_spm_pixels <- plyr::ldply(OLCI_A_2023_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_A_2024_spm_pixels <- plyr::ldply(OLCI_A_2024_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

OLCI_B_2018_spm_pixels <- plyr::ldply(OLCI_B_2018_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2019_spm_pixels <- plyr::ldply(OLCI_B_2019_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2020_spm_pixels <- plyr::ldply(OLCI_B_2020_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2021_spm_pixels <- plyr::ldply(OLCI_B_2021_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2022_spm_pixels <- plyr::ldply(OLCI_B_2022_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2023_spm_pixels <- plyr::ldply(OLCI_B_2023_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
OLCI_B_2024_spm_pixels <- plyr::ldply(OLCI_B_2024_dir, load_OLCI_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

OLCI_A_2016_2024_spm_pixels <- rbind(OLCI_A_2016_spm_pixels, OLCI_A_2017_spm_pixels, OLCI_A_2018_spm_pixels,
                                     OLCI_A_2019_spm_pixels, OLCI_A_2020_spm_pixels, OLCI_A_2021_spm_pixels, 
                                     OLCI_A_2022_spm_pixels, OLCI_A_2023_spm_pixels, OLCI_A_2024_spm_pixels)

OLCI_B_2018_2024_spm_pixels <- rbind(OLCI_B_2018_spm_pixels,
                                     OLCI_B_2019_spm_pixels, OLCI_B_2020_spm_pixels, OLCI_B_2021_spm_pixels, 
                                     OLCI_B_2022_spm_pixels, OLCI_B_2023_spm_pixels, OLCI_B_2024_spm_pixels)

# créer un df avec :
# quand 1 seul passage : on garde tous les pixels, qu'importe s'ils sont remplis ou pas
# quand deux passages (OLCI A et B): faire la moyenne des valeurs des pixels (ou alors
# garder la valeur maximale calculée ?)
# quand pas de passage, on garde quand même les pixels vides

# Combiner les deux satellites
OLCI_2016_2024_spm_pixels <- rbind(
  OLCI_A_2016_2024_spm_pixels,
  OLCI_B_2018_2024_spm_pixels
)

fusionner_annee <- function(annee) {
  
  a <- get(paste0("OLCI_A_", annee, "_spm_pixels"))
  
  if (annee >= 2018) {
    b <- get(paste0("OLCI_B_", annee, "_spm_pixels"))
    df_raw <- rbind(a, b)
    rm(a, b)
  } else {
    df_raw <- a
    rm(a)
  }
  gc()
  
  df_fusionne <- df_raw |>
    group_by(lon, lat, date) |>
    summarise(
      `SPM-G-PO_mean` = case_when(
        sum(!is.na(`SPM-G-PO_mean`)) == 2 ~ mean(`SPM-G-PO_mean`, na.rm = TRUE),
        sum(!is.na(`SPM-G-PO_mean`)) == 1 ~ sum(`SPM-G-PO_mean`,  na.rm = TRUE),
        TRUE                               ~ NA_real_
      ),
      .groups = "drop"
    )
  
  rm(df_raw)
  gc()
  
  return(df_fusionne)
}

annees <- 2016:2024
OLCI_2016_2024_spm_pixels <- map_dfr(annees, fusionner_annee)

save(OLCI_2016_2024_spm_pixels,
     file = "data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")



# Combien de pixels par satellite ce jour-là ?
OLCI_2016_2024_spm_pixels |>
  filter(date == "2024-01-12") |>
  summarise(
    n_total        = n(),
    n_valides      = sum(!is.na(`SPM-G-PO_mean`)),
    couverture     = n_valides / n()
  )

# Avant déduplication, les deux satellites avaient-ils des pixels ?
# Recharge les fichiers bruts pour ce jour
files_12jan_A <- list.files(
  "~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/",
  pattern = "20240112",
  full.names = TRUE
)
files_12jan_B <- list.files(
  "~/Downloads/OLCI/SPM/2024/OLCI_B_ODATIS_MR_2024_SPM/",
  pattern = "20240112",
  full.names = TRUE
)

cat("Fichiers S3A :", length(files_12jan_A), "\n")
cat("Fichiers S3B :", length(files_12jan_B), "\n")

# Pixels valides dans chaque fichier
if (length(files_12jan_A) > 0) {
  nc <- nc_open(files_12jan_A[1])
  spm_a <- ncvar_get(nc, "SPM-G-PO_mean")
  nc_close(nc)
  cat("S3A pixels valides :", sum(!is.na(spm_a)), "/", length(spm_a), "\n")
}

if (length(files_12jan_B) > 0) {
  nc <- nc_open(files_12jan_B[1])
  spm_b <- ncvar_get(nc, "SPM-G-PO_mean")
  nc_close(nc)
  cat("S3B pixels valides :", sum(!is.na(spm_b)), "/", length(spm_b), "\n")
}
# Combien de lignes avant déduplication pour ce jour ?
OLCI_2016_2024_spm_pixels |>
  filter(date == "2024-01-12") |>
  count(lon, lat) |>
  count(n)  # combien de pixels ont 1 observation vs 2 observations

OLCI_2016_2024_spm_pixels <- OLCI_2016_2024_spm_pixels |>
  group_by(lon, lat, date) |>
  summarise(
      `SPM_G_PO_mean` = if (all(is.na(`SPM-G-PO_mean`))) NA_real_
      else mean(`SPM-G-PO_mean`, na.rm = TRUE),
      .groups = "drop"
      )

OLCI_2016_2024_spm_pixels <- OLCI_2016_2024_spm_pixels |>
  filter(!is.na(`SPM-G-PO_mean`)) |>   # on vire les NA avant de grouper
  group_by(lon, lat, date) |>
  summarise(SPM-G-PO_mean = mean(`SPM-G-PO_mean`), .groups = "drop")

# OLCI_2016_2024_spm_pixels <- OLCI_2016_2024_spm_pixels

save(OLCI_2016_2024_spm_pixels, 
     file = "data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

save(OLCI_2016_2024_spm_pixels, 
     file = "data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")


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
seuil_99 <- quantile(OLCI_2016_2024_spm_pixels$`SPM-G-PO_mean`, 0.99, na.rm = TRUE)
# cat("Seuil 95ème percentile :", seuil_95, "mg/m³\n")
# seuil_95 = 1.7

# Seuil 95ème percentile : 0.4883142 mg/m³

# Stats du panache par jour
OLCI_2016_2024_spm_99 <- OLCI_2016_2024_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean` >= seuil_99, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99], na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

# save(OLCI_2016_2024_spm_95, file = "data/OLCI/SPM/OLCI_2016_2024_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# mean/median spm or panache area plot

# en échelle normale

# model_OLCI_2016_95 <- lm(aire_panache_km2 ~ date, data = OLCI_2016_2024_spm_95)
model_OLCI_2016_95 <- lm(median_spm ~ date, data = OLCI_2016_2024_spm_95)
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
  labs(title = "Évolution de la concentration moyenne en MES dans les panaches de la baie des Anges vu par le produit OLCI (ODATIS-MR)",
       x = "Date",
       y = "Concentration médiane en MES (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# en échelle log

data_log_spm <- OLCI_2016_2024_spm_95 |> 
  filter(aire_panache_km2> 0)

model_OLCI_2016_95_log <- lm(log10(median_spm) ~ date, data = data_log_spm)
p_value_OLCI_2016_95_log <- summary(model_OLCI_2016_95_log)$coefficients[2, 4]  # p-value pour la pente
intercept_OLCI_2016_95_log <- coef(model_OLCI_2016_95_log)[1]
slope_OLCI_2016_95_log <- coef(model_OLCI_2016_95_log)[2]

ggplot(data = data_log_spm, aes(x = date, y = median_spm)) +  # utilise data_log_spm
  geom_point(color = "red3", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, 
              formula = y ~ x,
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$median_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_OLCI_2016_95_log, 3), " + ", 
      round(slope_OLCI_2016_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_OLCI_2016_95_log < 0.001, "< 0.001", 
                           format(p_value_OLCI_2016_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de la concentration médiane en MES dans les panaches de la baie des Anges vu par le produit OLCI (ODATIS-Mr)(échelle log)",
       x = "Date",
       y = "Concentration médiane en MES (en mg/m3)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

### Plume extension / liquid flow rate --------------------------------------------------

All_debit <- All_debit |> 
  filter(date >= as.Date("2016-04-26"), date <= as.Date("2024-12-31")) |> 
  drop_na(debit_cumule)

# ── Tests Mann-Kendall + Theil-Sen sur l'aire des panaches ──
OLCI_2016_2024_spm_99_clean <- OLCI_2016_2024_spm_99 |>
  drop_na(aire_panache_km2) |>
  mutate(date_num = as.numeric(date - min(date)))

mk_panache      <- mk.test(OLCI_2016_2024_spm_99_clean$aire_panache_km2)
sen_panache     <- sens.slope(OLCI_2016_2024_spm_99_clean$aire_panache_km2)
slope_kmjan_pan <- sen_panache$estimates * 365

intercept_pan <- median(
  OLCI_2016_2024_spm_99_clean$aire_panache_km2 - sen_panache$estimates * OLCI_2016_2024_spm_99_clean$date_num,
  na.rm = TRUE
)

OLCI_2016_2024_spm_99_clean <- OLCI_2016_2024_spm_99_clean |>
  mutate(theilsen_fit = intercept_pan + sen_panache$estimates * date_num)

cat("Panaches — Mann-Kendall p =", round(mk_panache$p.value, 4),
    "| Theil-Sen pente =", round(slope_kmjan_pan, 2), "km²/an\n")

# ── Tests Mann-Kendall + Theil-Sen sur le débit cumulé ──

All_debit <- All_debit |> 
  filter(date >= as.Date("2016-04-26"), date <= as.Date("2024-12-31")) |> 
  drop_na(debit_cumule) |> 
  mutate(date_num = as.numeric(date - min(date)))

mk_debit       <- mk.test(All_debit$debit_cumule)
sen_debit      <- sens.slope(All_debit$debit_cumule)
slope_debit_an <- sen_debit$estimates * 365

intercept_debit <- median(
  All_debit$debit_cumule - sen_debit$estimates * All_debit$date_num,
  na.rm = TRUE
)

All_debit <- All_debit |>
  mutate(theilsen_fit_debit = intercept_debit + sen_debit$estimates * date_num)

cat("Débit — Mann-Kendall p =", round(mk_debit$p.value, 4),
    "| Theil-Sen pente =", round(slope_debit_an, 3), "m³/s/an\n")

# ── Mise à l'échelle ──
adjust_factors <- sec_axis_adjustement_factors(
  OLCI_2016_2024_spm_99_clean$aire_panache_km2,
  All_debit$debit_cumule
)

OLCI_2016_2024_spm_99_clean <- OLCI_2016_2024_spm_99_clean |>
  mutate(
    scaled_aire_panache_km2 = aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust,
    theilsen_fit_scaled     = theilsen_fit     * adjust_factors$diff + adjust_factors$adjust
  )

# ── Corrélation Spearman ──
merged_data <- merge(
  All_debit |> select(date, debit_cumule),
  OLCI_2016_2024_spm_99_clean |> select(date, aire_panache_km2),
  by = "date", all = FALSE
) |> drop_na()

correlation <- cor(merged_data$debit_cumule, merged_data$aire_panache_km2,
                   method = "spearman", use = "complete.obs")
p_value     <- cor.test(merged_data$debit_cumule, merged_data$aire_panache_km2,
                        method = "spearman")$p.value

n_panache <- nrow(OLCI_2016_2024_spm_99_clean)
n_commun  <- nrow(merged_data)

ggplot() +
  # Aire des panaches
  geom_line(
    data = OLCI_2016_2024_spm_99_clean,
    aes(x = date, y = scaled_aire_panache_km2, color = "Aire des panaches"),
    linewidth = 0.4, alpha = 0.6
  ) +
  # Tendance Theil-Sen panaches
  geom_line(
    data = OLCI_2016_2024_spm_99_clean,
    aes(x = date, y = theilsen_fit_scaled, color = "Tendance panaches"),
    linewidth = 1.2
  ) +
  # Débit cumulé
  geom_line(
    data = All_debit,
    aes(x = date, y = debit_cumule, color = "Débit cumulé"),
    linewidth = 0.4, alpha = 0.6
  ) +
  # Tendance Theil-Sen débit
  geom_line(
    data = All_debit,
    aes(x = date, y = theilsen_fit_debit, color = "Tendance débit"),
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Aire des panaches" = "darkcyan",
      "Tendance panaches" = "darkcyan",
      "Débit cumulé"      = "darkolivegreen3",
      "Tendance débit"    = "darkolivegreen3"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    name     = "Débit cumulé (m³/s)",
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Aire des panaches (km²)")
    )
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  # Corrélation Spearman — haut droite
  annotate(
    "text",
    x = max(OLCI_2016_2024_spm_99_clean$date, na.rm = TRUE),
    y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.97,
    hjust = 1, vjust = 1, size = 8,
    color = "grey20", fontface = "italic", family = "serif",
    label = paste0(
      "R = ", round(correlation, 2),
      "\np ", ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3)),
      "\nn = ", n_commun
    )
  ) +
  # Tendance Theil-Sen panaches — haut gauche
  annotate(
    "text",
    x = as.Date("2016-04-26"),
    y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.97,
    hjust = 0, vjust = 1, size = 8,
    color = "darkcyan", fontface = "italic", family = "serif",
    label = paste0(
      "Tendance panaches : ", round(slope_kmjan_pan, 2), " km²/an",
      "\np ", ifelse(mk_panache$p.value < 0.001, "< 0.001",
                                  ifelse(mk_panache$p.value < 0.05, "< 0.05",
                                         paste0("= ", round(mk_panache$p.value, 3)))),
      "\nn = ", n_panache
    )
  ) +
  # Tendance Theil-Sen débit — milieu gauche
  annotate(
    "text",
    x = as.Date("2016-04-26"),
    y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.70,
    hjust = 0, vjust = 1, size = 8,
    color = "darkolivegreen4", fontface = "italic", family = "serif",
    label = paste0(
      "Tendance débit : ", round(slope_debit_an, 2), " m³/s/an",
      "\np ", ifelse(mk_debit$p.value < 0.001, "< 0.001",
                                  ifelse(mk_debit$p.value < 0.05, "< 0.05",
                                         paste0("= ", round(mk_debit$p.value, 3))))
    )
  ) +
  labs(
    title    = "Évolution de l'extension des panaches turbides et du débit cumulé",
    subtitle = "Produit OLCI — ODATIS-MR (2016–2024)",
    x        = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 18, hjust = 0.5, family = "serif"),
    plot.subtitle    = element_text(size = 15, hjust = 0.5, color = "grey30", family = "serif"),
    axis.title       = element_text(face = "bold", family = "serif", size = 16),
    axis.text        = element_text(color = "grey30", family = "serif", size = 14),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    legend.position  = "top",
    legend.text      = element_text(size = 15),
    plot.margin      = margin(1, 1.5, 1, 1, "cm")
  )

# comparison between liquid flow rate and mean SPM concentration

adjust_factors <- sec_axis_adjustement_factors(OLCI_2016_2024_spm_95$median_spm, Y6442010_2016_2024$débit)

OLCI_2016_2024_spm_95$scaled_median_spm <- OLCI_2016_2024_spm_95$median_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Y6442010_2016_2024, 
             aes(x = date, y = débit, color = "Débit du Var"), size = 0.5) +
  geom_point(data = OLCI_2016_2024_spm_95, 
             aes(x = date, y = scaled_median_spm, color = "Concentration en MES"), size = 0.5) +
  scale_color_manual(values = c("Débit du Var" = "blue", "Concentration en MES" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension moyenne (en g/m³)")
  ) +
  labs(title = "Évolution de la concentration médiane en MES et du débit du Var vu par le produit OLCI (ODATIS-MR)",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# runoff vs mean SPM concentration correlation ---------------------------------

Var_OLCI_panache <- inner_join(Y6442010_2016_2024, OLCI_2016_2024_spm_95, by = "date")

cor.test(Var_OLCI_panache$débit, Var_OLCI_panache$median_spm, method = "spearman")



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






# Gangloff SPM ------------------------------------------------------------

## zones emboîtées ----------------------------------------------------------

# # définir les coordonnées de l'embouchure du Var (comme SEXTANT OC5 chla)
# lon_embouchure <- 7.199082
# lat_embouchure <- 43.654709
# 
# # définir des zones emboîtées de tailles croissantes autour de l'embouchure
# rayons_km <- c(5, 10, 20, 40, 70, 100)  # en km
# 
# # Convertir en degrés
# rayons_deg <- rayons_km / 111
# 
# # Calculer le percentile 95 pour chaque zone emboîtée
# seuils_zones <- lapply(rayons_deg, function(r) {
#   
#   pixels_zone <- OLCI_2016_2024_spm_pixels |>
#     filter(
#       lon >= lon_embouchure - r & lon <= lon_embouchure + r,
#       lat >= lat_embouchure - r & lat <= lat_embouchure + r
#     )
#   
#   aire_km2 <- nrow(distinct(pixels_zone, lon, lat)) * aire_pixel_km2
#   seuil    <- quantile(pixels_zone$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
#   
#   data.frame(rayon_km = r * 111, aire_km2 = aire_km2, seuil_95 = seuil)
# })
# 
# seuils_zones_df <- bind_rows(seuils_zones)
# print(seuils_zones_df)
# 
# # Visualiser le plateau
# ggplot(seuils_zones_df, aes(x = aire_km2, y = seuil_95)) +
#   geom_point(size = 3, color = "steelblue") +
#   geom_line() +
#   geom_hline(yintercept = seuil_95, linetype = "dashed", color = "red") +
#   labs(
#     title = "Détermination du seuil de détection du panache turbide",
#     subtitle = "Percentile 95 par zone emboîtée autour de l'embouchure",
#     x = "Aire de la zone (km²)",
#     y = "Percentile 95 des MES (g/m³)"
#   ) +
#   theme_bw()
# 
# # Le seuil retenu est la valeur du plateau (zone > ~5000 km²)
# seuil_retenu <- seuils_zones_df |>
#   filter(aire_km2 > 1900) |>
#   summarise(seuil = mean(seuil_95)) |>
#   pull(seuil)
# 
# cat("Seuil retenu :", seuil_retenu, "g/m³\n", na.rm = TRUE)
# # 0.4844611 g/m³

seuil_retenu <- 0.4883142

## ROPP --------------------------------------------------------------------

# Pixels où le panache est présent dans au moins 5% des images
n_images_total <- n_distinct(OLCI_2016_2024_spm_pixels$date)
# 2851

ROPP <- OLCI_2016_2024_spm_pixels |>
  group_by(lon, lat) |>
  summarise(
    freq_above_seuil = sum(`SPM-G-PO_mean` >= seuil_retenu, na.rm = TRUE) / n_images_total,
    .groups = "drop"
  ) |>
  filter(freq_above_seuil >= 0.05)   # au moins 5% des images

cat("Nombre de clean dans la ROPP :", nrow(ROPP), "\n")
# 1740

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

ggplot(ROPP) +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = freq_above_seuil)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Fréquence au-dessus du seuil",
    labels = scales::percent_format(accuracy = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Région d'occurrence des panaches turbides (ROPP)",
    # subtitle = "Pixels où la concentration en MES dépasse le seuil dans au moins 5% des images MODIS",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim        = range(ROPP$lon),
    ylim        = range(ROPP$lat),
    expand      = FALSE,
    default_crs = sf::st_crs(4326)
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

## stat panache ------------------------------------------------------------

# Statistiques du panache par jour
OLCI_panache_metrics <- OLCI_2016_2024_spm_pixels |>
  semi_join(ROPP, by = c("lon", "lat")) |>
  group_by(date) |>
  summarise(
    # Aire d'extension
    pixel_count = sum(`SPM-G-PO_mean` >= seuil_retenu, na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    # Concentrations
    mean_spm         = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    max_spm          = max(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    median_spm       = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    # Points extrêmes
    lat_sud          = min(lat[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    lon_ouest        = min(lon[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    lon_est          = max(lon[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    # Centroïde
    centroid_lon     = mean(lon[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    centroid_lat     = mean(lat[`SPM-G-PO_mean` >= seuil_retenu], na.rm = TRUE),
    .groups = "drop"
  ) |> 
  filter(pixel_count > 0)

# plotting ----------------------------------------------------------------

# correlation river flow and plume extension ------------------------------

debit_lags <- All_debit_2016_2024 |>
  arrange(date) |>
  select(date, debit_cumule) |> 
  mutate(
    debit_j1     = lag(debit_cumule, 1),             # j-1 seulement
    debit_2j     = (lag(debit_cumule, 1) + lag(debit_cumule, 2)) / 2,       # moyenne j-1, j-2
    debit_3j     = (lag(debit_cumule, 1) + lag(debit_cumule, 2) + lag(debit_cumule, 3)) / 3
  )

# Jointure et nettoyage
OLCI_panache_metrics <- OLCI_panache_metrics |>
  inner_join(debit_lags |> select(date, debit_cumule, debit_j1, debit_2j, debit_3j), by = "date") |>
  filter(
    aire_panache_km2 > 0,
    debit_cumule > 0,
    is.finite(debit_cumule)
  ) |>
  mutate(
    panache_log  = log10(aire_panache_km2),
    debit_j0_log = log10(debit_cumule),
    debit_j1_log = log10(debit_j1),
    debit_2j_log = log10(debit_2j),
    debit_3j_log = log10(debit_3j)
  )

# # Tout est déjà dans SEXTANT_panache_metrics, pas besoin de rejoindre !
OLCI_panache_metrics |>
  filter(aire_panache_km2 > 0) |>
  summarise(
    r_j0 = cor(log10(aire_panache_km2), debit_j0_log, use = "complete.obs", method = "spearman"),
    r_j1 = cor(log10(aire_panache_km2), debit_j1_log, use = "complete.obs", method = "spearman"),
    r_2j = cor(log10(aire_panache_km2), debit_2j_log, use = "complete.obs", method = "spearman"),
    r_3j = cor(log10(aire_panache_km2), debit_3j_log, use = "complete.obs", method = "spearman")
  )

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ debit_j0_log, data = OLCI_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1_log, data = OLCI_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j_log, data = OLCI_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j_log, data = OLCI_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ debit_cumule, data = OLCI_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1, data = OLCI_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j, data = OLCI_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j, data = OLCI_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

OLCI_panache_metrics <- OLCI_panache_metrics |> 
  filter(max_spm <= 4000)

## modèle log - log --------------------------------------------------------

# Modèle linéaire en log-log → relation puissance
modele_log <- lm(panache_log ~ debit_j0_log, data = OLCI_panache_metrics)
r2     <- summary(modele_log)$r.squared
print(r2)
pente  <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

# Équation puissance : aire = 10^ordonnee * debit_cumule^pente
# à afficher dans le graphique
label_eq <- paste0(
  "Aire = ", round(10^ordonnee, 3), " × Q^", round(pente, 2),
  "\nR² = ", round(r2, 2)
)

# Graphique log-log avec droite de régression
ggplot(OLCI_panache_metrics, aes(x = debit_cumule, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,            # régression sur les axes log
              color = "black", se = FALSE, linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text", 
           x = min(OLCI_panache_metrics$debit_log, na.rm = TRUE) * 1.5,
           y = max(OLCI_panache_metrics$aire_panache_km2, na.rm = TRUE) * 0.7,
           label = label_eq, hjust = 1.1, vjust = 2, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit cumulé (m"^{3}*".s"^{-1}*")"),
    y = "Aire du panache (km²)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title = element_text(face = "bold", family = "serif"),
    axis.text = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")  # Plus de marge à droite pour l'annotation
  )

## Modèle semi-log : log10(aire) ~ debit_cumule --------------------------------------------------------

modele_semilog <- lm(panache_log ~ debit_cumule, data = OLCI_panache_metrics)
r2       <- summary(modele_semilog)$r.squared
print(r2)
pente    <- coef(modele_semilog)[2]
ordonnee <- coef(modele_semilog)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(OLCI_panache_metrics, aes(x = debit_cumule, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  # geom_smooth(method = "lm", formula = y ~ x,
  #             color = "black", se = FALSE, linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  annotate("text",
           x = max(OLCI_panache_metrics$debit_cumule, na.rm = TRUE) * 0.7,
           y = min(OLCI_panache_metrics$aire_panache_km2, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = "Aire du panache (km²)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "serif"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey50", family = "serif"),
    axis.title = element_text(face = "bold", family = "serif"),
    axis.text = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")  # Plus de marge à droite pour l'annotation
  )

# corrélation MES mean et debit_cumule -------------------------------------------

modele <- lm(debit_cumule ~ mean_spm, data = OLCI_panache_metrics)
r2       <- summary(modele)$r.squared
print(r2)
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "MES = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(OLCI_panache_metrics, aes(x = debit_cumule, y = mean_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x        = Inf,
           y        = Inf,
           label    = label_eq,
           hjust    = 1.1,   # ancré à droite
           vjust    = 6,   # ancré en haut
           size     = 8,
           color    = "grey20",
           family   = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit cumulé (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration moyenne en MES (g m"^{-3}*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    plot.margin      = margin(1, 1.5, 1, 1, "cm")
  )

# corrélation MES max et debit_cumule -------------------------------------------

modele <- lm(debit_cumule ~ max_spm, data = OLCI_panache_metrics)
r2       <- summary(modele)$r.squared
print(r2)
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "MES = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(OLCI_panache_metrics, aes(x = debit_cumule, y = max_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x        = Inf,
           y        = Inf,
           label    = label_eq,
           hjust    = 1.1,   # ancré à droite
           vjust    = 6,   # ancré en haut
           size     = 8,
           color    = "grey20",
           family   = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit cumulé (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration maximale en MES (g m"^{-3}*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    plot.margin      = margin(1, 1.5, 1, 1, "cm")
  )

# cartographie ------------------------------------------------------------

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

OLCI_03_05_2013 <- OLCI_B_2020_spm_pixels |> 
  filter(date == "2013-05-03")

OLCI_03_10_2020 <- OLCI_B_2020_spm_pixels |> 
  filter(date == "2020-10-03")

OLCI_16_04_2024 <- OLCI_A_2024_spm_pixels |> 
  filter(date == "2024-04-16")

OLCI_2024 <- OLCI_2016_2024_spm_pixels |> 
  filter(date >= as.Date("2024-04-12"), date <= as.Date("2024-04-18"))

# max_spm <- max(OLCI_03_10_2020$`SPM-G-PO_mean`, na.rm = TRUE)
# max_spm <- 100
# max_spm <- max(OLCI_03_05_2013$`SPM-G-PO_mean`, na.rm = TRUE)
max_spm <- max(OLCI_16_04_2024$`SPM-G-PO_mean`, na.rm = TRUE)


pl_map <- OLCI_16_04_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
  geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
  
  # Flèche nord
  annotation_north_arrow(
    location = "tr",          # top-right
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"),
    width  = unit(1.5, "cm")
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression("MES (g m"^{-3}*")"),  # ← écriture scientifique
    limits = c(0, max_spm)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 20,
    barheight      = 2,
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Concentration en matières en suspension avec OLCI B — 16 avril 2024",
    subtitle = "Concentration en MES (g/m³)",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(OLCI_16_04_2024$lon),
    ylim   = range(OLCI_16_04_2024$lat),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 14, face = "bold", margin = margin(b = 5)),
    plot.subtitle    = element_text(size = 12, color = "grey50", margin = margin(b = 10)),
    panel.border     = element_rect(colour = "black", fill = NA),
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12)
  )

# Save as desired
ggsave("~/Satellite_analysis/Graphiques/OLCI/carto_OLCI_A_16_04_2024.png", pl_map, height = 9, width = 14)

# on plotte la semaine entière pour montrer qu'on loupe de la dynamique

library(ggpubr)

# Créer un graphique par jour
dates_semaine <- seq(as.Date("2024-04-12"), as.Date("2024-04-18"), by = "day")
lettres <- c("a)", "b)", "c)", "d)", "e)", "f)", "g)")

plots <- map2(dates_semaine, lettres, function(d, lettre) {
  
  df_jour <- OLCI_2024 |> filter(date == d)
  
  ggplot() +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    geom_tile(data = df_jour, aes(x = lon, y = lat, fill = `SPM-G-PO_mean`)) +
    geom_sf(data = countries_giscoR, colour = "black", fill = "grey80", linewidth = 0.3) +
    scale_fill_viridis_c(
      option   = "plasma",
      name     = expression("MES (g m"^{-3}*")"),
      limits   = c(0, max_spm),
      na.value = "transparent"
    ) +
    labs(
      title = paste0(lettre, " ", format(d, "%d %B %Y")),
      x = NULL, y = NULL
    ) +
    coord_sf(
      xlim        = range(OLCI_2024$lon),
      ylim        = range(OLCI_2024$lat),
      expand      = FALSE,
      default_crs = sf::st_crs(4326)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      panel.border     = element_rect(colour = "black", fill = NA),
      panel.grid.minor = element_blank(),
      legend.position  = "none",  # légende commune en dessous
      axis.text        = element_text(size = 7)
    )
})

# Légende commune
legende <- get_legend(
  plots[[1]] +
    guides(fill = guide_colorbar(
      barwidth       = 15,
      barheight      = 1,
      title.position = "top",
      title.hjust    = 0.5
    )) +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = 11),
          legend.text     = element_text(size = 9))
)

# Assembler
figure <- ggarrange(
  plotlist = plots,
  ncol     = 4,
  nrow     = 2,
  legend   = "none"
)

# Ajouter titre général et légende
ggarrange(
  figure,
  legende,
  ncol    = 1,
  heights = c(10, 1)
) |>
  annotate_figure(
    top = text_grob(
      "Distribution spatiale des MES — semaine du 12 au 18 avril 2024",
      face = "bold", size = 13
    ),
    bottom = text_grob(
      "OLCI - Correction atmosphérique = Polymer",
      color = "grey50", size = 10
    )
  )

# climatologie de la concentration en MES ----------------------------------

coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

OLCI_2016_2024_spm_clean <- OLCI_2016_2024_spm_pixels |> 
  filter(`SPM-G-PO_mean` >= 0, `SPM-G-PO_mean` < 100)

sum(is.na(OLCI_2016_2024_spm_clean$`SPM-G-PO_mean`))

OLCI_2016_2024_spm_clean <- OLCI_2016_2024_spm_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# ── Climatologie spatiale mensuelle (moyenne par pixel et par mois) ──
clim_spatiale_spm_month_OLCI <- OLCI_2016_2024_spm_clean |>
  filter(`SPM-G-PO_mean` >= 0) |>
  mutate(month = month(date)) |>
  group_by(lon, lat, month) |>
  summarise(
    mean_spm  = mean(`SPM-G-PO_mean`, na.rm = TRUE),
    median_spm = median(`SPM-G-PO_mean`, na.rm = TRUE),
    sd_spm    = sd(`SPM-G-PO_mean`, na.rm = TRUE),
    .groups   = "drop"
  )

range(clim_spatiale_spm_month_OLCI$mean_spm, na.rm = TRUE)  # OLCI
range(clim_spatiale_spm_month_sextant$mean_spm, na.rm = TRUE)  # SEXTANT

# Combien de valeurs très élevées ?
sum(OLCI_2016_2024_spm_clean$`SPM-G-PO_mean` > 100, na.rm = TRUE)
sum(OLCI_2016_2024_spm_clean$`SPM-G-PO_mean` > 1000, na.rm = TRUE)

# Quels pixels sont concernés ?
OLCI_2016_2024_spm_clean |>
  filter(`SPM-G-PO_mean` > 100) |>
  select(lon, lat, date, `SPM-G-PO_mean`) |>
  arrange(desc(`SPM-G-PO_mean`))

# plot mensuel de la concentration en MES
ggplot(clim_spatiale_spm_month_OLCI, aes(x = lon, y = lat, fill = mean_spm)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_spm_month_OLCI$lon),
    ylim = range(clim_spatiale_spm_month_OLCI$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("MES (g. m"^{-3}*")"),
    option   = "turbo",
    na.value = "white",
    limits   = c(0.1, 50),
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10"),
    oob      = scales::squish
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de la concentration en matières en suspension",
    subtitle = "2016-2024 · OLCI (ODATIS-MR)",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,   # ← large comme SEXTANT
    barheight      = 0.8,  # ← fine
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )

# plot mensuel de la sd de MES
ggplot(clim_spatiale_spm_month_OLCI, aes(x = lon, y = lat, fill = sd_spm)) +
  geom_raster() +
  geom_sf(data = countries_giscoR, fill = "grey92", color = "grey40",
          inherit.aes = FALSE, linewidth = 0.25) +
  coord_sf(
    xlim = range(clim_spatiale_spm_month_OLCI$lon),
    ylim = range(clim_spatiale_spm_month_OLCI$lat),
    expand = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(6.8, 7.4, by = 0.3),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    breaks = seq(43.2, 43.8, by = 0.3),
    labels = function(y) paste0(y, "°N")
  ) +
  scale_fill_viridis_c(
    trans    = "log10",
    name     = expression("MES (g. m"^{-3}*")"),
    option   = "turbo",
    na.value = "white",
    limits   = c(0.1, 50),
    breaks   = c(0.01, 0.1, 1, 10),
    labels   = c("0.01", "0.1", "1", "10"),
    oob      = scales::squish
  ) +
  facet_wrap(~ month, ncol = 6,
             labeller = labeller(month = c(
               "1"  = "Janvier",  "2"  = "Février",   "3"  = "Mars",
               "4"  = "Avril",    "5"  = "Mai",        "6"  = "Juin",
               "7"  = "Juillet",  "8"  = "Août",       "9"  = "Septembre",
               "10" = "Octobre",  "11" = "Novembre",   "12" = "Décembre"
             ))) +
  labs(
    title    = "Climatologie spatiale mensuelle de l'erreur standard à la concentration en MES",
    subtitle = "2016-2024 · OLCI (ODATIS-MR)",
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 15,   # ← large comme SEXTANT
    barheight      = 0.8,  # ← fine
    ticks          = TRUE,
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal"
  )) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 9),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey60", linewidth = 0.3),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey60", linewidth = 0.4),
    panel.spacing    = unit(0.15, "lines"),
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 13, margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0),
    plot.margin      = margin(10, 10, 10, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )


nc <- nc_open("~/Downloads/nommmmm/L3m_20020619__FRANCE_03_MER_CDOM-PO_DAY_00.nc")
print(nc)
