# TS analysis MODIS monthly

# 08/04/2026

# pathway : "~/Satellite_analysis/TS analysis MODIS monthly"


# This script will load MODIS monthly satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the 
# Var river at the Napoleon bridge and the total runoff

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(stars)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)
library(ggtext)

# setup -------------------------------------------------------------------

# Get satellite download function
source("~/sat_access/sat_access_script.R")

# lon lat ranges
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)

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

## df treatment function --------------------------------------------------------

filename <- "~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2015_SPM/L3m_20150101-20150131__FRANCE_03_MOD_SPM-G-NS_MO_00.nc"

nc <- tidync(filename)
print(nc)

load_MODIS_spm_pixels <- function(file_name, lon_range, lat_range){
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  # The necessary code
  MODIS_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    tidync::hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = file_date) |> 
    dplyr::select(lon, lat, date, `SPM-G-NS_mean`)
  
  # Exit
  return(MODIS_one)
}

# load data ---------------------------------------------------------------
## Hydro France data ---------------------------------------------------------------------

load("data/Hydro France/All_debit_2015_mensuel.Rdata")

load("data/Hydro France/Var_mensuel.Rdata")

## Satellite data ---------------------------------------------------------------------

load("data/MODIS/mensuel/MODIS_month_spm_pixels.Rdata")

load("data/MODIS/mensuel/MODIS_month_spm_95.Rdata")

# loading data ------------------------------------------------------------
## SPM ---------------------------------------------------------------------

MODIS_2015_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2015_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2016_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2016_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2017_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2017_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2018_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2018_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2019_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2019_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2020_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2020_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2021_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2021_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2022_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2022_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2023_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2023_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MODIS_2024_month_spm_dir <- dir("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2024_SPM/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### to define threshold with percentile 95 --------------------------------------------------------

# Load and combine
MODIS_2015_month_spm_pixels <- plyr::ldply(MODIS_2015_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2016_month_spm_pixels <- plyr::ldply(MODIS_2016_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2017_month_spm_pixels <- plyr::ldply(MODIS_2017_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2018_month_spm_pixels <- plyr::ldply(MODIS_2018_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2019_month_spm_pixels <- plyr::ldply(MODIS_2019_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2020_month_spm_pixels <- plyr::ldply(MODIS_2020_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2021_month_spm_pixels <- plyr::ldply(MODIS_2021_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2022_month_spm_pixels <- plyr::ldply(MODIS_2022_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2023_month_spm_pixels <- plyr::ldply(MODIS_2023_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MODIS_2024_month_spm_pixels <- plyr::ldply(MODIS_2024_month_spm_dir, load_MODIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

# Combine and save
MODIS_month_spm_pixels <- rbind(MODIS_2015_month_spm_pixels, MODIS_2016_month_spm_pixels, 
                                 MODIS_2017_month_spm_pixels, MODIS_2018_month_spm_pixels, 
                                 MODIS_2019_month_spm_pixels, MODIS_2020_month_spm_pixels, 
                                 MODIS_2021_month_spm_pixels, MODIS_2022_month_spm_pixels, 
                                 MODIS_2023_month_spm_pixels,MODIS_2024_month_spm_pixels)

save(MODIS_month_spm_pixels, file = "data/MODIS/mensuel/MODIS_month_spm_pixels.Rdata")

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

# Ou inspecter les coordonnées lon/lat directement
nc <- tidync("~/Downloads/MODIS ODATIS MR/mensuel/SPM/MODIS_2015_SPM/L3m_20150101-20150131__FRANCE_03_MOD_SPM-G-NS_MO_00.nc")

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
seuil_95 <- quantile(MODIS_month_spm_pixels$`SPM-G-NS_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "mg/m³\n")

# Seuil 95ème percentile : 0.576409 mg/m³

# Stats du panache par jour
MODIS_month_spm_95 <- MODIS_month_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-NS_mean` >= seuil_95, na.rm = TRUE),
    mean_spm = mean(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95], na.rm = TRUE),
    median_spm = median(`SPM-G-NS_mean`[`SPM-G-NS_mean` >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MODIS_month_spm_95, file = "data/MODIS/mensuel/MODIS_month_spm_95.Rdata")

# plotting ----------------------------------------------------------------

# mean / median SPM or panache area plot

# en échelle normale

# model_MODIS_2002_95 <- lm(aire_panache_km2 ~ date, data = MODIS_2002_2024_spm_95)
model_MODIS_month <- lm(median_spm ~ month, data = MODIS_month_spm_95)
p_value_MODIS_month <- summary(model_MODIS_month)$coefficients[2, 4]  # p-value pour la pente
intercept_MODIS_month <- coef(model_MODIS_month)[1]
slope_MODIS_month <- coef(model_MODIS_month)[2]

ggplot(data = MODIS_month_spm_95, aes(x = month, y = median_spm)) +
  
  geom_point(
    color     = "red3",
    size      = 0.8,
    alpha     = 0.6
  ) +
  geom_smooth(
    method    = "lm",
    se        = TRUE,
    color     = "darkslateblue",
    fill      = "#AED6F1",
    linewidth = 0.8,
    alpha     = 0.3
  ) +
  
  annotate(
    "text",
    x     = min(MODIS_month_spm_95$month, na.rm = TRUE),
    y     = max(MODIS_month_spm_95$median_spm, na.rm = TRUE) * 0.99,
    label = paste0("y = ", round(intercept_MODIS_month, 3),
                   " ", round(slope_MODIS_month, 7), " × x"),
    hjust = 0, vjust = 1,
    size  = 8, color = "red3", fontface = "italic", family = "serif"
  ) +
  annotate(
    "text",
    x     = min(MODIS_month_spm_95$month, na.rm = TRUE),
    y     = max(MODIS_month_spm_95$median_spm, na.rm = TRUE) * 0.88,
    label = paste0("p = ", ifelse(p_value_MODIS_month < 0.001, "< 0.001",
                                  format(p_value_MODIS_month, scientific = TRUE, digits = 3))),
    hjust = 0, vjust = 1,
    size  = 8, color = "red3", fontface = "italic", family = "serif"
  ) +
  
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  
  labs(
    title    = "Évolution de la concentration médiane mensuelle en MES dans les panaches près de Nice (2015–2024)",
    subtitle = "Régression linéaire — produit MODIS mensuel (ODATIS-MR)",
    x        = NULL,
    y        = expression("Concentration médiane en MES (en mg/m³)"),
    caption  = "Source : ODATIS-MR"
  ) +
  
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey96", linewidth = 0.2),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )

# échelle log

# Modèle en log
model_MODIS_month_log <- lm(log10(mean_spm) ~ month, data = MODIS_month_spm_95)
p_value_MODIS_month_log <- summary(model_MODIS_month_log)$coefficients[2, 4]
intercept_MODIS_month_log <- coef(model_MODIS_month_log)[1]
slope_MODIS_month_log <- coef(model_MODIS_month_log)[2]

ggplot(data = MODIS_month_spm_95, aes(x = month, y = mean_spm)) +
  
  geom_point(
    color     = "red3",
    size      = 0.8,
    alpha     = 0.6
  ) +
  geom_smooth(
    method    = "lm",
    se        = TRUE,
    color     = "darkslateblue",
    fill      = "#AED6F1",
    linewidth = 0.8,
    alpha     = 0.3
  ) +
  
  annotate(
    "text",
    x     = min(MODIS_month_spm_95$month, na.rm = TRUE),
    y     = max(MODIS_month_spm_95$mean_spm, na.rm = TRUE) * 0.99,
    label = paste0("log10(y) = ", round(intercept_MODIS_month_log, 3),
                   " ", round(slope_MODIS_month_log, 7), " × x"),
    hjust = 0, vjust = 1,
    size  = 8, color = "red3", fontface = "italic", family = "serif"
  ) +
  annotate(
    "text",
    x     = min(MODIS_month_spm_95$month, na.rm = TRUE),
    y     = max(MODIS_month_spm_95$mean_spm, na.rm = TRUE) * 0.9,  # ajuste position en log
    label = paste0("p = ", ifelse(p_value_MODIS_month_log < 0.001, "< 0.001",
                                  format(p_value_MODIS_month_log, scientific = TRUE, digits = 3))),
    hjust = 0, vjust = 1,
    size  = 8, color = "red3", fontface = "italic", family = "serif"
  ) +
  
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  
  labs(
    title    = "Évolution de la concentration moyenne mensuelle en MES dans les panaches près de Nice (2015–2024) (en échelle log)",
    subtitle = "Régression linéaire — produit MODIS mensuel (ODATIS-MR)",
    x        = NULL,
    y        = expression("Concentration moyenne en MES (en mg/m³)"),
    caption  = "Source : ODATIS-MR"
  ) +
  
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 11, margin = margin(r = 10)),
    axis.text        = element_text(size = 10, color = "grey30"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey96", linewidth = 0.2),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5)
  )


# comparison between liquid flow rate and panache extension / mean / median SPM

adjust_factors <- sec_axis_adjustement_factors(MODIS_month_spm_95$median_spm, Var_mensuel$débit)

MODIS_month_spm_95$scaled_median_spm <- MODIS_month_spm_95$median_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  
  geom_line(
    data  = All_debit_2015_mensuel,
    aes(x = month, y = mean_runoff_monthly, color = "Débit cumulé"),
    linewidth = 0.8, alpha = 0.6
  ) +
  geom_line(
    data  = MODIS_month_spm_95,
    aes(x = month, y = scaled_median_spm, color = "Concentration en MES"),
    size  = 0.8, alpha = 0.6
  ) +
  
  scale_color_manual(
    values = c("Débit cumulé" = "darkolivegreen3", "Concentration en MES" = "red3"),
    name   = NULL
  ) +
  
  scale_y_continuous(
    name     = expression("Débit (m³ s"^-1*")"),
    sec.axis = sec_axis(
      ~ (. - adjust_factors$adjust) / adjust_factors$diff,
      name = expression("Concentration médiane en MES (mg/m³)")
    ),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  
  scale_x_date(
    date_breaks       = "1 year",
    date_labels       = "%Y",
    date_minor_breaks = "6 months",
    expand            = expansion(mult = 0.01)
  ) +
  
  labs(
    title   = "Évolution de la concentration médiane mensuelle en MES et du débit cumulé du Var, du Paillon et du Magnan (2015–2024)",
    subtitle = "Produit MODIS mensuel (ODATIS-MR)",
    x       = NULL,
    caption = "Source : ODATIS-MR / Hydro France / MNCA"
  ) +
  
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.y.left  = element_text(size = 11, color = "darkolivegreen3", margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 11, color = "red3",   margin = margin(l = 10)),
    axis.text.y.left   = element_text(color = "darkolivegreen3"),
    axis.text.y.right  = element_text(color = "red3"),
    axis.text.x      = element_text(size = 10, color = "grey30", angle = 45, hjust = 1),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey96", linewidth = 0.2),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    legend.key       = element_blank()
  )

# runoff vs SPM concentration correlation ---------------------------------

Var_MODIS_month <- All_debit_2015_mensuel %>%
  select(month, mean_runoff_monthly) %>%
  inner_join(
    MODIS_month_spm_95 %>% select(month, mean_spm, median_spm, aire_panache_km2),
    by = "month"
  )

cor.test(Var_MODIS_month$mean_runoff_monthly, Var_MODIS_month$median_spm, method = "spearman")

# combine runoff for monthly data -----------------------------------------

# All_debit_2015_mensuel <- All_debit_2015 |> 
#   mutate(month = floor_date(date, "month")) |>  # crée une colonne mois
#   group_by(month) |> 
#   summarise(mean_runoff_monthly = mean(debit_cumule, na.rm = TRUE))
# 
# save(All_debit_2015_mensuel, file = "data/All_debit_2015_mensuel.Rdata")

# Var_mensuel <- read.csv("~/River_runoff_analysis/data/Hydro France/Y6442010_QmM_mensuel_2015_2026.csv")
# 
# Var_mensuel <- Var_mensuel |>
# rename(
#   month = "Date..TU.",
#   débit = "Valeur..en.m..s."
#   )
# 
# Var_mensuel <- Var_mensuel |> 
#   mutate(month = as.Date(month))

# Complete missing dates in the date range
Var_mensuel <- Var_mensuel |> 
  filter(débit <= 133)

Var_mensuel <- Var_mensuel |> 
  complete(month = seq(min(month), max(month), by = "month")) 


# les mois manquants sont ajoutés avec NA dans débit

# save(Var_mensuel, file = "data/Hydro France/Var_mensuel.Rdata")

# MODIS_month_spm_95 <- MODIS_month_spm_95 |> 
#   rename(month = "date")
