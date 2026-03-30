# MERIS spatio-temporal analysis
# 27/03/2026

# pathway : "Satellite analysis/MERIS spatio-temporal analysis.Rdata

# This script will load MERIS satellite data
# Then define a threshold based on the 95 percentile of spm datas
# Then define the extension of the plume per day and the average SPM mean
# in this plume

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

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

filename <- "~/Downloads/MERIS/SPM/2002/L3m_20020619__FRANCE_03_MER_SPM-G-PO_DAY_00.nc"

load_MERIS_spm_pixels <- function(file_name, lon_range, lat_range){
  file_caracter <- substr(basename(file_name), start = 5, stop = 12)
  file_date <- as.Date(file_caracter, format = "%Y%m%d")
  
  # The necessary code
   MERIS_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = file_date) |> 
    dplyr::select(lon, lat, date, `SPM-G-PO_mean`)  # tous les pixels, sans filtre
  
  # Exit
  return(MERIS_one)
}

# load data ---------------------------------------------------------------
## Hydro France data ---------------------------------------------------------------------

load("data/Hydro France/Y6442010_2002_2012.Rdata")

load("data/MERIS/MERIS_2002_2012_spm_pixels.Rdata")

load("data/MERIS/MERIS_2002_2012_spm_95.Rdata")

## SPM ---------------------------------------------------------------------
### direction ---------------------------------------------------------------

MERIS_2002_dir <- dir("~/Downloads/MERIS/SPM/2002/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2003_dir <- dir("~/Downloads/MERIS/SPM/2003/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2004_dir <- dir("~/Downloads/MERIS/SPM/2004/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2005_dir <- dir("~/Downloads/MERIS/SPM/2005/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2006_dir <- dir("~/Downloads/MERIS/SPM/2006/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2007_dir <- dir("~/Downloads/MERIS/SPM/2007/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2008_dir <- dir("~/Downloads/MERIS/SPM/2008/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2009_dir <- dir("~/Downloads/MERIS/SPM/2009/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2010_dir <- dir("~/Downloads/MERIS/SPM/2010/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2011_dir <- dir("~/Downloads/MERIS/SPM/2011/", pattern = ".nc", recursive = TRUE, full.names = TRUE)
MERIS_2012_dir <- dir("~/Downloads/MERIS/SPM/2012/", pattern = ".nc", recursive = TRUE, full.names = TRUE)

### define threshold with percentile 95 --------------------------------------------------------

MERIS_2002_spm_pixels <- plyr::ldply(MERIS_2002_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2003_spm_pixels <- plyr::ldply(MERIS_2003_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2004_spm_pixels <- plyr::ldply(MERIS_2004_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2005_spm_pixels <- plyr::ldply(MERIS_2005_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2006_spm_pixels <- plyr::ldply(MERIS_2006_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2007_spm_pixels <- plyr::ldply(MERIS_2007_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2008_spm_pixels <- plyr::ldply(MERIS_2008_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2009_spm_pixels <- plyr::ldply(MERIS_2009_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2010_spm_pixels <- plyr::ldply(MERIS_2010_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2011_spm_pixels <- plyr::ldply(MERIS_2011_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)
MERIS_2012_spm_pixels <- plyr::ldply(MERIS_2012_dir, load_MERIS_spm_pixels, .parallel = TRUE, lon_range = lon_range, lat_range = lat_range)

MERIS_2002_2012_spm_pixels <- rbind(MERIS_2002_spm_pixels, MERIS_2003_spm_pixels, MERIS_2004_spm_pixels,
                                    MERIS_2005_spm_pixels, MERIS_2006_spm_pixels, MERIS_2007_spm_pixels,
                                    MERIS_2008_spm_pixels, MERIS_2009_spm_pixels, MERIS_2010_spm_pixels,
                                    MERIS_2011_spm_pixels, MERIS_2012_spm_pixels)

save(MERIS_2002_2012_spm_pixels, file = "data/MERIS/MERIS_2002_2012_spm_pixels.Rdata")

load("data/MERIS/MERIS_2002_2012_spm_pixels.Rdata")

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

nc <- tidync("~/Downloads/MERIS/SPM/2002/L3m_20020619__FRANCE_03_MER_SPM-G-PO_DAY_00.nc")

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
seuil_95 <- quantile(MERIS_2002_2012_spm_pixels$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "mg/m³\n")

# Seuil 95ème percentile : 0.4883142 mg/m³

# Stats du panache par jour
MERIS_2002_2012_spm_95 <- MERIS_2002_2012_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(`SPM-G-PO_mean` >= seuil_95, na.rm = TRUE),
    mean_spm = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MERIS_2002_2012_spm_95, file = "data/MERIS/MERIS_2002_2012_spm_95.Rdata")

# plotting ----------------------------------------------------------------

MERIS_2002_2012_spm_95_propre <- MERIS_2002_2012_spm_95 |> 
  filter(mean_spm < 2000) 

# mean spm or panache area plot

# en échelle normale

model_MERIS_2002_95 <- lm(aire_panache_km2 ~ date, data = MERIS_2002_2012_spm_95)
model_MERIS_2002_95 <- lm(mean_spm ~ date, data = MERIS_2002_2012_spm_95_propre)
p_value_MERIS_2002_95 <- summary(model_MERIS_2002_95)$coefficients[2, 4]  # p-value pour la pente
intercept_MERIS_2002_95 <- coef(model_MERIS_2002_95)[1]
slope_MERIS_2002_95 <- coef(model_MERIS_2002_95)[2]

# ggplot(data = MERIS_2002_2012_spm_95_propre, aes(x = date, y = mean_spm)) +
#   geom_point(color = "red3", size = 0.5) +
#   # geom_point(data = MERIS_2002_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   annotate(
#     "text",
#     x = max(MERIS_2002_2012_spm_95_propre$date, na.rm = TRUE),
#     y = max(MERIS_2002_2012_spm_95_propre$mean_spm, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "y = ", round(intercept_MERIS_2002_95, 3), " + ", round(slope_MERIS_2002_95, 7), " * x",
#       "\n", "p = ", ifelse(p_value_MERIS_2002_95 < 0.001, "< 0.001", format(p_value_MERIS_2002_95, digits = 3))
#     ),
#     hjust = 1,  # Alignement à droite
#     vjust = 1,  # Alignement en haut
#     size = 6
#   ) +
#   labs(title = "Évolution de la concentration moyenne en MES dans les panaches de la baie des Anges vu par le produit MERIS (ODATIS-MR)",
#        x = "Date",
#        y = "Concentration moyenne en MES (en mg/m³)") +
#   theme_minimal() +
#   scale_x_date(
#     date_breaks = "1 year",
#     date_labels = "%Y"
#   )

ggplot(data = MERIS_2002_2012_spm_95, aes(x = date, y = aire_panache_km2)) +
  geom_point(color = "darkcyan", size = 0.5) +
  # geom_point(data = MERIS_2002_2024_spm_95, aes(x = date, y = mean_spm), color = "red", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  annotate(
    "text",
    x = max(MERIS_2002_2012_spm_95_propre$date, na.rm = TRUE),
    y = max(MERIS_2002_2012_spm_95_propre$aire_panache_km2, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_MERIS_2002_95, 3), " + ", round(slope_MERIS_2002_95, 7), " * x",
      "\n", "p = ", ifelse(p_value_MERIS_2002_95 < 0.001, "< 0.001", format(p_value_MERIS_2002_95, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la taille des panaches de la baie des Anges vu par le produit MERIS (ODATIS-MR)",
       x = "Date",
       y = "Aire des panaches (en km²)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# en échelle log

data_log_spm <- MERIS_2002_2012_spm_95 |>
  filter(mean_spm > 0)

# data_log_spm <- MERIS_2002_2012_spm_95 |>
#   filter(aire_panache_km2 > 0)

model_MODIS_2002_95_log <- lm(log10(mean_spm) ~ date, data = data_log_spm)
# model_MERIS_2002_95_log <- lm(log10(aire_panache_km2) ~ date, data = data_log_spm)
p_value_MERIS_2002_95_log <- summary(model_MERIS_2002_95_log)$coefficients[2, 4]  # p-value pour la pente
intercept_MERIS_2002_95_log <- coef(model_MERIS_2002_95_log)[1]
slope_MERIS_2002_95_log <- coef(model_MERIS_2002_95_log)[2]

# ggplot(data = data_log_spm, aes(x = date, y = aire_panache_km2)) +  # utilise data_log_spm
#   geom_point(color = "darkcyan", size = 0.5) +
#   geom_smooth(method = "lm", se = TRUE, 
#               formula = y ~ x,
#               color = "darkslateblue", fill = "pink", alpha = 0.2) +
#   scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
#   annotate(
#     "text",
#     x = max(data_log_spm$date, na.rm = TRUE),
#     y = max(data_log_spm$aire_panache_km2, na.rm = TRUE) * 0.9,
#     label = paste0(
#       "log(y) = ", round(intercept_MERIS_2002_95_log, 3), " + ", 
#       round(slope_MERIS_2002_95_log, 7), " * x",
#       "\n", "p = ", ifelse(p_value_MERIS_2002_95_log < 0.001, "< 0.001", 
#                            format(p_value_MERIS_2002_95_log, digits = 3))
#     ),
#     hjust = 1, vjust = 1, size = 6
#   ) +
#   labs(title = "Évolution de la taille des panaches de la baie des Anges vu par le produit MERIS (ODATIS-MR) (échelle log)",
#        x = "Date",
#        y = "Aire des panaches (en km²)") +
#   theme_minimal() +
#   scale_x_date(date_breaks = "1 year", date_labels = "%Y")

ggplot(data = data_log_spm, aes(x = date, y = mean_spm)) +  # utilise data_log_spm
  geom_point(color = "red3", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE,
              formula = y ~ x,
              color = "darkslateblue", fill = "pink", alpha = 0.2) +
  scale_y_log10(labels = scales::label_comma()) +  # force l'échelle log sur smooth aussi
  annotate(
    "text",
    x = max(data_log_spm$date, na.rm = TRUE),
    y = max(data_log_spm$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "log(y) = ", round(intercept_MERIS_2002_95_log, 3), " + ",
      round(slope_MERIS_2002_95_log, 7), " * x",
      "\n", "p = ", ifelse(p_value_MERIS_2002_95_log < 0.001, "< 0.001",
                           format(p_value_MERIS_2002_95_log, digits = 3))
    ),
    hjust = 1, vjust = 1, size = 6
  ) +
  labs(title = "Évolution de la concentration en MES dans les panaches de la baie des Anges vu par le produit MERIS (ODATIS-MR) (échelle log)",
       x = "Date",
       y = "Concentration moyenne en SPM (en mg/m³)") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")



# comparison between liquid flow rate and panache extension

adjust_factors <- sec_axis_adjustement_factors(MERIS_2002_2012_spm_95_propre$aire_panache_km2, Y6442010_2002_2012$débit)

MERIS_2002_2012_spm_95_propre$scaled_aire_panache_km2 <- MERIS_2002_2012_spm_95_propre$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(data = Y6442010_2002_2012, 
             aes(x = date, y = débit, color = "Débit"), size = 0.5) +
  geom_point(data = MERIS_2002_2012_spm_95_propre, 
             aes(x = date, y = scaled_aire_panache_km2, color = "SPM"), size = 0.5) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "limegreen")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(title = "Évolution de l'aire des panaches et du débit du Var vu par le produit MERIS (ODATIS-MR)",
       x = "Date") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# runoff vs SPM concentration correlation ---------------------------------

Var_MERIS <- Y6442010_2002_2012 %>%
  select(date, débit) %>%
  inner_join(
    MERIS_2002_2012_spm_95_propre %>% select(date, mean_spm, aire_panache_km2),
    by = "date"
  )

cor.test(Var_MERIS$débit, Var_MERIS$aire_panache_km2, method = "spearman")





# scatter plot

# Fusionner les données
Var_MERIS <- Y6442010_2002_2012 %>% 
  select(date, débit) %>% 
  left_join(
    MERIS_2002_2012_spm_95 %>% select(date, aire_panache_km2, mean_spm),
    by = "date"
  )

ggplot(data = Var_MERIS, aes(x = débit, y = mean_spm)) +
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
  labs(x = "Débit (m³/s)", y = "Concentration moyenne en MES (en mg/m³)", title = "Débit liquide du Var contre la concentration moyenne en MES dans les panaches vue par MERIS (ODATIS-MR)") +
  theme_minimal()


