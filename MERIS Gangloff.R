# MERIS Gangloff
# 15/05/2026

# pathway : "Satellite analysis/MERIS Gangloff.Rdata

# This script will load MERIS satellite data
# It will compute a threshold to the 95 percentile as Gangloff et al., 2017
# Then, it will identify the turbid plume using this threshold
# It will calculate statistics between river discharge and turbid plume area

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)
library(slider)

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

load("data/MERIS/MERIS_2002_2012_spm_pixels.Rdata")
load("data/Hydro France/Y6442010_2002_2012.Rdata")

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

# turbid plume identification ---------------------------------------------

## zones emboîtées ----------------------------------------------------------

# définir les coordonnées de l'embouchure du Var (comme SEXTANT OC5 chla)
lon_embouchure <- 7.199082
lat_embouchure <- 43.654709

# définir des zones emboîtées de tailles croissantes autour de l'embouchure
rayons_km <- c(5, 10, 20, 40, 70, 100)  # en km

# Convertir en degrés
rayons_deg <- rayons_km / 111

# Calculer le percentile 95 pour chaque zone emboîtée
seuils_zones <- lapply(rayons_deg, function(r) {
  
  pixels_zone <- MERIS_2002_2012_spm_pixels |>
    filter(
      lon >= lon_embouchure - r & lon <= lon_embouchure + r,
      lat >= lat_embouchure - r & lat <= lat_embouchure + r
    )
  
  aire_km2 <- nrow(distinct(pixels_zone, lon, lat)) * aire_pixel_km2
  seuil    <- quantile(pixels_zone$`SPM-G-PO_mean`, 0.95, na.rm = TRUE)
  
  data.frame(rayon_km = r * 111, aire_km2 = aire_km2, seuil_95 = seuil)
})

seuils_zones_df <- bind_rows(seuils_zones)
print(seuils_zones_df)

# Visualiser le plateau
ggplot(seuils_zones_df, aes(x = aire_km2, y = seuil_95)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line() +
  geom_hline(yintercept = seuil_95, linetype = "dashed", color = "red") +
  labs(
    title = "Détermination du seuil de détection du panache turbide",
    subtitle = "Percentile 95 par zone emboîtée autour de l'embouchure",
    x = "Aire de la zone (km²)",
    y = "Percentile 95 des MES (g/m³)"
  ) +
  theme_bw()

# Le seuil retenu est la valeur du plateau (zone > ~5000 km²)
seuil_retenu <- seuils_zones_df |>
  filter(aire_km2 > 2465) |>
  summarise(seuil = mean(seuil_95)) |>
  pull(seuil)

cat("Seuil retenu :", seuil_retenu, "g/m³\n", na.rm = TRUE)
# 0.4910305

## ROPP --------------------------------------------------------------------

# Pixels où le panache est présent dans au moins 5% des images
n_images_total <- n_distinct(MERIS_2002_2012_spm_pixels$date)

ROPP <- MERIS_2002_2012_spm_pixels |>
  group_by(lon, lat) |>
  summarise(
    freq_above_seuil = sum(`SPM-G-PO_mean` >= seuil_retenu, na.rm = TRUE) / n_images_total,
    .groups = "drop"
  ) |>
  filter(freq_above_seuil >= 0.05)   # au moins 5% des images

cat("Nombre de pixels dans la ROPP :", nrow(ROPP), "\n")

## filtrer les images avec trop de pixels manquants ------------------------

# Garder seulement les images avec > 80% de pixels valides sur la ROPP
pixels_ROPP <- MERIS_2002_2012_spm_pixels |>
  semi_join(ROPP, by = c("lon", "lat"))

images_valides <- pixels_ROPP |>
  group_by(date) |>
  summarise(
    n_pixels_valides = sum(!is.na(`SPM-G-PO_mean`)),
    n_pixels_total   = n(),
    pct_valide       = n_pixels_valides / n_pixels_total,
    .groups = "drop"
  ) |>
  filter(pct_valide >= 0.80)

cat("Images valides :", nrow(images_valides), "/", n_distinct(MERIS_2002_2012_spm_pixels$date), "\n")

## stat panache ------------------------------------------------------------

# Statistiques du panache par jour
MERIS_panache_metrics <- MERIS_2002_2012_spm_pixels |>
  filter(date %in% images_valides$date) |>
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
  )

# climatologie ------------------------------------------------------------

# période de 10 ans
MERIS_panache_metrics_stat <- MERIS_panache_metrics |>
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

# climatologie annuelle
MERIS_panache_metrics_year <- MERIS_panache_metrics_stat |> 
  group_by(year) |> 
  summarise(mean_spm_year_clim = mean(mean_spm, na.rm = TRUE), 
            median_spm_year_clim = median(mean_spm, na.rm = TRUE),
            sd_spm_year_clim = sd(mean_spm, na.rm = TRUE),
            max_spm_year_clim = max(mean_spm, na.rm = TRUE),
            mean_panache_extension_year_clim =  mean(aire_panache_km2, na.rm = TRUE),
            median_panache_extension_year_clim = median(aire_panache_km2, na.rm = TRUE),
            max_panache_extension_year_clim = max(aire_panache_km2, na.rm = TRUE))

# climatologie mensuelle
MERIS_panache_metrics_month <- MERIS_panache_metrics_stat |> 
  group_by(month) |>
  summarise(mean_spm_month_clim = mean(mean_spm, na.rm = TRUE), 
            median_spm_month_clim = median(mean_spm, na.rm = TRUE),
            sd_spm_month_clim = sd(mean_spm, na.rm = TRUE),
            max_spm_month_clim = max(mean_spm, na.rm = TRUE),
            mean_panache_extension_month_clim =  mean(aire_panache_km2, na.rm = TRUE),
            median_panache_extension_month_clim = median(aire_panache_km2, na.rm = TRUE),
            max_panache_extension_month_clim = max(aire_panache_km2, na.rm = TRUE))

# climatologie journanlière
MERIS_panache_metrics_doy <- MERIS_panache_metrics_stat |>
  group_by(doy) |> 
  summarise(mean_spm_doy_clim = mean(mean_spm, na.rm = TRUE), 
            median_spm_doy_clim = median(mean_spm, na.rm = TRUE),
            sd_spm_doy_clim = sd(mean_spm, na.rm = TRUE),
            max_spm_doy_clim = max(mean_spm, na.rm = TRUE),
            mean_panache_extension_doy_clim =  mean(aire_panache_km2, na.rm = TRUE),
            median_panache_extension_doy_clim = median(aire_panache_km2, na.rm = TRUE),
            max_panache_extension_doy_clim = max(aire_panache_km2, na.rm = TRUE))

# plotting ----------------------------------------------------------------

## climatologie annuelle ---------------------------------------------------

# create a line plot of the yearly climatology of spm
ggplot(MERIS_panache_metrics_year, aes(x = year, y = mean_spm_year_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_spm_year_clim - sd_spm_year_clim,
      ymax = mean_spm_year_clim + sd_spm_year_clim
    ),
    fill = "red3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie annuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie annuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie annuelle" = "red3")
    # ) +
    # scale_x_continuous(
    #   breaks = 1:12,
    #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
    #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie annuelle de la concentration moyenne en MES (2002-2012)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (g.m"^{-3}*")"),
    color   = NULL,
    caption = "Source : ODATIS-MR | MERIS"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the monthly climatology of spm
ggplot(MERIS_panache_metrics_month, aes(x = month, y = mean_spm_month_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_spm_month_clim - sd_spm_month_clim,
      ymax = mean_spm_month_clim + sd_spm_month_clim
    ),
    fill = "red3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie mensuelle"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie mensuelle"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie mensuelle" = "red3")
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
               "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  ) +
  labs(
    title   = "Climatologie mensuelle de la concentration moyenne en MES (2002-2012)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (g.m"^{-3}*")"),
    color   = NULL,
    caption = "Source : ODATIS-MR | MERIS"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# create a line plot of the daily climatology of spm
ggplot(MERIS_panache_metrics_doy, aes(x = doy, y = mean_spm_doy_clim)) +
  geom_ribbon(
    aes(
      ymin = mean_spm_doy_clim - sd_spm_doy_clim,
      ymax = mean_spm_doy_clim + sd_spm_doy_clim
    ),
    fill = "red3", alpha = 0.2
  ) +
  geom_line(aes(color = "Climatologie journalière"), linewidth = 0.8) +
  geom_point(aes(color = "Climatologie journalière"), size = 2.5) +
  scale_color_manual(
    values = c("Climatologie journalière" = "red3")
  ) +
  # scale_x_continuous(
  #   breaks = 1:12,
  #   labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
  #              "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
  # ) +
  labs(
    title   = "Climatologie journalière de la concentration en MES (2002-2012)",
    x       = NULL,
    y       = expression("Concentration moyenne en MES (g.m"^{-3}*")"),
    color   = NULL,
    caption = "Source : ODATIS-MR | MERIS"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption     = element_text(size = 13, color = "grey50", hjust = 0),
    axis.title.y     = element_text(size = 13, margin = margin(r = 10)),
    axis.text        = element_text(size = 12, color = "grey30"),
    axis.ticks       = element_line(color = "grey70"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", linewidth = 0.5),
    legend.position  = "top",
    legend.text      = element_text(size = 11)
  )

# correlation river flow and plume extension ------------------------------

# on doit moyenner le débit sur les 3 jours précédents l'acquisition satellitaire

debit_3j <- Y6442010_2002_2012 |>
  arrange(date) |>
  mutate(
    # Moyenne glissante sur les 3 jours PRÉCÉDENTS (j-3, j-2, j-1)
    # align = "right" + window = 3 donne j-2, j-1, j → on veut décaler d'un jour
    debit_3j_mean = slider::slide_dbl(
      débit,
      mean,
      .before = 3,   # les 3 jours avant
      .after  = -1,  # exclut le jour j lui-même
      na.rm   = TRUE
    )
  )

# Jointure et nettoyage
MERIS_panache_metrics <- MERIS_panache_metrics |>
  inner_join(debit_3j |> select(date, débit, debit_3j_mean), by = "date") |>
  filter(aire_panache_km2 > 0) |>                        # filtre AVANT le log
  mutate(
    panache_log = log10(aire_panache_km2),
    debit_log   = log10(debit_3j_mean)                            # log-log comme Fig. 14
  )

## modèle log - log --------------------------------------------------------

# Modèle linéaire en log-log → relation puissance
modele_log <- lm(panache_log ~ debit_log, data = MERIS_panache_metrics)
r2     <- summary(modele_log)$r.squared
pente  <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

# Équation puissance : aire = 10^ordonnee * débit^pente
# à afficher dans le graphique
label_eq <- paste0(
  "Aire = ", round(10^ordonnee, 3), " × Q^", round(pente, 2),
  "\nR² = ", round(r2, 2)
)

# Graphique log-log avec droite de régression
ggplot(MERIS_panache_metrics, aes(x = debit_log, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,            # régression sur les axes log
              color = "black", se = FALSE, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  annotate("text", 
           x = min(MERIS_panache_metrics$debit_log, na.rm = TRUE) * 1.5,
           y = max(MERIS_panache_metrics$aire_panache_km2, na.rm = TRUE) * 0.7,
           label = label_eq, hjust = 0, size = 4) +
  labs(
    title = "Relation entre débit et extension du panache turbide",
    x = "Débit (m³/s)",
    y = "Aire du panache (km²)"
  ) +
  theme_minimal()

## Modèle semi-log : log10(aire) ~ débit --------------------------------------------------------

modele_semilog <- lm(panache_log ~ debit_3j_mean, data = MERIS_panache_metrics)
r2       <- summary(modele_semilog)$r.squared
pente    <- coef(modele_semilog)[2]
ordonnee <- coef(modele_semilog)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(MERIS_panache_metrics, aes(x = debit_3j_mean, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  # geom_smooth(method = "lm", formula = y ~ x,
  #             color = "black", se = FALSE, linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  annotate("text",
           x = max(MERIS_panache_metrics$debit_3j_mean, na.rm = TRUE) * 0.7,
           y = min(MERIS_panache_metrics$aire_panache_km2, na.rm = TRUE) * 3,
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


# corrélation MES mean et débit -------------------------------------------

modele <- lm(debit_3j_mean ~ mean_spm, data = MERIS_panache_metrics)
r2       <- summary(modele)$r.squared
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(MERIS_panache_metrics, aes(x = debit_3j_mean, y = mean_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(MERIS_panache_metrics$debit_3j_mean, na.rm = TRUE) * 0.7,
           y = min(MERIS_panache_metrics$mean_spm, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration en MES (g m"^{-3}*")")
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

# corrélation MES max et débit -------------------------------------------

modele <- lm(debit_3j_mean ~ max_spm, data = MERIS_panache_metrics)
r2       <- summary(modele)$r.squared
pente    <- coef(modele)[2]
ordonnee <- coef(modele)[1]

label_eq <- paste0(
  "log10(Aire) = ", round(ordonnee, 3), " + ", round(pente, 5), " × Q",
  "\nR² = ", round(r2, 2)
)

ggplot(MERIS_panache_metrics, aes(x = debit_3j_mean, y = max_spm)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(MERIS_panache_metrics$debit_3j_mean, na.rm = TRUE) * 0.7,
           y = min(MERIS_panache_metrics$mean_spm, na.rm = TRUE) * 3,
           label = label_eq, hjust = 0.5, size = 8, color = "grey20",
           family = "serif",
           fontface = "italic") +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = expression("Concentration maximale en MES (g m"^{-3}*")")
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




