# relation ship wind / sextant data
# 18/05/2026

# Setup ------------------------------------------------------------------

# Install rnaturalearthhires as necessry
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(sf)
library(rnaturalearth)
library(giscoR) # Hi-res coastlines
library(ggpmisc)
library(doParallel); registerDoParallel(cores = parallel::detectCores()-2)
library(ggpubr)  # Pour stat_cor()
library(scales)
library(ggspatial)
library(effectsize)
library(moments)
library(rstatix)
library(dplyr)

load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

# SEXTANT -----------------------------------------------------------------

# combien de valeurs négatives
sum(SEXTANT_1998_2025_spm_pixels$analysed_spim < 0, na.rm = TRUE)
# [1] 21080

# supprimer seulement les valeurs négatives
SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_pixels |>
  filter(analysed_spim >= 0 | is.na(analysed_spim))

SEXTANT_1998_2025_spm_clean <- SEXTANT_1998_2025_spm_clean |> 
  mutate(
    date = as.Date(date),  
    year = year(date),     
    month = month(date),
    doy = yday(date)         
  )

## pixel area --------------------------------------------------------------

### extraction des valeurs en degré -----------------------------------------

# Lire les attributs du fichier pour trouver la résolution
tidync("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")[["attribute"]]

# Ou inspecter les coordonnées lon/lat directement
nc <- tidync("~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/1998/01/01/19980101-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")

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

### calcul de l'aire --------------------------------------------------------

# Conversion en km (pour ~43°N, zone Méditerranée/Atlantique Sud de France)
lat_ref <- 43  

res_lon_km <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km <- res_lat * 111

cat("Résolution lon :", round(res_lon_km, 3), "km\n")
cat("Résolution lat :", round(res_lat_km, 3), "km\n")

# Aire d'un pixel
aire_pixel_km2 <- res_lon_km * res_lat_km
cat("Aire d'un pixel :", round(aire_pixel_km2, 4), "km²\n")

## SPM ---------------------------------------------------------------------

### define 95ème percentile -------------------------------------------------

# Calculer le 95ème percentile
seuil_95 <- quantile(SEXTANT_1998_2025_spm_pixels$analysed_spim, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.94 g/m³

# Stats du panache par jour
SEXTANT_1998_2025_spm_95 <- SEXTANT_1998_2025_spm_pixels |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(analysed_spim >= seuil_95, na.rm = TRUE),
    mean_spm = mean(analysed_spim[analysed_spim >= seuil_95], na.rm = TRUE),
    median_spm = median(analysed_spim[analysed_spim >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

## Gangloff SPM ------------------------------------------------------------

### zones emboîtées ----------------------------------------------------------

# définir les coordonnées de l'embouchure du Var (comme SEXTANT OC5 chla)
lon_embouchure <- 7.199082
lat_embouchure <- 43.654709

# définir des zones emboîtées de tailles croissantes autour de l'embouchure
rayons_km <- c(5, 10, 20, 40, 70, 100)  # en km

# Convertir en degrés
rayons_deg <- rayons_km / 111

# Calculer le percentile 95 pour chaque zone emboîtée
seuils_zones <- lapply(rayons_deg, function(r) {
  
  pixels_zone <- SEXTANT_1998_2025_spm_clean |>
    filter(
      lon >= lon_embouchure - r & lon <= lon_embouchure + r,
      lat >= lat_embouchure - r & lat <= lat_embouchure + r
    )
  
  aire_km2 <- nrow(distinct(pixels_zone, lon, lat)) * aire_pixel_km2
  seuil    <- quantile(pixels_zone$analysed_spim, 0.95, na.rm = TRUE)
  
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
  filter(aire_km2 > 2459) |>
  summarise(seuil = mean(seuil_95)) |>
  pull(seuil)

cat("Seuil retenu :", seuil_retenu, "g/m³\n", na.rm = TRUE)
# 0.94 g/m³

### ROPP --------------------------------------------------------------------

# Pixels où le panache est présent dans au moins 5% des images
n_images_total <- n_distinct(SEXTANT_1998_2025_spm_clean$date)

ROPP <- SEXTANT_1998_2025_spm_clean |>
  group_by(lon, lat) |>
  summarise(
    freq_above_seuil = sum(analysed_spim >= seuil_retenu, na.rm = TRUE) / n_images_total,
    .groups = "drop"
  ) |>
  filter(freq_above_seuil >= 0.05)   # au moins 5% des images

cat("Nombre de clean dans la ROPP :", nrow(ROPP), "\n")

### filtrer les images avec trop de clean manquants ------------------------

# Garder seulement les images avec > 80% de clean valides sur la ROPP
clean_ROPP <- SEXTANT_1998_2025_spm_clean |>
  semi_join(ROPP, by = c("lon", "lat"))

images_valides <- clean_ROPP |>
  group_by(date) |>
  summarise(
    n_clean_valides = sum(!is.na(analysed_spim)),
    n_clean_total   = n(),
    pct_valide       = n_clean_valides / n_clean_total,
    .groups = "drop"
  ) |>
  filter(pct_valide >= 0.80)

cat("Images valides :", nrow(images_valides), "/", n_distinct(SEXTANT_1998_2025_spm_clean$date), "\n")


# Diagnostics à faire absolument
cat("Seuil retenu :", seuil_retenu, "\n")
cat("Nombre de pixels dans la ROPP :", nrow(ROPP), "\n")
cat("Images valides :", nrow(images_valides), "\n")
cat("Jours avec panache > 0 :", sum(SEXTANT_panache_metrics$aire_panache_km2 > 0), "\n")

# Distribution des aires de panache
summary(SEXTANT_panache_metrics$aire_panache_km2)
hist(SEXTANT_panache_metrics$aire_panache_km2, breaks = 50)

# Distribution du débit 3j associé
summary(SEXTANT_panache_metrics$debit_3j_mean)

### stat panache ------------------------------------------------------------

# Statistiques du panache par jour
SEXTANT_panache_metrics <- SEXTANT_1998_2025_spm_clean |>
  filter(date %in% images_valides$date) |>
  semi_join(ROPP, by = c("lon", "lat")) |>
  group_by(date) |>
  summarise(
    # Aire d'extension
    pixel_count = sum(analysed_spim >= seuil_retenu, na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    # Concentrations
    mean_spm         = mean(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    max_spm          = max(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    median_spm       = median(analysed_spim[analysed_spim >= seuil_retenu], na.rm = TRUE),
    # Points extrêmes
    lat_sud          = min(lat[analysed_spim >= seuil_retenu], na.rm = TRUE),
    lon_ouest        = min(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    lon_est          = max(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    # Centroïde
    centroid_lon     = mean(lon[analysed_spim >= seuil_retenu], na.rm = TRUE),
    centroid_lat     = mean(lat[analysed_spim >= seuil_retenu], na.rm = TRUE),
    .groups = "drop"
  )

## plotting ----------------------------------------------------------------

### correlation river flow and plume extension ------------------------------

debit_lags <- Y6442010_depuis_2000 |>
  arrange(date) |>
  select(date, débit) |> 
  mutate(
    debit_j1     = lag(débit, 1),             # j-1 seulement
    debit_2j     = (lag(débit, 1) + lag(débit, 2)) / 2,       # moyenne j-1, j-2
    debit_3j     = (lag(débit, 1) + lag(débit, 2) + lag(débit, 3)) / 3
  )

# Jointure et nettoyage
SEXTANT_panache_metrics <- SEXTANT_panache_metrics |>
  inner_join(debit_lags |> select(date, débit, debit_j1, debit_2j, debit_3j), by = "date") |>
  filter(
    aire_panache_km2 > 0,
    débit > 0,
    is.finite(débit)
  ) |>
  mutate(
    panache_log  = log10(aire_panache_km2),
    debit_j0_log = log10(débit),
    debit_j1_log = log10(debit_j1),
    debit_2j_log = log10(debit_2j),
    debit_3j_log = log10(debit_3j)
  )

# # Tout est déjà dans SEXTANT_panache_metrics, pas besoin de rejoindre !
SEXTANT_panache_metrics |>
  filter(aire_panache_km2 > 0) |>
  summarise(
    r_j0 = cor(log10(aire_panache_km2), debit_j0_log, use = "complete.obs", method = "spearman"),
    r_j1 = cor(log10(aire_panache_km2), debit_j1_log, use = "complete.obs", method = "spearman"),
    r_2j = cor(log10(aire_panache_km2), debit_2j_log, use = "complete.obs", method = "spearman"),
    r_3j = cor(log10(aire_panache_km2), debit_3j_log, use = "complete.obs", method = "spearman")
  )

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ debit_j0_log, data = SEXTANT_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1_log, data = SEXTANT_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j_log, data = SEXTANT_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j_log, data = SEXTANT_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

# Ensuite compare les trois modèles
modele_j0 <- lm(panache_log ~ débit, data = SEXTANT_panache_metrics)
modele_j1 <- lm(panache_log ~ debit_j1, data = SEXTANT_panache_metrics)
modele_2j <- lm(panache_log ~ debit_2j, data = SEXTANT_panache_metrics)
modele_3j <- lm(panache_log ~ debit_3j, data = SEXTANT_panache_metrics)

tibble(
  lag       = c("j-1", "j-1 à j-2", "j-1 à j-3"),
  R2        = c(summary(modele_j1)$r.squared,
                summary(modele_2j)$r.squared,
                summary(modele_3j)$r.squared)
)

### modèle log - log --------------------------------------------------------

# Modèle linéaire en log-log → relation puissance
modele_log <- lm(panache_log ~ debit_j0_log, data = SEXTANT_panache_metrics)
r2     <- summary(modele_log)$r.squared
print(r2)
pente  <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

# Équation puissance : aire = 10^ordonnee * débit^pente
# à afficher dans le graphique
label_eq <- paste0(
  "Aire = ", round(10^ordonnee, 3), " × Q^", round(pente, 2),
  "\nR² = ", round(r2, 2)
)

# Graphique log-log avec droite de régression
ggplot(SEXTANT_panache_metrics, aes(x = débit, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,            # régression sur les axes log
              color = "black", se = FALSE, linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text", 
           x = min(SEXTANT_panache_metrics$debit_log, na.rm = TRUE) * 1.5,
           y = max(SEXTANT_panache_metrics$aire_panache_km2, na.rm = TRUE) * 0.7,
           label = label_eq, hjust = 1.1, vjust = 2, size = 8, color = "grey20",
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

# WIND --------------------------------------------------------------------

Wind_T <- read.csv("~/Vent/data/Q_06_previous-1950-2024_RR-T-Vent.csv", 
                   header = TRUE, sep = ";")

Wind_T <- Wind_T |> 
  filter(NUM_POSTE == "6088001") |> 
  select("LAT", "LON", "NUM_POSTE", "FFM", "DXY", "HXI", "RR", "TM")

Wind_T <- Wind_T |> 
  mutate(date = seq(as.Date("1950-01-01"), as.Date("2024-12-31"), by = "day"))

Wind_T <- Wind_T |> 
  filter(date >= "1998-01-01")

Wind_2015_2024 <- Wind_T |> 
  filter(date >= as.Date("2015-01-01"), date <= as.Date("2024-12-31"))

## separate wind -----------------------------------------------------------

# West
West <- Wind_T |> 
  filter(DXY >= 255, DXY <= 285)
West <- West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# East
East <- Wind_T |> 
  filter(DXY >= 75, DXY <= 105)
East <- East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# North
North <- Wind_T |> 
  filter(DXY >= 345 | DXY <= 15)
North <- North %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Nord Est
North_East <- Wind_T |> 
  filter(DXY >= 30, DXY <= 60)
North_East <- North_East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# North West
North_West <- Wind_T |> 
  filter(DXY >= 300, DXY <= 330)
North_West <- North_West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# South
South <- Wind_T |> 
  filter(DXY >= 165, DXY <= 195)
South <- South %>%
  complete(date = seq(min(date), max(date), by = "day"))

# Sud Ouest
South_West <- Wind_T |> 
  filter(DXY >= 210, DXY <= 240)
South_West <- South_West %>%
  complete(date = seq(min(date), max(date), by = "day"))

# South East
South_East <- Wind_T |> 
  filter(DXY >= 120, DXY <= 150)
South_East <- South_East %>%
  complete(date = seq(min(date), max(date), by = "day"))

# combien chaque part de vent fait de mon data set de base 

n_total <- nrow(Wind_T)

vent_summary <- Wind_T |>
  mutate(
    secteur = case_when(
      DXY >= 345 | DXY <= 15             ~ "North",
      DXY >= 30  & DXY <= 60             ~ "North_East",
      DXY >= 75  & DXY <= 105            ~ "East",
      DXY >= 120 & DXY <= 150            ~ "South_East",
      DXY >= 165 & DXY <= 195            ~ "South",
      DXY >= 210 & DXY <= 240            ~ "South_West",
      DXY >= 255 & DXY <= 285            ~ "West",
      DXY >= 300 & DXY <= 330            ~ "North_West",
      TRUE                               ~ "Non classifié"
    )
  ) |>
  count(secteur) |>
  mutate(
    pct        = round(n / n_total * 100, 1),
    cumul_pct  = round(cumsum(pct), 1)
  ) |>
  arrange(desc(n))

print(vent_summary)
cat("\nTotal jours Wind_T :", n_total, "\n")
cat("Jours classifiés   :", sum(vent_summary$n[vent_summary$secteur != "Non classifié"]),
    "(", round(sum(vent_summary$pct[vent_summary$secteur != "Non classifié"]), 1), "%)\n")
cat("Jours non classifiés (zones grises entre secteurs) :",
    vent_summary$n[vent_summary$secteur == "Non classifié"], "\n")

# relationship ------------------------------------------------------------

Wind_classified <- Wind_T |>
  mutate(
    wind_type = case_when(
      DXY >= 300 & DXY <= 330 ~ "Nord-Ouest",  
      DXY >= 75  & DXY <= 105 ~ "Est",          
      TRUE ~ NA_character_                       
    )
  ) |>
  filter(!is.na(wind_type))                      

# Vérifier la répartition
Wind_classified |> 
  count(wind_type) |> 
  mutate(
    pct_secteur  = round(n / sum(n) * 100, 1),         
    pct_total    = round(n / nrow(Wind_T) * 100, 1)   
  )

# ── 2. Jointure avec les métriques du panache ──────────────────────────────

panache_vent <- SEXTANT_panache_metrics |>
  inner_join(
    Wind_classified |> select(date, FFM, DXY, wind_type),
    by = "date"
  ) |>
  filter(wind_type != "Autre")   # garder seulement les deux directions principales

# ── 3. Figure 8b : distance sud vs vitesse du vent ────────────────────────

# Il te faut la distance sud = distance entre l'embouchure et le point le plus au sud
lon_embouchure <- 7.199082
lat_embouchure <- 43.654709

lat_ref <- 43

panache_vent <- panache_vent |>
  mutate(
    dist_sud_km   = (lat_embouchure - lat_sud)   * 111,
    dist_ouest_km = (lon_embouchure - lon_ouest)  * 111 * cos(lat_embouchure * pi / 180),
    dist_est_km   = (lon_est - lon_embouchure)    * 111 * cos(lat_embouchure * pi / 180)
  )

# Sanity check
panache_vent |>
  group_by(wind_type) |>
  summarise(
    across(c(dist_sud_km, dist_ouest_km, dist_est_km),
           list(med = ~round(median(., na.rm=TRUE), 1),
                max = ~round(max(.,    na.rm=TRUE), 1)),
           .names = "{.col}_{.fn}")
  ) |> 
  print(width =  Inf)

# graphique propre
# Calculer les médianes par groupe
med_sud <- panache_vent |>
  group_by(wind_type) |>
  summarise(med = median(dist_sud_km, na.rm = TRUE))

med_ouest <- panache_vent |>
  group_by(wind_type) |>
  summarise(med = median(dist_ouest_km, na.rm = TRUE))

# Wilcoxon
w_sud   <- wilcox.test(dist_sud_km   ~ wind_type, data = panache_vent)
w_ouest <- wilcox.test(dist_ouest_km ~ wind_type, data = panache_vent)

# d de Cohen
d_sud   <- cohens_d(dist_sud_km   ~ wind_type, data = panache_vent)
d_ouest <- cohens_d(dist_ouest_km ~ wind_type, data = panache_vent)

p1 <- ggplot(panache_vent |> filter(!is.na(dist_sud_km)),
             aes(x = FFM, y = dist_sud_km, color = wind_type, shape = wind_type)) +
  geom_point(alpha = 0.5, size = 2) +
  # Lignes de médiane
  geom_hline(data = med_sud,
             aes(yintercept = med, color = wind_type),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Nord-Ouest" = "steelblue", "Est" = "tomato"),
                     name = "Direction du vent") +
  scale_shape_manual(values = c("Nord-Ouest" = 16, "Est" = 17),
                     name = "Direction du vent") +
  coord_cartesian(xlim = c(0, 10)) +
  # Annotation Wilcoxon + Cohen
  annotate("text",
           x = 9.5, y = max(panache_vent$dist_sud_km, na.rm = TRUE) * 0.95,
           hjust = 1, vjust = 1, size = 8, fontface = "italic", family = "serif",
           color = "grey20",
           label = paste0("Wilcoxon : p < 2.2×10⁻¹⁶",
                          "\nd de Cohen = ", round(d_sud$effsize, 3), " (petit)",
                          "\nMéd. Est = ", round(med_sud$med[med_sud$wind_type == "Est"], 1), " km",
                          "\nMéd. NO = ",  round(med_sud$med[med_sud$wind_type == "Nord-Ouest"], 1), " km")) +
  labs(
    x = expression("Vitesse du vent (m.s"^{-1}*")"),
    y = "Distance d'extension sud (km)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title     = element_text(face = "bold", family = "serif"),
    legend.text      = element_text(family = "serif"),
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(panache_vent |> filter(!is.na(dist_ouest_km)),
             aes(x = FFM, y = dist_ouest_km, color = wind_type, shape = wind_type)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(data = med_ouest,
             aes(yintercept = med, color = wind_type),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Nord-Ouest" = "steelblue", "Est" = "tomato"),
                     name = "Direction du vent") +
  scale_shape_manual(values = c("Nord-Ouest" = 16, "Est" = 17),
                     name = "Direction du vent") +
  coord_cartesian(xlim = c(0, 10)) +
  annotate("text",
           x = 9.5, y = max(panache_vent$dist_ouest_km, na.rm = TRUE) * 0.95,
           hjust = 1, vjust = 1, size = 8, fontface = "italic", family = "serif",
           color = "grey20",
           label = paste0("Wilcoxon : p < 2.2×10⁻¹⁶",
                          "\nd de Cohen = ", round(d_ouest$effsize, 3), " (petit)",
                          "\nMéd. Est = ", round(med_ouest$med[med_ouest$wind_type == "Est"], 1), " km",
                          "\nMéd. NO = ",  round(med_ouest$med[med_ouest$wind_type == "Nord-Ouest"], 1), " km")) +
  labs(
    x = expression("Vitesse du vent (m.s"^{-1}*")"),
    y = "Distance d'extension ouest (km)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title     = element_text(face = "bold", family = "serif"),
    legend.text      = element_text(family = "serif"),
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank()
  )

(p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Source : Sextant OC5",
    theme   = theme(
      plot.caption = element_text(size = 14, color = "grey50", hjust = 0)
    )
  ) &
  theme(legend.position = "top")







# Calculer les médianes par groupe
med_sud <- panache_vent |>
  group_by(wind_type) |>
  summarise(med = median(dist_sud_km, na.rm = TRUE))

# Ajouter sur p1
p1 + geom_hline(data = med_sud, 
                aes(yintercept = med, color = wind_type),
                linetype = "dashed", linewidth = 0.8)

p1 <- ggplot(panache_vent |> filter(!is.na(dist_sud_km)),
             aes(x = FFM, y = dist_sud_km, color = wind_type, shape = wind_type)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(
    values = c("Nord-Ouest" = "steelblue", "Est" = "tomato"),
    name = "Direction du vent"
  ) +
  scale_shape_manual(
    values = c("Nord-Ouest" = 16, "Est" = 17),
    name = "Direction du vent"
  ) +
  coord_cartesian(xlim = c(0, 10)) +   # ← coupe l'axe x à 10
  labs(
    x = expression("Vitesse du vent (m.s"^{-1}*")"),
    y = "Distance d'extension sud (km)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title      = element_text(face = "bold", family = "serif"),
    legend.text       = element_text(family = "serif"),
    axis.title        = element_text(face = "bold", family = "serif"),
    axis.text         = element_text(color = "grey30", family = "serif"),
    panel.grid.minor  = element_blank()
  )

p2 <- ggplot(panache_vent |> filter(!is.na(dist_ouest_km)),
             aes(x = FFM, y = dist_ouest_km, color = wind_type, shape = wind_type)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(
    values = c("Nord-Ouest" = "steelblue", "Est" = "tomato"),
    name = "Direction du vent"
  ) +
  scale_shape_manual(
    values = c("Nord-Ouest" = 16, "Est" = 17),
    name = "Direction du vent"
  ) +
  coord_cartesian(xlim = c(0, 10)) +   # ← coupe l'axe x à 10
  labs(
    x = expression("Vitesse du vent (m.s"^{-1}*")"),
    y = "Distance d'extension ouest (km)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title      = element_text(face = "bold", family = "serif"),
    legend.text       = element_text(family = "serif"),
    axis.title        = element_text(face = "bold", family = "serif"),
    axis.text         = element_text(color = "grey30", family = "serif"),
    panel.grid.minor  = element_blank()
  )

# patchwork côte à côte
(p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Source : Sextant OC5",
    theme   = theme(
      plot.caption = element_text(size = 14, color = "grey50", hjust = 0)
    )
  ) &
  theme(legend.position = "top")

Wind_new <- panache_vent |> 
  dplyr::select(FFM, wind_type, dist_sud_km, dist_ouest_km, dist_est_km) |> 
  pivot_longer(cols = c(dist_sud_km, dist_ouest_km, dist_est_km), values_to = "dist_km", 
              names_to = "panache_extention") |> 
  filter(!is.infinite(dist_km))

panache_vent <- panache_vent |> 
  filter(!is.infinite(dist_sud_km))

# ANOVA
aov(formula = dist_km ~ wind_type * panache_extention, data = Wind_new)
summary(aov(formula = dist_km ~ wind_type * panache_extention, data = Wind_new))
TukeyHSD(aov(formula = dist_km ~ wind_type * panache_extention, data = Wind_new)
)

sum(is.na(Wind_new))
# 0

# linear model
summary(lm(dist_sud_km ~ FFM, filter(panache_vent, wind_type == "Nord-Ouest")))
summary(lm(dist_est_km ~ FFM, filter(panache_vent, wind_type == "Nord-Ouest")))      
summary(lm(dist_ouest_km ~ FFM, filter(panache_vent, wind_type == "Nord-Ouest")))

summary(lm(dist_sud_km ~ FFM, filter(panache_vent, wind_type == "Est")))
summary(lm(dist_est_km ~ FFM, filter(panache_vent, wind_type == "Est")))      
summary(lm(dist_ouest_km ~ FFM, filter(panache_vent, wind_type == "Est")))

# ── 4. (Bonus) Figure 8a adaptée : aire du panache vs débit coloré par vent ─

ggplot(panache_vent,
       aes(x = débit, y = aire_panache_km2, color = wind_type, shape = wind_type)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(
    values = c("Nord-Ouest" = "steelblue", "Est" = "tomato"),
    name = "Direction du vent"
  ) +
  scale_shape_manual(
    values = c("Nord-Ouest" = 16, "Est" = 17),
    name = "Direction du vent"
  ) +
  labs(
    x = expression("Débit (m"^{3}*".s"^{-1}*")"),
    y = "Aire du panache (km²)",
    title = "Panache turbide du Var — réponse au débit et au vent"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "top",
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    panel.grid.minor = element_blank()
  )

# quantification ----------------------------------------------------------

# ── 1. Test de Wilcoxon (non-paramétrique) ────────────────────────────────
# Est-ce que les deux distributions sont significativement différentes ?

wilcox.test(dist_sud_km ~ wind_type, data = panache_vent)
# oui
wilcox.test(dist_ouest_km ~ wind_type, data = panache_vent)
# oui

# ── 2. Stats descriptives par groupe ─────────────────────────────────────

panache_vent |>
  group_by(wind_type) |>
  summarise(
    n           = n(),
    med_km_sud      = round(median(dist_sud_km, na.rm = TRUE), 1),
    mean_km_sud     = round(mean(dist_sud_km,   na.rm = TRUE), 1),
    sd_km_sud       = round(sd(dist_sud_km,     na.rm = TRUE), 1),
    q25_km_sud      = round(quantile(dist_sud_km, 0.25, na.rm = TRUE), 1),
    q75_km_sud      = round(quantile(dist_sud_km, 0.75, na.rm = TRUE), 1),
    max_km_sud      = round(max(dist_sud_km,    na.rm = TRUE), 1),
    med_km_ouest      = round(median(dist_ouest_km, na.rm = TRUE), 1),
    mean_km_ouest     = round(mean(dist_ouest_km,   na.rm = TRUE), 1),
    sd_km_ouest       = round(sd(dist_ouest_km,     na.rm = TRUE), 1),
    q25_km_ouest      = round(quantile(dist_ouest_km, 0.25, na.rm = TRUE), 1),
    q75_km_ouest      = round(quantile(dist_ouest_km, 0.75, na.rm = TRUE), 1),
    max_km_ouest      = round(max(dist_ouest_km,    na.rm = TRUE), 1)
  ) |> 
  print(width = Inf)

# les vents d'Est semblent étendrent plus le panache au sud mais pas de beaucoup
# Est, médiane = 13.8 et Nord Ouest, médiane = 11.6
# Pour l'extension à l'Ouest, les vents d'Est semblent étendrent de manière plus
# importante vers l'Ouest
# Est, médiane = 13.3 et Nord Ouest, médiane = 6

# ── 3. Taille d'effet : d de Cohen ───────────────────────────────────────
# Quantifie l'ampleur réelle de la différence, indépendamment de n
# < 0.2 = négligeable, 0.2-0.5 = petit, 0.5-0.8 = moyen, > 0.8 = grand

cohens_d(dist_sud_km ~ wind_type, data = panache_vent)
cohens_d(dist_ouest_km ~ wind_type, data = panache_vent)
# les différences sont petites

# ── 4. Régression avec wind_type + FFM ───────────────────────────────────
# La direction ajoute-t-elle de l'information au-delà de la vitesse seule ?

mod_vitesse_seule  <- lm(dist_sud_km ~ FFM,                   data = panache_vent)
mod_direction      <- lm(dist_sud_km ~ FFM + wind_type,       data = panache_vent)
mod_interaction    <- lm(dist_sud_km ~ FFM * wind_type,       data = panache_vent)

# Comparer les modèles
tibble(
  modele   = c("vitesse seule", "vitesse + direction", "vitesse × direction"),
  R2       = c(summary(mod_vitesse_seule)$r.squared,
               summary(mod_direction)$r.squared,
               summary(mod_interaction)$r.squared) |> round(3),
  AIC      = c(AIC(mod_vitesse_seule),
               AIC(mod_direction),
               AIC(mod_interaction))              |> round(1)
) |> print()

# Test ANOVA pour savoir si direction améliore significativement le modèle
anova(mod_vitesse_seule, mod_direction)
anova(mod_direction,     mod_interaction)

# ── 5. Visualisation avec boxplot pour mieux voir la différence ───────────

ggplot(panache_vent, aes(x = wind_type, y = dist_ouest_km,
                         fill = wind_type, color = wind_type)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.15, size = 1.5) +
  scale_fill_manual(values  = c("Nord-Ouest" = "steelblue", "Est" = "tomato")) +
  scale_color_manual(values = c("Nord-Ouest" = "steelblue", "Est" = "tomato")) +
  # Afficher la p-value du Wilcoxon directement sur le graph
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y = max(panache_vent$dist_ouest_km, na.rm=TRUE)) +
  labs(
    x        = NULL,
    y        = "Distance d'extension ouest du panache (km)",
    title    = "Extension du panache selon la direction du vent",
    subtitle = "Test de Wilcoxon — d de Cohen"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    axis.title       = element_text(face = "bold", family = "serif"),
    axis.text        = element_text(color = "grey30", family = "serif"),
    plot.title       = element_text(face = "bold", size = 16,
                                    hjust = 0.5, family = "serif"),
    plot.subtitle    = element_text(size = 13, hjust = 0.5,
                                    color = "grey50", family = "serif"),
    panel.grid.minor = element_blank()
  )


# on test la distribution des médianes

# Visualiser les deux distributions
ggplot(panache_vent, aes(x = dist_ouest_km, fill = wind_type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("Nord-Ouest" = "steelblue", "Est" = "tomato")) +
  theme_bw()

# Comparer forme et dispersion
panache_vent |>
  group_by(wind_type) |>
  summarise(
    sd       = round(sd(dist_ouest_km,                    na.rm = TRUE), 1),
    skewness = round(moments::skewness(dist_ouest_km,     na.rm = TRUE), 2),
    IQR      = round(IQR(dist_ouest_km,                   na.rm = TRUE), 1)
  )




