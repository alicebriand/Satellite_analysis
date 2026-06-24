# Comparaison satellites pixels

# 18/05/2026

# On veut comparer le nombre de pixels valides selon les différents produits 
# satellitaires

library(bit64)
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)

# load functions ----------------------------------------------------------

# to load sextant data
load_SEXTANT_spm_pixels <- function(file_name, lon_range, lat_range){
  
  # Find the date
  sextant_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
  
  # The necessary code
  sextant_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = sextant_one_date) |> 
    dplyr::select(lon, lat, date, analysed_spim)  # tous les pixels, sans filtre
  
  # Exit
  return(sextant_one)
}

# to load sextant data
load_SEXTANT_chl_pixels <- function(file_name, lon_range, lat_range){
  
  # Find the date
  sextant_one_date <- as.Date(tidync(file_name)[["attribute"]][["value"]][["start_date"]])
  
  # The necessary code
  sextant_one <- tidync(file_name) |> 
    hyper_filter(lon = lon >= lon_range[1] & lon <= lon_range[2],
                 lat = lat >= lat_range[1] & lat <= lat_range[2]) |> 
    hyper_tibble() |> 
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat),
           date = sextant_one_date) |> 
    dplyr::select(lon, lat, date, analysed_chl_a)  # tous les pixels, sans filtre
  
  # Exit
  return(sextant_one)
}

# directions --------------------------------------------------------------

# on choisit de travailler avec l'année 2024
# Pour cela on va compter le nombre de NA pour cette année là

# On charge tous les fichiers Sextant
load("data/SEXTANT/SPM/SEXTANT_1998_2025_spm_pixels.RData")

# On charge les fichiers MODIS
load("data/MODIS/SPM/MODIS_2002_2024_spm_pixels.RData")

# On charge les fichiers OLCI
load("data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

# NA count ----------------------------------------------------------------

detach("package:bit64", unload = TRUE)
detach("package:bit", unload = TRUE)

dates_2024 <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")

# Sextant
SEXTANT_2024 <- SEXTANT_1998_2025_spm_pixels |>
  mutate(annee = as.integer(format(date, "%Y"))) |>
  filter(annee == 2024) |>
  select(-annee)

sextant_valides <- sum(!is.na(SEXTANT_2024$analysed_spim))
SEXTANT_non_valides <- sum(is.na(SEXTANT_2024$analysed_spim))

# MODIS
MODIS_2024 <- MODIS_2002_2024_spm_pixels |> 
  filter(date >= ("2024-01-01"), date <= ("2024-12-31"))

MODIS_valides <- sum(!is.na(MODIS_2024$`SPM-G-NS_mean`))
MODIS_non_valides <- sum(is.na(MODIS_2024$`SPM-G-NS_mean`))

# OLCI
OLCI_2024 <- OLCI_2016_2024_spm_pixels |> 
  filter(date >= ("2024-01-01"), date <= ("2024-12-31"))

OLCI_valides <- sum(!is.na(OLCI_2024$`SPM-G-PO_mean`))
OLCI_non_valides <- sum(is.na(OLCI_2024$`SPM-G-PO_mean`))

# Sextant
SEXTANT_2023 <- SEXTANT_1998_2025_spm_pixels |> 
  filter(date >= ("2023-01-01"), date <= ("2023-12-31"))

sextant_valides <- sum(!is.na(SEXTANT_2023$analysed_spim))
SEXTANT_non_valides <- sum(is.na(SEXTANT_2023$analysed_spim))

# MODIS
MODIS_2023 <- MODIS_2002_2024_spm_pixels |> 
  filter(date >= ("2023-01-01"), date <= ("2023-12-31"))

MODIS_valides <- sum(!is.na(MODIS_2023$`SPM-G-NS_mean`))
MODIS_non_valides <- sum(is.na(MODIS_2023$`SPM-G-NS_mean`))

# OLCI
OLCI_2023 <- OLCI_2016_2024_spm_pixels |> 
  filter(date >= ("2023-01-01"), date <= ("2023-12-31"))

OLCI_valides <- sum(!is.na(OLCI_2023$`SPM-G-PO_mean`))
OLCI_non_valides <- sum(is.na(OLCI_2023$`SPM-G-PO_mean`))

SEXTANT_stats_complet <- SEXTANT_stats |>
  complete(date = dates_2024, fill = list(pixels_valides = 0, pixels_totaux = 1820, couverture = NA, produit = "SEXTANT")) |>
  mutate(produit = "SEXTANT")

MODIS_stats_complet <- MODIS_stats |>
  complete(date = dates_2024, fill = list(pixels_valides = 0, pixels_totaux = 25976, couverture = NA, produit = "MODIS")) |>
  mutate(produit = "MODIS")

OLCI_stats_complet <- OLCI_stats |>
  complete(date = dates_2024, fill = list(pixels_valides = 0, couverture = NA, produit = "OLCI")) |>
  mutate(produit = "OLCI")

# visualisation -----------------------------------------------------------

# chat --------------------------------------------------------------------

SEXTANT_stats <- SEXTANT_1998_2025_spm_pixels |>
  filter(date >= "2024-01-01",
         date <= "2024-12-31") |>
  group_by(date) |>
  summarise(
    pixels_valides = sum(!is.na(analysed_spim)),
    pixels_totaux = 1820,
    couverture = pixels_valides / pixels_totaux,
    produit = "SEXTANT"
  )

MODIS_stats <- MODIS_2002_2024_spm_pixels |>
  filter(date >= "2024-01-01",
         date <= "2024-12-31") |>
  group_by(date) |>
  summarise(
    pixels_valides = sum(!is.na(`SPM-G-NS_mean`)),
    pixels_totaux = 25976,
    couverture = pixels_valides / pixels_totaux,
    produit = "MODIS"
  )

OLCI_stats <- OLCI_2016_2024_spm_pixels |>
  filter(date >= "2024-01-01",
         date <= "2024-12-31") |>
  group_by(date) |>
  summarise(
    pixels_valides = sum(!is.na(`SPM-G-PO_mean`)),
    pixels_totaux = n(),
    couverture = pixels_valides / pixels_totaux,
    produit = "OLCI"
  )

comparaison <- bind_rows(
  SEXTANT_stats,
  MODIS_stats,
  OLCI_stats
)

# boxplot
ggplot(comparaison,
       aes(x = produit,
           y = couverture,
           fill = produit)) +
  
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_fill_manual(values = c(
    "SEXTANT" = "#1b9e77",
    "MODIS"   = "#d95f02",
    "OLCI"    = "#7570b3"
  )) +
  
  labs(
    title = "Couverture valide des produits satellitaires en 2024",
    x = "Produit satellitaire",
    y = "Couverture valide (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggplot(comparaison,
       aes(x = date,
           y = couverture,
           color = produit)) +
  
  geom_line(alpha = 0.7) +
  
  geom_point(size = 1.5,
             alpha = 0.7) +
  
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  
  scale_color_manual(values = c(
    "SEXTANT" = "#1b9e77",
    "MODIS"   = "#d95f02",
    "OLCI"    = "#7570b3"
  )) +
  
  labs(
    title = "Évolution temporelle de la couverture valide en 2024",
    x = "Date",
    y = "Couverture valide (%)",
    color = "Produit"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold")
  )

ggplot(comparaison,
       aes(date, couverture, color = produit)) +
  
  geom_point(alpha = 0.3) +
  
  geom_smooth(se = FALSE, linewidth = 1.2) +
  
  scale_y_continuous(labels = scales::percent) +
  
  theme_minimal()

stats_resume <- comparaison |>
  group_by(produit) |>
  summarise(
    mediane = median(couverture),
    n_jours_avec_data = sum(couverture > 0)
  )

stats_resume <- comparaison |>
  group_by(produit) |>
  summarise(
    mediane = median(couverture),
    n_jours_avec_data = sum(couverture > 0)
  )

ggplot(comparaison,
       aes(x = produit,
           y = couverture,
           fill = produit)) +
  
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  geom_text(
    data = stats_resume,
    aes(x = produit,
        y = -0.05,
        label = paste0("n = ", n_jours_avec_data, " jours")),
    inherit.aes = FALSE,
    size = 3.5,
    color = "grey30"
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(-0.08, 1)  # laisser de la place pour le texte en bas
  ) +
  
  scale_fill_manual(values = c(
    "SEXTANT OC5" = "#1b9e77",
    "MODIS"   = "#d95f02",
    "OLCI"    = "#7570b3"
  )) +
  
  labs(
    title = "Couverture valide des produits satellitaires en 2024",
    subtitle = "Chaque point représente un jour de l'année",
    x = "Produit satellitaire",
    y = "Couverture valide (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )



# Assembler
comparaison_complet <- bind_rows(
  SEXTANT_stats_complet,
  MODIS_stats_complet,
  OLCI_stats_complet
)

# Corriger SEXTANT et remplacer 0 par NA
comparaison_plot <- comparaison_complet |>
  mutate(
    couverture = ifelse(produit == "SEXTANT" & !is.na(couverture) & couverture > 0, 1, couverture),
    couverture = ifelse(!is.na(couverture) & couverture == 0, NA, couverture)
  )

comparaison_plot |>
  mutate(
    mois = lubridate::month(date, label = TRUE),
    jour = lubridate::day(date)
  ) |>
  ggplot(aes(x = jour, y = mois, fill = couverture)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~produit, ncol = 1) +
  scale_fill_gradientn(
    colors   = c("#ffffcc", "#41b6c4", "#0c2c84"),
    labels   = scales::percent,
    na.value = "grey85",
    limits   = c(0, 1)
  ) +
  labs(
    title    = "Couverture journalière par produit en 2024",
    fill     = "Couverture",
    x        = "Jour",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey50", size = 10),
    strip.text    = element_text(face = "bold")
  )

# comparaison_plot <- comparaison |>
#   mutate(
#     couverture = ifelse(produit == "SEXTANT" & couverture > 0, 1, couverture)
#   )
# 
# comparaison_plot |>
#   mutate(
#     mois = lubridate::month(date, label = TRUE),
#     jour = lubridate::day(date)
#   ) |>
#   ggplot(aes(x = jour, y = mois, fill = couverture)) +
#   geom_tile(color = "white", linewidth = 0.3) +
#   facet_wrap(~produit, ncol = 1) +
#   scale_fill_gradientn(
#     colors = c("grey90", "#ffffcc", "#41b6c4", "#0c2c84"),
#     labels = scales::percent,
#     na.value = "grey95",
#     limits = c(0, 1)           # ← forcer l'échelle de 0 à 1
#   ) +
#   labs(
#     title    = "Couverture journalière par produit en 2024",
#     subtitle = "SEXTANT OC5 : couverture corrigée du masque terrestre",
#     fill     = "Couverture",
#     x        = "Jour",
#     y        = NULL
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     plot.title    = element_text(face = "bold"),
#     plot.subtitle = element_text(color = "grey50", size = 10),
#     strip.text    = element_text(face = "bold")
#   )

# comparaison |>
#   mutate(
#     mois = lubridate::month(date, label = TRUE),
#     jour = lubridate::day(date)
#   ) |>
#   ggplot(aes(x = jour, y = mois, fill = couverture)) +
#   geom_tile(color = "white", linewidth = 0.3) +
#   facet_wrap(~produit, ncol = 1) +
#   scale_fill_gradientn(
#     colors = c("grey90", "#ffffcc", "#41b6c4", "#0c2c84"),
#     labels = scales::percent,
#     na.value = "grey95"
#   ) +
#   labs(title = "Couverture journalière par produit en 2024",
#        fill = "Couverture")
# 
# ggplot(comparaison, aes(date, couverture, color = produit, fill = produit)) +
#   geom_ribbon(stat = "smooth", alpha = 0.15, color = NA) +
#   geom_smooth(se = FALSE, linewidth = 1.2, method = "loess", span = 0.2) +
#   geom_point(data = ~filter(.x, couverture == 0),  # marquer les jours sans données
#              aes(y = 0), shape = 124, size = 3, alpha = 0.5) +
#   scale_y_continuous(labels = scales::percent) +
#   facet_wrap(~produit, ncol = 1)

comparaison |>
  mutate(mois = lubridate::floor_date(date, "month"),
         classe = case_when(
           couverture == 0   ~ "Aucune donnée",
           couverture < 0.25 ~ "<25%",
           couverture < 0.50 ~ "25-50%",
           couverture >= 0.50 ~ ">50%"
         )) |>
  count(produit, mois, classe) |>
  ggplot(aes(x = mois, y = n, fill = classe)) +
  geom_col() +
  facet_wrap(~produit, ncol = 1) +
  scale_fill_brewer(palette = "RdYlGn", direction = 1)

couleurs <- c("SEXTANT" = "#1b9e77", "MODIS" = "#d95f02", "OLCI" = "#7570b3")

theme_propre <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey50", size = 11, margin = margin(b = 10)),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    strip.text       = element_text(face = "bold", size = 12),
    axis.title       = element_text(size = 11, color = "grey40"),
    plot.margin      = margin(12, 16, 12, 16)
  )

# 1. Boxplot ---------------------------------------------------------------
stats_resume <- comparaison |>
  group_by(produit) |>
  summarise(
    med = median(couverture),
    n   = sum(couverture > 0)
  )

ggplot(comparaison, aes(x = produit, y = couverture, fill = produit, color = produit)) +
  geom_jitter(width = 0.18, alpha = 0.5, size = 2, show.legend = FALSE) +
  geom_boxplot(width = 0.45, alpha = 0.75, outlier.shape = NA, linewidth = 0.6) +
  geom_text(
    data = stats_resume,
    aes(x = produit, y = -0.06, label = paste0(n, " jours\navec données")),
    inherit.aes = FALSE, size = 5, color = "grey45", lineheight = 0.9
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(-0.1, 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values  = couleurs) +
  scale_color_manual(values = couleurs) +
  labs(
    title    = "Couverture valide par produit satellitaire — 2024",
    # subtitle = "Chaque point représente un jour de l'année",
    x = NULL, y = "Couverture valide"
  ) +
  theme_propre

# 2. Heatmap calendrier ----------------------------------------------------
comparaison |>
  mutate(
    mois = month(date, label = TRUE, abbr = TRUE),
    jour = day(date)
  ) |>
  ggplot(aes(x = jour, y = fct_rev(mois), fill = couverture)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~produit, ncol = 1) +
  scale_fill_gradientn(
    colors   = c("grey92", "#c7e9b4", "#41b6c4", "#225ea8"),
    na.value = "grey96",
    labels   = percent_format(accuracy = 1),
    name     = "Couverture"
  ) +
  scale_x_continuous(breaks = c(1, 10, 20, 31), expand = c(0, 0)) +
  labs(
    title    = "Couverture journalière par produit — 2024",
    # subtitle = "Gris = aucune donnée ce jour-là",
    x = "Jour du mois", y = NULL
  ) +
  theme_propre +
  theme(
    legend.position  = "right",
    panel.grid       = element_blank(),
    axis.text.y      = element_text(size = 9),
    strip.text       = element_text(face = "bold")
  )

# 3. Barres mensuelles -----------------------------------------------------
comparaison |>
  mutate(
    mois   = floor_date(date, "month"),
    classe = case_when(
      couverture == 0    ~ "Aucune donnée",
      couverture < 0.25  ~ "< 25 %",
      couverture < 0.50  ~ "25 – 50 %",
      couverture >= 0.50 ~ "> 50 %"
    ),
    classe = factor(classe, levels = c("> 50 %", "25 – 50 %", "< 25 %", "Aucune donnée"))
  ) |>
  count(produit, mois, classe) |>
  ggplot(aes(x = mois, y = n, fill = classe)) +
  geom_col(width = 25, color = "white", linewidth = 0.3) +
  facet_wrap(~produit, ncol = 1) +
  scale_fill_manual(values = c(
    "> 50 %"       = "#2c7bb6",
    "25 – 50 %"    = "#abd9e9",
    "< 25 %"       = "#fdae61",
    "Aucune donnée" = "grey88"
  )) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = c(0.02, 0)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30), expand = c(0, 0)) +
  labs(
    title    = "Jours utilisables par mois — 2022",
    # subtitle = "Un jour > 50 % de couverture est considéré exploitable",
    x = NULL, y = "Nombre de jours", fill = NULL
  ) +
  theme_propre +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text     = element_text(size = 10)
  )

# extraction des bandes qualité -------------------------------------------

## OLCI --------------------------------------------------------------------

# Inspecter la structure complète du fichier
tidync("~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/L3m_20240722__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc")

# Voir toutes les variables disponibles
nc <- ncdf4::nc_open("~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/L3m_20240722__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc")
print(nc)
ncdf4::nc_close(nc)


nc_file <- "~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/L3m_20240722__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc"

nc <- nc_open(nc_file)

# Extraire les deux variables en matrices brutes
spm   <- ncvar_get(nc, "SPM-G-PO_mean")
flags <- ncvar_get(nc, "SPM-G-PO_flags")
lon   <- ncvar_get(nc, "lon")
lat   <- ncvar_get(nc, "lat")

nc_close(nc)

# Voir la distribution brute des flags (incluant 0)
table(flags, useNA = "always")

# Construire le dataframe complet
olci_df <- expand.grid(lon = lon, lat = lat) |>
  mutate(
    spm   = as.vector(spm),
    flags = as.vector(flags)
  )

# Pixels valides selon SPM (FillValue = -999)
olci_df |>
  filter(spm != -999) |>
  count(flags) |>
  arrange(flags)

# Bit 0 = valide (commun à 65, 129, 2177, 2241)
# Exclure flag = 9 (suspect, pas de SPM associé de toute façon)
# Garder uniquement les pixels avec SPM valide ET flag avec bit 0 actif

olci_df_clean <- olci_df |>
  filter(
    !is.na(spm),
    spm != -999,
    bitwAnd(flags, 1L) == 1L  # bit 0 activé = pixel valide
  )

nrow(olci_df_clean)  # doit donner ~9972

# Vérification : quels flags restent ?
count(olci_df_clean, flags)

load_OLCI_flags <- function(file_name, lon_range, lat_range) {
  
  nc <- nc_open(file_name)
  spm   <- ncvar_get(nc, "SPM-G-PO_mean")
  flags <- ncvar_get(nc, "SPM-G-PO_flags")
  lon   <- ncvar_get(nc, "lon")
  lat   <- ncvar_get(nc, "lat")
  nc_close(nc)
  
  # Date depuis le nom de fichier
  date_str <- stringr::str_extract(basename(file_name), "\\d{8}")
  date <- as.Date(date_str, format = "%Y%m%d")
  
  expand.grid(lon = lon, lat = lat) |>
    mutate(
      spm   = as.vector(spm),
      flags = as.vector(flags),
      date  = date
    ) |>
    filter(
      lon >= lon_range[1], lon <= lon_range[2],
      lat >= lat_range[1], lat <= lat_range[2],
      !is.na(spm),
      spm != -999,
      bitwAnd(flags, 1L) == 1L
    ) |>
    dplyr::select(lon, lat, date, spm, flags)
}

library(ncdf4)
library(dplyr)
library(stringr)
library(purrr)

# Définis ta zone autour de Nice
lon_range <- c(6.8925000, 7.4200000)
lat_range <- c(43.2136389, 43.7300000)

# Liste tous les fichiers OLCI 2024
files_olci_a <- list.files(
  "~/Downloads/OLCI/SPM/2024/OLCI_A_ODATIS_MR_2024_SPM/",
  pattern = "\\.nc$",
  full.names = TRUE
)

files_olci_b <- list.files(
  "~/Downloads/OLCI/SPM/2024/OLCI_B_ODATIS_MR_2024_SPM/",
  pattern = "\\.nc$",
  full.names = TRUE
)

# Combiner les deux listes de fichiers
files_olci <- c(files_olci_a, files_olci_b)

# Applique la fonction sur tous les fichiers
OLCI_2024_flags <- files_olci |>
  map_dfr(
    \(f) {
      tryCatch(
        load_OLCI_flags(f, lon_range, lat_range),
        error = function(e) {
          message("Erreur sur : ", basename(f), " — ", e$message)
          NULL
        }
      )
    }
  )

# Résumé
glimpse(OLCI_2024_flags)
nrow(OLCI_2024_flags)
n_distinct(OLCI_2024_flags$date)

OLCI_2024_flags |>
  group_by(date) |>
  summarise(n = n()) |>
  arrange(desc(n)) |>
  head(10)

count(OLCI_2024_flags, flags) |> arrange(flags)

# Pixels totaux réels de ta grille (lon x lat de ton subset Nice)
# À calculer une fois depuis un fichier
nc <- nc_open(files_olci[1])
lon_all <- ncvar_get(nc, "lon")
lat_all <- ncvar_get(nc, "lat")
nc_close(nc)

pixels_grille <- sum(lon_all >= lon_range[1] & lon_all <= lon_range[2]) *
  sum(lat_all >= lat_range[1] & lat_all <= lat_range[2])

cat("Pixels dans la grille Nice :", pixels_grille, "\n")

OLCI_stats_flags <- OLCI_2024_flags |>
  # Garder un seul pixel par position par jour (priorité au plus fort SPM ou au premier)
  group_by(date, lon, lat) |>
  slice_max(spm, n = 1, with_ties = FALSE) |>  # ou slice_head(n=1) si tu veux juste dédoublonner
  ungroup() |>
  # Ensuite calculer la couverture
  group_by(date) |>
  summarise(
    pixels_valides = n(),
    pixels_totaux  = pixels_grille,
    couverture     = pixels_valides / pixels_totaux,
    produit        = "OLCI"
  )

summary(OLCI_stats_flags$couverture)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000385 0.1819083 0.4920311 0.4420101 0.7129369 0.7694025 

## MODIS --------------------------------------------------------------------

# Inspecter la structure complète du fichier
tidync("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2024_SPM/L3m_20240723__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc")

# Voir toutes les variables disponibles
nc <- ncdf4::nc_open("~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2024_SPM/L3m_20240723__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc")
print(nc)
ncdf4::nc_close(nc)


nc_file <- "~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2024_SPM/L3m_20240723__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc"

nc <- nc_open(nc_file)

# Extraire les deux variables en matrices brutes
spm   <- ncvar_get(nc, "SPM-G-NS_mean")
flags <- ncvar_get(nc, "SPM-G-NS_flags")
lon   <- ncvar_get(nc, "lon")
lat   <- ncvar_get(nc, "lat")

nc_close(nc)

# Voir la distribution brute des flags (incluant 0)
table(flags, useNA = "always")

# Construire le dataframe complet
modis_df <- expand.grid(lon = lon, lat = lat) |>
  mutate(
    spm   = as.vector(spm),
    flags = as.vector(flags)
  )

# Pixels valides selon SPM (FillValue = -999)
modis_df |>
  filter(spm != -999) |>
  count(flags) |>
  arrange(flags)

load_MODIS_flags <- function(file_name, lon_range, lat_range) {
  
  nc <- nc_open(file_name)
  spm   <- ncvar_get(nc, "SPM-G-NS_mean")
  flags <- ncvar_get(nc, "SPM-G-NS_flags")
  lon   <- ncvar_get(nc, "lon")
  lat   <- ncvar_get(nc, "lat")
  nc_close(nc)
  
  date_str <- stringr::str_extract(basename(file_name), "\\d{8}")
  date <- as.Date(date_str, format = "%Y%m%d")
  
  expand.grid(lon = lon, lat = lat) |>
    mutate(
      spm   = as.vector(spm),
      flags = as.vector(flags),
      date  = date
    ) |>
    filter(
      lon >= lon_range[1], lon <= lon_range[2],
      lat >= lat_range[1], lat <= lat_range[2],
      !is.na(spm),
      bitwAnd(flags, 16384L) == 16384L  # bit 14 = valide pour MODIS NS
    ) |>
    dplyr::select(lon, lat, date, spm, flags)
}

# Appliquer sur tous les fichiers MODIS 2024
files_modis <- list.files(
  "~/Downloads/MODIS ODATIS MR/SPM/MODIS_ODATIS_MR_2024_SPM/",
  pattern = "\\.nc$",
  full.names = TRUE
)

MODIS_2024_flags <- files_modis |>
  map_dfr(
    \(f) tryCatch(
      load_MODIS_flags(f, lon_range, lat_range),
      error = function(e) {
        message("Erreur : ", basename(f))
        NULL
      }
    )
  )

# Résumé
cat("Pixels valides :", nrow(MODIS_2024_flags), "\n")
cat("Jours couverts :", n_distinct(MODIS_2024_flags$date), "\n")

# Stats couverture
MODIS_stats_flags <- MODIS_2024_flags |>
  group_by(date) |>
  summarise(
    pixels_valides = n(),
    pixels_totaux  = pixels_grille,
    couverture     = pixels_valides / pixels_totaux,
    produit        = "MODIS"
  )

summary(MODIS_stats_flags$couverture)


## SEXTANT -----------------------------------------------------------------

nc_file_sextant <- "~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2024/07/22/20240722-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"

nc <- nc_open(nc_file_sextant)
print(nc)
nc_close(nc)

load_SEXTANT_flags <- function(file_name, lon_range, lat_range) {
  
  nc <- nc_open(file_name)
  spim <- ncvar_get(nc, "analysed_spim")
  mask <- ncvar_get(nc, "mask")
  lon  <- ncvar_get(nc, "lon")
  lat  <- ncvar_get(nc, "lat")
  nc_close(nc)
  
  date_str <- stringr::str_extract(basename(file_name), "\\d{8}")
  date <- as.Date(date_str, format = "%Y%m%d")
  
  expand.grid(lon = lon, lat = lat) |>
    mutate(
      spim = as.vector(spim),
      mask = as.vector(mask),
      date = date
    ) |>
    filter(
      lon >= lon_range[1], lon <= lon_range[2],
      lat >= lat_range[1], lat <= lat_range[2]
    ) |>
    dplyr::select(lon, lat, date, spim, mask)
}

# Test sur un fichier
test <- load_SEXTANT_flags(nc_file_sextant, lon_range, lat_range)

# Distribution du masque
count(test, mask)

files_sextant <- list.files(
  "~/pCloudDrive/Stage/SEXTANT/SPM/merged/Standard/DAILY/2024/",
  pattern = "\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)

length(files_sextant)  # doit être ~366

cat("Pixels SEXTANT sur Nice :", nrow(test), "\n")
# 1820
# vs pixels_grille OLCI/MODIS = 25976


# cartographie ------------------------------------------------------------

# on observe une différence pour certains jours dans la cartographie entre OLCI
# et MODIS, en effet, pour MODIS il mnaque des données à la côte

# Exemples sur le jour du 12/01/2024 : 
coastline_giscoR <- gisco_get_coastallines(resolution = "01")
countries_giscoR  <- gisco_get_countries(region = "Europe", resolution = "01")

# OLCI

OLCI_12_01_2024 <- OLCI_2016_2024_spm_pixels |> 
  filter(date == "2024-01-12")

max_spm <- max(OLCI_12_01_2024$`SPM-G-PO_mean`, na.rm = TRUE)

pl_map <- OLCI_12_01_2024 %>%
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
    title    = "Concentration en matières en suspension — 12 janvier 2024",
    subtitle = "OLCI - Correction atmosphérique = Polymer",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(OLCI_12_01_2024$lon),
    ylim   = range(OLCI_12_01_2024$lat),
    expand = FALSE,
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

ggsave("~/Satellite_analysis/Graphiques/OLCI/spatial/OLCI_12_01_2024.png", pl_map, height = 9, width = 14)

# MODIS
MODIS_12_01_2024 <- MODIS_2002_2024_spm_pixels |> 
  filter(date == "2024-01-12")

max_spm <- max(MODIS_12_01_2024$`SPM-G-NS_mean`, na.rm = TRUE)

pl_map <- MODIS_12_01_2024 %>%
  ggplot() +
  annotation_borders(fill = "grey80") +
  geom_tile(aes(x = lon, y = lat, fill = `SPM-G-NS_mean`)) +
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
    title    = "Concentration en matières en suspension — 12 janvier 2024",
    subtitle = "MODIS - Correction atmosphérique = NirSwir",
    x        = "Longitude (°E)",
    y        = "Latitude (°N)"
  ) +
  coord_sf(
    xlim   = range(MODIS_12_01_2024$lon),
    ylim   = range(MODIS_12_01_2024$lat),
    expand = FALSE,
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

ggsave("~/Satellite_analysis/Graphiques/MODIS/spatial/MODIS_12_01_2024.png", pl_map, height = 9, width = 14)

# On se demande si c'est sytématique 

library(sf)
library(dplyr)

# 1. Créer un masque côtier -----------------------------------------------
# Récupérer la côte
coastline_giscoR <- gisco_get_coastallines(resolution = "01")

# Créer une grille de référence depuis OLCI (tous les pixels uniques)
grille_ref <- OLCI_2016_2024_spm_pixels |>
  distinct(lon, lat) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Calculer la distance à la côte pour chaque pixel (en mètres)
# Attention : opération lente, à faire une seule fois
dist_cote <- st_distance(grille_ref, coastline_giscoR |> st_union())

grille_ref$dist_cote_km <- as.numeric(dist_cote) / 1000

# Définir les zones
grille_ref <- grille_ref |>
  mutate(
    zone = case_when(
      dist_cote_km <= 5  ~ "0-5 km",
      dist_cote_km <= 15 ~ "5-15 km",
      dist_cote_km <= 30 ~ "15-30 km",
      TRUE               ~ "> 30 km"
    )
  )

# Récupérer les coordonnées avec la zone
grille_zones <- grille_ref |>
  st_drop_geometry() |>
  bind_cols(grille_ref |> st_coordinates() |> as.data.frame() |> rename(lon = X, lat = Y))

# Sauvegarder pour ne pas recalculer
save(grille_zones, file = "data/grille_zones_cote.RData")

# 2. Joindre la zone à chaque produit -------------------------------------

OLCI_zones <- OLCI_2016_2024_spm_pixels |>
  filter(date >= "2024-01-01", date <= "2024-12-31") |>
  left_join(grille_zones, by = c("lon", "lat"))

MODIS_zones <- MODIS_2002_2024_spm_pixels |>
  filter(date >= "2024-01-01", date <= "2024-12-31") |>
  left_join(grille_zones, by = c("lon", "lat"))

# 3. Couverture par zone et par produit -----------------------------------

couverture_zone <- bind_rows(
  OLCI_zones |>
    group_by(zone, date) |>
    summarise(
      couverture = sum(!is.na(`SPM-G-PO_mean`)) / n(),
      produit = "OLCI",
      .groups = "drop"
    ),
  MODIS_zones |>
    group_by(zone, date) |>
    summarise(
      couverture = sum(!is.na(`SPM-G-NS_mean`)) / n(),
      produit = "MODIS",
      .groups = "drop"
    )
) |>
  mutate(zone = factor(zone, levels = c("0-5 km", "5-15 km", "15-30 km", "> 30 km")))

# ou

couverture_zone_modis <- MODIS_zones |>
  filter(!is.na(zone)) |>  # exclure terre et pixels sans zone
  group_by(zone, date) |>
  summarise(
    couverture = sum(!is.na(`SPM-G-NS_mean`)) / n(),
    .groups = "drop"
  ) |>
  mutate(zone = factor(zone, levels = c("0-5 km", "5-15 km", "15-30 km", "> 30 km")))

# 4. Visualisation --------------------------------------------------------

p2 <- ggplot(couverture_zone_modis,
       aes(x = zone, y = couverture)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA,
               fill = "#d95f02", color = "#993c1d") +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8, color = "#d95f02") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    # title    = "Couverture valide MODIS par distance à la côte — 2024",
    # subtitle = "Un biais côtier se manifeste si la couverture décroît vers la zone 0–5 km",
    x        = "Distance à la côte",
    y        = "Couverture valide"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey50", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )

# Visualiser la distribution spatiale des zones
p1 <- grille_zones |>
  ggplot(aes(x = lon, y = lat, color = zone)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "RdYlGn") +
  coord_fixed() +
  theme_minimal()

(p1 | p2) +
  plot_annotation(
    title   = "Couverture valide de MODIS par rapport à la distance à la côte",
    caption = "Source : MODIS",
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
    )
  )

# Combien de pixels MODIS sans zone assignée ?
MODIS_zones |>
  filter(date == "2024-07-22") |>
  count(is.na(zone))

# Visualiser où sont ces pixels sans zone
MODIS_zones |>
  filter(is.na(zone)) |>
  distinct(lon, lat) |>
  ggplot(aes(x = lon, y = lat)) +
  geom_point(size = 0.5, color = "red") +
  theme_minimal()

# Ces pixels sans zone ont-ils des valeurs SPM valides ?
MODIS_zones |>
  filter(date == "2024-07-22") |>
  group_by(is.na(zone)) |>
  summarise(n_valides = sum(!is.na(`SPM-G-NS_mean`)))
