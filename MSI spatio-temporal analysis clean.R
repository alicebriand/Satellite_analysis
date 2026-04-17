# MSI spatio-temporal analysis clean

# 02/03/2026

# pathway : "~/Satellite_analysis/MSI spatio-temporal analysis clean"

# This script will load MSI df product and perform temporal analysis

# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

# loading data ------------------------------------------------------------

load("data/MSI/SPM/MSI df/MSI_2020_SPM_df.Rdata")
load("data/MSI/SPM/MSI df/MSI_2021_SPM_df.Rdata")
load("data/MSI/SPM/MSI df/MSI_2022_SPM_df.Rdata")
load("data/MSI/SPM/MSI df/MSI_2023_SPM_df.Rdata")
load("data/MSI/SPM/MSI df/MSI_2024_SPM_df.Rdata")
load("data/MSI/SPM/MSI df/MSI_2025_SPM_df.Rdata")

load("data/MSI/SPM/95ème percentile/MSI_2020_SPM_95.Rdata")
load("data/MSI/SPM/95ème percentile/MSI_2021_SPM_95.Rdata")
load("data/MSI/SPM/95ème percentile/MSI_2022_SPM_95.Rdata")
load("data/MSI/SPM/95ème percentile/MSI_2023_SPM_95.Rdata")
load("data/MSI/SPM/95ème percentile/MSI_2024_SPM_95.Rdata")
load("data/MSI/SPM/95ème percentile/MSI_2025_SPM_95.Rdata")

load("data/MSI/CHL/MSI df/MSI_2020_CHL_df.Rdata")
load("data/MSI/CHL/MSI df/MSI_2021_CHL_df.Rdata")
load("data/MSI/CHL/MSI df/MSI_2022_CHL_df.Rdata")
load("data/MSI/CHL/MSI df/MSI_2023_CHL_df.Rdata")
load("data/MSI/CHL/MSI df/MSI_2024_CHL_df.Rdata")
load("data/MSI/CHL/MSI df/MSI_2025_CHL_df.Rdata")

load("data/Hydro France/All_debit_2020.Rdata")
load("data/Hydro France/All_debit_2021.Rdata")
load("data/Hydro France/All_debit_2022.Rdata")
load("data/Hydro France/All_debit_2023.Rdata")
load("data/Hydro France/All_debit_2024.Rdata")
load("data/Hydro France/All_debit_2025.Rdata")
load("data/Hydro France/All_debit_2020_2025.Rdata")

# Functions ---------------------------------------------------------------

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

# pixel area --------------------------------------------------------------

## extraction des valeurs en degré -----------------------------------------

# Résolution en longitude
res_lon <- diff(sort(unique(MSI_2020_SPM$longitude)))[1]

# Résolution en latitude
res_lat <- diff(sort(unique(MSI_2020_SPM$latitude)))[1]

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

# define 95ème percentile -------------------------------------------------

# 2020

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2020_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.5759109 g/m³

# 2021

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2020_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.4928144 g/m³

# 2022

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2022_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.4302868 g/m³

# 2023

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2023_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.6213881 g/m³

# 2024

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2024_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.7444398 g/m³

# 2025

# Calculer le 95ème percentile
seuil_95 <- quantile(MSI_2025_SPM$SPM, 0.95, na.rm = TRUE)
cat("Seuil 95ème percentile :", seuil_95, "g/m³\n")

# seuil = 0.5575231 g/m³

# Comme on a des seuils différents pour chaque année, on moyenne les seuils sur les 
# 5 ans pour faire un seuil moyen

seuil_global <- 0.57039385

# 2020 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2020_SPM_95 <- MSI_2020_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_95, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_95], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2020_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2020_SPM_95.Rdata")

# 2021 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2021_SPM_95 <- MSI_2021_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_global, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_global], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_global], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2021_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2021_SPM_95.Rdata")

# 2022 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2022_SPM_95 <- MSI_2022_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_95, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_95], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2022_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2022_SPM_95.Rdata")

# 2023 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2023_SPM_95 <- MSI_2023_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_95, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_95], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2023_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2023_SPM_95.Rdata")

# 2024 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2024_SPM_95 <- MSI_2024_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_95, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_95], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2024_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2024_SPM_95.Rdata")

# 2025 --------------------------------------------------------------------

# Stats du panache par jour
MSI_2025_SPM_95 <- MSI_2025_SPM |> 
  group_by(date) |> 
  summarise(
    pixel_count = sum(SPM >= seuil_95, na.rm = TRUE),
    mean_spm = mean(SPM[SPM >= seuil_95], na.rm = TRUE),
    median_spm = median(SPM[SPM >= seuil_95], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2  # si tu as déjà calculé aire_pixel_km2
  )

save(MSI_2025_SPM_95, file = "data/MSI/SPM/95ème percentile/MSI_2025_SPM_95.Rdata")

# plotting ----------------------------------------------------------------

# 2020 --------------------------------------------------------------------

MSI_2020_SPM_95_deb <- MSI_2020_SPM_95 |> 
  left_join(All_debit_2020, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2020_SPM_95_deb$aire_panache_km2, MSI_2020_SPM_95_deb$debit_cumule)

MSI_2020_SPM_95_deb$scaled_aire_panache <- MSI_2020_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2020_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2020_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2020_SPM_95_deb$debit_cumule, MSI_2020_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2020_SPM_95_deb$debit_cumule, MSI_2020_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2020_SPM_95_deb$debit_cumule, 
                       MSI_2020_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2020_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.4) +
  geom_point(data = MSI_2020_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2020_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2020_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2020 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# 2021 --------------------------------------------------------------------

MSI_2021_SPM_95_deb <- MSI_2021_SPM_95 |> 
  left_join(All_debit_2021, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2021_SPM_95_deb$aire_panache_km2, MSI_2021_SPM_95_deb$debit_cumule)

MSI_2021_SPM_95_deb$scaled_aire_panache <- MSI_2021_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2021_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2021_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2021_SPM_95_deb$debit_cumule, MSI_2021_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2021_SPM_95_deb$debit_cumule, MSI_2021_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2021_SPM_95_deb$debit_cumule, 
                       MSI_2021_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2021_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 1.5, alpha = 0.4) +
  geom_point(data = MSI_2021_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 1.5, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2021_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2021_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2021 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# 2022 --------------------------------------------------------------------

MSI_2022_SPM_95_deb <- MSI_2022_SPM_95 |> 
  left_join(All_debit_2022, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2022_SPM_95_deb$aire_panache_km2, MSI_2022_SPM_95_deb$debit_cumule)

MSI_2022_SPM_95_deb$scaled_aire_panache <- MSI_2022_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2022_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2022_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2022_SPM_95_deb$debit_cumule, MSI_2022_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2022_SPM_95_deb$debit_cumule, MSI_2022_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2022_SPM_95_deb$debit_cumule, 
                       MSI_2022_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2022_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.4) +
  geom_point(data = MSI_2022_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2022_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2022_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2022 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# 2023 --------------------------------------------------------------------

MSI_2023_SPM_95_deb <- MSI_2023_SPM_95 |> 
  left_join(All_debit_2023, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2023_SPM_95_deb$aire_panache_km2, MSI_2023_SPM_95_deb$debit_cumule)

MSI_2023_SPM_95_deb$scaled_aire_panache <- MSI_2023_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2023_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2023_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2023_SPM_95_deb$debit_cumule, MSI_2023_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2023_SPM_95_deb$debit_cumule, MSI_2023_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2023_SPM_95_deb$debit_cumule, 
                       MSI_2023_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2023_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.4) +
  geom_point(data = MSI_2023_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2023_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2023_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2023 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# 2024 --------------------------------------------------------------------

MSI_2024_SPM_95_deb <- MSI_2024_SPM_95 |> 
  left_join(All_debit_2024, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2024_SPM_95_deb$aire_panache_km2, MSI_2024_SPM_95_deb$debit_cumule)

MSI_2024_SPM_95_deb$scaled_aire_panache <- MSI_2024_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2024_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2024_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2024_SPM_95_deb$debit_cumule, MSI_2024_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2024_SPM_95_deb$debit_cumule, MSI_2024_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2024_SPM_95_deb$debit_cumule, 
                       MSI_2024_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2024_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.4) +
  geom_point(data = MSI_2024_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2024_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2024_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np = ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2024 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )

# 2025 --------------------------------------------------------------------

MSI_2025_SPM_95_deb <- MSI_2025_SPM_95 |> 
  left_join(All_debit_2025, by = "date")

adjust_factors <- sec_axis_adjustement_factors(MSI_2025_SPM_95_deb$aire_panache_km2, MSI_2025_SPM_95_deb$debit_cumule)

MSI_2025_SPM_95_deb$scaled_aire_panache <- MSI_2025_SPM_95_deb$aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust

# 1. Tester la normalité
shapiro.test(MSI_2025_SPM_95_deb$debit_cumule)
shapiro.test(MSI_2025_SPM_95_deb$aire_panache_km2)
# Si p-value < 0.05 → pas normal → Spearman

# 2. Visualiser la relation
plot(MSI_2025_SPM_95_deb$debit_cumule, MSI_2025_SPM_95_deb$aire_panache_km2)
# Si la relation est courbe → Spearman
cor.test(MSI_2025_SPM_95_deb$debit_cumule, MSI_2025_SPM_95_deb$aire_panache_km2, method = "spearman")

# 1. Stocker le résultat du cor.test
cor_result <- cor.test(MSI_2025_SPM_95_deb$debit_cumule, 
                       MSI_2025_SPM_95_deb$aire_panache_km2, 
                       method = "spearman")

# 2. Extraire les valeurs
rho <- round(as.numeric(cor_result$estimate), 3)
p_value <- cor_result$p.value

# 3. Plotting
ggplot() +
  geom_point(data = MSI_2025_SPM_95_deb,
             aes(x = date, y = debit_cumule, color = "Débit cumulé"),
             size = 2, alpha = 0.4) +
  geom_point(data = MSI_2025_SPM_95_deb,
             aes(x = date, y = scaled_aire_panache, color = "Aire des panaches"),
             size = 2, alpha = 0.4) +
  annotate(
    "text",
    x = min(MSI_2025_SPM_95_deb$date, na.rm = TRUE),
    y = max(MSI_2025_SPM_95_deb$scaled_aire_panache, na.rm = TRUE) * 0.95,
    label = paste0(
      "ρ = ", rho,
      "\np ", ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))
    ),
    hjust = 0, vjust = 1,
    size = 8,
    color = "grey20",
    fontface = "italic"
  ) +   # ← le + manquait ici !
  scale_color_manual(values = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan")) +
  scale_fill_manual(values  = c("Débit cumulé" = "darkolivegreen3", "Aire des panaches" = "darkcyan"),
                    guide = "none") +
  scale_y_continuous(
    name = "Débit (m³/s)",
    expand = expansion(mult = c(0.02, 0.08)),
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff,
                        name = "Aire des panaches (en km²)")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title   = "Aire des panaches et débit cumulé journalier du Var, Paillon et Magnan en 2025 — MSI (Copernicus Marine)",
    x       = NULL,
    color   = NULL,
    caption = "Source : Mediterranean Sea, Bio-Geo-Chemical, L3, daily observation | product layer : cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m "
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", margin = margin(b = 10)),
    plot.caption      = element_text(size = 12, color = "grey50", hjust = 0),
    axis.title.y      = element_text(size = 14, margin = margin(r = 10)),  # axe Y gauche
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)), # axe Y droite
    axis.title.x      = element_text(size = 14, margin = margin(t = 10)),  # axe X
    axis.text         = element_text(size = 13, color = "grey30"),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.ticks        = element_line(color = "grey70"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "grey70", linewidth = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 10)
  )


