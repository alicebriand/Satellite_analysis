# OLCI propre

library(dplyr)
library(tidync)
library(seasonal)
library(trend)

# ── Corrélation débit cumulé ~ aire panaches OLCI (seuil P99) ───────────────
# Méthode propre : seulement les vrais jours d'observation

load("data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")
load("data/Hydro France/All_debit.Rdata")

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

# 1. Recalculer OLCI_pixels_par_jour si pas déjà fait
OLCI_pixels_par_jour <- OLCI_2016_2024_spm_pixels |>
  group_by(date) |>
  summarise(
    n_pixels_total   = n(),
    n_pixels_valides = sum(!is.na(`SPM-G-PO_mean`)),
    pct_couverture   = n_pixels_valides / n_pixels_total * 100,
    .groups = "drop"
  )

# 2. Recalculer les métriques de panache avec le seuil P99
seuil_99 <- quantile(OLCI_2016_2024_spm_pixels$`SPM-G-PO_mean`,
                     0.99, na.rm = TRUE)
cat("Seuil P99 :", round(seuil_99, 3), "g/m³\n")

OLCI_2016_2024_spm_99 <- OLCI_2016_2024_spm_pixels |>
  group_by(date) |>
  summarise(
    pixel_count      = sum(`SPM-G-PO_mean` >= seuil_99, na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    mean_spm         = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99],
                            na.rm = TRUE),
    median_spm       = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99],
                              na.rm = TRUE),
    .groups = "drop"
  )

# 3. Joindre avec la couverture et filtrer les faux zéros
OLCI_spm_99_clean <- OLCI_2016_2024_spm_99 |>
  left_join(
    OLCI_pixels_par_jour |> select(date, pct_couverture, n_pixels_valides),
    by = "date"
  ) |>
  mutate(
    statut = case_when(
      n_pixels_valides == 0              ~ "aucune observation",
      aire_panache_km2 == 0              ~ "observation sans panache",
      TRUE                               ~ "panache détecté"
    )
  )

# Résumé
cat("\n── Répartition des jours ──\n")
OLCI_spm_99_clean |> count(statut) |> print()

# 4. Garder seulement les vrais jours d'observation avec panache
OLCI_spm_99_valide <- OLCI_spm_99_clean |>
  filter(
    n_pixels_valides > 0,       # ← au moins un pixel observé
    pct_couverture   >= 10,     # ← au moins 10% de couverture
    aire_panache_km2  > 0       # ← panache détecté
  )

cat("\nNombre de jours valides avec panache :", nrow(OLCI_spm_99_valide), "\n")

# All_debit_2016_2024 <- All_debit |> 
#   filter(date >= "2016-01-01", date <= "2024-12-31")

load("data/Hydro France/Y6442010_2016_2024.Rdata")

# 5. Joindre avec le débit cumulé
merged_corr <- OLCI_spm_99_valide |>
  inner_join(
    Y6442010_2016_2024 |> select(date, débit),
    by = "date"
  ) |>
  drop_na(aire_panache_km2, débit) |>
  filter(
    débit     > 0,
    aire_panache_km2 > 0
  ) |>
  mutate(
    log_aire  = log10(aire_panache_km2),
    log_debit = log10(débit)
  )

cat("Nombre de jours avec débit ET panache :", nrow(merged_corr), "\n")

# 6. Corrélation de Spearman
cor_result <- cor.test(
  merged_corr$aire_panache_km2,
  merged_corr$débit,
  method = "spearman"
)

# Corrélation log-log
cor_log <- cor.test(
  merged_corr$log_aire,
  merged_corr$débit,
  method = "spearman"
)

cat("\n── Résultats ──\n")
cat("r Spearman (échelle normale) =", round(cor_result$estimate, 3), "\n")
cat("p-value =", format(cor_result$p.value, scientific = TRUE, digits = 3), "\n")
cat("r Spearman (log-log) =", round(cor_log$estimate, 3), "\n")
cat("p-value =", format(cor_log$p.value, scientific = TRUE, digits = 3), "\n")

# 7. Modèle log-log
modele_log <- lm(log_aire ~ log_debit, data = merged_corr)
r2       <- summary(modele_log)$r.squared
pente    <- coef(modele_log)[2]
ordonnee <- coef(modele_log)[1]

cat("R² =", round(r2, 3), "\n")

# 8. Graphique log-log
label_eq <- paste0(
  "r = ", round(cor_log$estimate, 2),
  "\nR² = ", round(r2, 2),
  "\np ", ifelse(cor_log$p.value < 0.001, "< 0.001",
                 paste0("= ", format(cor_log$p.value, digits = 2, scientific = TRUE)))
)

ggplot(merged_corr, aes(x = débit, y = aire_panache_km2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", se = TRUE,
              linewidth = 0.8, fill = "grey80") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  annotate("text",
           x        = max(merged_corr$débit, na.rm = TRUE),
           y        = min(merged_corr$aire_panache_km2, na.rm = TRUE) * 2,
           label    = label_eq,
           hjust    = 1, vjust = 0,
           size     = 8, color  = "grey20",
           fontface = "italic", family = "serif") +
  labs(
    x        = expression("Débit cumulé (m"^{3}*" s"^{-1}*")"),
    y        = "Aire du panache (km²)",
    title    = "Relation débit du Var et aire du panache turbide",
    subtitle = paste0("OLCI 2016–2024 — seuil P99 = ",
                      round(seuil_99, 2), " g/m³")
                      # nrow(merged_corr), " jours")
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(color = "grey40", hjust = 0.5, size = 11)
  )

# ── 1. Préparer les panaches sans faux zéros ─────────────────────────────────
OLCI_spm_99_valide <- OLCI_2016_2024_spm_99 |>
  left_join(
    OLCI_pixels_par_jour |> select(date, pct_couverture, n_pixels_valides),
    by = "date"
  ) |>
  filter(
    n_pixels_valides > 0,
    pct_couverture   >= 10,
    aire_panache_km2  > 0
  ) |>
  drop_na(aire_panache_km2) |>
  mutate(date_num = as.numeric(date - min(date)))

# ── 2. Mann-Kendall + Theil-Sen sur les panaches ─────────────────────────────
mk_panache      <- mk.test(OLCI_spm_99_valide$aire_panache_km2)
sen_panache     <- sens.slope(OLCI_spm_99_valide$aire_panache_km2)
slope_kmjan_pan <- sen_panache$estimates * 365

intercept_pan <- median(
  OLCI_spm_99_valide$aire_panache_km2 - sen_panache$estimates * OLCI_spm_99_valide$date_num,
  na.rm = TRUE
)

OLCI_spm_99_valide <- OLCI_spm_99_valide |>
  mutate(theilsen_fit = intercept_pan + sen_panache$estimates * date_num)

cat("Panaches — Mann-Kendall p =", round(mk_panache$p.value, 4),
    "| Theil-Sen pente =", round(slope_kmjan_pan, 2), "km²/an\n")

# ── 3. Mann-Kendall + Theil-Sen sur le débit ─────────────────────────────────
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

# ── 4. Corrélation Spearman sur les vrais jours ───────────────────────────────
merged_data <- inner_join(
  All_debit        |> select(date, debit_cumule),
  OLCI_spm_99_valide |> select(date, aire_panache_km2),
  by = "date"
) |> drop_na()

correlation <- cor(merged_data$debit_cumule, merged_data$aire_panache_km2,
                   method = "spearman", use = "complete.obs")
p_value     <- cor.test(merged_data$debit_cumule, merged_data$aire_panache_km2,
                        method = "spearman")$p.value

n_panache <- nrow(OLCI_spm_99_valide)
n_commun  <- nrow(merged_data)

cat("r Spearman =", round(correlation, 3), "\n")
cat("p-value =", format(p_value, scientific = TRUE, digits = 3), "\n")
cat("n commun =", n_commun, "\n")

# ── 5. Mise à l'échelle pour double axe ──────────────────────────────────────
adjust_factors <- sec_axis_adjustement_factors(
  OLCI_spm_99_valide$aire_panache_km2,
  All_debit$debit_cumule
)

OLCI_spm_99_valide <- OLCI_spm_99_valide |>
  mutate(
    scaled_aire      = aire_panache_km2 * adjust_factors$diff + adjust_factors$adjust,
    theilsen_scaled  = theilsen_fit     * adjust_factors$diff + adjust_factors$adjust
  )

# ── 6. Graphique ─────────────────────────────────────────────────────────────
ggplot() +
  # Aire des panaches
  geom_line(
    data = OLCI_spm_99_valide,
    aes(x = date, y = scaled_aire, color = "Aire des panaches"),
    linewidth = 0.4, alpha = 0.6
  ) +
  # Tendance Theil-Sen panaches
  geom_line(
    data = OLCI_spm_99_valide,
    aes(x = date, y = theilsen_scaled, color = "Tendance panaches"),
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
      "Tendance panaches" = "cyan4",
      "Débit cumulé"      = "darkolivegreen3",
      "Tendance débit"    = "darkolivegreen4"
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
  # Corrélation — haut droite
  annotate("text",
           x = max(OLCI_spm_99_valide$date, na.rm = TRUE),
           y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.97,
           hjust = 1, vjust = 1, size = 8,
           color = "grey20", fontface = "italic", family = "serif",
           label = paste0(
             "R = ", round(correlation, 2),
             "\np ", ifelse(p_value < 0.001, "< 0.001",
                            format(p_value, digits = 3)),
             "\nn = ", n_commun
           )) +
  # Tendance panaches — haut gauche
  annotate("text",
           x = as.Date("2016-04-26"),
           y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.97,
           hjust = 0, vjust = 1, size = 8,
           color = "cyan4", fontface = "italic", family = "serif",
           label = paste0(
             "Tendance panaches : ", round(slope_kmjan_pan, 2), " km²/an",
             "\np ", ifelse(mk_panache$p.value < 0.001, "< 0.001",
                            ifelse(mk_panache$p.value < 0.05, "< 0.05",
                                   paste0("= ", round(mk_panache$p.value, 3)))),
             "\nn = ", n_panache
           )) +
  # Tendance débit — milieu gauche
  annotate("text",
           x = as.Date("2016-04-26"),
           y = max(All_debit$debit_cumule, na.rm = TRUE) * 0.70,
           hjust = 0, vjust = 1, size = 8,
           color = "darkolivegreen4", fontface = "italic", family = "serif",
           label = paste0(
             "Tendance débit : ", round(slope_debit_an, 2), " m³/s/an",
             "\np ", ifelse(mk_debit$p.value < 0.001, "< 0.001",
                            ifelse(mk_debit$p.value < 0.05, "< 0.05",
                                   paste0("= ", round(mk_debit$p.value, 3))))
           )) +
  labs(
    title    = "Évolution de l'extension des panaches turbides et du débit cumulé",
    subtitle = "Produit OLCI — ODATIS-MR (2016–2024) — jours avec observations valides uniquement",
    x        = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 18, hjust = 0.5, family = "serif"),
    plot.subtitle    = element_text(size = 13, hjust = 0.5, color = "grey30", family = "serif"),
    axis.title       = element_text(face = "bold", family = "serif", size = 16),
    axis.text        = element_text(color = "grey30", family = "serif", size = 14),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    legend.position  = "top",
    legend.text      = element_text(size = 15),
    plot.margin      = margin(1, 1.5, 1, 1, "cm")
  )



# ── 1. Préparer le débit sans lag ─────────────────────────────────────────────
debit_j0 <- All_debit_2016_2024 |>
  arrange(date) |>
  select(date, debit_cumule) |>
  drop_na(debit_cumule)

# ── 2. Préparer les métriques de panache (sans faux zéros) ───────────────────
OLCI_panache_gangloff <- OLCI_panache_metrics |>
  left_join(
    OLCI_pixels_par_jour |> select(date, pct_couverture, n_pixels_valides),
    by = "date"
  ) |>
  filter(
    n_pixels_valides > 0,
    pct_couverture   >= 10,
    pixel_count      > 0
  ) |>
  inner_join(
    debit_j0 |> select(date, debit_cumule),
    by = "date"
  ) |>
  drop_na(aire_panache_km2, mean_spm, max_spm, debit_cumule) |>
  filter(debit_cumule > 0)

cat("Nombre d'images valides :", nrow(OLCI_panache_gangloff), "\n")

# ── 3. Corrélations ───────────────────────────────────────────────────────────
cor_aire <- cor.test(OLCI_panache_gangloff$debit_cumule,
                     OLCI_panache_gangloff$aire_panache_km2,
                     method = "pearson")

cor_mean <- cor.test(OLCI_panache_gangloff$debit_cumule,
                     OLCI_panache_gangloff$mean_spm,
                     method = "pearson")

cor_max  <- cor.test(OLCI_panache_gangloff$debit_cumule,
                     OLCI_panache_gangloff$max_spm,
                     method = "pearson")

cat("\n── Résultats ──\n")
cat("Aire ~ Q  : R² =", round(cor_aire$estimate^2, 2),
    "| r =", round(cor_aire$estimate, 2),
    "| p =", format(cor_aire$p.value, scientific = TRUE, digits = 2), "\n")
cat("Cmean ~ Q : R² =", round(cor_mean$estimate^2, 2),
    "| r =", round(cor_mean$estimate, 2),
    "| p =", format(cor_mean$p.value, scientific = TRUE, digits = 2), "\n")
cat("Cmax ~ Q  : R² =", round(cor_max$estimate^2, 2),
    "| r =", round(cor_max$estimate, 2),
    "| p =", format(cor_max$p.value, scientific = TRUE, digits = 2), "\n")

# ── 4. Modèles linéaires ──────────────────────────────────────────────────────
modele_aire <- lm(aire_panache_km2 ~ debit_cumule, data = OLCI_panache_gangloff)
modele_mean <- lm(mean_spm         ~ debit_cumule, data = OLCI_panache_gangloff)
modele_max  <- lm(max_spm          ~ debit_cumule, data = OLCI_panache_gangloff)

# ── 5. Graphiques ─────────────────────────────────────────────────────────────
plot_corr_gangloff <- function(df, y_var, y_label, modele, cor_result, titre) {
  
  r2    <- round(cor_result$estimate^2, 2)
  pente <- coef(modele)[2]
  ord   <- coef(modele)[1]
  
  label_eq <- paste0(
    y_label, " = ", round(ord, 2), " + ", round(pente, 5), " × Q",
    "\nR² = ", r2,
    "\np ", ifelse(cor_result$p.value < 0.001, "< 0.001",
                   format(cor_result$p.value, digits = 2))
  )
  
  ggplot(df, aes(x = debit_cumule, y = .data[[y_var]])) +
    geom_point(alpha = 0.5, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", se = TRUE,
                linewidth = 0.8, fill = "grey80") +
    annotate("text",
             x        = max(df$debit_cumule, na.rm = TRUE),
             y        = max(df[[y_var]], na.rm = TRUE) * 0.95,
             label    = label_eq,
             hjust    = 1, vjust = 1,
             size     = 5, color = "grey20",
             fontface = "italic") +
    labs(
      title = titre,
      x     = expression("Débit cumulé (m"^{3}*" s"^{-1}*")"),
      y     = y_label
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5)
    )
}

p_aire <- plot_corr_gangloff(
  OLCI_panache_gangloff, "aire_panache_km2",
  "Aire du panache (km²)",
  modele_aire, cor_aire,
  "Extension du panache ~ débit"
)

p_mean <- plot_corr_gangloff(
  OLCI_panache_gangloff, "mean_spm",
  expression("C"["mean"]*" (g m"^{-3}*")"),
  modele_mean, cor_mean,
  "Concentration moyenne ~ débit"
)

p_max <- plot_corr_gangloff(
  OLCI_panache_gangloff, "max_spm",
  expression("C"["max"]*" (g m"^{-3}*")"),
  modele_max, cor_max,
  "Concentration maximale ~ débit"
)

p_aire / p_mean / p_max



