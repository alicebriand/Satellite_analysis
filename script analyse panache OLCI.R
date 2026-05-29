# script analyse panache OLCI

library(dunn.test)

# load data ---------------------------------------------------------------

load("data/Hydro France/All_debit.Rdata")
load("data/OLCI/SPM/OLCI_2016_2024_spm_pixels.Rdata")

# nettoyage SPM -----------------------------------------------------------

OLCI_2016_2024_spm_clean <- OLCI_2016_2024_spm_pixels |> 
  filter(`SPM-G-PO_mean` >= 0, `SPM-G-PO_mean` < 100) |>
  mutate(date = as.Date(date))

# pixel area --------------------------------------------------------------

nc <- tidync("~/Downloads/OLCI/SPM/2016/OLCI_A_ODATIS_MR_2016_SPM/L3m_20160426__FRANCE_03_OLA_SPM-G-PO_DAY_00.nc")

tmp <- hyper_tibble(nc) |>
  mutate(lon = as.numeric(lon), lat = as.numeric(lat))

res_lon <- diff(sort(unique(tmp$lon)))[1]
res_lat <- diff(sort(unique(tmp$lat)))[1]

cat("Résolution lon :", res_lon, "°\n")
cat("Résolution lat :", res_lat, "°\n")

lat_ref        <- 43
res_lon_km     <- res_lon * 111 * cos(lat_ref * pi / 180)
res_lat_km     <- res_lat * 111
aire_pixel_km2 <- res_lon_km * res_lat_km

cat("Aire d'un pixel OLCI :", round(aire_pixel_km2, 6), "km²\n")

# seuil 99ème percentile --------------------------------------------------

seuil_99 <- quantile(OLCI_2016_2024_spm_clean$`SPM-G-PO_mean`, 0.99, na.rm = TRUE)
cat("Seuil 99ème percentile OLCI :", seuil_99, "g/m³\n")

# stats du panache par jour -----------------------------------------------

panache_journalier_olci <- OLCI_2016_2024_spm_clean |>
  group_by(date) |>
  summarise(
    pixel_count      = sum(`SPM-G-PO_mean` >= seuil_99, na.rm = TRUE),
    mean_spm         = mean(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99], na.rm = TRUE),
    sd_spm           = sd(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99],   na.rm = TRUE),
    median_spm       = median(`SPM-G-PO_mean`[`SPM-G-PO_mean` >= seuil_99], na.rm = TRUE),
    aire_panache_km2 = pixel_count * aire_pixel_km2,
    .groups          = "drop"
  )

# identification des crues ------------------------------------------------

seuil_var     <- 121
seuil_paillon <- 85
seuil_magnan  <- 10

All_debit_olci <- All_debit |>
  filter(date >= as.Date("2016-04-26"), date <= as.Date("2024-12-31")) |>
  drop_na(debit_cumule)

crues_olci <- All_debit_olci |>
  mutate(
    en_crue_var     = débit          >= seuil_var,
    en_crue_paillon = ABA_debit_mean >= seuil_paillon,
    en_crue_magnan  = AAM_debit_mean >= seuil_magnan,
    en_crue_any     = en_crue_var | en_crue_paillon | en_crue_magnan
  )

dates_crue_olci <- crues_olci |>
  filter(en_crue_any) |>
  pull(date)

# jointure débit + panache ------------------------------------------------

data_crue_panache_olci <- merge(
  crues_olci,
  panache_journalier_olci,
  by  = "date",
  all = FALSE
)

# calcul jours_depuis_crue + période --------------------------------------

data_crue_panache_olci <- data_crue_panache_olci |>
  rowwise() |>
  mutate(
    jours_depuis_crue = {
      diffs <- as.numeric(date - dates_crue_olci)
      diffs[abs(diffs) == min(abs(diffs))][1]
    }
  ) |>
  ungroup() |>
  mutate(
    periode = factor(case_when(
      en_crue_any                                        ~ "Crue",
      jours_depuis_crue >= -3 & jours_depuis_crue < 0   ~ "J-1 à J-3",
      jours_depuis_crue >  0  & jours_depuis_crue <= 7  ~ "J+1 à J+7",
      TRUE                                               ~ "Hors crue"
    ), levels = c("J-1 à J-3", "Crue", "J+1 à J+7", "Hors crue"))
  )

# statistiques par période ------------------------------------------------

data_crue_panache_olci |>
  group_by(periode) |>
  summarise(
    mean_panache   = mean(aire_panache_km2,   na.rm = TRUE),
    median_panache = median(aire_panache_km2, na.rm = TRUE),
    sd_panache     = sd(aire_panache_km2,     na.rm = TRUE),
    n              = n()
  )

# tests statistiques ------------------------------------------------------

kruskal.test(aire_panache_km2 ~ periode, data = data_crue_panache_olci)

dunn.test(data_crue_panache_olci$aire_panache_km2,
          data_crue_panache_olci$periode,
          method = "bonferroni")

# graphiques --------------------------------------------------------------

## courbe de réponse jour par jour ----------------------------------------

courbe_reponse_olci <- data_crue_panache_olci |>
  filter(jours_depuis_crue >= -5, jours_depuis_crue <= 10) |>
  group_by(jours_depuis_crue) |>
  summarise(
    mean_panache = mean(aire_panache_km2, na.rm = TRUE),
    sd_panache   = sd(aire_panache_km2,   na.rm = TRUE),
    .groups      = "drop"
  )

ggplot(courbe_reponse_olci, aes(x = jours_depuis_crue, y = mean_panache)) +
  # Zones colorées avant/après
  annotate("rect",
           xmin = -5, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "steelblue", alpha = 0.05) +
  annotate("rect",
           xmin = 0, xmax = 10,
           ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.05) +
  # Enveloppe écart-type
  geom_ribbon(aes(ymin = pmax(mean_panache - sd_panache, 0),
                  ymax = mean_panache + sd_panache),
              fill = "steelblue", alpha = 0.2) +
  # Ligne et points
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_point(color = "white",    size = 4) +
  geom_point(color = "steelblue", size = 2.8) +
  # Ligne verticale début crue
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "firebrick", linewidth = 0.9) +
  # Labels zones
  annotate("text",
           x = -2.5,
           y = max(courbe_reponse_olci$mean_panache, na.rm = TRUE) * 0.98,
           label = "Avant crue", color = "steelblue",
           fontface = "italic", size = 3.5, hjust = 0.5) +
  annotate("text",
           x = 5,
           y = max(courbe_reponse_olci$mean_panache, na.rm = TRUE) * 0.98,
           label = "Après crue", color = "firebrick",
           fontface = "italic", size = 3.5, hjust = 0.5) +
  annotate("text",
           x = 0.3,
           y = max(courbe_reponse_olci$mean_panache, na.rm = TRUE) * 0.88,
           label = "J = 0\n(début crue)",
           color = "firebrick", hjust = 0, size = 3.2, fontface = "italic") +
  scale_x_continuous(breaks = seq(-5, 10, by = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title    = "Réponse des panaches turbides aux épisodes de crue",
    subtitle = "Produit OLCI — ODATIS-MR (2016–2024) · enveloppe = ± 1 écart-type",
    x        = "Jours par rapport au début de la crue",
    y        = expression("Aire moyenne des panaches (km²)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 18, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 15, margin = margin(b = 8)),
    axis.text        = element_text(color = "grey30", size = 15),
    axis.title       = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    panel.border     = element_rect(color = "grey60", linewidth = 0.5)
  )

## boxplot par période -----------------------------------------------------

medianes_olci <- data_crue_panache_olci |>
  filter(!is.na(periode), aire_panache_km2 > 0) |>
  group_by(periode) |>
  summarise(
    mediane = median(aire_panache_km2, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  )

ggplot(data_crue_panache_olci |> filter(!is.na(periode), aire_panache_km2 > 0),
       aes(x = periode, y = aire_panache_km2, fill = periode)) +
  geom_boxplot(
    alpha         = 0.75,
    outlier.alpha = 0.2,
    outlier.size  = 1,
    linewidth     = 0.4,
    color         = "grey30"
  ) +
  geom_text(
    data    = medianes_olci,
    aes(x = periode, y = mediane,
        label = paste0("n = ", n)),
    vjust       = -0.6,
    size        = 3.2,
    color       = "grey20",
    fontface    = "italic",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(
    "J-1 à J-3" = "steelblue",
    "Crue"       = "firebrick",
    "J+1 à J+7" = "darkorange",
    "Hors crue"  = "grey70"
  )) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100, 1000),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    title    = "Extension des panaches turbides selon la période de crue",
    subtitle = "Produit OLCI — ODATIS-MR (2016–2024) · échelle logarithmique · valeurs nulles exclues",
    x        = NULL,
    y        = expression("Aire des panaches (km²)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 18, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 15, margin = margin(b = 8)),
    axis.text.x      = element_text(size = 15, color = "grey20", face = "bold"),
    axis.text.y      = element_text(size = 15, color = "grey30"),
    axis.title.y     = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    panel.border     = element_rect(color = "grey60", linewidth = 0.5),
    legend.position  = "none"
  )
