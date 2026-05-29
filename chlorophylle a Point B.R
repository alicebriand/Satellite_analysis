# chlorophylle a Point B

library(tidyverse)
library(trend)
library(seasonal)

# ── Chargement ──
somlit <- read_delim(
  "~/Downloads/Somlit_Extraction_ctd_20260529_145836_7ea626b4e98a8859.csv",
  delim     = ";",
  skip      = 3,          # sauter les 3 lignes d'en-tête
  col_names = c("ID_SITE", "DATE", "HEURE", "nomSite",
                "gpsLong", "gpsLat", "FLUORESCENCE", "PROFONDEUR"),
  col_types = cols(
    DATE          = col_date(format = "%Y-%m-%d"),
    FLUORESCENCE  = col_double(),
    PROFONDEUR    = col_double()
  ),
  locale = locale(decimal_mark = ".")
) |>
  filter(row_number() > 1)   # supprimer la ligne des unités (//)

# ── Nettoyage ──
somlit_clean <- somlit |>
  filter(
    nomSite      == "Point B",
    FLUORESCENCE != 999999,        # valeurs manquantes codées 999999
    FLUORESCENCE >= 0,
    !is.na(FLUORESCENCE)
  ) |>
  mutate(
    date       = as.Date(DATE),
    year       = year(date),
    month      = month(date),
    profondeur = as.numeric(PROFONDEUR)
  )

cat("Période :", as.character(range(somlit_clean$date)), "\n")
cat("Profondeurs :", range(somlit_clean$profondeur), "m\n")
cat("N observations :", nrow(somlit_clean), "\n")

# ── Couche de surface (0-10m) ──
somlit_surface <- somlit_clean |>
  filter(profondeur <= 10)

# ── Moyenne journalière surface ──
somlit_surface_daily <- somlit_surface |>
  group_by(date, year, month) |>
  summarise(
    mean_fluo   = mean(FLUORESCENCE,   na.rm = TRUE),
    max_fluo    = max(FLUORESCENCE,    na.rm = TRUE),
    median_fluo = median(FLUORESCENCE, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  drop_na(mean_fluo)

# ── Climatologie mensuelle surface ──
somlit_clim_surface <- somlit_surface |>
  group_by(month) |>
  summarise(
    mean_fluo   = mean(FLUORESCENCE,   na.rm = TRUE),
    median_fluo = median(FLUORESCENCE, na.rm = TRUE),
    sd_fluo     = sd(FLUORESCENCE,     na.rm = TRUE),
    .groups     = "drop"
  )

# ── Mann-Kendall surface ──
mk_fluo <- mk.test(somlit_surface_daily$mean_fluo)
sen_fluo <- sens.slope(somlit_surface_daily$mean_fluo)
cat("Mann-Kendall p =", round(mk_fluo$p.value, 4), "\n")
cat("Theil-Sen pente =", round(sen_fluo$estimates * 365, 4), "u.a./an\n")

# ════════════════════════════════════════════════════════
# GRAPHIQUES
# ════════════════════════════════════════════════════════

mois_fr <- c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
             "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")

## 1. Série temporelle surface ────────────────────────────
ggplot(somlit_surface_daily, aes(x = date, y = mean_fluo)) +
  geom_line(color = "chartreuse4", linewidth = 0.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "firebrick",
              fill = "firebrick", alpha = 0.15, linewidth = 1) +
  annotate("text",
           x     = max(somlit_surface_daily$date),
           y     = max(somlit_surface_daily$mean_fluo, na.rm = TRUE) * 0.97,
           hjust = 1, vjust = 1, size = 3.8, fontface = "italic",
           color = "grey20",
           label = paste0(
             "Mann-Kendall p ",
             ifelse(mk_fluo$p.value < 0.001, "< 0.001",
                    paste0("= ", round(mk_fluo$p.value, 3)))
           )) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "Série temporelle de la fluorescence au Point B",
    subtitle = "Couche de surface (0–10 m) · moyenne journalière · SOMLIT",
    x        = NULL,
    y        = "Fluorescence (u.a.)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )

## 2. Climatologie mensuelle surface ──────────────────────
somlit_clim_surface |>
  mutate(month_label = factor(mois_fr[month], levels = mois_fr)) |>
  ggplot(aes(x = month_label, y = mean_fluo, group = 1)) +
  geom_ribbon(aes(ymin = pmax(mean_fluo - sd_fluo, 0),
                  ymax = mean_fluo + sd_fluo),
              fill = "chartreuse4", alpha = 0.2) +
  geom_line(color  = "chartreuse4", linewidth = 1.1) +
  geom_point(color = "white",       size = 4) +
  geom_point(color = "chartreuse4", size = 2.8) +
  labs(
    title    = "Climatologie mensuelle de la fluorescence au Point B",
    subtitle = "Couche de surface (0–10 m) · SOMLIT · enveloppe = ± 1 écart-type",
    x        = NULL,
    y        = "Fluorescence moyenne (u.a.)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )


# X11 ---------------------------------------------------------------------

# ── Agréger en mensuel ──
somlit_mensuel <- somlit_surface_daily |>
  mutate(mois = floor_date(date, "month")) |>
  group_by(mois) |>
  summarise(fluo_mois = mean(mean_fluo, na.rm = TRUE), .groups = "drop") |>
  drop_na(fluo_mois)

cat("Période :", as.character(range(somlit_mensuel$mois)), "\n")
cat("N mois   :", nrow(somlit_mensuel), "\n")

# ── Série temporelle mensuelle ──
fluo_ts <- ts(
  data      = somlit_mensuel$fluo_mois,
  start     = c(year(min(somlit_mensuel$mois)),
                month(min(somlit_mensuel$mois))),
  frequency = 12
)

# ── Décomposition X11 ──
x11_fluo <- seasonal::seas(fluo_ts, x11 = "")
summary(x11_fluo)

# ── Extraire les composantes ──
composantes_fluo <- data.frame(
  date         = somlit_mensuel$mois,
  observed     = as.numeric(original(x11_fluo)),
  tendance     = as.numeric(trend(x11_fluo)),
  saisonnalite = as.numeric(series(x11_fluo, "d10")),
  residus      = as.numeric(irregular(x11_fluo))
)

# ── Climatologie mensuelle de la saisonnalité ──
mois_fr <- c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
             "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")

saisonnalite_clim_fluo <- composantes_fluo |>
  mutate(
    month_num   = month(date),
    month_label = factor(mois_fr[month_num], levels = mois_fr)
  ) |>
  group_by(month_label) |>
  summarise(
    mean_sais = mean(saisonnalite, na.rm = TRUE),
    min_sais  = min(saisonnalite,  na.rm = TRUE),
    max_sais  = max(saisonnalite,  na.rm = TRUE),
    .groups   = "drop"
  )

# ════════════════════════════════════════════════════════
# GRAPHIQUES
# ════════════════════════════════════════════════════════

## 1. Décomposition X11 complète ──────────────────────────

p1 <- ggplot(composantes_fluo, aes(x = date)) +
  geom_line(aes(y = observed, color = "Signal brut"),
            linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = tendance, color = "Tendance"),
            linewidth = 1.1) +
  scale_color_manual(values = c("Signal brut" = "chartreuse4",
                                "Tendance"    = "firebrick")) +
  labs(title = "a) Signal observé et tendance",
       x = NULL, y = "Fluorescence (u.a.)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        legend.position  = "top",
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

p2 <- ggplot(composantes_fluo, aes(x = date, y = saisonnalite)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(saisonnalite, 1), ymax = 1),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(saisonnalite, 1)),
              fill = "chartreuse4", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  labs(title = "b) Composante saisonnière",
       x = NULL, y = "Facteur saisonnier") +
  theme_bw(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

p3 <- ggplot(composantes_fluo, aes(x = date, y = residus)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = pmin(residus, 1), ymax = 1),
              fill = "steelblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(residus, 1)),
              fill = "tomato", alpha = 0.3) +
  geom_line(color = "grey30", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "c) Résidus (irrégulier)",
       x = NULL, y = "Facteur") +
  theme_bw(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1))

(p1 / p2 / p3) +
  plot_annotation(
    title    = "Décomposition X11 de la fluorescence au Point B",
    subtitle = "Couche de surface (0–10 m) · SOMLIT",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey50")
    )
  )

## 2. Saisonnalité moyenne annuelle avec min/max ──────────

ggplot(saisonnalite_clim_fluo,
       aes(x = month_label, y = mean_sais, group = 1)) +
  geom_ribbon(aes(ymin = min_sais, ymax = max_sais),
              fill = "chartreuse4", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_line(color  = "chartreuse4", linewidth = 1.1) +
  geom_point(color = "white",       size = 4) +
  geom_point(color = "chartreuse4", size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(
    title    = "Saisonnalité X11 de la fluorescence au Point B",
    subtitle = "Couche de surface (0–10 m) · enveloppe = min/max interannuel · SOMLIT",
    x        = NULL,
    y        = "Facteur saisonnier moyen"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.text.x      = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )
