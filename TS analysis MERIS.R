# TS analysis MERIS

# 19/03/2026

# pathway : "~/Satellite_analysis/TS analysis MERIS.R


# This script will load MERIS satellite data
# Then perform a temporal analysis
# and then compare it to the runoff of the Var river from 2002 to 2012


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = detectCores()-2) # Detects cores automagically
library(heatwaveR)
library(ggpmisc)

# load function ---------------------------------------------------------------

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

# load data ---------------------------------------------------------------

load("data/MERIS/MERIS_2002_2012.Rdata")
load("data/Hydro France/Y6442010_depuis_2000.Rdata")

# data-treatment ----------------------------------------------------------

MERIS_2002_2012 <- MERIS_2002_2012 |> 
  filter(mean_spm <= 1000)

# Complete missing dates in the date range
MERIS_2002_2012 <- MERIS_2002_2012 %>%
  complete(date = seq(min(date), max(date), by = "day"))

Y6442010_2002_2012 <- Y6442010_depuis_2000 |> 
  filter(date >= as.Date("2002-07-21"), date <= as.Date("2012-04-07"))
  
# plotting ----------------------------------------------------------------

## MERIS only ----------------------------------------------------------------

    # en échelle normale

model_MERIS_2002_mean <- lm(mean_spm ~ date, data = MERIS_2002_2012)
p_value_MERIS_2002_mean <- summary(model_MERIS_2002_mean)$coefficients[2, 4]  # p-value pour la pente
intercept_MERIS_2002_mean <- coef(model_MERIS_2002_mean)[1]
slope_MERIS_2002_mean <- coef(model_MERIS_2002_mean)[2]

ggplot(data = MERIS_2002_2012, aes(x = date, y = mean_spm)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkslateblue", fill = "pink", alpha = 0.2) +
  geom_point(color = "red3") +
  annotate(
    "text",
    x = max(MERIS_2002_2012$date, na.rm = TRUE),
    y = max(MERIS_2002_2012$mean_spm, na.rm = TRUE) * 0.9,
    label = paste0(
      "y = ", round(intercept_MERIS_2002_mean, 3), " + ", round(slope_MERIS_2002_mean, 7), " * x",
      "\n", "p = ", ifelse(p_value_MERIS_2002_mean < 0.001, "< 0.001", format(p_value_MERIS_2002_mean, digits = 3))
    ),
    hjust = 1,  # Alignement à droite
    vjust = 1,  # Alignement en haut
    size = 6
  ) +
  labs(title = "Évolution de la concentration en matière particulaire en suspension moyenne entre 2002 et 2012 avec le produit MERIS issu de ODATIS-MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )



    # en échelle log

MERIS_filtered <- MERIS_2002_2012 %>%
  filter(mean_spm > 0 & !is.na(mean_spm))

model_log_MERIS_2002_mean <- lm(log10(mean_spm) ~ date, data = MERIS_filtered)
p_value_log_MERIS_2002_mean <- summary(model_log_MERIS_2002_mean)$coefficients[2, 4]

intercept_log_MERIS_2002_mean <- coef(model_log_MERIS_2002_mean)[1]
slope_log_MERIS_2002_mean <- coef(model_log_MERIS_2002_mean)[2]
p_value_log_MERIS_2002_mean <- summary(model_log_MERIS_2002_mean)$coefficients[2, 4]

# Formater l'équation
equation_text_log <- paste0(
  "log10(y) = ", round(intercept_log_MERIS_2002_mean, 4),
  ifelse(sign(slope_log_MERIS_2002_mean) == 1, " + ", " - "),
  abs(round(slope_log_MERIS_2002_mean, 4)), " * x",
  "\n",  # Saut de ligne
  "p-value = ", format.pval(p_value_log_MERIS_2002_mean, digits = 3)
)


# Graphique
ggplot(data = MERIS_filtered, aes(x = date, y = mean_spm)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "darkslateblue") +
  geom_point(color = "red3") +
  labs(
    title = "Évolution de la concentration en matière particulaire en suspension médiane entre 2002 et 2025 (échelle log) avec MERIS",
    x = "Date",
    y = "Concentration médiane (g/m³, échelle log)"
  ) +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_log10() +
  annotate(
    "text",
    x = as.Date("2010-01-01"),
    y = max(MERIS_filtered$mean_spm, na.rm = TRUE) * 0.8,
    label = equation_text_log,
    hjust = 0,
    vjust = 1,
    size = 5,
    color = "black"
  )

## comparison liquid flow vs SPM concentration -----------------------------

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(MERIS_2002_2012$mean_spm, Y6442010_2002_2012$débit)

MERIS_2002_2012$scaled_mean_spm <- MERIS_2002_2012$mean_spm * adjust_factors$diff + adjust_factors$adjust

# en échelle normale
ggplot() +
  geom_point(
    data = Y6442010_2002_2012,
    aes(x = date, y = débit, color = "Débit")
  ) +
  geom_point(
    data = MERIS_2002_2012,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "deeppink")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension moyenne (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension entre 2002 et 2012 avec le produit MERIS issu de ODATIS-MR",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# scatter plot ------------------------------------------------------------

ggplot(Var_MERIS_SPM, aes(x = débit, y = mean_spm)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "red3", fill = "pink", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Relation entre débit et concentration en SPM en échelle log",
    x = "Débit (m³/s)",
    y = "Concentration moyenne en SPM (g/m³)"
  ) +
  theme_minimal()


# runoff vs SPM concentration correlation ---------------------------------

Var_MERIS_SPM <- inner_join(Y6442010_2002_2012, MERIS_2002_2012, by = "date")

cor.test(Var_MERIS_SPM$débit, Var_MERIS_SPM$mean_spm, method = "spearman")
