# TS analysis OLCI

# 03/03/2026

# pathway : "~/Documents/Alice/Code/TS analysis OLCI


# This script will load OLCI satellite data
# Then perform a temporal analysis analysis
# and then compare it to the runoff of the 
# Var river at the Napoleon bridge


# Setup ------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(tidync)
library(gganimate)
library(doParallel); registerDoParallel(cores = 14)

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

# loading data ------------------------------------------------------------

load("data/OLCI/SPM/OLCI_A_2016_2024_SPM.Rdata")
load("data/OLCI/SPM/OLCI_B_2018_2024_SPM.Rdata")
load("data/OLCI/SPM/OLCI_A_B_2016_2024_SPM.Rdata")

load("data/OLCI/CHL/OLCI_A_2016_2024_CHL.Rdata")

load("data/Hydro France/Y6442010_2016_2024.Rdata")

# plotting ----------------------------------------------------------------
## SPM ----------------------------------------------------------------

### OLCI A ------------------------------------------------------------------

# Générer les dates de début de mois
dates_mois_OLCI_A_spm <- seq(
  from = as.Date("2016-04-26"),
  to = as.Date("2024-12-30"),
  by = "month"
)

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI
ggplot(data = OLCI_A_2016_2024_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  geom_vline(
    xintercept = dates_mois_OLCI_A_spm,
    color = "gray80",
    alpha = 0.5,
    linetype = "dashed",
    size = 0.2
  ) +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 2016 et 2024 avec le produit OLCI A issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

### OLCI B ------------------------------------------------------------------

# Générer les dates de début de mois
dates_mois_OLCI_B_spm <- seq(
  from = as.Date("2018-05-15"),
  to = as.Date("2024-12-28"),
  by = "month"
)

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI
ggplot(data = OLCI_B_2018_2024_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "red3") +
  geom_vline(
    xintercept = dates_mois_OLCI_B_spm,
    color = "gray80",
    alpha = 0.5,
    linetype = "dashed",
    size = 0.2
  ) +
  labs(title = "Evolution de la concentration en matière particulaire en suspension moyenne entre 2018 et 2024 avec le produit OLCI B issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

### OLCI global ------------------------------------------------------------------

# Générer les dates de début de mois
dates_mois_OLCI_A_B_spm <- seq(
  from = as.Date("2016-04-26"),
  to = as.Date("2024-12-30"),
  by = "month"
)

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI
ggplot(data = OLCI_A_B_2016_2024_SPM, aes(x = date, y = mean_spm)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_point(color = "red3", size = 0.5) +
  geom_vline(
    xintercept = dates_mois_OLCI_A_B_spm,
    color = "gray80",
    alpha = 0.5,
    linetype = "dashed",
    size = 0.2
  ) +
  labs(title = "Évolution de la concentration en matière particulaire en suspension moyenne entre 2016 et 2024 avec le produit OLCI issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en matière particulaire en suspension (en g/m³)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(OLCI_A_B_2016_2024_SPM$mean_spm, Y6442010_2016_2024$débit)

OLCI_A_B_2016_2024_SPM$scaled_mean_spm <- OLCI_A_B_2016_2024_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_point(
    data = Y6442010_2016_2024,
    aes(x = date, y = débit, color = "Débit"), size = 0.5
  ) +
  geom_point(
    data = OLCI_A_B_2016_2024_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM"), size = 0.5
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension entre 2016 et 2024 avec le produit OLCI issu de ODATIS MR",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# zoom sur une année, ex : 2023

OLCI_A_B_2023_SPM <- OLCI_A_B_2016_2024_SPM %>%
  dplyr::filter(date >= as.Date("2023-01-01"), date <= as.Date("2023-12-31"))

Y6442010_Hydro_complete_2023 <- Y6442010_Hydro_complete %>% 
  dplyr::filter(Date >= as.Date("2023-01-01"), Date <= as.Date("2023-12-31"))

adjust_factors <- sec_axis_adjustement_factors(OLCI_A_B_2023_SPM$mean_spm, Y6442010_Hydro_complete_2023$débit)

OLCI_A_B_2023_SPM$scaled_mean_spm <- OLCI_A_B_2023_SPM$mean_spm * adjust_factors$diff + adjust_factors$adjust

ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete_2023,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = OLCI_A_B_2023_SPM,
    aes(x = date, y = scaled_mean_spm, color = "SPM")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "SPM" = "red3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Matière particulaire en suspension (en g/m³)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en matière en suspension en 2023 avec le produit OLCI issu de ODATIS MR",
    x = "Date"
  ) +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

# scatter plot ------------------------------------------------------------

ggplot(Var_OLCI_SPM, aes(x = débit, y = mean_spm)) +
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

Var_OLCI_SPM <- inner_join(Y6442010_2016_2024, OLCI_A_B_2016_2024_SPM, by = "date")

cor.test(Var_OLCI_SPM$débit, Var_OLCI_SPM$mean_spm, method = "spearman")


## CHL ----------------------------------------------------------------
### OLCI A ------------------------------------------------------------------

# Générer les dates de début de mois
dates_mois_OLCI_A_chl <- seq(
  from = as.Date("2016-04-26"),
  to = as.Date("2024-12-30"),
  by = "month"
)

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI
ggplot(data = OLCI_A_2016_2024_CHL, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "chartreuse3") +
  labs(title = "Evolution de la concentration en chlorophylle moyenne entre 2016 et 2024 avec le produit OLCI A issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en chlorophylle (en µg/L)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )

### OLCI B ------------------------------------------------------------------

# Générer les dates de début de mois
dates_mois_OLCI_B_chl <- seq(
  from = as.Date("2018-04-26"),
  to = as.Date("2024-12-30"),
  by = "month"
)

# on plotte seulement la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI
ggplot(data = OLCI_B_2018_2023_CHL, aes(x = date, y = mean_chl)) +
  # geom_ribbon(aes(ymin = mean_spm - std_spm, ymax = mean_spm + std_spm,
  #                 alpha = 0.2, fill = "blue")) +
  geom_smooth(method = "lm", se = FALSE, color = "darkslateblue") +
  geom_line(color = "chartreuse3") +
  geom_vline(
    xintercept = dates_mois_OLCI_B_chl,
    color = "gray80",
    alpha = 0.5,
    linetype = "dashed",
    size = 0.2
  ) +
  labs(title = "Evolution de la concentration en chlorophylle moyenne entre 2018 et 2023 avec le produit OLCI B issu de ODATIS MR",
       x = "Date",
       y = "Concentration moyenne en chlorophylle (en µg/L)") +
  theme_minimal() +
  scale_x_date(
    date_breaks = "1 year",  
    date_labels = "%Y"       
  )


### OLCI global ------------------------------------------------------------------

# on plotte la série temporelle de la concentration moyenne en SPM entre
# 2015 et 2025 avec OLCI contre le débit liquide du Var 

# pour cela on a besoin de facteur d'ajustement : 

# adjusting scale
adjust_factors <- sec_axis_adjustement_factors(OLCI_A_B_2016_2024_CHL$mean_chl, Y6442010_Hydro_complete$débit)

OLCI_A_B_2016_2024_CHL$scaled_mean_chl <- OLCI_A_B_2016_2024_CHL$mean_chl * adjust_factors$diff + adjust_factors$adjust


ggplot() +
  geom_line(
    data = Y6442010_Hydro_complete,
    aes(x = Date, y = débit, color = "Débit")
  ) +
  geom_line(
    data = OLCI_A_B_2016_2024_CHL,
    aes(x = date, y = scaled_mean_chl, color = "CHL")
  ) +
  scale_color_manual(values = c("Débit" = "blue", "CHL" = "chartreuse3")) +
  scale_y_continuous(
    name = "Débit (m³/s)",
    sec.axis = sec_axis(~ (. - adjust_factors$adjust) / adjust_factors$diff, name = "Chlorophylle (en µg/L)")
  ) +
  labs(
    title = "Débit du Var au pont Napoléon et concentration en chlorophylle entre 2016 et 2024 avec le produit OLCI issu de ODATIS MR",
    x = "Date"
  ) +
  theme_minimal()
