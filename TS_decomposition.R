# TS_decomposition.R

# This script provides a workflow for decomposing a time series
# using the STL and X11 methods

# STL provides a decomposed time series with a more-or-less fixed seasonal signal
# While the X11 method allows for a changing seasonal signal
# Note however that X11 does not work with daily data

# A very nice explanation, with code, of the differences between the two methods can be found here:
# https://rpubs.com/oddskid/ts_decomposition


# Setup -------------------------------------------------------------------

# Note, if there are any errors loading these packages, make sure they have been installed first

library(tidyverse) # For basic data manipulation and plotting
library(zoo) # For various TS utilities
library(fpp3) # For more advanced tidyverse-like TS modelling workflows # NB: This is a big package, so it may take a while to install and load.
library(GGally) # For ggplot2 functionality with advanced TS models
library(seasonal) # For X11


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


# Data --------------------------------------------------------------------

# TS of Rhone débit
rhone_debit <- read_csv("data/grand_rhone_debit.csv") |> 
  mutate(debit = zoo::na.approx(debit))

# Plot
ggplot(rhone_debit, aes(x = date, y = debit)) + geom_line()

# TS of Rhone panache size
rhone_panache <- read_csv("data/grand_rhone_plume.csv") |> 
  dplyr::select(date, area_of_the_plume_mask_in_km2) |> 
  dplyr::rename(area = area_of_the_plume_mask_in_km2) |> 
  # Remove two anomolous days of data
  filter(area <= 10000) |> 
  # Ensure that all date values are present and fill in the NA values
  complete(date = seq(min(date), max(date), by = "day"), fill = list(value = NA)) |> 
  # NB: Be careful here with how this interpolates the missing panache data
  mutate(area = zoo::na.approx(area))
  # If it interpolates between two peaks this will be incorrect (check the plot to see if this is true)
  # It would rather be better to replace missing values with 0
  # mutate(area = ifelse(is.na(area), 0, area))

# Plot
ggplot(rhone_panache, aes(x = date, y = area)) + geom_line()


# fpp3::us_employment

# tsibble -----------------------------------------------------------------

# The following models require that the time series be 'tsibble' objects
# This also requires that there are no NA values, which was addressed above

# For STL these may be left as daily objects
rhone_debit_ts <- rhone_debit |> tsibble(index = date)
rhone_panache_ts <- rhone_panache |> tsibble(index = date)

# For X11 this requires that they be monthly data
rhone_debit_monthly <- rhone_debit |> 
  mutate(month = yearmonth(date)) |> 
  summarise(debit = mean(debit, na.rm = TRUE), .by = "month")
rhone_debit_monthly_ts <- rhone_debit_monthly |> tsibble(index = month)
rhone_panache_monthly <- rhone_panache |> 
  mutate(month = yearmonth(date)) |> 
  summarise(area = mean(area, na.rm = TRUE), .by = "month")
rhone_panache_monthly_ts <- rhone_panache_monthly |> tsibble(index = month)


# STL ---------------------------------------------------------------------

# Run the STL analysis
rhone_debit_stl <- rhone_debit_ts |> model(stl = STL(debit))
rhone_panache_stl <- rhone_panache_ts |> model(stl = STL(area))

# Basic plot of output
components(rhone_debit_stl) |> autoplot()
components(rhone_panache_stl) |> autoplot()

# Extract components for comparison/plotting
rhone_debit_wide <- rhone_debit |> 
  mutate(STL_seas = rhone_debit_stl[[1]][[1]]$fit$decomposition$season_year,
         STL_inter = rhone_debit_stl[[1]][[1]]$fit$decomposition$trend)
rhone_panache_wide <- rhone_panache |> 
  mutate(STL_seas = rhone_panache_stl[[1]][[1]]$fit$decomposition$season_year,
         STL_inter = rhone_panache_stl[[1]][[1]]$fit$decomposition$trend)

# Calculate linear models and prepare labels for plotting
rhone_debit_trend <- rhone_debit_wide |> 
  summarise(slope = coef(lm(STL_inter ~ date))["date"] * 365.25) |> 
  # NB: These are the x and y coordinates at which the label will be plotted
  mutate(x = as.Date("2015-01-01"), y = 1300)
rhone_panache_trend <- rhone_panache_wide |> 
  summarise(slope = coef(lm(STL_inter ~ date))["date"] * 365.25) |> 
  # NB: These are the x and y coordinates at which the label will be plotted
  mutate(x = as.Date("1998-01-01"), y = 1300)

# Scale the columns for better plotting
rhone_STL_seas_scaling_factor <- sec_axis_adjustement_factors(var_to_scale = rhone_debit_wide$STL_seas, 
                                                               var_ref = rhone_panache_wide$STL_seas)
rhone_STL_inter_scaling_factor <- sec_axis_adjustement_factors(var_to_scale = rhone_debit_wide$STL_inter, 
                                                               var_ref = rhone_panache_wide$STL_inter)

# Make the conversion
rhone_debit_wide <- rhone_debit_wide |> 
  mutate(STL_inter_scaled = STL_inter * rhone_STL_inter_scaling_factor$diff + rhone_STL_inter_scaling_factor$adjust,
         STL_seas_scaled = STL_seas * rhone_STL_seas_scaling_factor$diff + rhone_STL_seas_scaling_factor$adjust)

# Time series comparison plot
ggplot() + 
  
  # Dot and line plot for Panache
  # NB: Switch STL_seas to STL_trend to see the different variables
  geom_point(data = rhone_panache_wide, aes(x = date, y = STL_inter), color = "brown", alpha = 0.7) + 
  geom_path(data = rhone_panache_wide, aes(x = date, y = STL_inter), color = "brown", alpha = 0.7) + 
  
  # Dot and line plot for débit
  # NB: Switch STL_seas_scaled to STL_trend_scaled to see the different variables
  geom_point(data = rhone_debit_wide, aes(x = date, y = STL_inter_scaled), color = "blue", alpha = 0.7) + 
  geom_path(data = rhone_debit_wide, aes(x = date, y = STL_inter_scaled), color = "blue", alpha = 0.7) + 
  
  # Labels for linear model stats
  geom_label(data = rhone_panache_trend,  show.legend = FALSE, colour = "brown", size = 6, hjust = 0,
             aes(x = x, y = y, label = paste0("Panache area trend = ",round(slope, 2), " km-2 yr-1", sep = ""))) +
  geom_label(data = rhone_debit_trend,  show.legend = FALSE, colour = "blue", size = 6, hjust = 0,
             aes(x = x, y = y, label = paste0("River flow trend = ",round(slope, 2), " m-3 s-1 yr-1", sep = ""))) +
  
  scale_y_continuous(name = "Plume area (km²)",
                     sec.axis = sec_axis(transform = ~ {. - rhone_STL_inter_scaling_factor$adjust} / rhone_STL_inter_scaling_factor$diff, 
                                         name = "River flow (m³/s)")) +
  
  labs(title = "Interannual STL decomposition of Rhône plume size (Y1; left) and river flow (Y2; right)") +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle = element_text(hjust = 0.5),
        plot.title = element_text(size = 30, colour = "black"),
        text = element_text(size = 25, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 30, colour = "black")) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y.left = element_text(color = "brown"), 
        axis.ticks.y.left = element_line(color = "brown"),
        axis.line.y.left = element_line(color = "brown"),
        axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
        
        axis.text.y.right = element_text(color = "blue"), 
        axis.ticks.y.right = element_line(color = "blue"),
        axis.line.y.right = element_line(color = "blue"),
        axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
        
        panel.border = element_rect(linetype = "solid", fill = NA))

# Scatterplot of river flow versus panache area - this doesn't look incredible...
# Merge both data.frames
left_join(rhone_debit_wide |> select(date, STL_inter), rhone_panache_wide |> select(date, STL_inter), by = "date") |>
  dplyr::rename(stl_inter_river_flow = STL_inter.x, stl_inter_panache_area = STL_inter.y) |>
  
  # Start plot
  ggplot(aes(x = stl_inter_river_flow, y = stl_inter_panache_area)) + 
  
  # Add points
  geom_point() + 
  
  # Add linear model
  geom_smooth(method = "lm") + 
  
  # Force equal coords
  coord_equal() +
   
  # Add labels
  labs(x = "River flow (m³/s)", y = "Plume area (km²)", title = "Scatterplot of interannual trends in river flow and plume size") +
  
  # Theme tweaks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y.left = element_text(color = "brown"), 
        axis.ticks.y.left = element_line(color = "brown"),
        axis.line.y.left = element_line(color = "brown"),
        axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
        
        axis.text.x.bottom = element_text(color = "blue"), 
        axis.ticks.x.bottom = element_line(color = "blue"),
        axis.line.x.bottom = element_line(color = "blue"),
        axis.title.x.bottom = element_text(color = "blue", margin = unit(c(7.5, 0, 0, 0), "mm")),
        
        panel.border = element_rect(linetype = "solid", fill = NA))


# X11 ---------------------------------------------------------------------

# Perform X11 on the 'tsibble' objects
rhone_debit_x11 <- rhone_debit_monthly_ts |>
  model(x11 = X_13ARIMA_SEATS(debit ~ x11()))
rhone_panache_x11 <- rhone_panache_monthly_ts |>
  model(x11 = X_13ARIMA_SEATS(area ~ x11()))

# Basic plot
components(rhone_debit_x11) |> autoplot()
components(rhone_panache_x11) |> autoplot()

# Extract components for comparison/plotting
rhone_debit_monthly_wide <- rhone_debit_monthly |> 
  mutate(X11_seas = rhone_debit_x11[[1]][[1]][["fit"]][["fit"]][["data"]][,2],
         X11_inter = rhone_debit_x11[[1]][[1]][["fit"]][["fit"]][["data"]][,4])
rhone_panache_monthly_wide <- rhone_panache_monthly |> 
  mutate(X11_seas = rhone_panache_x11[[1]][[1]][["fit"]][["fit"]][["data"]][,2],
         X11_inter = rhone_panache_x11[[1]][[1]][["fit"]][["fit"]][["data"]][,4])

# Calculate linear models and prepare labels for plotting
rhone_debit_trend_X11 <- rhone_debit_monthly_wide |> 
  summarise(slope = coef(lm(X11_inter ~ month))["month"] * 365.25) |>
  # NB: These are the x and y coordinates at which the label will be plotted
  mutate(x = as.Date("2015-01-01"), y = 1300)
rhone_panache_trend_X11 <- rhone_panache_monthly_wide |> 
  summarise(slope = coef(lm(X11_inter ~ month))["month"] * 365.25) |> 
  # NB: These are the x and y coordinates at which the label will be plotted
  mutate(x = as.Date("1998-01-01"), y = 1300)

# Scale the columns for better plotting
rhone_X11_seas_scaling_factor <- sec_axis_adjustement_factors(var_to_scale = rhone_debit_monthly_wide$X11_seas, 
                                                              var_ref = rhone_panache_monthly_wide$X11_seas)
rhone_X11_inter_scaling_factor <- sec_axis_adjustement_factors(var_to_scale = rhone_debit_monthly_wide$X11_inter, 
                                                               var_ref = rhone_panache_monthly_wide$X11_inter)

# Make the conversion
rhone_debit_monthly_wide <- rhone_debit_monthly_wide |> 
  mutate(X11_inter_scaled = X11_inter * rhone_X11_inter_scaling_factor$diff + rhone_X11_inter_scaling_factor$adjust,
         X11_seas_scaled = X11_seas * rhone_X11_seas_scaling_factor$diff + rhone_X11_seas_scaling_factor$adjust)

# Time series comparison plot
ggplot() + 
  
  # Dot and line plot for Panache
  # NB: Switch STL_seas to STL_trend to see the different variables
  geom_point(data = rhone_panache_monthly_wide, aes(x = month, y = X11_inter), color = "brown", alpha = 0.7) + 
  geom_path(data = rhone_panache_monthly_wide, aes(x = month, y = X11_inter), color = "brown", alpha = 0.7) + 
  
  # Dot and line plot for débit
  # NB: Switch STL_seas_scaled to STL_trend_scaled to see the different variables
  geom_point(data = rhone_debit_monthly_wide, aes(x = month, y = X11_inter_scaled), color = "blue", alpha = 0.7) + 
  geom_path(data = rhone_debit_monthly_wide, aes(x = month, y = X11_inter_scaled), color = "blue", alpha = 0.7) + 
  
  # Labels for linear model stats
  geom_label(data = rhone_panache_trend_X11,  show.legend = FALSE, colour = "brown", size = 6, hjust = 0,
             aes(x = x, y = y, label = paste0("Panache area trend = ",round(slope, 2), " km-2 yr-1", sep = ""))) +
  geom_label(data = rhone_debit_trend_X11,  show.legend = FALSE, colour = "blue", size = 6, hjust = 0,
             aes(x = x, y = y, label = paste0("River flow trend = ",round(slope, 2), " m-3 s-1 yr-1", sep = ""))) +
  
  scale_y_continuous(name = "Plume area (km²)",
                     sec.axis = sec_axis(transform = ~ {. - rhone_X11_inter_scaling_factor$adjust} / rhone_X11_inter_scaling_factor$diff, 
                                         name = "River flow (m³/s)")) +
  
  labs(title = "Interannual STL decomposition of Rhône plume size (Y1; left) and river flow (Y2; right)") +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.subtitle = element_text(hjust = 0.5),
        plot.title = element_text(size = 30, colour = "black"),
        text = element_text(size = 25, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 30, colour = "black")) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y.left = element_text(color = "brown"), 
        axis.ticks.y.left = element_line(color = "brown"),
        axis.line.y.left = element_line(color = "brown"),
        axis.title.y.left = element_text(color = "brown", margin = unit(c(0, 7.5, 0, 0), "mm")),
        
        axis.text.y.right = element_text(color = "blue"), 
        axis.ticks.y.right = element_line(color = "blue"),
        axis.line.y.right = element_line(color = "blue"),
        axis.title.y.right = element_text(color = "blue", margin = unit(c(0, 0, 0, 7.5), "mm")),
        
        panel.border = element_rect(linetype = "solid", fill = NA))

