library(tidyverse)
library(sf)

topofire = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv") %>%
  mutate(model = 'Topofire')

soilwat2 = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv") %>%
  mutate(model = 'SOILWAT2')

all_results = bind_rows(topofire, soilwat2)

# plot rneonSoilFlux# plot results
states = read_sf("~/mco-drought-indicators/processing/base-data/raw/states.shp") %>%
  filter(!STATE_NAME %in% c('Hawaii', 'Alaska', 'Puerto Rico', 'Virgin Islands'))

station_data = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

plot_data = all_results %>%
  left_join(., station_data) %>%
  drop_na(latitude, longitude) %>%
  st_as_sf(., coords = c('longitude', 'latitude'), crs = 4326) %>%
  mutate(name_new = ifelse(generalized_depth == 'Deep', 'Deep (>50cm)', 
                           ifelse(generalized_depth == 'Shallow', 'Shallow (<= 10cm)', 
                                  ifelse(generalized_depth == 'Middle', 'Middle (>10cm - <=50cm)', 'Depth Averaged'))),
         name_new = factor(name_new, levels = c('Depth Averaged', 'Shallow', 
                                                'Shallow (<= 10cm)', 'Middle (>10cm - <=50cm)', 
                                                'Deep (>50cm)')))

make_plot = function(data, variable, states, 
                     model = NULL,
                     clamp = FALSE, 
                     export_path = NULL, 
                     height = 7, 
                     width = 8,
                     depth_avg_only = FALSE) {
  
  variable <- rlang::ensym(variable)  # Convert to symbol
  
  # Optional filter by model name
  if (!is.null(model)) {
    data <- data %>% filter(model == !!model)
  }
  
  # Optional filter to include only "Depth Averaged"
  if (depth_avg_only) {
    data <- data %>% filter(name_new == "Depth Averaged")
  }
  
  # Apply clamping to avoid negative values
  if (clamp) {
    data <- data %>% mutate(!!variable := pmax(!!variable, 0))
  }
  
  # Reproject geometries to EPSG:5070
  data <- sf::st_transform(data, 5070)
  states <- sf::st_transform(states, 5070)
  
  # Dynamic title
  plot_title <- if (!is.null(model)) paste(model, "Correlations") else "Model Correlations"
  subtitle_text <- if (depth_avg_only) "Depth Averaged Only (without NEON stations)" else "(without NEON stations)"
  
  # Create the plot
  cor_plot <- ggplot(data, aes(fill = !!variable)) +
    geom_sf(data = states, fill = 'transparent') +
    geom_sf(shape = 21, color = 'black') +
    scale_fill_gradientn(colours = viridis::turbo(5) %>% rev()) +
    facet_wrap(~name_new) + 
    theme_bw(base_size = 16) +
    guides(fill = guide_colourbar(title.position = "bottom", title.hjust = 0.5)) +
    theme(
      legend.key = element_blank(),
      strip.background = element_rect(colour = "transparent", fill = "transparent"),
      legend.position = 'bottom', 
      legend.key.width = unit(2, "cm"),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    ggtitle(plot_title, subtitle = subtitle_text)
  
  # Save the plot if export path is given
  if (!is.null(export_path)) {
    ggsave(export_path, plot = cor_plot, height = height, width = width)
    message("Plot saved to: ", export_path)
  }
  
  return(cor_plot)
}

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `Pearson's r`,
  states = states,
  model = "Topofire",
  clamp = TRUE,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/topofire-pearson-r-clamped.png'
)

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `Pearson's r`,
  states = states,
  model = "SOILWAT2",
  clamp = TRUE,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-pearson-r-clamped.png'
)

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `Bias`,
  states = states,
  model = "Topofire",
  clamp = F,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/topofire-bias.png'
)


# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `Bias`,
  states = states,
  model = "SOILWAT2",
  clamp = F,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-bias.png'
)

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `KGE`,
  states = states,
  model = "Topofire",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/topofire-kge.png'
)


# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `KGE`,
  states = states,
  model = "SOILWAT2",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-kge.png'
)

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `NSE`,
  states = states,
  model = "Topofire",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/topofire-nse.png'
)


# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `NSE`,
  states = states,
  model = "SOILWAT2",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-nse.png'
)

# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `RMSE`,
  states = states,
  model = "Topofire",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/topofire-rmse.png'
)


# Create a plot for the Topofire model, clamped, only depth-averaged data
make_plot(
  data = plot_data,
  variable = `RMSE`,
  states = states,
  model = "SOILWAT2",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-rmse.png'
)

model_diff_plot = function(data, variable, states,
                           model_1, model_2,
                           clamp_to_percentile = FALSE,
                           export_path = NULL,
                           height = 7, width = 8,
                           depth_avg_only = FALSE) {
  
  variable <- rlang::ensym(variable)
  var_label <- rlang::as_label(variable)
  
  # Define metrics where lower values are better
  lower_is_better <- c("rmse", "mae", "bias")
  metric_direction <- ifelse(tolower(var_label) %in% lower_is_better, "lower", "higher")
  
  # Filter data for selected models
  data_filtered <- data %>%
    filter(model %in% c(model_1, model_2))
  
  if (depth_avg_only) {
    data_filtered <- data_filtered %>%
      filter(name_new == "Depth Averaged")
  }
  
  # Extract selected variable
  summarized <- data_filtered %>%
    st_drop_geometry() %>%
    select(site_id, model, name_new, value = !!variable)
  
  # Pivot to wide and compute differences
  data_diff <- summarized %>%
    pivot_wider(names_from = model, values_from = value)
  
  if (tolower(var_label) == "bias") {
    data_diff <- data_diff %>%
      mutate(difference = abs(!!sym(model_2)) - abs(!!sym(model_1)))
    direction_label <- paste0("Δ|", var_label, "| (closer to 0 is better)")
  } else if (tolower(var_label) %in% lower_is_better) {
    data_diff <- data_diff %>%
      mutate(difference = !!sym(model_2) - !!sym(model_1))
    direction_label <- paste0("Δ", var_label, "")
  } else {
    data_diff <- data_diff %>%
      mutate(difference = !!sym(model_1) - !!sym(model_2))
    direction_label <- paste0("Δ", var_label, "")
  }
  
  # Remove NAs
  data_diff <- data_diff %>%
    filter(!is.na(difference))
  
  # Reattach geometry and depth labels
  site_meta <- data_filtered %>%
    select(site_id, name_new, geometry) %>%
    distinct()
  
  data_diff <- data_diff %>%
    left_join(site_meta, by = c("site_id", "name_new")) %>%
    st_as_sf()
  
  # Reproject
  data_diff <- st_transform(data_diff, 5070)
  states <- st_transform(states, 5070)
  
  # Clamp to 10th/90th percentiles (symmetric around 0)
  if (clamp_to_percentile) {
    p10 <- quantile(data_diff$difference, 0.10, na.rm = TRUE)
    p90 <- quantile(data_diff$difference, 0.90, na.rm = TRUE)
    clamp_val <- max(abs(c(p10, p90)))
    
    data_diff <- data_diff %>%
      mutate(difference_clamped = case_when(
        difference < -clamp_val ~ -clamp_val,
        difference >  clamp_val ~  clamp_val,
        TRUE ~ difference
      ))
    
    plot_var <- sym("difference_clamped")
    color_limits <- c(-clamp_val, clamp_val)
    legend_breaks <- c(-clamp_val, 0, clamp_val)
    max_label <- round(clamp_val, 2)
  } else {
    plot_var <- sym("difference")
    max_abs_diff <- max(abs(data_diff$difference), na.rm = TRUE)
    max_lim <- ceiling(max_abs_diff * 10) / 10
    color_limits <- c(-max_lim, max_lim)
    legend_breaks <- c(-max_lim, 0, max_lim)
    max_label <- max_lim
  }
  
  # Define fill labels (consistent color ramp)
  fill_labels <- c(
    paste0("-", max_label, "\n", model_2, " is better"),
    "0\nNo difference",
    paste0(max_label, "\n", model_1, " is better")
  )
  
  # Title logic
  plot_title <- paste(model_1, "vs", model_2, var_label, "Difference")
  subtitle_text <- if (depth_avg_only) "Depth Averaged Only (without NEON stations)" else "(without NEON stations)"
  
  # Plot
  diff_plot <- ggplot(data_diff, aes(fill = !!plot_var)) +
    geom_sf(data = states, fill = 'transparent') +
    geom_sf(shape = 21, color = 'black') +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue", midpoint = 0,
      name = direction_label,
      limits = color_limits,
      breaks = legend_breaks,
      labels = fill_labels
    ) +
    facet_wrap(~name_new) +
    theme_bw(base_size = 16) +
    guides(fill = guide_colourbar(title.position = "bottom", title.hjust = 0.5)) +
    theme(
      legend.key = element_blank(),
      strip.background = element_rect(colour = "transparent", fill = "transparent"),
      legend.position = 'bottom',
      legend.key.width = unit(2, "cm"),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    ggtitle(plot_title, subtitle = subtitle_text)
  
  # Export
  if (!is.null(export_path)) {
    ggsave(export_path, plot = diff_plot, height = height, width = width)
    message("Plot saved to: ", export_path)
  }
  
  return(diff_plot)
}

model_diff_plot(
  data = plot_data,
  variable = 'KGE',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "Topofire",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-topofire-vs-soilwat2-kge.png"
)

model_diff_plot(
  data = plot_data,
  variable = "Pearson's r",
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "Topofire",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-topofire-vs-soilwat2-r.png"
)

model_diff_plot(
  data = plot_data,
  variable = 'Bias',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "Topofire",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-topofire-vs-soilwat2-bias.png"
)

model_diff_plot(
  data = plot_data,
  variable = 'RMSE',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "Topofire",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-topofire-vs-soilwat2-RMSE.png"
)

model_diff_plot(
  data = plot_data,
  variable = 'NSE',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "Topofire",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-topofire-vs-soilwat2-RMSE.png"
)


# hitsogram plot

histogram_plot = function(data, variable, network, 
                          model = NULL,
                          clamp = FALSE, 
                          export_path = NULL, 
                          height = 7, 
                          width = 8,
                          depth_avg_only = FALSE) {
  
  variable = rlang::ensym(variable)  # Convert variable to a symbol
  network = rlang::ensym(network)    # Convert network column to a symbol
  
  # Optional filter by model
  if (!is.null(model)) {
    data <- data %>% filter(model == !!model)
  }
  
  # Optional filter to include only "Depth Averaged"
  if (depth_avg_only) {
    data <- data %>% filter(name_new == "Depth Averaged")
  }
  
  # Apply clamping
  if (clamp) {
    data <- data %>% mutate(!!variable := pmax(!!variable, 0))
  }
  
  # Dynamic title
  plot_title <- if (!is.null(model)) paste(model, "Correlations") else "Model Correlations"
  subtitle_text <- if (depth_avg_only) "Depth Averaged Only (without NEON stations)" else "(without NEON stations)"
  
  # Generate the histogram
  hist_plot = ggplot(data, aes(x = !!variable, color = !!network)) +
  geom_density(alpha = 0.6) +
    facet_wrap(~name_new) + 
    theme_bw(base_size = 16) +
    theme(
      legend.key = element_blank(),
      strip.background = element_rect(colour = "transparent", fill = "transparent"),
      legend.position = 'bottom', 
      legend.key.width = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    ggtitle(plot_title, subtitle = subtitle_text) +
    labs(y = 'Station Density', fill = NULL, color = NULL)
  
  if (!is.null(export_path)) {
    ggsave(export_path, plot = hist_plot, height = height, width = width)
    message("Plot saved to: ", export_path)
  }
  
  return(hist_plot)
}

histogram_plot(
  plot_data, `Pearson's r`, network,
  model = "Topofire",
  clamp = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/topofire-r-hist-clamped.png"
)

histogram_plot(
  plot_data, `Pearson's r`, network,
  model = "SOILWAT2",
  clamp = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/soilwat2-r-hist-clamped.png"
)

histogram_plot(
  plot_data, `KGE`, network,
  model = "Topofire",
  clamp = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/topofire-kge-hist-clamped.png"
)

histogram_plot(
  plot_data, `KGE`, network,
  model = "SOILWAT2",
  clamp = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/soilwat2-kge-hist-clamped.png"
)
