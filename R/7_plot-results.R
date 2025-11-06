library(tidyverse)
library(sf)

sport = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/SPoRT-LIS-corelations-generalized-depth.csv") |>
  mutate(model = 'SPoRT-LIS') |>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup()

soilwat2 = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv") |>
  mutate(model = 'SOILWAT2')|>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup()

all_results = bind_rows(sport, soilwat2)

# plot rneonSoilFlux# plot results
states = read_sf("~/mco-drought-indicators/processing/base-data/raw/states.shp") %>%
  filter(!STATE_NAME %in% c('Hawaii', 'Alaska', 'Puerto Rico', 'Virgin Islands'))

station_data = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

plot_data = all_results %>%
  left_join(., station_data) %>%
  drop_na(latitude, longitude) %>%
  st_as_sf(., coords = c('longitude', 'latitude'), crs = 4326) %>%
  mutate(name_new = ifelse(generalized_depth == 'Deep', 'Deep (>40cm)', 
                           ifelse(generalized_depth == 'Shallow', 'Shallow (<= 10cm)', 
                                  ifelse(generalized_depth == 'Middle', 'Middle (>10cm - <=40cm)', 'Depth Averaged'))),
         name_new = factor(name_new, levels = c('Depth Averaged', 'Shallow', 
                                                'Shallow (<= 10cm)', 'Middle (>10cm - <=40cm)', 
                                                'Deep (>40cm)')))


make_plot = function(data, variable, states, 
                     model = NULL,
                     clamp = FALSE, 
                     export_path = NULL, 
                     height = 7, 
                     width = 8,
                     depth_avg_only = FALSE) {
  
  variable <- rlang::ensym(variable)
  
  # Optional model filter
  if (!is.null(model)) data <- dplyr::filter(data, model == !!model)
  
  # Optional depth filter
  if (depth_avg_only) data <- dplyr::filter(data, name_new == "Depth Averaged")
  
  # --- Normalize any legacy labels (safety) and set desired facet order ---
  desired_levels <- c(
    "Shallow (<= 10cm)",
    "Middle (>10cm - <=40cm)",
    "Deep (>40cm)",
    "Depth Averaged"
  )
  
  data <- data %>%
    dplyr::mutate(
      # If any plain labels got through, upgrade them to the new detailed ones
      name_new = dplyr::case_when(
        name_new == "Shallow" ~ "Shallow (<= 10cm)",
        name_new == "Middle"  ~ "Middle (>10cm - <=40cm)",
        name_new == "Deep"    ~ "Deep (>40cm)",
        TRUE                  ~ as.character(name_new)
      ),
      # Reorder facets: Shallow, Middle, Deep, Depth Averaged
      name_new = forcats::fct_relevel(name_new, desired_levels)
    )
  
  # Apply clamping
  if (clamp) data <- dplyr::mutate(data, !!variable := pmax(!!variable, 0))
  
  # Reproject to Albers Equal Area
  data   <- sf::st_transform(data, 5070)
  states <- sf::st_transform(states, 5070)
  
  # Legend title
  var_label <- rlang::as_name(variable)
  legend_title <- dplyr::case_when(
    var_label == "KGE"         ~ "Kling-Gupta Efficiency (KGE)",
    var_label == "NSE"         ~ "Nash-Sutcliffe Efficiency (NSE)",
    var_label == "RMSE"        ~ "Root Mean Square Error (RMSE)",
    var_label == "Bias"        ~ "Bias",
    var_label == "Pearson's r" ~ "Pearson's r",
    TRUE                       ~ var_label
  )
  
  plot_title   <- if (!is.null(model)) paste(model, legend_title) else paste("Model", legend_title)
  subtitle_txt <- if (depth_avg_only) "Depth Averaged Only" else ""
  
  label_func <- function(x) {
    formatted <- format(round(x, 2), nsmall = 2)
    if (clamp && length(x) > 0) formatted[which.min(x)] <- paste0("<", formatted[which.min(x)])
    formatted
  }
  
  cor_plot <- ggplot2::ggplot(data, ggplot2::aes(fill = !!variable)) +
    ggplot2::geom_sf(data = states, fill = 'transparent') +
    ggplot2::geom_sf(shape = 21, color = 'black') +
    ggplot2::scale_fill_gradientn(
      colours = rev(viridis::turbo(5)),
      labels  = label_func
    ) +
    ggplot2::facet_wrap(~ name_new) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title = legend_title, title.position = "bottom", title.hjust = 0.5)) +
    ggplot2::theme(
      legend.key       = ggplot2::element_blank(),
      legend.position  = 'bottom',
      legend.key.width = grid::unit(2, "cm"),
      strip.background = ggplot2::element_rect(colour = "transparent", fill = "transparent"),
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 18),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 12)
    ) +
    ggplot2::ggtitle(plot_title, subtitle = subtitle_txt)
  
  if (!is.null(export_path)) {
    ggplot2::ggsave(export_path, plot = cor_plot, height = height, width = width)
    message("Plot saved to: ", export_path)
  }
  cor_plot
}

make_plot(
  data = plot_data,
  variable = `Pearson's r`,
  states = states,
  model = "SPoRT-LIS",
  clamp = TRUE,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/sport-pearson-r-clamped.png'
)

make_plot(
  data = plot_data,
  variable = `Pearson's r`,
  states = states,
  model = "SOILWAT2",
  clamp = TRUE,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-pearson-r-clamped.png'
)

make_plot(
  data = plot_data,
  variable = `Bias`,
  states = states,
  model = "SPoRT-LIS",
  clamp = F,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/sport-bias.png'
)

make_plot(
  data = plot_data,
  variable = `Bias`,
  states = states,
  model = "SOILWAT2",
  clamp = F,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-bias.png'
)

make_plot(
  data = plot_data,
  variable = `KGE`,
  states = states,
  model = "SPoRT-LIS",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/sport-kge.png'
)

make_plot(
  data = plot_data,
  variable = `KGE`,
  states = states,
  model = "SOILWAT2",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/soilwat2-kge.png'
)

make_plot(
  data = plot_data,
  variable = `RMSE`,
  states = states,
  model = "SPoRT-LIS",
  clamp = T,
  depth_avg_only = FALSE,
  export_path = '~/nrcs-scd-soil-moisture-eval/figs/sport-rmse.png'
)

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
    select(network, site_id, model, name_new, value = !!variable)
  
  # Pivot to wide and compute differences
  data_diff <- summarized |>
    pivot_wider(
      names_from = model,
      values_from = value,
      id_cols = c(network, site_id, name_new)
    )
  
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
  
  t.test(data_diff$`SPoRT-LIS`, data_diff$SOILWAT2, paired = TRUE) 
  wilcox.test(data_diff$`SPoRT-LIS`, data_diff$SOILWAT2, paired = TRUE)
  
  # Reattach geometry and depth labels
  site_meta <- data_filtered %>%
    select(site_id, name_new, geometry) %>%
    distinct()
  
  data_diff <- data_diff %>%
    left_join(site_meta, by = c("site_id", "name_new")) %>%
    st_as_sf()
  
  desired_levels <- c(
    "Shallow (<= 10cm)",
    "Middle (>10cm - <=40cm)",
    "Deep (>40cm)",
    "Depth Averaged"
  )
  
  data_diff <- data_diff %>%
    dplyr::mutate(
      name_new = as.character(name_new),
      name_new = dplyr::case_when(
        name_new == "Shallow" ~ "Shallow (<= 10cm)",
        name_new == "Middle"  ~ "Middle (>10cm - <=40cm)",
        name_new == "Deep"    ~ "Deep (>40cm)",
        TRUE                  ~ name_new
      ),
      name_new = forcats::fct_relevel(name_new, desired_levels)
    )
  
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
  subtitle_text <- if (depth_avg_only) "Depth Averaged Only" else ""
  
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
  model_2 = "SPoRT-LIS",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-sport-vs-soilwat2-kge.png"
)

model_diff_plot(
  data = plot_data,
  variable = "Pearson's r",
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "SPoRT-LIS",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-sport-vs-soilwat2-r.png"
)

model_diff_plot(
  data = plot_data,
  variable = 'Bias',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "SPoRT-LIS",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-sport-vs-soilwat2-bias.png"
)

model_diff_plot(
  data = plot_data,
  variable = 'RMSE',
  states = states,
  model_1 = "SOILWAT2",
  model_2 = "SPoRT-LIS",
  clamp_to_percentile = TRUE,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/diff-sport-vs-soilwat2-RMSE.png"
)
