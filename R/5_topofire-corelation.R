library(tidyverse)
library(magrittr)
library(sf)
library(khroma)

#porosity from https://data.usgs.gov/datacatalog/data/USGS:5fd7c19cd34e30b9123cb51f

`%notin%` = Negate(`%in%`)

staion_data = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv") %>%
  pivot_wider(names_from = generalized_depth, values_from = soil_moisture)

valid_time = staion_data %$%
  date %>%
  unique()

topofire_data = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/topofire/raw", full.names = T) %>%
  purrr::map(read_csv, show_col_types = F) %>%
  bind_rows() 

topofire_dates = colnames(topofire_data)[2:length(colnames(topofire_data))] %>%
  as.Date()

common_dates = topofire_dates[topofire_dates %in% valid_time]

#memory management
gc();gc()

topofire_data_long = topofire_data %>% 
  dplyr::select(site_id, common_dates %>% as.character()) %>%
  pivot_longer(cols = -c(site_id)) %>%
  dplyr::select(site_id, date = name, Topofire = value) %>%
  drop_na(Topofire) %>%
  mutate(date = date %>% as.Date())

rm(topofire_data); gc(); gc()

#compute VWC using available water holding capacity and porosity

params = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv") %>%
  dplyr::select(site_id, topofire_aws, topofire_porosity)

topofire_data_final = topofire_data_long %>% 
  left_join(params, 'site_id')  %>%
  mutate(Topofire_VWC = (Topofire / topofire_aws) * (topofire_porosity / 1000)) %>%
  dplyr::select(site_id, date, Topofire = Topofire_VWC)

rm(topofire_data_long); gc(); gc()

data_joined = staion_data %>% 
  left_join(., topofire_data_final, by = c('site_id', 'date'))

compute_bias = function(actual, predicted) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors must be of the same length.")
  }
  return(mean(predicted - actual))
}

compute_rmse = function(actual, predicted) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors must be of the same length.")
  }
  return(sqrt(mean((predicted - actual)^2)))
}

compute_kge = function(actual, predicted) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors must be of the same length.")
  }
  
  actual_mean = mean(actual)
  predicted_mean = mean(predicted)
  
  r = cor(actual, predicted)
  alpha = sd(predicted) / sd(actual)
  beta = predicted_mean / actual_mean
  
  kge = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  return(kge)
}

compute_nse = function(actual, predicted) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors must be of the same length.")
  }
  
  numerator = sum((actual - predicted)^2)
  denominator = sum((actual - mean(actual))^2)
  
  nse = 1 - (numerator / denominator)
  return(nse)
}

#run correlation by site and depth
site_depth_corelations = data_joined %>%
  pivot_longer(cols = -c(network, site_id, date, Topofire), names_to = 'generalized_depth', values_to = 'soil_moisture') %>%
  drop_na(Topofire, soil_moisture) %>%
  group_by(network, site_id, generalized_depth) %>%
  summarise(
    `Pearson's r` = cor(soil_moisture, Topofire),
    `Bias` = compute_bias(soil_moisture, Topofire),
    `RMSE` = compute_rmse(soil_moisture, Topofire),
    `NSE` = compute_nse(soil_moisture, Topofire),
    `KGE` = compute_kge(soil_moisture, Topofire)
  ) %>%
  ungroup() %>%
  mutate(
    NSE = ifelse(NSE<0, 0, NSE),
    KGE = ifelse(KGE<0, 0, KGE)
  )

write_csv(site_depth_corelations, "~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv")

scatter_model_vs_obs = function(data, sims, obs, network, depth, 
                                model = NULL,
                                clamp = FALSE, 
                                export_path = NULL, 
                                height = 10, 
                                width = 12,
                                depth_avg_only = FALSE) {
  
  sims <- rlang::ensym(sims)
  obs <- rlang::ensym(obs)
  network <- rlang::ensym(network)
  depth <- rlang::ensym(depth)
  
  # Optional filter by model
  if (!is.null(model)) {
    data <- data %>% filter(model == !!model)
  }
  
  # Optional filter to only depth-averaged data
  if (depth_avg_only) {
    data <- data %>% filter(!!depth == "Depth Averaged")
  }
  
  # Apply clamping
  if (clamp) {
    data <- data %>%
      mutate(!!obs := pmax(!!obs, 0), !!sims := pmax(!!sims, 0))
  }
  
  # Set depth order
  data <- data %>%
    mutate(!!depth := factor(!!depth, levels = c("Depth Averaged", "Shallow", "Middle", "Deep")))
  
  color_scale <- khroma::color("roma")
  
  # Dynamic title and subtitle
  plot_title <- if (!is.null(model)) paste(model, "vs Observations") else "Model vs Observations"
  subtitle_text <- if (depth_avg_only) "Depth Averaged Only Across Networks" else "(Across Depths & Networks)"
  
  scatter_plot <- ggplot(data, aes(y = !!obs, x = !!sims)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(ndensity)),
      contour = FALSE,
      bins = 50
    ) +
    geom_smooth(method = lm) +
    facet_grid(rows = vars(!!network), cols = vars(!!depth)) +
    scale_fill_gradientn(
      colours = color_scale(100),
      name = "Normalized Density",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
      limits = c(0, 0.6),
      breaks = c(0, 0.2, 0.4, 0.6),
      labels = c('0', '0.2', '0.4', '> 0.6'),
      na.value = color_scale(100)[100]
    ) +
    theme_bw(base_size = 16) +
    xlim(c(0, 0.6)) +
    ylim(c(0, 0.6)) +
    theme(
      legend.position = 'bottom',
      legend.key.width = unit(4, "cm"),
      legend.title = element_text(size = 14, hjust = 0.5),
      strip.background = element_rect(fill = "transparent", color = "transparent"),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    ggtitle(plot_title, subtitle = subtitle_text) +
    labs(y = 'Observed Soil Moisture (m³/m³)', x = 'Predicted Soil Moisture (m³/m³)')
  
  if (!is.null(export_path)) {
    ggsave(export_path, plot = scatter_plot, height = height, width = width)
    message("Plot saved to: ", export_path)
  }
  
  return(scatter_plot)
}

scatter_model_vs_obs(
  data = data_joined %>%
    pivot_longer(cols = -c(network, site_id, date, Topofire), names_to = 'generalized_depth', values_to = 'soil_moisture') %>%
    drop_na(Topofire, soil_moisture),
  sims = Topofire,
  obs = soil_moisture,
  network = network,
  depth = generalized_depth,
  model = "Topofire",
  clamp = F,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/topofire-scatter.png"
)
