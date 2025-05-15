library(tidyverse)
library(sf)

read_soilwat = function(x){
  tryCatch({
    name = stringr::str_extract(x, "(?<=__)[^/]+(?=_sc1\\.rds)") %>%
      stringr::str_split(., "-", simplify = TRUE) 
    
    if (!name[1] %in% c('NEON', 'OKMesonet', 'MTMesonet')) {
      network = name[1]
      site_id = name[3]
    } else {
      network = name[1]
      site_id = name[2]
    }
    
    network = ifelse(network == 'SNOTEL', 'SNTL',
                     ifelse(network == 'MTMesonet', 'MT Mesonet',
                            ifelse(network == 'OKMesonet', 'OK Mesonet', network)))
    
    site_id = ifelse(network %in% c('SNTL', 'SCAN'), paste0(network, ':', site_id), site_id)
    
    df = readRDS(x) %>% tibble::as_tibble()
    
    deep_cols = names(df)[stringr::str_detect(names(df), "Sim_VWC_0[5-9][0-9]?to\\d+_cm")]
    
    df = df %>%
      dplyr::mutate(
        network = network,
        site_id = site_id,
        Shallow = rowMeans(dplyr::across(c("Sim_VWC_000to005_cm", "Sim_VWC_005to010_cm")), na.rm = TRUE),
        Middle  = rowMeans(dplyr::across(c("Sim_VWC_010to020_cm", "Sim_VWC_020to030_cm", 
                                           "Sim_VWC_030to040_cm", "Sim_VWC_040to050_cm")), na.rm = TRUE),
        Deep    = rowMeans(dplyr::across(dplyr::all_of(deep_cols)), na.rm = TRUE),
        `Depth Averaged` = rowMeans(dplyr::across(c('Shallow', 'Middle', 'Deep'))),
        date    = lubridate::make_date(Year, Month, Day)
      ) %>%
      dplyr::select(network, site_id, date, `Depth Averaged`, Shallow, Middle, Deep)
    
    return(df)
  }, error = function(e) {
    message("Skipping file due to error: ", x)
    return(NULL)
  })
}

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

station_meta = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

obs = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv")

files = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/20250409_SOILWAT2_OutputVWCShared__20250320_SoilClimateDynamics_SOILWAT2_simulations", full.names = T)

soilwat_raw = files %>%
  purrr::map(read_soilwat) %>%
  bind_rows()

soilwat = soilwat_raw %>%
  pivot_longer(cols = -c(network, site_id, date), names_to = 'generalized_depth', values_to = 'SOILWAT2')

data_joined = obs %>%
  left_join(soilwat, by = c('network', 'site_id', 'date', 'generalized_depth')) %>%
  drop_na()

#run correlation by site and depth
site_depth_corelations = data_joined %>%
  group_by(network, site_id, generalized_depth) %>%
  summarise(
    `Pearson's r` = cor(soil_moisture, SOILWAT2),
    `Bias` = compute_bias(soil_moisture, SOILWAT2),
    `RMSE` = compute_rmse(soil_moisture, SOILWAT2),
    `NSE` = compute_nse(soil_moisture, SOILWAT2),
    `KGE` = compute_kge(soil_moisture, SOILWAT2)
  ) %>%
  ungroup() %>%
  mutate(
    NSE = ifelse(NSE<0, 0, NSE),
    KGE = ifelse(KGE<0, 0, KGE)
  )

plot_data = site_depth_corelations %>%
  left_join(., station_meta) %>%
  drop_na(c('latitude', 'longitude')) %>%
  st_as_sf(., coords = c('longitude', 'latitude'), crs = 4326) %>%
  mutate(name_new = ifelse(generalized_depth == 'Deep', 'Deep (>50cm)', 
                           ifelse(generalized_depth == 'Shallow', 'Shallow (<= 10cm)', 
                                  ifelse(generalized_depth == 'Middle', 'Middle (>10cm - <=50cm)', 'Depth Averaged'))),
         name_new = factor(name_new, levels = c('Depth Averaged', 'Shallow', 
                                                'Shallow (<= 10cm)', 'Middle (>10cm - <=50cm)', 
                                                'Deep (>50cm)')))

write_csv(site_depth_corelations, "~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv")

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
  data = data_joined,
  sims = SOILWAT2,
  obs = soil_moisture,
  network = network,
  depth = generalized_depth,
  model = "SOILWAT2",
  clamp = F,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/soilwat2-scatter.png"
)
