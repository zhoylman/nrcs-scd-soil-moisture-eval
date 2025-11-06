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
    
    # 40-deepest cms
    deep_cols = names(df)[stringr::str_detect(names(df), "Sim_VWC_0[4-9][0-9]?to\\d+_cm")]
    
    df = df %>%
      dplyr::mutate(
        network = network,
        site_id = site_id,
        Shallow = rowMeans(dplyr::across(c("Sim_VWC_000to005_cm", "Sim_VWC_005to010_cm")), na.rm = TRUE),
        Middle  = rowMeans(dplyr::across(c("Sim_VWC_010to020_cm", "Sim_VWC_020to030_cm", 
                                           "Sim_VWC_030to040_cm")), na.rm = TRUE),
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


station_meta = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

obs = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv") |>
  mutate(site_id_long = site_id,
         date = as.Date(date),
         site_id = case_when(
           network == "NEON" ~ str_remove(site_id, "_.*"),
           TRUE              ~ site_id
         )
  )


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
  group_by(network, site_id, site_id_long, generalized_depth) %>%
  summarise(
    `Pearson's r` = cor(soil_moisture, SOILWAT2),
    `Bias` = compute_bias(soil_moisture, SOILWAT2),
    `RMSE` = compute_rmse(soil_moisture, SOILWAT2),
    `KGE` = compute_kge(soil_moisture, SOILWAT2)
  ) %>%
  ungroup() %>%
  mutate(
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

scatter_model_vs_obs <- function(data, sims, obs, network, depth, 
                                 model = NULL,
                                 clamp = FALSE, 
                                 export_path = NULL, 
                                 height = 14, 
                                 width = 14,
                                 depth_avg_only = FALSE) {
  
  sims    <- rlang::ensym(sims)
  obs     <- rlang::ensym(obs)
  network <- rlang::ensym(network)
  depth   <- rlang::ensym(depth)
  
  if (!is.null(model))  data <- dplyr::filter(data, model == !!model)
  if (depth_avg_only)   data <- dplyr::filter(data, !!depth == "Depth Averaged")
  if (clamp) {
    data <- dplyr::mutate(data,
                          !!obs  := pmax(!!obs, 0),
                          !!sims := pmax(!!sims, 0)
    )
  }
  
  data <- dplyr::mutate(
    data,
    !!depth := factor(!!depth, levels = c("Shallow","Middle","Deep","Depth Averaged"))
  )
  net_vals   <- unique(dplyr::pull(data, !!network))
  priority   <- c("MT Mesonet","OK Mesonet")
  net_levels <- c(priority[priority %in% net_vals], setdiff(net_vals, priority))
  data <- dplyr::mutate(data, !!network := factor(!!network, levels = net_levels))
  
  plot_df <- dplyr::filter(
    data,
    dplyr::between(!!sims, 0, 0.6),
    dplyr::between(!!obs,  0, 0.6)
  )
  
  color_scale <- khroma::color("roma")
  plot_title  <- if (!is.null(model)) paste(model, "vs Observations") else "Model vs Observations"
  subtitle_txt <- if (depth_avg_only) "Depth Averaged Only Across Networks" else "(Across Depths & Networks)"
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = !!sims, y = !!obs)) +
    ggplot2::stat_density_2d(
      data = plot_df,  # ensure density also matches the filtered domain
      geom = "raster",
      ggplot2::aes(fill = ggplot2::after_stat(ndensity)),
      contour = FALSE,
      bins = 50,
      na.rm = TRUE
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "black", linewidth = 0.8) +
    ggplot2::geom_smooth(
      data = plot_df,
      mapping = ggplot2::aes(group = interaction(!!network, !!depth)),
      method = stats::lm, formula = y ~ x, se = TRUE, color = "black", na.rm = TRUE
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(!!network),
                        cols  = ggplot2::vars(!!depth)) +
    ggplot2::scale_fill_gradientn(
      colours = color_scale(100),
      name = "Normalized Density",
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5),
      limits = c(0, 0.6), breaks = c(0, 0.2, 0.4, 0.6),
      labels = c("0", "0.2", "0.4", "> 0.6"),
      na.value = color_scale(100)[100]
    ) +
    ggplot2::labs(
      x = "Predicted Soil Moisture (m³/m³)",
      y = "Observed Soil Moisture (m³/m³)"
    ) +
    ggplot2::ggtitle(plot_title, subtitle = subtitle_txt) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::coord_equal(xlim = c(0, 0.6), ylim = c(0, 0.6), expand = FALSE) +
    scale_x_continuous(
      limits = c(0, 0.6),
      breaks = c(0, 0.3, 0.6),
      labels = c("0", "0.3", "0.6"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 0.6),
      breaks = c(0, 0.3, 0.6),
      labels = c("0", "0.3", "0.6"),
      expand = c(0, 0)
    ) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.key.width = grid::unit(4, "cm"),
      legend.title     = ggplot2::element_text(size = 14, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "transparent", color = "transparent"),
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 18),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 12),
      panel.spacing.x  = grid::unit(1.2, "lines"),
      panel.spacing.y  = grid::unit(1.2, "lines"),
      axis.title.x     = element_text(margin = margin(t = 10)),
      axis.title.y     = element_text(margin = margin(r = 10))
    )
  
  if (!is.null(export_path)) ggplot2::ggsave(export_path, p, height = height, width = width)
  p
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
