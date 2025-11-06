library(tidyverse)
library(magrittr)
library(sf)
library(khroma)

`%notin%` = Negate(`%in%`)

#import station data
staion_data = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv") %>%
  pivot_wider(names_from = generalized_depth, values_from = soil_moisture) %>%
  mutate(site_id_long = site_id,
    site_id = case_when(
      network == "NEON" ~ str_remove(site_id, "_.*"),
      TRUE              ~ site_id
    )
  ) |> 
  select(-site_id_long) |>
  pivot_longer(cols = -c(network, site_id, date), values_to = 'obs', names_to = 'generalized_depth')

#import raw data
SPoRT_files = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/SPoRT-LIS", full.names = T) |>
  lapply(read_csv) |>
  bind_rows() 

#clean, pivot and compute site_id and network
SPoRT_data = SPoRT_files |>
  rename(date = time) |>
  mutate(date = as.Date(date)) |>
  pivot_longer(cols = -c('var', 'date')) |>
  separate_wider_delim(
    name,
    delim = "_",
    names = c("site_id", "network"),
    too_few = "align_start"
  ) |>
  mutate(generalized_depth = ifelse(var == 'SPoRT_raw_0-10cm', 'Shallow',
                        ifelse(var == 'SPoRT_raw_10-40cm', 'Middle','Deep'))) |>
  select(network, site_id, date, generalized_depth, SPoRT = value) 
  
#compute depth averaged data and merge all together
SPoRT_data_final = SPoRT_data |>
  group_by(network, site_id, date) |>
  summarize(
    generalized_depth = "Depth Averaged",
    SPoRT = mean(SPoRT, na.rm = TRUE),
    .groups = "drop"
  ) |>
  bind_rows(SPoRT_data) |>
  arrange(network, site_id, date, generalized_depth)

#join with the obs
data_joined = staion_data %>% 
  left_join(., SPoRT_data_final, by = c('network','site_id', 'date', 'generalized_depth'))

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

#run correlation by site and depth
site_depth_corelations = data_joined %>%
  drop_na(SPoRT, obs) %>%
  group_by(network, site_id, generalized_depth) %>%
  summarise(
    `Pearson's r` = cor(obs, SPoRT),
    `Bias` = compute_bias(obs, SPoRT),
    `RMSE` = compute_rmse(obs, SPoRT),
    `KGE` = compute_kge(obs, SPoRT)
  ) %>%
  ungroup() %>%
  mutate(
    KGE = ifelse(KGE<0, 0, KGE)
  )

write_csv(site_depth_corelations, "~/nrcs-scd-soil-moisture-eval-data/processed/SPoRT-LIS-corelations-generalized-depth.csv")

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
  sims = SPoRT,
  obs = obs,
  network = network,
  depth = generalized_depth,
  model = "SPoRT-LIS",
  clamp = F,
  depth_avg_only = F,
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/SPoRT-LIS-scatter.png"
)


