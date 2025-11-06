# =============================================================
# Title: Stability-Selected LASSO for Soil Moisture Model Comparison
# Author: Zachary H. Hoylman
# Date: 8/10/2025
#
# Description:
# This script implements LASSO regression with stability selection 
# using the `stabs` R package to identify robust predictors of 
# Topofire and SOILWAT2 model performance across soil moisture sites.
#
# Why stability-selected LASSO?
# - Traditional AIC-based model selection (e.g., stepwise regression)
#   does not handle collinearity well and lacks regularization.
# - Standard LASSO (via `glmnet`) improves on this by shrinking 
#   coefficients and selecting features, but its selection is 
#   unstable under collinearity — small changes in data can lead 
#   to different predictors being chosen.
# - Stability selection wraps LASSO with a resampling-based procedure 
#   that identifies predictors consistently selected across many 
#   subsamples of the data, offering greater interpretability and 
#   control over false positive rates.
#
# This approach is especially important here due to the high likelihood 
# of collinearity among climate and soil covariates, and the need for 
# scientifically interpretable predictor sets in manuscript figures 
# and tables.
# =============================================================

# ========== Load Libraries ==========
library(tidyverse)
library(glmnet)
library(gridExtra)
library(grid)
library(ragg)
library(stabs)

# ========== Load Data ==========
station_covariates = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv") #%>%
  #dplyr::select(network, site_id, ppt, soltotal, tdmean, tmean, vpdmax, sand, silt, clay, depth_to_restrictive, topofire_aws, porosity_gNATSGO)

station_meta = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

topofire_corelations = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv") %>%
  mutate(model = "Topofire") |>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            NSE = median(NSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup() |>
  rename(r = `Pearson's r`)

soilwat2_corelations = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv") %>%
  mutate(model = "SOILWAT2")|>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            NSE = median(NSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup() |>
  rename(r = `Pearson's r`)

# Combine
all_corelations = bind_rows(topofire_corelations, soilwat2_corelations) %>%
  filter(generalized_depth == 'Depth Averaged')

#diffrerence data
# 1. Pivot to wide format by model
all_corelations_difference = all_corelations %>%
  filter(generalized_depth == "Depth Averaged") %>%  # optional if you want just one depth
  pivot_wider(
    names_from = model,
    values_from = c(r, Bias, RMSE, NSE, KGE),
    names_sep = "_"
  ) %>%
  mutate(
    r_diff    = r_Topofire - r_SOILWAT2,
    bias_diff = Bias_Topofire - Bias_SOILWAT2,
    rmse_diff = RMSE_Topofire - RMSE_SOILWAT2,
    nse_diff  = NSE_Topofire - NSE_SOILWAT2,
    kge_diff  = KGE_Topofire - KGE_SOILWAT2
  ) %>%
  dplyr::select(network,site_id,generalized_depth, r_diff, bias_diff, rmse_diff, nse_diff, kge_diff)

# ========== Prepare Full Merged Dataset ==========
merged_all_models = all_corelations %>%
  left_join(station_covariates, by = c("network", "site_id"))

# For difference model
merged_all_models_differnece = all_corelations_difference %>%
  left_join(station_covariates, by = c("network", "site_id")) 

# ========== Variables ==========
predictors = setdiff(names(station_covariates), c("network", "site_id"))

target_sets = list(
  Topofire = c("Pearson's r", "KGE"),
  SOILWAT2 = c("Pearson's r", "KGE"),
  Difference = c("r_diff", "kge_diff")
)

# ========== LASSO + OLS Function ==========
fit_lasso_and_refit_ols = function(df, response) {
  df = df %>% select(all_of(c(response, predictors))) %>% drop_na()
  if (nrow(df) == 0) return(NULL)
  
  X = model.matrix(as.formula(paste(response, "~", paste(predictors, collapse = "+"))), data = df)[, -1]
  y = df[[response]]
  cv_fit = cv.glmnet(X, y, alpha = 1)
  best_lambda = cv_fit$lambda.min
  final_lasso = glmnet(X, y, alpha = 1, lambda = best_lambda)
  coefs = coef(final_lasso)
  nonzero = rownames(coefs)[coefs[, 1] != 0]
  final_feats = setdiff(nonzero, "(Intercept)")
  if (length(final_feats) == 0) return(lm(as.formula(paste(response, "~ 1")), data = df))
  lm(as.formula(paste(response, "~", paste(final_feats, collapse = "+"))), data = df)
}

fit_stability_selected_ols = function(df, response, cutoff = 0.75, PFER = 1) {
  df = df %>% select(all_of(c(response, predictors))) %>% drop_na()
  if (nrow(df) == 0) return(NULL)
  
  X = model.matrix(as.formula(paste(response, "~", paste(predictors, collapse = "+"))), data = df)[, -1]
  y = df[[response]]
  
  stabsel_fit = stabsel(
    x = X,
    y = y,
    fitfun = glmnet.lasso,
    cutoff = cutoff,
    PFER = PFER
  )
  
  selected_vars = colnames(X)[stabsel_fit$selected]
  
  if (length(selected_vars) == 0) {
    return(lm(as.formula(paste(response, "~ 1")), data = df))
  } else {
    return(lm(as.formula(paste(response, "~", paste(selected_vars, collapse = "+"))), data = df))
  }
}

# ========== Variable Name Map ==========
variable_name_map = c(
  ppt = "Precipitation (mm/month)",
  tmean = "Mean Temperature (°C)",
  tmin = "Min Temperature (°C)",
  tmax = "Max Temperature (°C)",
  tdmean = "Dew Point Temperature (°C)",
  vpdmin = "Min Vapor Pressure Deficit (hPa)",
  vpdmax = "Max Vapor Pressure Deficit (hPa)",
  solclear = "Clear Sky Solar Radiation (MJ/m²/day)",
  solslope = "Sloped Surface Solar Radiation (MJ/m²/day)",
  soltotal = "Total Solar Radiation (MJ/m²/day)",
  soltrans = "Atmospheric Transmittance (Fraction)",
  elevation = "Elevation (m)",
  sand = "Sand (%)",
  silt = "Silt (%)",
  clay = "Clay (%)",
  som = "Soil Organic Matter (%)",
  bulk_density = "Bulk Density (g/cm³)",
  awc = "Available Water Holding Capacity (cm)",
  depth_to_restrictive = "Depth to Restrictive Layer (cm)",
  topofire_aws = "gNATSGO AWS (mm)",
  porosity_gNATSGO = "gNATSGO Porosity (mm/m)"
)

# ========== Table Export Function ==========
format_regression_table_png <- function(model_summary, soil_depth, response_variable, export_path, prefix = NULL) {
  if (is.null(model_summary) || !inherits(model_summary, "summary.lm")) return()
  
  # Label map
  response_label <- dplyr::case_when(
    response_variable == "r" ~ "Pearson's r",
    response_variable == "KGE" ~ "Kling-Gupta Efficiency (KGE)",
    response_variable == "r_diff" ~ "Pearson's r Difference (Topofire - SOILWAT2)",
    response_variable == "kge_diff" ~ "KGE Difference (Topofire - SOILWAT2)",
    TRUE ~ response_variable
  )
  
  title_line1 <- "Stability-Selected LASSO Regression Results"
  title_line2 <- paste(prefix, soil_depth, response_label, sep = " ")
  
  # Extract and format coefficients
  coefs <- model_summary$coefficients
  r2 <- round(model_summary$r.squared, 3)
  adj_r2 <- round(model_summary$adj.r.squared, 3)
  
  reg_table <- as_tibble(coefs, rownames = "Predictor") %>%
    dplyr::mutate(
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `Estimate ± SE` = paste0(Estimate, " ± ", `Std. Error`),
      Predictor = dplyr::recode(Predictor, !!!variable_name_map)
    ) %>%
    dplyr::select(Predictor, `Estimate ± SE`)
  
  # Adjust height dynamically
  n_rows <- nrow(reg_table)
  row_height <- 35
  base_height <- 200  # Increased to allow for 2-line title
  img_height <- base_height + (n_rows * row_height)
  
  # Table grob
  gt_grob <- gridExtra::tableGrob(reg_table, rows = NULL)
  gt_grob$theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 14), bg_params = list(fill = c("white", "gray95"))),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  # Plot
  ragg::agg_png(export_path, width = 1000, height = img_height, res = 150)
  grid.newpage()
  
  # Use a more compact layout
  layout <- grid.layout(
    nrow = 3,
    heights = unit(c(2.5, 1, 1), c("lines", "null", "lines"))  # reduce top and bottom padding
  )
  pushViewport(viewport(layout = layout))
  
  # Title lines
  pushViewport(viewport(layout.pos.row = 1))
  grid.text(title_line1, y = unit(0.7, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text(title_line2, y = unit(0.3, "npc"), gp = gpar(fontsize = 12))
  popViewport()
  
  # Table
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(gt_grob)
  popViewport()
  
  # Footer (R²)
  pushViewport(viewport(layout.pos.row = 3))
  grid.text(
    paste("R² =", r2, "| Adjusted R² =", adj_r2),
    y = unit(0.5, "npc"),
    gp = gpar(fontsize = 11, fontface = "italic")
  )
  popViewport()
  
  dev.off()
  message("Saved: ", export_path)
}

# ========== MODEL FITTING AND EXPORT LOOP ==========
fit_and_export_models = function(data, response_vars, prefix) {
  map(response_vars, function(resp) {
    #old method, original LASSO
    #model = fit_lasso_and_refit_ols(data, resp)
    #new method with stabs (stability selection)
    model = fit_stability_selected_ols(data, resp)
    if (!is.null(model)) {
      summary_obj = summary(model)
      out_path = glue::glue("~/nrcs-scd-soil-moisture-eval/tables/lasso-{prefix}-{tolower(resp)}.png")
      format_regression_table_png(summary_obj, "Depth Averaged", resp, export_path = out_path, prefix = prefix)
    }
  })
}

fit_and_export_models(merged_all_models %>% filter(model == 'Topofire'), c("r", "KGE"), "Topofire")
fit_and_export_models(merged_all_models %>% filter(model == 'SOILWAT2'), c("r", "KGE"), "SOILWAT2")

# quartile box plots
plot_box_by_quartile <- function(data, predictor, target, n_bins = 4, 
                                 predictor_label = NULL, target_label = NULL, 
                                 title = NULL, fill_color = "#2c7fb8",
                                 export_path = NULL,
                                 width = 8, height = 6, dpi = 300) {
  
  predictor_sym <- rlang::ensym(predictor)
  target_sym <- rlang::ensym(target)
  
  # Drop NA and bin predictor
  plot_data <- data %>%
    drop_na(!!predictor_sym, !!target_sym) %>%
    mutate(bin = ntile(!!predictor_sym, n_bins))
  
  # Compute mean of predictor for each bin
  bin_means <- plot_data %>%
    group_by(bin) %>%
    summarise(
      mean_predictor = round(mean(!!predictor_sym, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(bin_label = paste0("Mean: ", mean_predictor))
  
  # Join bin labels back to plot data
  plot_data <- plot_data %>%
    left_join(bin_means, by = "bin")
  
  # Auto-generate axis labels if not provided
  if (is.null(predictor_label)) predictor_label <- rlang::as_label(predictor_sym)
  if (is.null(target_label)) target_label <- rlang::as_label(target_sym)
  if (is.null(title)) title <- paste(target_label, "by", predictor_label)
  
  # Recalculate range using IQR-based limits (no outliers)
  y_vals <- plot_data[[rlang::as_label(target_sym)]]
  q1 <- quantile(y_vals, 0.25, na.rm = TRUE)
  q3 <- quantile(y_vals, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  y_max <- max(abs(q3 + 1.5 * iqr), abs(q1 - 1.5 * iqr))
  y_lim <- c(-y_max, y_max)
  zero_offset <- 0.03 * (2 * y_max)
  
  # Plot
  p <- ggplot(plot_data, aes(x = factor(bin_label, levels = bin_means$bin_label), y = !!target_sym)) +
    geom_boxplot(fill = fill_color, color = "black", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    annotate("text", x = 0.5, y = 0 + zero_offset, label = "Topofire is better →", angle = 90, hjust = 0, size = 4.5) +
    annotate("text", x = 0.5, y = 0 - zero_offset, label = "← SOILWAT2 is better", angle = 90, hjust = 1, size = 4.5) +
    scale_y_continuous(limits = y_lim) +
    labs(
      title = title,
      x = predictor_label,
      y = target_label
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.x = element_blank()
    )
  
  # Save if export path is specified
  if (!is.null(export_path)) {
    ggsave(export_path, plot = p, width = width, height = height, dpi = dpi)
    message("Plot saved to: ", export_path)
  }
  
  return(p)
}

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = ppt,
  target = kge_diff,
  n_bins = 4,
  predictor_label = "Precipitation (mm/month)",
  target_label = "ΔKGE (Topofire - SOILWAT2)",
  fill_color = "blue",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/kge_diff_by_ppt.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = tmean,
  target = kge_diff,
  n_bins = 4,
  predictor_label = "Mean Temperature (°C)",
  target_label = "ΔKGE (Topofire - SOILWAT2)",
  fill_color = "red",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/kge_diff_by_tmean.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = ppt,
  target = r_diff,
  n_bins = 4,
  predictor_label = "Precipitation (mm/month)",
  target_label = "ΔPearson's R (Topofire - SOILWAT2)",
  fill_color = "blue",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/r_diff_by_ppt.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = tmean,
  target = r_diff,
  n_bins = 4,
  predictor_label = "Mean Temperature (°C)",
  target_label = "ΔPearson's R (Topofire - SOILWAT2)",
  fill_color = "red",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/r_diff_by_tmean.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = topofire_aws,
  target = kge_diff,
  n_bins = 4,
  predictor_label = "SSURGO AWS (mm)",
  target_label = "ΔKGE (Topofire - SOILWAT2)",
  fill_color = "darkgrey",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/kge_diff_by_ssurgo_aws.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = topofire_aws,
  target = r_diff,
  n_bins = 4,
  predictor_label = "SSURGO AWS (mm)",
  target_label = "ΔPearson's R (Topofire - SOILWAT2)",
  fill_color = "darkgrey",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/r_diff_by_ssurgo_aws.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = sand,
  target = kge_diff,
  n_bins = 4,
  predictor_label = "Sand (%)",
  target_label = "ΔKGE (Topofire - SOILWAT2)",
  fill_color = "#836539",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/kge_diff_by_ssurgo_sand.png"
)

plot_box_by_quartile(
  data = merged_all_models_differnece,
  predictor = sand,
  target = r_diff,
  n_bins = 4,
  predictor_label = "Sand (%)",
  target_label = "ΔPearson's R (Topofire - SOILWAT2)",
  fill_color = "#836539",
  export_path = "~/nrcs-scd-soil-moisture-eval/figs/r_diff_by_ssurgo_sand.png"
)