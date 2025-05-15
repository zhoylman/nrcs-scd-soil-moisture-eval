library(tidyverse)
library(glmnet)       # For LASSO
library(gridExtra)
library(grid)
library(ragg)

# ============== LOAD AND PREP DATA (same as before) =========================
station_meta <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")
station_covariates <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv")
topofire_corelations <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv")

merged_data <- topofire_corelations %>%
  dplyr::select(-c("Pearson's r (Normalized)")) %>%
  left_join(., station_covariates, by = c('network', 'site_id')) %>%
  dplyr::select(-c('network', 'site_id')) %>%
  rename('r' = `Pearson's r`, 'bias' = `Bias (Normalized)`, 'rmse' = `RMSE (Normalized)`)

target_vars <- c("r", "bias", "rmse")
predictors <- setdiff(names(merged_data), c("name", target_vars))

# ============== LASSO FIT FUNCTION (returns final OLS model) ================
# 1) Fit LASSO with cross-validation
# 2) Extract nonzero features (excluding intercept)
# 3) Refit OLS on that subset => returns an lm model for easy summary/p-values
fit_lasso_and_refit_ols <- function(df, response) {
  # Subset relevant columns and drop NA
  df <- df %>%
    dplyr::select(all_of(c(response, predictors))) %>%
    drop_na()
  
  if (nrow(df) == 0) {
    warning(paste("Skipping model for", response, "- No valid data after removing NAs"))
    return(NULL)
  }
  
  # Create model matrix for LASSO
  X <- model.matrix(as.formula(paste(response, "~", paste(predictors, collapse = " + "))), data = df)[, -1]
  y <- df[[response]]
  
  # LASSO with cross-validation
  cv_fit <- cv.glmnet(X, y, alpha = 1)
  
  # Best lambda
  best_lambda <- cv_fit$lambda.min
  
  # Full LASSO model at best lambda
  final_lasso <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  
  # Extract nonzero coefficients (including intercept)
  coefs <- coef(final_lasso)
  coefs_df <- as.data.frame(as.matrix(coefs))
  coefs_df <- tibble::rownames_to_column(coefs_df, "term")
  names(coefs_df)[2] <- "Estimate"
  
  # Identify nonzero features (excluding intercept if Estimate != 0)
  # The intercept is "term == '(Intercept)'"
  nonzero_features <- coefs_df %>%
    filter(Estimate != 0) %>%
    pull(term)
  
  # If only intercept is nonzero, there's no point refitting OLS
  if (length(nonzero_features) == 1 && "(Intercept)" %in% nonzero_features) {
    warning(paste("LASSO selected only the intercept for", response))
    # We'll just do a trivial OLS with intercept
    # or return NULL to skip
    return(lm(as.formula(paste(response, "~ 1")), data = df))
  }
  
  # Build formula for OLS with the selected features
  # remove the intercept from feature list for formula
  final_feats <- setdiff(nonzero_features, "(Intercept)")
  formula_ols <- as.formula(
    paste(response, "~", paste(final_feats, collapse = " + "))
  )
  
  # Refit OLS
  final_ols <- lm(formula_ols, data = df)
  return(final_ols)
}

# ============== GROUP DATA AND FIT ================
results <- merged_data %>%
  group_by(name) %>%
  group_split() %>%
  set_names(unique(merged_data$name)) %>%
  map(~ {
    df <- as.data.frame(.x)
    map(setNames(target_vars, target_vars), ~ fit_lasso_and_refit_ols(df, .))
  })

# Summaries
results_summaries <- map(results, function(model_list) { 
  map(model_list, ~ if (!is.null(.)) summary(.) else NULL) 
}) 

# ============== VARIABLE NAME MAPPING ================
variable_name_map <- c(
  # Climate
  ppt = "Precipitation (mm)",
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
  
  # Soil
  elevation = "Elevation (m)",
  sand = "Sand (%)",
  silt = "Silt (%)",
  clay = "Clay (%)",
  som = "Soil Organic Matter (%)",
  bulk_density = "Bulk Density (g/cm³)",
  awc = "Available Water Holding Capacity (cm)",
  depth_to_restrictive = "Depth to Restrictive Layer (cm)"
)

format_regression_table_png <- function(model_summary,
                                        soil_depth = "Depth Averaged",
                                        response_variable = "r",
                                        export_path = "regression_results.png") {
  # Check summary
  if (is.null(model_summary) || !inherits(model_summary, "summary.lm")) {
    stop("Invalid model: Model is NULL or not a valid summary.lm object.")
  }
  
  # Map response variable name
  response_variable_name <- dplyr::case_when(
    response_variable == "r"    ~ "Pearson's r",
    response_variable == "bias" ~ "Bias",
    response_variable == "rmse" ~ "RMSE",
    TRUE                        ~ response_variable
  )
  
  # Title
  table_title <- paste("LASSO Regression Results:", soil_depth, "Soil Moisture (", response_variable_name, ")", sep = " ")
  
  # Extract coefficients from OLS summary
  coef_summary <- model_summary$coefficients
  
  # R² and Adjusted R²
  r_squared <- round(model_summary$r.squared, 3)
  adj_r_squared <- round(model_summary$adj.r.squared, 3)
  
  # Convert to tibble
  reg_table <- tibble::as_tibble(coef_summary, rownames = "Predictor") %>%
    dplyr::select(-`t value`) %>%  # remove t-value if you want
    dplyr::rename(
      Estimate = `Estimate`,
      `Std. Error` = `Std. Error`,
      `p value` = `Pr(>|t|)`
    ) %>%
    mutate(
      # Round numeric columns
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `p value` = dplyr::case_when(
        `p value` < 0.001 ~ "<0.001",
        `p value` < 0.01  ~ "<0.01",
        `p value` < 0.05  ~ "<0.05",
        `p value` < 0.1   ~ "<0.1",
        TRUE             ~ as.character(round(`p value`, 4))
      ),
      # Recode predictor names if you have a mapping vector called variable_name_map
      Predictor = dplyr::recode(Predictor, !!!variable_name_map)
    )
  
  # Build tableGrob
  gt_grob <- gridExtra::tableGrob(reg_table, rows = NULL)
  
  # Styling
  gt_grob$theme <- gridExtra::ttheme_default(
    core = list(
      fg_params = list(fontsize = 14, fontface = "plain"),
      bg_params = list(fill = c("white", "gray95"), col = NA)
    ),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  # Prepare final PNG
  ragg::agg_png(export_path, width = 1200, height = 800, res = 150)
  grid.newpage()
  
  # Define a 3-row layout: row1 for title, row2 for table, row3 for R2 note
  layout <- grid.layout(
    nrow = 3,
    ncol = 1,
    heights = unit(c(1.5, 1, 1.5), c("lines","null","lines"))  # lines, null, lines
  )
  pushViewport(viewport(layout = layout))
  
  # 1) Title in row=1
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(table_title, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  
  # 2) Table in row=2 (auto-scaling space)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(gt_grob)
  popViewport()
  
  # 3) R2 note in row=3
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  r2_note <- paste0("R² = ", r_squared, " | Adjusted R² = ", adj_r_squared)
  grid.text(r2_note, gp = gpar(fontsize = 12, fontface = "italic"))
  popViewport()
  
  dev.off()
  message(paste("Table saved as", export_path))
}

# ============== EXAMPLE USAGE ================
# For "Depth Averaged" group, response = "r"
if(!is.null(results_summaries$`Depth Averaged`$r)) {
  format_regression_table_png(
    model_summary = results_summaries$`Depth Averaged`$r,
    soil_depth = "Depth Averaged",
    response_variable = "r",
    export_path = "~/nrcs-scd-soil-moisture-eval/tables/lasso-regression-results-depth-averaged-r.png"
  )
}

# Similarly for "bias" and "rmse"
if(!is.null(results_summaries$`Depth Averaged`$bias)) {
  format_regression_table_png(
    model_summary = results_summaries$`Depth Averaged`$bias,
    soil_depth = "Depth Averaged",
    response_variable = "bias",
    export_path = "~/nrcs-scd-soil-moisture-eval/tables/lasso-regression-results-depth-averaged-bias.png"
  )
}
if(!is.null(results_summaries$`Depth Averaged`$rmse)) {
  format_regression_table_png(
    model_summary = results_summaries$`Depth Averaged`$rmse,
    soil_depth = "Depth Averaged",
    response_variable = "rmse",
    export_path = "~/nrcs-scd-soil-moisture-eval/tables/lasso-regression-results-depth-averaged-rmse.png"
  )
}
