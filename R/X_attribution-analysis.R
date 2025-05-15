library(tidyverse)
library(MASS)
library(gt)
library(dplyr)
library(gridExtra)
library(grid)
library(ragg)

# Load data
station_meta = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")
station_covariates = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv")
topofire_corelations = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv")

# Merge and clean data
merged_data = topofire_corelations %>%
  dplyr::select(-c("Pearson's r (Normalized)")) %>%
  left_join(., station_covariates, by = c('network', 'site_id')) %>%
  dplyr::select(-c('network', 'site_id')) %>%
  rename('r' = `Pearson's r`, 'bias' = `Bias (Normalized)`, 'rmse' = `RMSE (Normalized)`)

# Define response variables
target_vars = c("r", "bias", "rmse")

# Define predictor variables (excluding 'name' and response variables)
predictors = setdiff(names(merged_data), c("name", target_vars))

# Function to fit models with AIC selection (removing NAs)
fit_mlr_aic = function(df, response) {
  df = df %>%
    dplyr::select(all_of(c(response, predictors))) %>%  # Select only relevant columns
    drop_na()  # Remove rows with NA values
  
  if (nrow(df) == 0) {
    warning(paste("Skipping model for", response, "- No valid data after removing NAs"))
    return(NULL)  # Return NULL if no data left
  }
  
  formula = as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  full_model = lm(formula, data = df)
  
  # Check if model has enough degrees of freedom for AIC selection
  if (length(coef(full_model)) >= nrow(df)) {
    warning(paste("Skipping AIC selection for", response, "- Too few observations"))
    return(full_model)  # Return full model without AIC selection
  }
  
  best_model = stepAIC(full_model, direction = "both", trace = FALSE)
  return(best_model)
}

# Split data by 'name' and fit models for each target variable
results = merged_data %>%
  group_by(name) %>%
  group_split() %>%
  set_names(unique(merged_data$name)) %>%
  map(~ {
    df = as.data.frame(.x)  # Convert tibble to data frame
    map(setNames(target_vars, target_vars), ~ fit_mlr_aic(df, .))
  })


# Print summary of one model
results_summaries = map(results, function(model_list) { 
  map(model_list, ~ if (!is.null(.)) summary(.) else NULL) 
})

# Example: View results for 'Depth Averaged' and 'r'
results_summaries$`Depth Averaged`$r

# Define a named list to map variable names to real names
variable_name_map = c(
  # Climate variables
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
  
  # Soil variables
  elevation = "Elevation (m)",
  sand = "Sand (%)",
  silt = "Silt (%)",
  clay = "Clay (%)",
  som = "Soil Organic Matter (%)",
  bulk_density = "Bulk Density (g/cm³)",
  awc = "Available Water Holding Capacity (cm)",
  depth_to_restrictive = "Depth to Restrictive Layer (cm)"
)

# Function to create and export a regression table PNG with R² as a footnote
format_regression_table_png = function(model_summary, soil_depth = "Depth Averaged", response_variable = "r", export_path = "regression_results.png") {
  # Ensure the input is a valid summary.lm object
  if (is.null(model_summary) || !inherits(model_summary, "summary.lm")) {
    stop("Invalid model: Model is NULL or not a valid summary.lm object. Please select a valid model.")
  }
  
  # Properly format the response variable name
  response_variable_name = case_when(
    response_variable == "r" ~ "Pearson's r",
    response_variable == "bias" ~ "Bias",
    response_variable == "rmse" ~ "RMSE",
    TRUE ~ response_variable  # Default fallback (if needed)
  )
  
  # Generate dynamic title
  table_title = paste("Regression Results:", soil_depth, "Soil Moisture (", response_variable_name, ")", sep = " ")
  
  # Extract coefficients
  coef_summary = model_summary$coefficients
  
  # Extract R² and Adjusted R²
  r_squared = round(model_summary$r.squared, 4)
  adj_r_squared = round(model_summary$adj.r.squared, 4)
  
  # Convert to a tibble, rename columns, round values, and apply significance formatting in `p value`
  reg_table = as_tibble(coef_summary, rownames = "Predictor") %>%
    dplyr::select(-`t value`) %>%  # Remove t-value column
    rename(
      Estimate = `Estimate`,
      `Std. Error` = `Std. Error`,
      `p value` = `Pr(>|t|)`
    ) %>%
    mutate(
      Estimate = round(Estimate, 4),         
      `Std. Error` = round(`Std. Error`, 4),
      `p value` = case_when(
        `p value` < 0.001 ~ "<0.001", 
        `p value` < 0.01  ~ "<0.01",
        `p value` < 0.05  ~ "<0.05",
        `p value` < 0.1   ~ "<0.1",
        TRUE              ~ as.character(round(`p value`, 4))
      ),
      # Apply variable name mapping using dplyr::recode()
      Predictor = recode(Predictor, !!!variable_name_map)
    ) %>%
    dplyr::select(Predictor, Estimate, `Std. Error`, `p value`)  # Ensure correct column order
  
  # Convert to `gridExtra::tableGrob()` for PNG export (Removes row numbers)
  gt_grob = gridExtra::tableGrob(reg_table, rows = NULL)  # rows = NULL removes row numbers
  
  # Apply styling for a `kableExtra`-like appearance
  gt_grob$theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 14, fontface = "plain"),
      bg_params = list(fill = c("white", "gray95"), col = NA)  # Alternating row colors
    ),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  # **Save as PNG using ragg (No Chrome Required!)**
  ragg::agg_png(export_path, width = 1200, height = 800, res = 150)
  grid.newpage()
  
  # Add title manually (consistent spacing)
  grid.text(table_title, y = 0.88, gp = gpar(fontsize = 14, fontface = "bold"))
  
  # Draw the table closer to the title
  grid.draw(gt_grob)
  
  # Add R² note at the bottom
  grid.text(
    paste("R² =", r_squared, "| Adjusted R² =", adj_r_squared),
    y = 0.1, gp = gpar(fontsize = 12, fontface = "italic")
  )
  
  dev.off()
  
  message(paste("Table saved as", export_path))
}

format_regression_table_png(
  model_summary = results_summaries$`Depth Averaged`$r,
  soil_depth = "Depth Averaged",
  response_variable = "r",
  export_path = "~/nrcs-scd-soil-moisture-eval/tables/regression-results-depth-averaged-r.png"
)

 format_regression_table_png(
  model_summary = results_summaries$`Depth Averaged`$bias,
  soil_depth = "Depth Averaged",
  response_variable = "bias",
  export_path = "~/nrcs-scd-soil-moisture-eval/tables/regression-results-depth-averaged-bias.png"
)

format_regression_table_png(
  model_summary = results_summaries$`Depth Averaged`$rmse,
  soil_depth = "Depth Averaged",
  response_variable = "rmse",
  export_path = "~/nrcs-scd-soil-moisture-eval/tables/regression-results-depth-averaged-rmse.png"
)
