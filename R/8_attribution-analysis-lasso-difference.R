# ========== Load Libraries ==========
library(tidyverse)
library(glmnet)
library(gridExtra)
library(grid)
library(ragg)

# ========== Load Data ==========
station_covariates <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv")
station_meta <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv")

topofire_corelations <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv") %>%
  rename(KGE = `KGE`, r = `Pearson's r`) %>%
  mutate(model = "Topofire")

soilwat2_corelations <- read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv") %>%
  rename(KGE = `KGE`, r = `Pearson's r`) %>%
  mutate(model = "SOILWAT2")

# Combine
all_corelations <- bind_rows(topofire_corelations, soilwat2_corelations) %>%
  filter(generalized_depth == 'Depth Averaged')

# ========== Prepare Full Merged Dataset ==========
merged_all_models <- all_corelations %>%
  left_join(station_covariates, by = c("network", "site_id")) %>%
  filter(!is.na(r) & !is.na(KGE))

# For difference model
accuracy_diff_data <- merged_all_models %>%
  group_by(network, site_id, model) %>%
  summarise(
    r = mean(r, na.rm = TRUE),
    KGE = mean(KGE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = model, values_from = c(r, KGE)) %>%
  left_join(station_covariates, by = c("network", "site_id")) %>%
  mutate(
    r_diff = r_Topofire - r_SOILWAT2,
    KGE_diff = KGE_Topofire - KGE_SOILWAT2
  )
# ========== Variables ==========
predictors <- setdiff(names(station_covariates), c("network", "site_id"))
target_sets <- list(
  Topofire = c("r", "KGE"),
  SOILWAT2 = c("r", "KGE"),
  Difference = c("r_diff", "KGE_diff")
)

# ========== LASSO + OLS Function ==========
fit_lasso_and_refit_ols <- function(df, response) {
  df <- df %>% select(all_of(c(response, predictors))) %>% drop_na()
  if (nrow(df) == 0) return(NULL)
  
  X <- model.matrix(as.formula(paste(response, "~", paste(predictors, collapse = "+"))), data = df)[, -1]
  y <- df[[response]]
  cv_fit <- cv.glmnet(X, y, alpha = 1)
  best_lambda <- cv_fit$lambda.min
  final_lasso <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  coefs <- coef(final_lasso)
  nonzero <- rownames(coefs)[coefs[, 1] != 0]
  final_feats <- setdiff(nonzero, "(Intercept)")
  if (length(final_feats) == 0) return(lm(as.formula(paste(response, "~ 1")), data = df))
  lm(as.formula(paste(response, "~", paste(final_feats, collapse = "+"))), data = df)
}

# ========== Variable Name Map ==========
variable_name_map <- c(
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
  elevation = "Elevation (m)",
  sand = "Sand (%)",
  silt = "Silt (%)",
  clay = "Clay (%)",
  som = "Soil Organic Matter (%)",
  bulk_density = "Bulk Density (g/cm³)",
  awc = "Available Water Holding Capacity (cm)",
  depth_to_restrictive = "Depth to Restrictive Layer (cm)",
  topofire_aws = "SSURGO AWS (mm)",
  topofire_porosity = "SSURGO Porosity (mm/m)"
)

# ========== Table Export Function ==========
format_regression_table_png <- function(model_summary, soil_depth, response_variable, export_path) {
  if (is.null(model_summary) || !inherits(model_summary, "summary.lm")) return()
  
  response_label <- case_when(
    response_variable == "r" ~ "Pearson's r",
    response_variable == "KGE" ~ "Kling-Gupta Efficiency (KGE)",
    response_variable == "r_diff" ~ "Pearson's r Difference",
    response_variable == "KGE_diff" ~ "KGE Difference",
    TRUE ~ response_variable
  )
  
  title <- paste("LASSO Regression Results:", soil_depth, response_label)
  
  coefs <- model_summary$coefficients
  r2 <- round(model_summary$r.squared, 3)
  adj_r2 <- round(model_summary$adj.r.squared, 3)
  
  reg_table <- as_tibble(coefs, rownames = "Predictor") %>%
    select(-`t value`) %>%
    rename(`Std. Error` = `Std. Error`, `p value` = `Pr(>|t|)`) %>%
    mutate(
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `p value` = case_when(
        `p value` < 0.001 ~ "<0.001",
        `p value` < 0.01 ~ "<0.01",
        `p value` < 0.05 ~ "<0.05",
        `p value` < 0.1 ~ "<0.1",
        TRUE ~ as.character(round(`p value`, 4))
      ),
      Predictor = recode(Predictor, !!!variable_name_map)
    )
  
  gt_grob <- tableGrob(reg_table, rows = NULL)
  gt_grob$theme <- ttheme_default(
    core = list(fg_params = list(fontsize = 14), bg_params = list(fill = c("white", "gray95"))),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  ragg::agg_png(export_path, width = 1200, height = 800, res = 150)
  grid.newpage()
  layout <- grid.layout(nrow = 3, heights = unit(c(1.5, 1, 1.5), c("lines", "null", "lines")))
  pushViewport(viewport(layout = layout))
  pushViewport(viewport(layout.pos.row = 1)); grid.text(title, gp = gpar(fontsize = 14, fontface = "bold")); popViewport()
  pushViewport(viewport(layout.pos.row = 2)); grid.draw(gt_grob); popViewport()
  pushViewport(viewport(layout.pos.row = 3)); grid.text(paste("R² =", r2, "| Adjusted R² =", adj_r2), gp = gpar(fontsize = 12, fontface = "italic")); popViewport()
  dev.off()
  message("Saved:", export_path)
}

# ========== MODEL FITTING AND EXPORT LOOP ==========
fit_and_export_models <- function(data, response_vars, prefix) {
  map(response_vars, function(resp) {
    model <- fit_lasso_and_refit_ols(data, resp)
    if (!is.null(model)) {
      summary_obj <- summary(model)
      out_path <- glue::glue("~/nrcs-scd-soil-moisture-eval/tables/lasso-{prefix}-{tolower(resp)}.png")
      format_regression_table_png(summary_obj, "Depth Averaged", resp, out_path)
    }
  })
}

fit_and_export_models(all_corelations %>% filter(model == 'Topofire'), c("r", "KGE"), "Topofire")
fit_and_export_models(soilwat2_data, c("r", "KGE"), "soilwat2")
fit_and_export_models(accuracy_diff_data,     c("r_diff", "KGE_diff"), "difference")
