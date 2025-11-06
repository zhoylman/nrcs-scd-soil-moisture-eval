# install.packages(c("gt","webshot2"))  # if needed
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
library(gt)

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

# helper: "mean ± sd"
summarise_metrics <- function(df) {
  df %>%
    summarise(
      r    = sprintf("%.2f ± %.2f", mean(`Pearson's r`, na.rm = TRUE), sd(`Pearson's r`, na.rm = TRUE)),
      KGE  = sprintf("%.2f ± %.2f", mean(KGE,            na.rm = TRUE), sd(KGE,            na.rm = TRUE)),
      Bias = sprintf("%.2f ± %.2f", mean(Bias,           na.rm = TRUE), sd(Bias,           na.rm = TRUE)),
      RMSE = sprintf("%.2f ± %.2f", mean(RMSE,           na.rm = TRUE), sd(RMSE,           na.rm = TRUE)),
      .groups = "drop"
    )
}

make_eval_table <- function(all_results, outfile_png,
                            row_order = c("MT Mesonet","OK Mesonet","NEON","SCAN","SNTL","USCRN")) {
  
  depth <- "Depth Averaged"
  
  # Per-network
  per_net <-
    all_results %>%
    filter(generalized_depth == depth) %>%
    group_by(network, model) %>%
    summarise_metrics() %>%
    ungroup() %>%
    pivot_longer(c(r, KGE, Bias, RMSE), names_to = "Metric", values_to = "val") %>%
    mutate(colname = paste(model, Metric, sep = "_")) %>%
    select(network, colname, val) %>%
    pivot_wider(names_from = colname, values_from = val) %>%
    mutate(network = factor(network, levels = row_order)) %>%
    arrange(network) %>%
    mutate(network = as.character(network))
  
  # Nationwide (pooling all networks)
  nat <-
    all_results %>%
    filter(generalized_depth == depth) %>%
    group_by(model) %>%
    summarise_metrics() %>%
    ungroup() %>%
    pivot_longer(c(r, KGE, Bias, RMSE), names_to = "Metric", values_to = "val") %>%
    mutate(colname = paste(model, Metric, sep = "_")) %>%
    select(colname, val) %>%
    pivot_wider(names_from = colname, values_from = val) %>%
    mutate(network = "Nationwide")
  
  tbl <- bind_rows(per_net, nat) %>% rename(Network = network)
  
  gt_tbl <-
    tbl %>%
    gt(rowname_col = "Network") %>%
    # ---- ADD TITLE HERE ----
  tab_header(
    title = md("**Depth Averaged Results**"),
    subtitle = md("Comparison of SOILWAT2 and SPoRT-LIS models across networks")
  ) %>%
    # -------------------------
  tab_spanner("SOILWAT2", columns = c(SOILWAT2_r, SOILWAT2_KGE, SOILWAT2_Bias, SOILWAT2_RMSE)) %>%
    tab_spanner("SPoRT-LIS", columns = c(`SPoRT-LIS_r`, `SPoRT-LIS_KGE`, `SPoRT-LIS_Bias`, `SPoRT-LIS_RMSE`)) %>%
    cols_label(
      SOILWAT2_r = html("<b>r</b>"), SOILWAT2_KGE = html("<b>KGE</b>"),
      SOILWAT2_Bias = html("<b>Bias</b>"), SOILWAT2_RMSE = html("<b>RMSE</b>"),
      `SPoRT-LIS_r` = html("<b>r</b>"), `SPoRT-LIS_KGE` = html("<b>KGE</b>"),
      `SPoRT-LIS_Bias` = html("<b>Bias</b>"), `SPoRT-LIS_RMSE` = html("<b>RMSE</b>")
    ) %>%
    tab_options(
      table.font.size = px(14),
      data_row.padding = px(6),
      column_labels.font.weight = "bold",
      heading.align = "center"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_stub(rows = Network == "Nationwide")
    )
  
  gtsave(gt_tbl, outfile_png, vwidth = 1200, vheight = 420, expand = 5)
  invisible(gt_tbl)
}

# Run it
make_eval_table(
  all_results,
  outfile_png = "~/nrcs-scd-soil-moisture-eval/tables/summary_table_by_network.png"
)

make_eval_table_all_depths <- function(
    all_results,
    outfile_png,
    row_order  = c("MT Mesonet","OK Mesonet","NEON","SCAN","SNTL","USCRN"),
    depth_order = c("Shallow","Middle","Deep","Depth Averaged")
) {
  row_order_full <- c(row_order, "Nationwide")
  
  # ---- Per-network + depth ----
  per_net <-
    all_results %>%
    dplyr::filter(generalized_depth %in% depth_order) %>%
    dplyr::group_by(generalized_depth, network, model) %>%
    summarise_metrics() %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(c(r, KGE, Bias, RMSE), names_to = "Metric", values_to = "val") %>%
    dplyr::mutate(colname = paste(model, Metric, sep = "_")) %>%
    dplyr::select(generalized_depth, network, colname, val) %>%
    tidyr::pivot_wider(names_from = colname, values_from = val) %>%
    dplyr::mutate(
      generalized_depth = factor(generalized_depth, levels = depth_order)
    )
  
  # ---- Nationwide (pooled) for each depth ----
  nat <-
    all_results %>%
    dplyr::filter(generalized_depth %in% depth_order) %>%
    dplyr::group_by(generalized_depth, model) %>%
    summarise_metrics() %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(c(r, KGE, Bias, RMSE), names_to = "Metric", values_to = "val") %>%
    dplyr::mutate(colname = paste(model, Metric, sep = "_")) %>%
    dplyr::select(generalized_depth, colname, val) %>%
    tidyr::pivot_wider(names_from = colname, values_from = val) %>%
    dplyr::mutate(
      network = "Nationwide",
      generalized_depth = factor(generalized_depth, levels = depth_order)
    )
  
  tbl <- dplyr::bind_rows(per_net, nat) %>%
    dplyr::rename(Network = network, Depth = generalized_depth) %>%
    dplyr::mutate(
      # factor AFTER binding so "Nationwide" is included in levels
      Network = factor(as.character(Network), levels = row_order_full)
    ) %>%
    dplyr::arrange(Depth, Network)
  
  gt_tbl <-
    tbl %>%
    gt::gt(rowname_col = "Network", groupname_col = "Depth") %>%
    gt::tab_header(
      title = gt::md("**Model Evaluation Results by Depth**"),
      subtitle = gt::md("Comparison of SOILWAT2 and SPoRT-LIS across networks and soil depths")
    ) %>%
    gt::tab_spanner("SOILWAT2",
                    columns = c(SOILWAT2_r, SOILWAT2_KGE, SOILWAT2_Bias, SOILWAT2_RMSE)) %>%
    gt::tab_spanner("SPoRT-LIS",
                    columns = c(`SPoRT-LIS_r`, `SPoRT-LIS_KGE`, `SPoRT-LIS_Bias`, `SPoRT-LIS_RMSE`)) %>%
    gt::cols_label(
      SOILWAT2_r    = gt::html("<b>r</b>"),
      SOILWAT2_KGE  = gt::html("<b>KGE</b>"),
      SOILWAT2_Bias = gt::html("<b>Bias</b>"),
      SOILWAT2_RMSE = gt::html("<b>RMSE</b>"),
      `SPoRT-LIS_r`    = gt::html("<b>r</b>"),
      `SPoRT-LIS_KGE`  = gt::html("<b>KGE</b>"),
      `SPoRT-LIS_Bias` = gt::html("<b>Bias</b>"),
      `SPoRT-LIS_RMSE` = gt::html("<b>RMSE</b>")
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(14),
      data_row.padding = gt::px(6),
      column_labels.font.weight = "bold",
      heading.align = "center"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_stub(rows = Network == "Nationwide")
    )
  
  gt::gtsave(gt_tbl, outfile_png, vwidth = 1200, vheight = 800, expand = 5)
  message("Saved: ", outfile_png)
  invisible(gt_tbl)
}

make_eval_table_all_depths(
  all_results,
  outfile_png = "~/nrcs-scd-soil-moisture-eval/tables/summary_table_all_depths.png"
)
