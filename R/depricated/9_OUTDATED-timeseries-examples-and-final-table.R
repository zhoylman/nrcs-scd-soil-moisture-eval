library(tidyverse)

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

topofire = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/topofire-corelations-generalized-depth.csv") |>
  mutate(model = 'Topofire') |>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            NSE = median(NSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup()

soilwat2 = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/soilwat2-corelations-generalized-depth.csv") |>
  mutate(model = 'SOILWAT2')|>
  #average results across horizontal positions in neon
  group_by(network, site_id,  generalized_depth, model) |>
  summarise(`Pearson's r` = median(`Pearson's r`, na.rm = T),
            Bias = median(Bias, na.rm = T),
            RMSE = median(RMSE, na.rm = T),
            NSE = median(NSE, na.rm = T),
            KGE = median(KGE, na.rm = T)) |>
  ungroup()

all_results = bind_rows(topofire, soilwat2)

#pull in topofire data
staion_data = {
  read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv") %>%
    pivot_wider(names_from = generalized_depth, values_from = soil_moisture) %>%
    mutate(
      site_id_long = site_id,
      site_id = case_when(
        network == "NEON" ~ str_remove(site_id, "_.*"),
        TRUE              ~ site_id
      )
    )
} %>% 
  left_join(
    {
      read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv") %>%
        pivot_wider(names_from = generalized_depth, values_from = soil_moisture) %>%
        mutate(
          site_id = case_when(
            network == "NEON" ~ str_remove(site_id, "_.*"),
            TRUE              ~ site_id
          )
        ) %>%
        filter(network == "NEON") %>%
        group_by(site_id, date) %>%
        summarise(
          Depth_Averaged_median = median(`Depth Averaged`, na.rm = TRUE),
          Depth_Averaged_median = `Depth Averaged`[2],
          Middle_median         = median(Middle, na.rm = TRUE),
          Shallow_median        = median(Shallow, na.rm = TRUE),
          .groups = "drop"
        )
    },
    by = c("site_id", "date")
  ) %>%
  mutate(
    `Depth Averaged` = if_else(network == "NEON", Depth_Averaged_median, `Depth Averaged`),
    Middle           = if_else(network == "NEON", Middle_median, Middle),
    Shallow          = if_else(network == "NEON", Shallow_median, Shallow)
  ) %>%
  select(-Depth_Averaged_median, -Middle_median, -Shallow_median)

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

topofire_infill = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/topofire_infill/raw", full.names = T) |>
  map(read_csv, show_col_types = F) |>
  bind_rows() |>
  select(-nc_id) |>
  pivot_longer(cols = -c(time)) |>
  # need to account for 'WOOD' duplicates here. 
  separate(name, into = c("network", "site_id"), sep = "_", remove = FALSE) |>
  select(site_id, date = time, Topofire = value)

topofire_data_long = topofire_data %>% 
  dplyr::select(site_id, common_dates %>% as.character()) %>%
  pivot_longer(cols = -c(site_id)) %>%
  dplyr::select(site_id, date = name, Topofire = value) %>%
  drop_na(Topofire) %>%
  mutate(date = date %>% as.Date()) |>
  bind_rows(topofire_infill)
  
#rm(topofire_data); gc(); gc()

#compute VWC using available water holding capacity and porosity

max_topofire = topofire_data_long |>
  group_by(site_id) |>
  summarise(max_topofire = max(Topofire, na.rm = T))

params = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv") |>
  dplyr::select(site_id, topofire_aws, awc_gNATSGO, porosity_gNATSGO, fc_gNATSGO, aws0_999) |>
  left_join(max_topofire)

topofire_data_final = topofire_data_long %>% 
  left_join(params, 'site_id')  %>%
  mutate(Topofire_VWC_old = (Topofire / max_topofire) * (porosity_gNATSGO / 1000)) %>%
  mutate(Topofire_VWC_porosity = (Topofire / porosity_gNATSGO) * (porosity_gNATSGO / 1000)) %>%
  mutate(Topofire_VWC = ((Topofire / aws0_999) * (awc_gNATSGO/1000)) + ((awc_gNATSGO - fc_gNATSGO) / 1000)) %>%
  dplyr::select(site_id, date, Topofire = Topofire_VWC_old)

#rm(topofire_data_long); gc(); gc()

#soilwat2

files = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/20250409_SOILWAT2_OutputVWCShared__20250320_SoilClimateDynamics_SOILWAT2_simulations", full.names = T)

soilwat_raw = files %>%
  purrr::map(read_soilwat) %>%
  bind_rows()

soilwat = soilwat_raw %>%
  pivot_longer(cols = -c(network, site_id, date), names_to = 'generalized_depth', values_to = 'SOILWAT2') |>
  filter(generalized_depth == 'Depth Averaged') |>
  select(-c(generalized_depth, network))


############

sites_of_interest = all_results |>
  select(network, site_id, generalized_depth,
         model, KGE) %>%
  filter(generalized_depth == 'Depth Averaged') |>
  pivot_wider(names_from = model, values_from = KGE) |>
  mutate(Topofire = ifelse(Topofire == 0, NA, Topofire),
         SOILWAT2 = ifelse(SOILWAT2 == 0, NA, SOILWAT2),
         diff = Topofire - SOILWAT2) |>
  filter(Topofire > 0.2 & SOILWAT2 > 0.3,
         diff < 0 & diff >-0.8) |>
  drop_na() |>
  mutate(avg_KGE = (Topofire + SOILWAT2) / 2,
         abs_diff = abs(diff)) |>
  group_by(network) |>
  arrange(desc(avg_KGE)) |>
  filter(!site_id %in% c('KONZ', 'BLAN', 'SCBI', 'MLBS', 'SCAN:2227', 'SCAN:2223', 'SNTL:523', 'SNTL:767', 'PORT', 'BIXB')) |>
  slice(1) |>
  ungroup() |>
  select(-avg_KGE, -abs_diff) 

year_of_interest = staion_data %>% 
  left_join(., topofire_data_final, by = c('site_id', 'date')) |>
  filter(site_id %in% sites_of_interest$site_id) |>
  drop_na() |>
  mutate(year = lubridate::year(date)) |>
  group_by(site_id, year) |>
  summarise(n = length(year)) |>
  group_by(site_id) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() 

topofire_data_of_interest = staion_data %>% 
  left_join(., topofire_data_final, by = c('site_id', 'date')) |>
  filter(site_id %in% sites_of_interest$site_id) |>
  mutate(year = lubridate::year(date)) %>%
  semi_join(year_of_interest, by = c("site_id", "year"))

soilwat2_data_of_interest = staion_data %>% 
  left_join(., soilwat, by = c('site_id', 'date')) |>
  filter(site_id %in% sites_of_interest$site_id) |>
  drop_na() |>
  mutate(year = lubridate::year(date)) %>%
  semi_join(year_of_interest, by = c("site_id", "year"))

station_data_of_interest = staion_data %>% 
  filter(site_id %in% sites_of_interest$site_id) |>
  mutate(year = lubridate::year(date)) %>%
  semi_join(year_of_interest, by = c("site_id", "year")) 

custom_labels = c(
  "mdamoies"    = "MT Mesonet: Moiese N",
  "DELA"        = "NEON: DELA",
  "STIG"        = "OK Mesonet: STIG",
  "SCAN:2030"   = "SCAN: 2030",
  "SNTL:1056"    = "SNOTEL: 1056",
  "1037"        = "USCRN: 1037"
)

plot = ggplot(station_data_of_interest, aes(x = date)) +
  geom_line(
    aes(y = `Depth Averaged`, color = "Observed (Sensor)"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  geom_line(
    data = topofire_data_of_interest,
    aes(y = Topofire, color = "Topofire Model"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  geom_line(
    data = soilwat2_data_of_interest,
    aes(x = date, y = SOILWAT2, color = "SOILWAT2 Model"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  facet_wrap(
    ~site_id,
    scales = 'free_x',
    ncol = 2,
    labeller = labeller(site_id = custom_labels)
  ) +
  scale_x_date(
    date_labels = "%b",
    date_breaks = "1 month",
    expand = c(0.01, 0.01)
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed (Sensor)" = "black",
      "Topofire Model"    = "#d95f02",
      "SOILWAT2 Model"    = "#7570b3"
    )
  ) +
  labs(
    title = "Soil Moisture Time Series Comparison",
    subtitle = "Observed vs. Topofire and SOILWAT2 Model Estimates",
    x = "Time",
    y = "Volumetric Water Content (m³/m³)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.spacing = unit(1.5, "lines")
  )

plot

ggsave('~/nrcs-scd-soil-moisture-eval/figs/timeseries_example.png', plot, width = 10, height = 8, dpi = 300)


# make table that describes model performance

summary_table_by_network = all_results |>
  pivot_longer(-c(network, site_id, generalized_depth, model),
               names_to = 'Metric') |>
  filter(generalized_depth == 'Depth Averaged') |>
  group_by(network, model, Metric) |>
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) |>
  mutate(summary = glue::glue('{round(mean, 2)} ± {round(sd, 2)}')) |>
  select( network, model, Metric, summary) |>
  pivot_wider(
              names_from = model, values_from = summary
              ) |>
  select(Network = network,
         Metric, SOILWAT2, Topofire) |>
  ungroup()

format_summary_table_png <- function(summary_table, export_path, prefix = NULL) {
  if (is.null(summary_table) || !is.data.frame(summary_table)) return()
  
  title_line1 <- "Model Evaluation Summary"
  title_line2 <- if (!is.null(prefix)) prefix else NULL
  
  # Replace repeated Network values with blanks (for visual grouping)
  display_table <- summary_table %>%
    dplyr::mutate(Network = dplyr::case_when(
      dplyr::row_number() == 1 ~ Network,
      Network != dplyr::lag(Network) ~ Network,
      TRUE ~ ""
    )) %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric,
                             "KGE" = "Kling-Gupta Efficiency (KGE)",
                             "NSE" = "Nash-Sutcliffe Efficiency (NSE)",
                             "RMSE" = "Root Mean Square Error (RMSE)",
                             "Bias" = "Bias",
                             "Pearson's r" = "Pearson's r"
      )
    ) %>%
    dplyr::rename(
      `SOILWAT2` = SOILWAT2,
      `Topofire` = Topofire
    )
  
  # Adjust plot height
  n_rows <- nrow(display_table)
  row_height <- 35
  base_height <- 400
  img_height <- base_height + (n_rows * row_height)
  
  # Create table grob
  gt_grob <- gridExtra::tableGrob(display_table, rows = NULL)
  gt_grob$theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 13), bg_params = list(fill = c("white", "gray95"))),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  # Identify row indices where the Network label starts (i.e., first row of each group)
  group_start_rows <- which(summary_table$Network != dplyr::lag(summary_table$Network, default = ""))
  
  # Add lines *above* each new group (except the first one)
  border_rows <- group_start_rows[group_start_rows != 1] - 1
  
  for (r in border_rows) {
    gt_grob <- gtable::gtable_add_grob(
      gt_grob,
      grobs = segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
        gp = gpar(col = "black", lwd = 2)
      ),
      t = r + 1, b = r + 1, l = 1, r = ncol(gt_grob)
    )
  }
  
  # Render to PNG
  ragg::agg_png(export_path, width = 1000, height = img_height, res = 150)
  grid.newpage()
  layout <- grid.layout(nrow = 3, heights = unit(c(3.5, 1, 1.5), c("lines", "null", "lines")))
  pushViewport(viewport(layout = layout))
  
  # Title
  pushViewport(viewport(layout.pos.row = 1))
  grid.text(title_line1, y = unit(0.75, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))
  if (!is.null(title_line2)) {
    grid.text(title_line2, y = unit(0.35, "npc"), gp = gpar(fontsize = 13))
  }
  popViewport()
  
  # Table
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(gt_grob)
  popViewport()
  
  # Footer
  pushViewport(viewport(layout.pos.row = 3))
  grid.text("Values shown as: mean ± sd", gp = gpar(fontsize = 11, fontface = "italic"))
  popViewport()
  
  dev.off()
  message("Saved: ", export_path)
}

format_summary_table_png <- function(summary_table, export_path, prefix = NULL) {
  if (is.null(summary_table) || !is.data.frame(summary_table)) return()
  
  title_line1 <- "Model Evaluation Summary"
  title_line2 <- if (!is.null(prefix)) prefix else NULL
  
  display_table <- summary_table %>%
    dplyr::mutate(Network = dplyr::case_when(
      dplyr::row_number() == 1 ~ Network,
      Network != dplyr::lag(Network) ~ Network,
      TRUE ~ ""
    )) %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric,
                             "KGE" = "Kling-Gupta Efficiency (KGE)",
                             "NSE" = "Nash-Sutcliffe Efficiency (NSE)",
                             "RMSE" = "Root Mean Square Error (RMSE)",
                             "Bias" = "Bias",
                             "Pearson's r" = "Pearson's r"
      )
    ) %>%
    dplyr::rename(`SOILWAT2` = SOILWAT2, `Topofire` = Topofire)
  
  n_rows <- nrow(display_table)
  row_height <- 35
  base_height <- 400
  img_height <- base_height + (n_rows * row_height)
  
  gt_grob <- gridExtra::tableGrob(display_table, rows = NULL)
  gt_grob$theme <- gridExtra::ttheme_default(
    core   = list(fg_params = list(fontsize = 13),
                  bg_params = list(fill = c("white", "gray95"))),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
  
  # add bold horizontal rule above each new Network group (except first)
  group_start_rows <- which(summary_table$Network !=
                              dplyr::lag(summary_table$Network, default = ""))
  border_rows <- group_start_rows[group_start_rows != 1] - 1
  for (r in border_rows) {
    gt_grob <- gtable::gtable_add_grob(
      gt_grob,
      grobs = grid::segmentsGrob(
        x0 = grid::unit(0, "npc"), x1 = grid::unit(1, "npc"),
        y0 = grid::unit(0, "npc"), y1 = grid::unit(0, "npc"),
        gp = grid::gpar(col = "black", lwd = 2)
      ),
      t = r + 1, b = r + 1, l = 1, r = ncol(gt_grob)
    )
  }
  
  ragg::agg_png(export_path, width = 1000, height = img_height, res = 150)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 3,
                           heights = grid::unit(c(3.5, 1, 1.5), c("lines","null","lines")))
  grid::pushViewport(grid::viewport(layout = lay))
  
  # Title
  grid::pushViewport(grid::viewport(layout.pos.row = 1))
  grid::grid.text(title_line1, y = grid::unit(0.75, "npc"),
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))
  if (!is.null(title_line2)) {
    grid::grid.text(title_line2, y = grid::unit(0.35, "npc"),
                    gp = grid::gpar(fontsize = 13))
  }
  grid::popViewport()
  
  # Table
  grid::pushViewport(grid::viewport(layout.pos.row = 2))
  grid::grid.draw(gt_grob)
  grid::popViewport()
  
  # Footer
  grid::pushViewport(grid::viewport(layout.pos.row = 3))
  grid::grid.text("Values shown as: mean ± sd",
                  gp = grid::gpar(fontsize = 11, fontface = "italic"))
  grid::popViewport()
  
  grDevices::dev.off()
  message("Saved: ", export_path)
}

summary_table = all_results |>
  pivot_longer(-c(network, site_id, generalized_depth, model),
               names_to = 'Metric') |>
  filter(generalized_depth == 'Depth Averaged') |>
  group_by(model, Metric) |>
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) |>
  mutate(summary = glue::glue('{round(mean, 2)} ± {round(sd, 2)}')) |>
  select(model, Metric, summary) |>
  pivot_wider(
    names_from = model, values_from = summary
  ) 

format_summary_table_national_png <- function(summary_table, export_path, prefix = NULL) {
  if (is.null(summary_table) || !is.data.frame(summary_table)) return()
  
  title_line1 <- "National Model Evaluation Summary"
  title_line2 <- if (!is.null(prefix)) prefix else NULL
  
  # Format Metric names
  display_table <- summary_table %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric,
                             "KGE" = "Kling-Gupta Efficiency (KGE)",
                             "NSE" = "Nash-Sutcliffe Efficiency (NSE)",
                             "RMSE" = "Root Mean Square Error (RMSE)",
                             "Bias" = "Bias",
                             "Pearson's r" = "Pearson's r"
      )
    ) %>%
    dplyr::rename(
      `SOILWAT2` = SOILWAT2,
      `Topofire` = Topofire
    )
  
  # Image height
  n_rows <- nrow(display_table)
  row_height <- 45
  base_height <- 300
  img_height <- base_height + (n_rows * row_height)
  
  # Table grob
  gt_grob <- gridExtra::tableGrob(display_table, rows = NULL)
  gt_grob$theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 14), bg_params = list(fill = c("white", "gray95"))),
    colhead = list(fg_params = list(fontsize = 15, fontface = "bold"))
  )
  
  # Render PNG
  ragg::agg_png(export_path, width = 800, height = img_height, res = 150)
  grid.newpage()
  layout <- grid.layout(nrow = 3, heights = unit(c(4, 1, 1.5), c("lines", "null", "lines")))
  pushViewport(viewport(layout = layout))
  
  # Titles
  pushViewport(viewport(layout.pos.row = 1))
  grid.text(title_line1, y = unit(0.75, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
  if (!is.null(title_line2)) {
    grid.text(title_line2, y = unit(0.35, "npc"), gp = gpar(fontsize = 13))
  }
  popViewport()
  
  # Table
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(gt_grob)
  popViewport()
  
  # Footer
  pushViewport(viewport(layout.pos.row = 3))
  grid.text("Values shown as: mean ± sd", gp = gpar(fontsize = 11, fontface = "italic"))
  popViewport()
  
  dev.off()
  message("Saved: ", export_path)
}

library(grid)

format_summary_table_national_png(
  summary_table = summary_table,
  export_path = "~/nrcs-scd-soil-moisture-eval/tables/summary_table_all.png",
  prefix = "Depth Averaged Soil Moisture Results"
)

format_summary_table_png(
  summary_table = summary_table_by_network,
  export_path = "~/nrcs-scd-soil-moisture-eval/tables/summary_table_by_network.png",
  prefix = "Depth Averaged Soil Moisture Results by Network"
)

