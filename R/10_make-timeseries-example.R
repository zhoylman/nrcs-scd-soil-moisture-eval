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

SPoRT = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/SPoRT-LIS-corelations-generalized-depth.csv") |>
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

all_results = bind_rows(SPoRT, soilwat2)

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

#sport
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
  arrange(network, site_id, date, generalized_depth) |>
  filter(generalized_depth == 'Depth Averaged') |>
  select(site_id, date, `SPoRT-LIS` = SPoRT)

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
  mutate(`SPoRT-LIS` = ifelse(`SPoRT-LIS` == 0, NA, `SPoRT-LIS`),
         SOILWAT2 = ifelse(SOILWAT2 == 0, NA, SOILWAT2),
         diff = `SPoRT-LIS` - SOILWAT2) |>
  filter(`SPoRT-LIS` > 0.2 & SOILWAT2 > 0.3,
         diff < 0 & diff >-0.8) |>
  drop_na() |>
  mutate(avg_KGE = (`SPoRT-LIS` + SOILWAT2) / 2,
         abs_diff = abs(diff)) |>
  group_by(network) |>
  arrange(desc(avg_KGE)) |>
  #filter(!site_id %in% c('KONZ', 'BLAN', 'SCBI', 'MLBS', 'SCAN:2227', 'SCAN:2223', 'SNTL:523', 'SNTL:767', 'PORT', 'BIXB')) |>
  filter(!site_id %in% c('ALTU', 'reedpoin')) |>
  slice(1) |>
  ungroup() |>
  select(-avg_KGE, -abs_diff) 

year_of_interest = staion_data %>% 
  left_join(., SPoRT_data_final, by = c('site_id', 'date')) |>
  filter(site_id %in% sites_of_interest$site_id) |>
  drop_na() |>
  mutate(year = lubridate::year(date)) |>
  group_by(site_id, year) |>
  summarise(n = length(year)) |>
  group_by(site_id) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() 

sport_data_of_interest = staion_data %>% 
  left_join(., SPoRT_data_final, by = c('site_id', 'date')) |>
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
  "namlower"    = "MT Mesonet: Hot Springs NE Lower",
  "SJER"        = "NEON: SJER",
  "STIG"        = "OK Mesonet: STIG",
  "SCAN:2223"   = "SCAN: 2223",
  "SNTL:1084"    = "SNOTEL: 1084",
  "1510"        = "USCRN: 1510"
)

#order
order_levels = c("namlower", "STIG", "SJER", "SCAN:2223", "SNTL:1084", "1510")

plot = ggplot(station_data_of_interest, aes(x = date)) +
  geom_line(
    aes(y = `Depth Averaged`, color = "Observed (Sensor)"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  geom_line(
    data = sport_data_of_interest,
    aes(y = `SPoRT-LIS`, color = "SPoRT-LIS Model"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  geom_line(
    data = soilwat2_data_of_interest,
    aes(x = date, y = SOILWAT2, color = "SOILWAT2 Model"),
    size = 0.7, alpha = 0.8, na.rm = FALSE
  ) +
  facet_wrap(
    vars(factor(site_id, levels = order_levels)),
    labeller = labeller(`factor(site_id, levels = order_levels)` = custom_labels),
    scales = 'free_x',
    ncol = 2
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
      "SPoRT-LIS Model"    = "#d95f02",
      "SOILWAT2 Model"    = "#7570b3"
    )
  ) +
  labs(
    title = "Soil Moisture Time Series Comparison",
    subtitle = "Observed vs. SPoRT-LIS and SOILWAT2 Model Estimates",
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
