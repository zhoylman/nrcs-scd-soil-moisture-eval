library(tidyverse)
library(data.table)
library(sf)

#read in data 
mt_mesonet = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/mt-mesonet-raw.csv") %>%
  pivot_longer(cols = -c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         var = str_sub(name, 6,6),
         var = ifelse(var == 't', 'temperature', 'moisture'),
         network = 'MT Mesonet') %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = var, values_from = value) %>%
  rowwise() %>%
  #filter for temperature threshold
  mutate(moisture_corrected = ifelse(temperature < 1, NA, moisture)/100) %>%
  dplyr::select(network, site_id, date, depth, moisture_corrected)
  
ok_mesonet = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/ok-mesonet-raw.csv") %>%
  pivot_longer(cols = -c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         network = 'OK Mesonet') %>%
  dplyr::select(network, site_id, date, depth, moisture_corrected = value)

uscrn_nrcs = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/uscrn-sntl-scan-raw.csv") %>%
  pivot_longer(cols = -c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         var = str_sub(name, 6,6),
         var = ifelse(var == 't', 'temperature', 'moisture'),
         network = ifelse(str_detect(site_id, 'SCAN'), 'SCAN',
                          ifelse(str_detect(site_id, 'SNTL'), 'SNTL', 'USCRN'))) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = var, values_from = value) %>%
  rowwise() %>%
  #filter for temperature threshold 
  mutate(moisture_corrected = ifelse(temperature < 34, NA, moisture),
         moisture_corrected = ifelse(network != 'USCRN', moisture_corrected/100, moisture_corrected)) %>%
  dplyr::select(network, site_id, date, depth, moisture_corrected)


#NEON
#function to join depth id with data
join_soilmoisture_with_depth = function(x, depth_lookup) {
  # ------------------------------------------------------------------------------
  # Function: join_soilmoisture_with_depth
  # Author: Zachary H. Hoylman
  # Date: 6-18-2025
  #
  # Description:
  #   Joins daily soil moisture estimates with time-varying sensor depth metadata
  #   using a non-equi interval join. Ensures that each soil moisture observation
  #   is matched to the correct sensor depth based on its siteID, position, depth,
  #   and measurement date, accounting for changes in sensor deployment over time.
  #
  # Inputs:
  #   x - A data.frame or data.table containing daily soil moisture values with:
  #         - siteID
  #         - horizontalPosition
  #         - depth (as character)
  #         - date (as Date)
  #         - daily_soil_moisture
  #
  #   depth_lookup - A data.frame or data.table with sensor metadata including:
  #         - siteID
  #         - horizontalPosition
  #         - depth (as character)
  #         - startDateTime (POSIXct)
  #         - endDateTime (POSIXct or NA)
  #         - sensorDepth (numeric, meters below surface)
  #
  # Returns:
  #   A data.table with matched rows from `x` and corresponding sensorDepth
  #   based on interval-aware joining.
  # ------------------------------------------------------------------------------
  
  # Convert to data.table
  x_dt = as.data.table(x)
  dl_dt = as.data.table(depth_lookup)
  
  # Ensure proper date columns and NA handling
  dl_dt[, startDate := as.Date(startDateTime)]
  dl_dt[, endDate := as.Date(endDateTime)]
  dl_dt[is.na(endDate), endDate := as.Date("2100-01-01")]
  
  x_dt[, date := as.Date(date)]
  x_dt[, startDate := date]
  x_dt[, endDate := date]
  
  # Set keys for efficient overlap join
  setkey(dl_dt, siteID, horizontalPosition, depth, startDate, endDate)
  
  # Perform interval-based join
  result = foverlaps(
    x_dt[, .(siteID, horizontalPosition, depth, date, daily_soil_moisture, startDate, endDate)],
    dl_dt,
    by.x = c("siteID", "horizontalPosition", "depth", "startDate", "endDate"),
    type = "within",
    nomatch = NULL
  )
  
  # Return clean result
  return(result[, .(siteID, horizontalPosition, depth, date = date, daily_soil_moisture, sensorDepth)])
}

process_neon = function(path){
  # ------------------------------------------------------------------------------
  # Function: process_neon
  # Author: Zachary H. Hoylman
  # Date: 6-18-2025
  #
  # Description:
  #   Processes a NEON SWC (soil water content) CSV file to produce daily,
  #   quality-controlled soil moisture summaries joined with sensor depth metadata.
  #   Performs unit conversions and outputs a normalized, tidy dataset ready for
  #   cross-network or model comparison.
  #
  # Inputs:
  #   path - File path (character) to a NEON-formatted soil moisture CSV file
  #          containing correctedVSWC and correctedVSWCFinalQF variables for depths
  #          501–508 (though not all depths are guaranteed to be present).
  #
  # Assumptions:
  #   - A shared file "swc_depthsV2.csv" is available in:
  #       "~/nrcs-scd-soil-moisture-eval-data/raw/NEON_SWC_forZach/"
  #   - Data columns follow standard NEON variable naming conventions.
  #
  # Returns:
  #   A tibble with columns:
  #     - network (always "NEON")
  #     - site_id (formatted as "SITE_hHORIZONTALPOSITION")
  #     - date (Date)
  #     - depth (sensor depth in cm, absolute value)
  #     - moisture_corrected (daily mean VSWC, NA if QC failed)
  #
  # Dependencies:
  #   Requires: readr, dplyr, tidyr, stringr, lubridate, glue, join_soilmoisture_with_depth()
  # ------------------------------------------------------------------------------
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(glue)
  
  # Read input CSV (let readr infer types)
  x_raw = read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  # Ensure startDateTime is POSIXct (fix if read as character)
  if (inherits(x_raw$startDateTime, "character")) {
    x_raw$startDateTime = parse_date_time(
      x_raw$startDateTime,
      orders = c("Ymd HMS", "Ymd HM", "Ymd", "mdY HMS", "mdY HM", "mdY"),
      tz = "UTC"
    )
  }
  
  # Load and prepare sensor depth metadata
  depth_lookup = read_csv(
    "~/nrcs-scd-soil-moisture-eval-data/raw/NEON_SWC_forZach/swc_depthsV2.csv",
    show_col_types = FALSE, progress = FALSE
  ) |>
    mutate(
      horizontalPosition = as.numeric(horizontalPosition.HOR),
      depth = as.character(verticalPosition.VER)
    ) |>
    select(siteID, horizontalPosition, depth, startDateTime, endDateTime, sensorDepth)
  
  # Extract soil moisture values
  x_data = x_raw |>
    select(siteID, horizontalPosition, startDateTime, matches("^correctedVSWCMean_")) |>
    pivot_longer(
      cols = matches("^correctedVSWCMean_"),
      names_to = "name_soil_moisture",
      values_to = "soil_moisture"
    ) |>
    mutate(depth = str_extract(name_soil_moisture, "\\d+")) |>
    select(siteID, horizontalPosition, startDateTime, soil_moisture, depth)
  
  # Extract QC flags
  x_qc = x_raw |>
    select(siteID, horizontalPosition, startDateTime, matches("^correctedVSWCFinalQF_")) |>
    pivot_longer(
      cols = matches("^correctedVSWCFinalQF_"),
      names_to = "name_qc",
      values_to = "qc"
    ) |>
    mutate(depth = str_extract(name_qc, "\\d+")) |>
    select(siteID, horizontalPosition, startDateTime, qc, depth)
  
  # Merge values and QC, apply filtering
  x = x_data |>
    left_join(x_qc, by = c("siteID", "horizontalPosition", "startDateTime", "depth")) |>
    mutate(
      soil_moisture_qc = case_when(
        qc == 0 ~ soil_moisture,
        qc != 0 ~ NA_real_,
        TRUE    ~ soil_moisture
      )
    ) |>
    select(siteID, horizontalPosition, startDateTime, depth, soil_moisture = soil_moisture_qc) |>
    mutate(date = as_date(startDateTime)) |>
    group_by(siteID, horizontalPosition, depth, date) |>
    summarise(daily_soil_moisture = mean(soil_moisture, na.rm = TRUE), .groups = "drop") |>
    drop_na(daily_soil_moisture)
  
  # Join with sensor metadata and finalize output
  x_joined = join_soilmoisture_with_depth(x, depth_lookup) |>
    as_tibble() |>
    mutate(
      network = "NEON",
      depth = abs(sensorDepth) * 100,
      site_h_position = glue("{siteID}_h{horizontalPosition}")
    ) |>
    select(
      network,
      site_id = site_h_position,
      date,
      depth,
      moisture_corrected = daily_soil_moisture
    )
  
  return(x_joined)
}

#paths to the neon data from Ed Ayers (provided to us via NEON for authoratative data)
paths = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/NEON_SWC_forZach", full.names = T) |>
  (\(x) x[!grepl("swc_depthsV2.csv", x)])()

neon = paths |>
  purrr::map(process_neon) |>
  bind_rows()

#compute "generalized depth" data
all_data = list(mt_mesonet, ok_mesonet, uscrn_nrcs, neon) %>%
  bind_rows() %>%
  mutate(moisture_corrected = ifelse(moisture_corrected > 1, NA, moisture_corrected),
         moisture_corrected = ifelse(moisture_corrected < 0, NA, moisture_corrected))

generalized_depths = all_data %>%
  mutate(generalized_depth = ifelse(depth <= 10, 'Shallow',
                                    ifelse(depth > 10 & depth <= 40, 'Middle', 'Deep'))) %>%
  group_by(network, site_id, date, generalized_depth) %>%
  summarise(soil_moisture = median(moisture_corrected, na.rm = T))

full_profile = generalized_depths %>%
  group_by(network, site_id, date) %>%
  summarise(n = sum(!is.na(soil_moisture)),
            soil_moisture = median(soil_moisture, na.rm = T),
            soil_moisture = ifelse(n == 3, soil_moisture, NA)) %>%
  mutate(generalized_depth = 'Depth Averaged') %>%
  dplyr::select(-n)

final_data = bind_rows(generalized_depths, full_profile) %>%
  arrange(network, site_id, date, generalized_depth) %>%
  mutate(soil_moisture = ifelse(soil_moisture > 0.7, NA, soil_moisture))
  
#write out all data
write_csv(final_data, "~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv")

#plot of stations

all_sites_with_data = final_data |>
  ungroup() |>
  distinct(network, site_id) |>
  mutate(site_id = ifelse(network == 'NEON', substr(site_id,1,4), site_id)) |>
  distinct(network, site_id)


us_sf = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') 

station_meta_conus_with_data = read_csv('~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-raw.csv') %>% 
  filter(site_id %in% all_sites_with_data$site_id)

station_meta_conus_with_data_plot = station_meta_conus_with_data %>%
  left_join(., station_meta_conus_with_data %>% 
              group_by(network) %>%
              summarise(n = length(network)))

#make CONUS domain map with data final
conus_domain_with_data <- ggplot() +
  geom_sf(
    data = us_sf %>%
      filter(NAME %notin% c("Virgin Islands", "Alaska", "Hawaii", "Puerto Rico")) %>%
      st_transform(5070),
    fill = "gray95", color = "gray60"
  ) +
  geom_sf(
    data = station_meta_conus_with_data_plot %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_transform(5070) %>%
      mutate(network_label = glue::glue("{network} (n = {n})")),
    aes(fill = network_label),
    shape = 21, color = "black", size = 2, stroke = 0.3
  ) +
  theme_bw(base_size = 16) +
  ggtitle(
    "Soil Moisture Observations Across the U.S.",
    subtitle = glue::glue("(n sites = {nrow(station_meta_conus_with_data_plot)})")
  ) +
  labs(fill = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm"),
    legend.text = element_text(size = 12),
    strip.background = element_rect(colour = "transparent", fill = "transparent")
  )

write_csv(station_meta_conus_with_data, '~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv')
ggsave(conus_domain_with_data, file = '~/nrcs-scd-soil-moisture-eval/figs/conus_domain_with_data_final.png')