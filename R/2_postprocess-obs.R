library(tidyverse)

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

#still need to deal with neon data
neon = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/neon-raw.csv")

#compute "generalized depth" data
all_data = list(mt_mesonet, ok_mesonet, uscrn_nrcs) %>%
  bind_rows() %>%
  mutate(moisture_corrected = ifelse(moisture_corrected > 1, NA, moisture_corrected),
         moisture_corrected = ifelse(moisture_corrected < 0, NA, moisture_corrected))

generalized_depths = all_data %>%
  mutate(generalized_depth = ifelse(depth <= 10, 'Shallow',
                                    ifelse(depth > 10 & depth <= 50, 'Middle', 'Deep'))) %>%
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
  arrange(network, site_id, date, generalized_depth)

#write out all data
write_csv(final_data, "~/nrcs-scd-soil-moisture-eval-data/processed/final-soil-moisture-data-generalized.csv")
