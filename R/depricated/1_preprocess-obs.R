library(tidyverse)
library(sf)
library(magrittr)
library(neonSoilFlux)
library(furrr)
library(lubridate)
library(progressr)
library(soilDB)

`%notin%` = Negate(`%in%`)

read_ok_mesonet = function(x) {
  
  raw = readxl::read_excel(
    x,
    col_types = "text",                       # <-- do NOT let readxl guess
    na = c("", "NA", "NaN", "-9999")
  )
  
  out = raw |>
    dplyr::mutate(
      # site id from path
      site_id = stringr::str_extract(x, "(?<=Mesonet/)(.*)(?=\\.xls$)")
    ) |>
    # convert Date: Excel serial -> Date (keeps non-serials as NA)
    dplyr::mutate(
      Date = suppressWarnings(as.numeric(.data$Date))
    ) |>
    # everything else to numeric; turn TRUE/FALSE/T/F etc. into NA
    dplyr::mutate(
      dplyr::across(
        .cols = -dplyr::all_of(c("Date", "site_id")),
        .fns  = ~ readr::parse_double(
          .x,
          na = c("", "NA", "NaN", "-9999", "TRUE", "FALSE", "T", "F", "—", "–")
        )
      )
    ) |>
    dplyr::select(site_id, dplyr::everything())
  
  return(out)
}

data('SCAN_SNOTEL_metadata', package = 'soilDB')

us_sf = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') 

mt_mesonet_meta = read_csv("~/nrcs-scd-soil-moisture-eval-data/raw/MT_mesonet_precise_station_locations_filtered.csv") %>%
  mutate(network = 'MT Mesonet',
         site_id = station) %>%
  dplyr::select(c('network', 'site_id', 'latitude', 'longitude'))

station_meta_raw = bind_rows(mt_mesonet_meta, 
                         readxl::read_excel("~/nrcs-scd-soil-moisture-eval-data/raw/NRCS_SCD_Project_In_situ_station_meta_data.xlsx",
                                            sheet = "combined") %>%
                           mutate(network = cleaned_text <- gsub("'", "", Network),
                                  site_id = cleaned_text <- gsub("'", "", `Site ID`)) %>%
                           dplyr::select(c('network', 'site_id', latitude = 'Latitude', longitude = 'Longitude')) %>%
                           filter(network %notin% c('MT Mesonet')))

nrcs_meta = SCAN_SNOTEL_metadata %>% 
  as_tibble() %>%
  mutate(
    Network = ifelse(Network == 'SNTL', 'SNOTEL', Network),
    site_network = glue::glue('{Network}:{Site}')
    ) %>%
  filter(site_network %in% glue::glue('{station_meta_raw$network}:{station_meta_raw$site_id}')) %>%
  dplyr::select(
    network = Network, 
    site_id = Site,
    latitude = Latitude,
    longitude = Longitude) %>%
  mutate(
    network = ifelse(network == 'SNOTEL', 'SNTL', network),
    site_id = glue::glue('{network}:{site_id}')
  )

#lookup for differnet USCRN ids 
uscrn_lookup = read_csv('~/nrcs-scd-soil-moisture-eval-data/raw/uscrn_name_conversion.csv')

#update NRCS SNOTEL and SCAN metadata
station_meta = station_meta_raw %>%
  filter(network %notin% c('SCAN', 'SNOTEL')) %>%
  bind_rows(., nrcs_meta) %>%
  mutate(network = ifelse(network == 'neon', 'NEON', network))

station_meta_plot = station_meta %>%
  left_join(., station_meta %>% 
              group_by(network) %>%
              summarise(n = length(network)))

#make full domain map
full_domain = ggplot()+
  geom_sf(data = us_sf) +
  geom_point(data = station_meta_plot, aes(x = longitude, y = latitude, fill = glue::glue('{station_meta_plot$network} (n = {station_meta_plot$n})')), color = 'black',  shape = 21) +
  xlim(c(-170, -53)) +
  theme_bw() +
  ggtitle(glue::glue('Full Domain (n sites = {length(station_meta_plot$network)})'))+
  labs(fill = NULL)+
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave(full_domain, file = '~/nrcs-scd-soil-moisture-eval/figs/full_domain.png')

#crop to CONUS
station_meta_conus = st_as_sf(station_meta, coords = c('longitude', 'latitude')) %>%
  st_set_crs(st_crs('EPSG:4326')) %>%
  st_intersection(., 
                  us_sf %>% 
                    filter(NAME %notin% c('Vigin Islands', 'Alaska', 'Hawaii', 'Puerto Rico')) %>%
                    st_make_valid()
                  ) %>%
  mutate(longitude = st_coordinates(.)[,1],
         latitude = st_coordinates(.)[,2]) %>%
  dplyr::select(c('network', 'site_id', 'latitude', 'longitude')) %>%
  st_drop_geometry()

write_csv(station_meta_conus, '~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-raw.csv')

station_meta_conus_plot = station_meta_conus %>%
  left_join(., station_meta_conus %>% 
              group_by(network) %>%
              summarise(n = length(network)))

#make CONUS domain map
conus_domain = ggplot()+
  geom_sf(data =  us_sf %>% 
            filter(NAME %notin% c('Vigin Islands', 'Alaska', 'Hawaii', 'Puerto Rico'))
          ) +
  geom_point(data = station_meta_conus_plot, aes(x = longitude, y = latitude, fill = glue::glue('{station_meta_conus_plot$network} (n = {station_meta_conus_plot$n})')), color = 'black',  shape = 21) +
  theme_bw() +
  ggtitle(glue::glue('CONUS Domain (n sites = {length(station_meta_conus_plot$network)})'))+
  labs(fill = NULL)+
  theme(plot.title = element_text(hjust = 0.5))


ggsave(conus_domain, file = '~/nrcs-scd-soil-moisture-eval/figs/conus_domain.png')

##### preproces soil moisture data ######

# read in data MT MESONET
mt_mesonet_vwc_temp_data = httr::GET("https://mesonet.climate.umt.edu/api/v2/observations/daily/",
          query = list(
            premade = TRUE,
            rm_na = TRUE, #includes data even if there is missing values
            type = "csv",
            units = "si",
            end_time = "2099-01-01",
            start_time = "1900-01-01", #setting generic (way early) year to capture all start dates
            #retrieving both because of need to filter VWC by soil temp to exclude inaccuracies in sensor readings where soil temp ≤ 0ºC
            elements = "soil_vwc, soil_temp",
            stations = station_meta_conus %>% 
              filter(network == 'MT Mesonet') %$% 
              site_id %>% 
              paste(., collapse = ", ")
            )
          ) %>% 
  httr::content() 

mt_mesonet_vwc_temp_data_final = mt_mesonet_vwc_temp_data %>% 
  mutate(date = as.Date(datetime)) %>%
  rowwise() %>%
  filter(any(!is.na(c_across(-c(station, datetime, date))))) %>%
  dplyr::select(-c('datetime')) %>%
  dplyr::select(c('site_id' = 'station', 'date', everything())) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name) %>%
           abs,
         var = stringr::word(name, 2),
         var = ifelse(var == 'Temperature', 'soil_temperature', 'soil_moisture'),
         var_depth = glue::glue('{var}_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth)
  
write_csv(mt_mesonet_vwc_temp_data_final, '~/nrcs-scd-soil-moisture-eval-data/processed/mt-mesonet-raw.csv')
# read in data USCRN SNOTEL and SCAN

uscrn_sntl_scan = read_csv("~/nrcs-scd-soil-moisture-eval-data/raw/data from Hoylman/in_situ_soil_moisture_data Zachary Hoylman.csv") %>%
  filter(site_id %notin% unique(c(mt_mesonet_vwc_temp_data_final$site_id, "sidneymt", "wrsround", 'conradmt'))) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name) %>%
           abs) %>%
  mutate(depth = depth * 2.54 %>%
           round(1),
         depth = ifelse(depth == 90, 91, depth),
         varname = str_extract(name, "^[^_]+_[^_]+"),
         var_depth = glue::glue('{varname}_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth) %>%
  left_join(., uscrn_lookup %>% 
              dplyr::select(site_id = WBAN, 
                            new_name = STATION_ID)
            ) %>%
  mutate(site_id = ifelse(is.na(new_name), site_id, new_name)) %>%
  dplyr::select(-c(new_name))

write_csv(uscrn_sntl_scan, '~/nrcs-scd-soil-moisture-eval-data/processed/uscrn-sntl-scan-raw.csv')

# read in data OK MESONET - reading in as logical ()
ok_mesonet = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/Soil Moisture Data Oklahoma Mesonet", full.names = T) %>%
  lapply(., read_ok_mesonet) %>%
  bind_rows() %>%
  mutate(date = as.Date('0000-01-01') + Date-1,) %>%
  dplyr::select(site_id, date, VWC05, VWC25, VWC60) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         var_depth = glue::glue('soil_moisture_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth)

write_csv(ok_mesonet, '~/nrcs-scd-soil-moisture-eval-data/processed/ok-mesonet-raw.csv')

