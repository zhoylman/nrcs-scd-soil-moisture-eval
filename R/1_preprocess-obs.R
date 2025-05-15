library(tidyverse)
library(sf)
library(magrittr)
library(neonSoilFlux)
library(furrr)
library(lubridate)
library(progressr)
library(soilDB)

`%notin%` = Negate(`%in%`)

read_ok_mesonet = function(x){
  temp = readxl::read_excel(x) %>%
    mutate(site_id = str_extract(x, "(?<=Mesonet/)(.*)(?=.xls)")) %>%
    select('site_id', everything())
  
  return(temp)
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

# read in data OK MESONET
ok_mesonet = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/Soil Moisture Data Oklahoma Mesonet", full.names = T) %>%
  lapply(., read_ok_mesonet) %>%
  bind_rows() %>%
  mutate(date = as.Date('0000-01-01') + Date-1) %>%
  dplyr::select(site_id, date, VWC05, VWC25, VWC60) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         var_depth = glue::glue('soil_moisture_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth)

write_csv(ok_mesonet, '~/nrcs-scd-soil-moisture-eval-data/processed/ok-mesonet-raw.csv')

#neon!
sensor_depth_lookup = read_csv("~/nrcs-scd-soil-moisture-eval-data/raw/Soil Moisture Data Neon/neon/swc_depthsV2.csv") %>%
  mutate(verticalPosition.VER = as.character(verticalPosition.VER))

neon_sites = station_meta_conus |>
  filter(network == 'neon')

# Define the start and end dates
start_date <- ymd("2016-04-01")
end_date <- ymd("2024-10-01")  # Get the current month

# Generate a sequence of dates from start to end by month
month_year_vector <- seq(start_date, end_date, by = "month") %>%
  format(., "%Y-%m")

# soil moisture!

plan('multisession', workers = 8)

# Initialize a progress handler
handlers(global = TRUE)

# Set up progress bar
with_progress({
  # Define the total steps in the progress bar
  p <- progressor(along = month_year_vector)
  
  month_year_vector %>%
    purrr::map(function(.t){
      p(sprintf("Processing month-year: %s", .t))  # Update progress bar
      
      neon_sites %$%
        site_id %>%
        furrr::future_map(function(.x){
          tryCatch({
            neon = neonUtilities::loadByProduct( 
              dpID = "DP1.00094.001", # "DP1.00094.001" = soil moisture, "DP1.00041.001" = soil temperature
              site = .x, 
              startdate = .t,
              enddate = .t, 
              timeIndex = "30", 
              package = "expanded", 
              check.size = 8, 
              include.provisional = TRUE,
              nCores = 8
            )
            
            positions = neon$sensor_positions_00094
            
            # Then correct the swc
            site_swc <- swc_correct(neon, "WREF","2022-06")
            # Select the columns and the time frequency
            time_frequency <- "30_minute"
            column_selectors <- c("Mean")
            swc <- site_swc |>
              purrr::pluck(paste0("SWS_",time_frequency)) |>
              dplyr::select(tidyselect::all_of(c("domainID","siteID",
                                                 "horizontalPosition","verticalPosition","startDateTime","VSWCFinalQF")),
                            tidyselect::matches(stringr::str_c("VSWC",column_selectors)))
            # Determine a data frame of the different horizontal and vertical positions
            swc_positions <- site_swc |>
              purrr::pluck(paste0("sensor_positions_","00094"))
            
            # Add on the positions for swc
            swc <- determine_position(swc_positions,swc) 
            
            swc_min_date = swc %>% 
              group_by(verticalPosition) %>%
              summarise(minStartDateTime = min(startDateTime))
            
            average_depths = sensor_depth_lookup %>%
              group_by(siteID, 
                       horizontalPosition.HOR, 
                       verticalPosition.VER) %>%
              summarise(median_sensor_depth = median(sensorDepth),
                        iqr_sensor_depth = IQR(sensorDepth, na.rm = T))
            
            final = swc |>
              dplyr::select(-c(zOffset)) %>%
              left_join(average_depths, by = c('siteID',
                                               'verticalPosition' = 'verticalPosition.VER',
                                               'horizontalPosition' = 'horizontalPosition.HOR')) %>%
              mutate(label = paste0('Horizontal Pos. (',horizontalPosition,')'),
                     HOR.VER = paste0(horizontalPosition, '.', verticalPosition)) %>% 
              left_join(., positions %>%
                          select(zOffset, siteID, HOR.VER))
            
            # fac_depth = ggplot(data = final, aes(x = startDateTime, y = VSWCMean, color = verticalPosition))+
            #   geom_line() +
            #   facet_wrap(~label, nrow = 3)+
            #   labs(color = 'Sensor Depth') +
            #   theme_bw(base_size = 16)+
            #   theme(strip.background = element_rect(colour="black", fill='transparent'))
            # 
            abs_depth = ggplot(data = final, aes(x = startDateTime, y = VSWCMean, color = zOffset |> as_factor()))+
              geom_line() +
              facet_wrap(~label, nrow = 3)+
              labs(color = 'Sensor Depth') +
              theme_bw(base_size = 16) +
              theme(strip.background = element_rect(colour="black", fill="transparent"))
            
            write_csv(final %>%
                        dplyr::select(-label), glue::glue("~/nrcs-scd-soil-moisture-eval-data/raw/raw_neon_data/soil_moisture/{.x}_{.t}.csv"))
            # ggsave(fac_depth, file = glue::glue('~/nrcs-scd-soil-moisture-eval-data/figs/neon_figs/factor_depth/{.x}_{.t}_factored_plot.png'), width = 11, height = 8)
            ggsave(abs_depth, file = glue::glue('~/nrcs-scd-soil-moisture-eval-data/figs/neon_figs/absolute_depth/{.x}_{.t}_absolute_depth_plot.png'), width = 11, height = 8)
          }, error = function(e){
            
          })
          
        }
        )
    })
})

plan(sequential)


# temperature!
# Define the start and end dates
start_date <- ymd("2016-04-01")
end_date <- ymd("2024-10-01")  # Get the current month

# Generate a sequence of dates from start to end by month
month_year_vector <- seq(start_date, end_date, by = "month") %>%
  format(., "%Y-%m")

plan('multisession', workers = 8)

# Initialize a progress handler
handlers(global = TRUE)

# Set up progress bar
with_progress({
  # Define the total steps in the progress bar
  p <- progressor(along = month_year_vector)
  
  month_year_vector %>%
    purrr::map(function(.t){
      p(sprintf("Processing month-year: %s", .t))  # Update progress bar
      
      neon_sites %$%
        site_id %>%
        furrr::future_map(function(.x){
          tryCatch({
            neon = neonUtilities::loadByProduct( 
              dpID = "DP1.00041.001", # "DP1.00094.001" = soil moisture, "DP1.00041.001" = soil temperature
              site = .x, 
              startdate = .t,
              enddate = .t, 
              timeIndex = "30", 
              package = "expanded", 
              check.size = 8, 
              include.provisional = TRUE,
              nCores = 8
            )
            
            soil_temp_position = neon$sensor_positions_00041 %>%
              as_tibble()
            
            soil_temp = neon$ST_30_minute %>%
              as_tibble()
            
            average_depths = sensor_depth_lookup %>%
              group_by(siteID, 
                       horizontalPosition.HOR, 
                       verticalPosition.VER) %>%
              summarise(median_sensor_depth = median(sensorDepth),
                        iqr_sensor_depth = IQR(sensorDepth, na.rm = T))
            
            final = soil_temp |>
              left_join(average_depths, by = c('siteID',
                                               'verticalPosition' = 'verticalPosition.VER',
                                               'horizontalPosition' = 'horizontalPosition.HOR')) %>%
              mutate(
                     HOR.VER = paste0(horizontalPosition, '.', verticalPosition)
                     ) %>% 
              left_join(soil_temp_position %>%
                          select(zOffset, siteID, HOR.VER)) %>%
              dplyr::select(domainID, siteID, horizontalPosition, verticalPosition, startDateTime,
                            soilTempMean, soilTempMinimum, soilTempMaximum,
                            median_sensor_depth, iqr_sensor_depth,   
                            HOR.VER, zOffset)
            
            write_csv(final, glue::glue("~/nrcs-scd-soil-moisture-eval-data/raw/raw_neon_data/soil_temperature/{.x}_{.t}.csv"))
          }, error = function(e){
            
          })
          
        }
        )
    })
})

plan(sequential)

rm(list = ls()); gc()


read_neon_to_daily = function(.f, var_id = 'soil_moisture'){
  options(dplyr.summarise.inform = FALSE)
  data = read_csv(.f, show_col_types = FALSE) %>%
    mutate(date = as.Date(startDateTime)) %>%
    group_by(siteID, horizontalPosition, date, zOffset) %>% 
    summarise(soil_moisture = median(VSWCMean, na.rm = T)) %>%
    ungroup() %>%
    mutate(var = var_id,
           depth = abs(zOffset) * 100) %>%
    dplyr::select(site_id = siteID, horizontalPosition, date, depth, soil_moisture)
}

read_neon_to_daily_temp = function(.f, var_id = 'soil_temperature'){
  options(dplyr.summarise.inform = FALSE)
  data = read_csv(.f, show_col_types = FALSE) %>%
    mutate(date = as.Date(startDateTime)) %>%
    group_by(siteID, horizontalPosition, date, zOffset) %>% 
    summarise(soil_temperature = median(soilTempMean, na.rm = T)) %>%
    ungroup() %>%
    mutate(var = var_id,
           depth = abs(zOffset) * 100) %>%
    dplyr::select(site_id = siteID, horizontalPosition, date, depth, soil_temperature)
}

plan('multisession', workers = 8)
tictoc::tic()
neon_soil_moisture = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/raw_neon_data/soil_moisture", full.names = T) %>%
  furrr::future_map(., read_neon_to_daily) %>%
  bind_rows()

neon_soil_temperature = list.files("~/nrcs-scd-soil-moisture-eval-data/raw/raw_neon_data/soil_temperature", full.names = T) %>%
  furrr::future_map(., read_neon_to_daily_temp) %>%
  bind_rows()
tictoc::toc()
plan(sequential)

neon_merged = left_join(neon_soil_moisture, neon_soil_temperature)

write_csv(neon_merged, '~/nrcs-scd-soil-moisture-eval-data/processed/neon-raw.csv')

# if you dont want to re-run all the pre-processing:
# mt_mesonet_vwc_temp_data_final = read_csv('~/nrcs-scd-soil-moisture-eval-data/processed/mt-mesonet-raw.csv')
# uscrn_sntl_scan = read_csv('~/nrcs-scd-soil-moisture-eval-data/processed/uscrn-sntl-scan-raw.csv')
# ok_mesonet = read_csv('~/nrcs-scd-soil-moisture-eval-data/processed/ok-mesonet-raw.csv')
# neon_merged = read_csv('~/nrcs-scd-soil-moisture-eval-data/processed/neon-raw.csv')

all_sites_with_data = c(
  unique(neon_merged$site_id),
  unique(mt_mesonet_vwc_temp_data_final$site_id),
  unique(uscrn_sntl_scan$site_id),
  unique(ok_mesonet$site_id)
)

station_meta_conus_with_data = station_meta_conus %>% 
  filter(site_id %in% all_sites_with_data)

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
    "Sites With Data",
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
