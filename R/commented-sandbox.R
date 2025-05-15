# Load required libraries
library(tidyverse)      # Core tidy data manipulation and visualization functions
library(sf)             # Tools for working with spatial data (simple features)
library(magrittr)       # Enhanced piping functionality
library(neonSoilFlux)   # Functions for handling NEON soil flux data
library(furrr)          # Parallel processing with futures
library(lubridate)      # Tools for working with dates and times
library(progressr)      # Progress reporting for long-running operations
library(soilDB)         # Soil data retrieval from NRCS databases

# Define a custom negation operator for filtering
`%notin%` = Negate(`%in%`)

# Define a function to read Oklahoma Mesonet data
read_ok_mesonet = function(x) {
  # Read data from an Excel file and extract the site_id from the file path
  temp = readxl::read_excel(x) %>%
    mutate(site_id = str_extract(x, "(?<=Mesonet/)(.*)(?=.xls)")) %>%
    select('site_id', everything())  # Reorganize columns to place 'site_id' first
  return(temp)  # Return the processed data
}

# Load station metadata
# SCAN and SNOTEL metadata from the soilDB package
data('SCAN_SNOTEL_metadata', package = 'soilDB')

# Read US state boundaries as a spatial object
us_sf = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json')

# Retrieve metadata for Montana Mesonet stations
mt_mesonet_meta = read_csv('https://mesonet.climate.umt.edu/api/v2/stations/?public=true&active=true&type=csv') %>%
  mutate(network = 'MT Mesonet', site_id = station) %>%
  dplyr::select(c('network', 'site_id', 'latitude', 'longitude'))

# Combine metadata from Montana Mesonet and other station sources
station_meta_raw = bind_rows(
  mt_mesonet_meta,
  readxl::read_excel("~/nrcs-scd-soil-moisture-eval-data/raw/NRCS_SCD_Project_In_situ_station_meta_data.xlsx",
                     sheet = "combined") %>%
    mutate(
      network = gsub("'", "", Network),    # Clean up network names
      site_id = gsub("'", "", `Site ID`)  # Clean up site IDs
    ) %>%
    dplyr::select(c('network', 'site_id', latitude = 'Latitude', longitude = 'Longitude')) %>%
    filter(network %notin% c('MT Mesonet'))  # Exclude Montana Mesonet from this set
)

# Process metadata for NRCS stations (SCAN and SNOTEL)
nrcs_meta = SCAN_SNOTEL_metadata %>%
  as_tibble() %>%
  mutate(
    Network = ifelse(Network == 'SNTL', 'SNOTEL', Network),  # Standardize network names
    site_network = glue::glue('{Network}:{Site}')           # Create composite identifier
  ) %>%
  filter(site_network %in% glue::glue('{station_meta_raw$network}:{station_meta_raw$site_id}')) %>%
  dplyr::select(
    network = Network, 
    site_id = Site,
    latitude = Latitude,
    longitude = Longitude
  ) %>%
  mutate(
    network = ifelse(network == 'SNOTEL', 'SNTL', network), # Adjust SNOTEL naming
    site_id = glue::glue('{network}:{site_id}')             # Update site IDs
  )

# Load USCRN lookup table for station identifiers
uscrn_lookup = read_csv('~/nrcs-scd-soil-moisture-eval-data/raw/uscrn_name_conversion.csv')

# Update metadata to include processed NRCS data and standardize network names
station_meta = station_meta_raw %>%
  filter(network %notin% c('SCAN', 'SNOTEL')) %>%
  bind_rows(., nrcs_meta) %>%
  mutate(network = ifelse(network == 'neon', 'NEON', network))  # Standardize "neon" to "NEON"

# Add summary counts for stations by network
station_meta_plot = station_meta %>%
  left_join(
    station_meta %>% 
      group_by(network) %>%
      summarise(n = length(network))  # Count stations per network
  )

# Create a full domain map showing station locations
full_domain = ggplot() +
  geom_sf(data = us_sf) +  # Add US state boundaries
  geom_point(
    data = station_meta_plot, 
    aes(
      x = longitude, 
      y = latitude, 
      fill = glue::glue('{station_meta_plot$network} (n = {station_meta_plot$n})')
    ), 
    color = 'black', 
    shape = 21
  ) +
  xlim(c(-170, -53)) +  # Set longitude bounds
  theme_bw() +
  ggtitle(glue::glue('Full Domain (n sites = {length(station_meta_plot$network)})')) +
  labs(fill = NULL) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the full domain map to a file
ggsave(full_domain, file = '~/nrcs-scd-soil-moisture-eval/figs/full_domain.png')

# Filter metadata to include only stations within the CONUS region
station_meta_conus = st_as_sf(station_meta, coords = c('longitude', 'latitude')) %>%
  st_set_crs(st_crs('EPSG:4326')) %>%  # Set CRS to WGS84
  st_intersection(
    us_sf %>%
      filter(NAME %notin% c('Virgin Islands', 'Alaska', 'Hawaii', 'Puerto Rico')) %>%  # Exclude non-CONUS regions
      st_make_valid()
  ) %>%
  mutate(
    longitude = st_coordinates(.)[, 1],  # Extract longitude
    latitude = st_coordinates(.)[, 2]   # Extract latitude
  ) %>%
  dplyr::select(c('network', 'site_id', 'latitude', 'longitude')) %>%
  st_drop_geometry()

# Save CONUS station metadata to a CSV
write_csv(station_meta_conus, '~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-raw.csv')

# Add summary counts for CONUS stations by network
station_meta_conus_plot = station_meta_conus %>%
  left_join(
    station_meta_conus %>%
      group_by(network) %>%
      summarise(n = length(network))  # Count stations per network
  )

# Create a CONUS map showing station locations
conus_domain = ggplot() +
  geom_sf(
    data = us_sf %>%
      filter(NAME %notin% c('Virgin Islands', 'Alaska', 'Hawaii', 'Puerto Rico'))
  ) +
  geom_point(
    data = station_meta_conus_plot, 
    aes(
      x = longitude, 
      y = latitude, 
      fill = glue::glue('{station_meta_conus_plot$network} (n = {station_meta_conus_plot$n})')
    ), 
    color = 'black', 
    shape = 21
  ) +
  theme_bw() +
  ggtitle(glue::glue('CONUS Domain (n sites = {length(station_meta_conus_plot$network)})')) +
  labs(fill = NULL) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the CONUS domain map to a file
ggsave(conus_domain, file = '~/nrcs-scd-soil-moisture-eval/figs/conus_domain.png')

##### Preprocess Soil Moisture Data #####

# Load MT Mesonet data
# Query soil volumetric water content (VWC) and soil temperature (Temp) data.
# Filter data to exclude inaccurate sensor readings where soil temp ≤ 0ºC.
mt_mesonet_vwc_temp_data <- httr::GET("https://mesonet.climate.umt.edu/api/v2/observations/daily/",
                                      query = list(
                                        premade = TRUE,       # Premade query
                                        rm_na = TRUE,         # Include data even if missing values
                                        type = "csv",         # Data format: CSV
                                        units = "si",         # Use SI units
                                        end_time = "2099-01-01", # Generic end date
                                        start_time = "1900-01-01", # Generic early start date to capture all records
                                        elements = "soil_vwc, soil_temp", # Variables: soil moisture and temperature
                                        stations = station_meta_conus %>%
                                          filter(network == 'MT Mesonet') %$%
                                          site_id %>%
                                          paste(., collapse = ", ") # Filter by site IDs
                                      )
) %>% 
  httr::content() 

# Process MT Mesonet data
mt_mesonet_vwc_temp_data_final <- mt_mesonet_vwc_temp_data %>% 
  mutate(date = as.Date(datetime)) %>%      # Extract date
  rowwise() %>%
  filter(any(!is.na(c_across(-c(station, datetime, date))))) %>% # Remove rows with all NA values
  dplyr::select(-c('datetime')) %>%
  dplyr::select(c('site_id' = 'station', 'date', everything())) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name) %>% abs,      # Parse and take absolute depth
         var = stringr::word(name, 2),           # Extract variable name
         var = ifelse(var == 'Temperature', 'soil_temperature', 'soil_moisture'),
         var_depth = glue::glue('{var}_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth)

write_csv(mt_mesonet_vwc_temp_data_final, '~/nrcs-scd-soil-moisture-eval-data/processed/mt-mesonet-raw.csv')

# Load USCRN, SNOTEL, and SCAN data
# Combine data, process depths, and standardize site IDs
uscrn_sntl_scan <- read_csv("~/nrcs-scd-soil-moisture-eval-data/raw/data from Hoylman/in_situ_soil_moisture_data Zachary Hoylman.csv") %>%
  filter(site_id %notin% unique(c(mt_mesonet_vwc_temp_data_final$site_id, "sidneymt", "wrsround", 'conradmt'))) %>%
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name) %>% abs) %>%
  mutate(depth = depth * 2.54 %>% round(1),     # Convert depth to cm
         depth = ifelse(depth == 90, 91, depth),
         varname = str_extract(name, "^[^_]+_[^_]+"), # Extract variable name
         var_depth = glue::glue('{varname}_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth) %>%
  left_join(., uscrn_lookup %>% 
              dplyr::select(site_id = WBAN, 
                            new_name = STATION_ID)) %>%
  mutate(site_id = ifelse(is.na(new_name), site_id, new_name)) %>%
  dplyr::select(-c(new_name))

write_csv(uscrn_sntl_scan, '~/nrcs-scd-soil-moisture-eval-data/processed/uscrn-sntl-scan-raw.csv')

# Load Oklahoma Mesonet data
# Process soil moisture at multiple depths
ok_mesonet <- list.files("~/nrcs-scd-soil-moisture-eval-data/raw/Soil Moisture Data Oklahoma Mesonet", full.names = T) %>%
  lapply(., read_ok_mesonet) %>%
  bind_rows() %>%
  mutate(date = as.Date('0000-01-01') + Date - 1) %>% # Convert numeric date to Date object
  dplyr::select(site_id, date, VWC05, VWC25, VWC60) %>% # Select VWC at specific depths
  pivot_longer(-c(site_id, date)) %>%
  mutate(depth = parse_number(name),
         var_depth = glue::glue('soil_moisture_{depth}cm')) %>%
  dplyr::select(site_id, date, value, var_depth) %>%
  pivot_wider(values_from = value, names_from = var_depth)

write_csv(ok_mesonet, '~/nrcs-scd-soil-moisture-eval-data/processed/ok-mesonet-raw.csv')

# Load and preprocess NEON data
sensor_depth_lookup <- read_csv("~/nrcs-scd-soil-moisture-eval-data/raw/Soil Moisture Data Neon/neon/swc_depthsV2.csv") %>%
  mutate(verticalPosition.VER = as.character(verticalPosition.VER))

neon_sites <- station_meta_conus |>
  filter(network == 'neon')

# Define the start and end dates for NEON data processing
start_date <- ymd("2016-04-01")
end_date <- ymd("2024-10-01")

# Generate a sequence of monthly intervals
month_year_vector <- seq(start_date, end_date, by = "month") %>%
  format(., "%Y-%m")

# Parallel processing setup for soil moisture
plan('multisession', workers = 8)

# Initialize progress handlers for processing
handlers(global = TRUE)