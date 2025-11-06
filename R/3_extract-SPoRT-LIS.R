##############################################################
# Title: Extract SPoRT Time Series for ROIs
# Author: Dr. Zachary H. Hoylman
# Date: 11-6-2025
##############################################################

## MUST BE RUN ON MCO SERVER! DATA IS VERY LARGE!

library(terra)
library(exactextractr)
library(tidyverse)
library(sf)
library(tictoc)
library(glue)

# ---- Import ROIs ----
roi = read_csv("~/Desktop/station-meta-conus-w-data-final.csv") %>%
  st_as_sf(coords = c('longitude', 'latitude')) %>%
  st_set_crs('EPSG:4326') |>
  mutate(site_id_unique = paste0(site_id, '_', network))

# ---- NEW: Raw SPoRT Extraction Function (WIDE format, matching GridMET) ----
extract_sport_raw = function(roi, site_id_col, tag) {
  message(glue::glue("[{tag}] Extracting raw SPoRT time series (wide format)"))
  
  # File paths and corresponding depths
  files = list(
    # "0-10cm"   = "/data/ssd2/soil-moisture-models/nc/SPoRT_mean_soil_moisture_0-10cm.nc"
    # "10-40cm"  = "/data/ssd2/soil-moisture-models/nc/SPoRT_mean_soil_moisture_10-40cm.nc"
    "40-100cm" = "/data/ssd2/soil-moisture-models/nc/SPoRT_mean_soil_moisture_40-100cm.nc"
    #"0-100cm" = "/data/ssd2/soil-moisture-models/nc/SPoRT_mean_soil_moisture_0-100cm.nc"
  )
  
  # Loop over files and extract rasters
  results = lapply(names(files), function(depth_label) {
    sport_rast = terra::rast(files[[depth_label]])
    time_vals = terra::time(sport_rast)
    #time_vals = read_csv("/data/ssd2/soil-moisture-models/nc/SPoRT_mean_soil_moisture_0-100cm_time.csv")$time
    
    tictoc::tic(glue::glue("Extracting depth {depth_label}"))
    raw = terra::extract(sport_rast, terra::vect(roi)) %>%
      t() %>%
      tibble::as_tibble()
    tictoc::toc()
    
    colnames(raw) = roi[[site_id_col]]
    
    # Build wide tibble: var, time, site_1, site_2, ...
    wide = raw[-1, ] %>%
      dplyr::mutate(
        time = time_vals,
        var = paste0("SPoRT_raw_", depth_label)
      ) %>%
      dplyr::select(var, time, everything())
    
    return(wide)
  })
  
  # Combine all variable-depths into one tibble
  sport_wide = dplyr::bind_rows(results)
  return(sport_wide)
}

# ---- Run for Observational ROI ----
sport_raw_operational = extract_sport_raw(roi, "site_id_unique", "40-100cm")

write_csv(sport_raw_operational, "~/temp/SPoRT-data-40-100cm.csv")
