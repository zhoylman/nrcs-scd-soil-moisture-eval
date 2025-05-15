#load required libraries
library(reticulate) # allows for python interfacing
library(rgee) # R wrapper for the python GEE library
library(sf) # simple feature library - used for vectors
library(tidyverse) # package for tidy syntax etc
library(geojsonio) # package to send ROI SF objects to GEE

#set up the gee environment
use_condaenv("gee-base", conda = "auto",required = TRUE)
ee = import("ee")
ee_Initialize(drive = TRUE)

stations = read_csv("~/nrcs-scd-soil-moisture-eval-data/processed/station-meta-conus-w-data-final.csv") %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs('EPSG:4326')

stations_ee = stations %>%
  sf_as_ee()

#soils data from https://gee-community-catalog.org/projects/soilprop/
#units https://casoilresource.lawr.ucdavis.edu/soil-properties/download.php
# soils = percent by weight	
sand = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/physical/sand') 
silt = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/physical/silt') 
clay = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/physical/clay') 
som = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/chemical/som_max')
bulk_density = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/physical/bulk_density')
awc = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/physical/water_storage')
depth_to_restrictive = ee$Image('projects/earthengine-legacy/assets/projects/sat-io/open-datasets/CSRL_soil_properties/land_use/resdept') 
topofire_aws = ee$Image('users/zhoylman/aws_conus_240m_filled_rf_full')$rename('topofire_aws')
topofire_poro = ee$Image('projects/ee-zhoylman/assets/porosity_gNATSGO_US')$rename('topofire_porosity')


soils =  ee$ImageCollection(c(sand, silt, clay, som, bulk_density, awc, depth_to_restrictive, topofire_aws, topofire_poro))$
  toBands()

soils_extract = ee_extract(x = soils,
                           y = stations_ee,
                           scale = 800)

soils_extract_final = soils_extract %>%
  as_tibble() %>%
  set_names(c("network", "site_id", 'sand', 'silt', 'clay', 
              'som', 'bulk_density', 'awc', 'depth_to_restrictive'))


# Load PRISM 30-year normals dataset
prism = ee$ImageCollection("OREGONSTATE/PRISM/Norm91m")$
  mean()

prism_extract = ee_extract(x = prism,
                           y = stations_ee,
                           scale = 928)

prism_extract_final = prism_extract %>%
  as_tibble()

# elevation
dem = ee$Image('NASA/NASADEM_HGT/001')$
  select('elevation')

dem_extract = ee_extract(x = dem,
                         y = stations_ee,
                         scale = 30)

dem_extract_final = dem_extract %>%
  as_tibble()

# merge and export 
final = dem_extract_final %>%
  left_join(prism_extract_final, by = c('network', 'site_id')) %>%
  left_join(soils_extract_final, by = c('network', 'site_id'))

write_csv(final,'~/nrcs-scd-soil-moisture-eval-data/processed/station-covariates.csv')
