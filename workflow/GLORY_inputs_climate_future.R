################################################################################
# For Linear Programming Model
# Generate Precipitation and ET inputs
# For future climate
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-08-11
################################################################################

library(data.table)
library(minpack.lm)
library(dplyr)
library(tidyr)
library(rmap)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/workflow' # update correspondingly
setwd(work.dir)
data.dir <- file.path(work.dir, 'data')
output.dir <- file.path(work.dir, 'outputs')

xanthos.dir <- file.path(data.dir, 'hydrology_xanthos')

gcm <- 'MIROC-ESM-CHEM'
rcp <- 'rcp6p0'
scenario_name <- paste(gcm, rcp, sep = '_')
smooth_runoff_var_name <- 'runoff_impacts'


f_runoff_smooth <- file.path(data.dir,
                      paste0(paste(smooth_runoff_var_name, gcm, rcp, sep = '_'), 
                             '.csv'))
f_runoff <- file.path(xanthos.dir, 
                      paste0('Basin_runoff_km3permonth_', scenario_name, '_1950_2099.csv'))
f_pet <- paste0('pet_km3permonth_', scenario_name, '_1950_2099.csv')

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
# Read smoothed Xanthos runoff output in km3
runoff_smooth <- data.table::fread(f_runoff_smooth, header = TRUE, skip = 4)
runoff <- data.table::fread(f_runoff, header = TRUE)

# PET
pet <- data.table::fread(file.path(xanthos.dir, f_pet))

# Read Xanthos grid area (ha)
grid_area <- data.table::fread(file.path(data.dir, 'Grid_Areas_ID.csv'),
                               header = FALSE) %>% 
  dplyr::mutate(xanthos_id = row_number(),
                grid_area_KM2 = 0.01 * V1) %>%  # convert to km2
  dplyr::select(-V1)

# Read dam data and GCAM basins

basin_id <- data.frame(basin_id = rmap::mapGCAMBasins$subRegionAlt,
                       basin_name = rmap::mapGCAMBasins$subRegion)
gcam_basin_mapping <- rmap::mapping_gcambasins %>% 
  dplyr::rename(basin_name = subRegionMap,
                gcam_basin_name = subRegion) %>% 
  dplyr::mutate(gcam_basin_name = sub('Madasgacar', 'Madagascar', gcam_basin_name),
                basin_name = sub('Madasgacar', 'Madagascar', basin_name))
grid_info <- data.table::fread(file.path(data.dir, 'basin.csv'), 
                                header = TRUE) %>% 
  dplyr::mutate(xanthos_id = row_number()) %>% 
  dplyr::rename(basin_id = basin) %>% 
  dplyr::left_join(grid_area, by = 'xanthos_id') %>% 
  dplyr::left_join(basin_id, by = 'basin_id')


# Georeferenced GranD data
grand_georef <-
  data.table::fread(file.path(data.dir, 
                              'GranD_v1.3_remap_to_xanthos.csv'), 
                    header = TRUE) %>% 
  dplyr::left_join(basin_id, by = c('basins'='basin_id')) %>% 
  dplyr::rename(basin_id = basins)

# Get reservoir capacity, inflow, and evaporation from reservoir surface area
# Reservoir purpose is also labeled for further filtering
dam_grid_info <- grand_georef %>% 
  dplyr::left_join(world.dam %>%
                     dplyr::select(GRAND_ID, AREA_SKM),
                   by = c('GRAND_ID')) %>% 
  dplyr::mutate(CAP_KM3 = CAP_MCM / 1000,
                purpose = case_when(purpose == 0 ~ 'none',
                                    purpose == 1 ~ 'hydropower',
                                    purpose == 2 ~ 'irrigation',
                                    purpose == 3 ~ 'flood control')) %>% 
  dplyr::select(basin_id, basin_name, xanthos_id, CAP_KM3, AREA_SKM, purpose)

# selected time series
ts_select <- format(seq.Date(from = as.Date('2016/1/1'),
                             to = as.Date('2050/12/1'),
                             by = 'month'),
                    format = '%Y%m')

# ------------------------------------------------------------------------------
# Climate Inputs
# ------------------------------------------------------------------------------

# Inflow ---
# Use the same runoff with the smoothed runoff from xanthos (km3/year)
interval <- 5
start_period <- 2020
end_period <- 2050
missing_gcam_basin <- c('North_Marina_Islands_and_Guam', 
                        'Andaman_Nicobar_Islands',
                        'Micronesia', 
                        'Antarctica')
# use unsmoothed runoff to supplement the smoothed runoff because of those 4
# missing basins
runoff_supplement <- runoff %>% 
  dplyr::select(c('id', all_of(ts_select))) %>% 
  dplyr::rename(basin_id = id) %>% 
  dplyr::left_join(basin_id, by = c('basin_id')) %>%
  dplyr::filter(basin_name %in% missing_gcam_basin) %>% 
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(
    time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d'),
    period = cut(time, breaks = '5 years', labels = FALSE),
    period = 2020 + (period - 1) * 5
  ) %>%
  # get average value over each 5-year interval for each basin
  # unit: from km3/month to km3/year
  dplyr::group_by(basin_name, period) %>%
  dplyr::summarise(smooth_runoff_km3 = mean(value) * 12) %>%
  dplyr::ungroup()


runoff_basin <- runoff_smooth %>% 
  dplyr::filter(year.fillout %in% seq(start_period, end_period, interval)) %>% 
  tidyr::separate(col = renewresource, 
                  into = c('gcam_basin_name', 'sector'), sep = '_') %>% 
  dplyr::rename(period = year.fillout,
                smooth_runoff_km3 = maxSubResource) %>% 
  dplyr::select(gcam_basin_name, period, smooth_runoff_km3) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('gcam_basin_name')) %>% 
  dplyr::select(basin_name, period, smooth_runoff_km3) %>% 
  dplyr::bind_rows(runoff_supplement)


# PET ---
# Calculate average annual PET over each 5-year intervals in the future
# 2020 to 2050 (km3/month -> km/year), convert to depth
pet_ann <- pet %>% 
  dplyr::select(c('id', all_of(ts_select))) %>%
  dplyr::rename(xanthos_id = id) %>%  
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d'),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = 2020 + (period - 1) * 5) %>% 
  # get average value over each 5-year interval for each cell, to unit: km3/year
  dplyr::group_by(xanthos_id, period) %>% 
  dplyr::summarise(pet_km3 = mean(value) * 12) %>% 
  dplyr::ungroup() %>%  
  dplyr::left_join(grid_info, by = c('xanthos_id')) %>%
  # Convert pet to depth: km/year
  dplyr::mutate(pet_km = pet_km3 / grid_area_KM2) %>% 
  dplyr::select(xanthos_id, basin_id, basin_name, period, pet_km)

# get mean pet depth (km/year) from irrigation and flood control reservoirs cells
# get total pet volume (km3/year) from irrigation and fc reservoirs
pet_basin <- pet_ann %>% 
  dplyr::group_by(basin_id, basin_name, period) %>% 
  dplyr::summarise(pet_basin_km = mean(pet_km)) %>% 
  dplyr::ungroup()
  
  
pet_nonhydro <- pet_ann %>%
  dplyr::left_join(dam_grid_info, by = c('xanthos_id', 'basin_id', 'basin_name')) %>%
  dplyr::mutate(purpose = dplyr::if_else(
    purpose %in% c('irrigation', 'flood control', 'hydropower'),
    'select',
    'discard'
  )) %>%
  dplyr::group_by(basin_id, basin_name, period, purpose) %>%
  dplyr::summarise(pet_km = mean(pet_km)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(pet_km = dplyr::if_else(purpose == 'select', pet_km, 0)) %>% 
  dplyr::group_by(basin_id, basin_name, period) %>%
  dplyr::summarise(pet_km = sum(pet_km)) %>%
  dplyr::ungroup()

pet_output <- pet_basin %>% 
  dplyr::left_join(pet_nonhydro, by = c('basin_id', 'basin_name', 'period')) %>% 
  dplyr::mutate(pet_km = dplyr::if_else(pet_km == 0, pet_basin_km, pet_km)) %>% 
  dplyr::select(basin_id, basin_name, period, pet_km)
  
# merge runoff and pet depth
output_climate <- runoff_basin %>% 
  dplyr::left_join(pet_output, by = c('basin_name', 'period')) %>% 
  dplyr::select(basin_id, basin_name, period, smooth_runoff_km3, pet_km)

# write output
write.csv(x = output_climate,
          file = file.path(output.dir,
                           paste0('LP_climate_', gcm, '_', rcp, '.csv')),
          row.names = FALSE)