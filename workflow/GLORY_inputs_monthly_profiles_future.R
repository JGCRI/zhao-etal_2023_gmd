################################################################################
# For Linear Programming Model                                   
# Monthly Profiles for Demand, Evaporation, and Inflow using projected xanthos    
# Contact: Mengqi Zhao                                            
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-08-10
################################################################################

# Load packages
library(dplyr)
library(data.table)
library(lubridate)
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

f_avgchflow <- paste0('avgchflow_m3persec_', scenario_name, '_1950_2099.csv')
f_pet <- paste0('pet_km3permonth_', scenario_name, '_1950_2099.csv')


# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
# average channel flow
avgchflow <- data.table::fread(file.path(xanthos.dir, f_avgchflow))

# VIC WATCH observation
vic <- data.table::fread(file.path(xanthos.dir, 'vic_watch_basin_km3_1971_2001_monthly.csv'))

# PET
pet <- data.table::fread(file.path(xanthos.dir, f_pet))

# reservoir purpose
# (1-hydropower, 2-irrigation, 3-flood-control and others)
purpose <- data.table::fread(file.path(data.dir, 'reservoir_purpose.csv')) %>% 
  dplyr::rename(purpose = V1)

# xanthos coordinates (from xanthos input data)
coord <- data.table::fread(file.path(data.dir, 'coordinates.csv'), header = FALSE)
coord_basin <- data.table::fread(file.path(data.dir, 'basin.csv'), header = TRUE)
basin <- data.table::fread(file.path(data.dir, 'basin_ID.csv'), header = TRUE, skip = 6)
basin_rmap <- data.table(basin_id = as.numeric(as.character(rmap::mapGCAMBasins$subRegionAlt)),
                         basin_name = rmap::mapGCAMBasins$subRegion)

# merge coordinates and basin info
basin_grid <- coord_basin %>% 
  dplyr::mutate(xanthos_id = row_number()) %>% 
  dplyr::left_join(basin, by = c('basin' = 'basin_id')) %>% 
  dplyr::rename(basin_id = basin)
coord_grid <- coord %>% 
  dplyr::select(V1, V2, V3) %>% 
  dplyr::rename(xanthos_id = V1,
                lon = V2,
                lat = V3) %>% 
  dplyr::left_join(basin_grid, by = c('xanthos_id'))

# read contributing grids for the outlet cell in each basin
outlet_grid <- data.table::fread(file.path(data.dir, 'outlet.csv')) %>%
  dplyr::select(V1, V2, V3) %>%
  dplyr::mutate(basin_id = seq(1, 235, 1)) %>% 
  dplyr::left_join(basin, by = 'basin_id')

# selected time series
ts_select <- format(seq.Date(from = as.Date('2016/1/1'),
                             to = as.Date('2050/12/1'),
                             by = 'month'),
                    format = '%Y%m')

# ------------------------------------------------------------------------------
# Inflow
# ------------------------------------------------------------------------------

# Calculate average annual inflow in each grid at 5 year interval
inflow_reservoir <- avgchflow %>%
  dplyr::select(c('id', all_of(ts_select))) %>% 
  dplyr::rename(xanthos_id = id) %>%
  dplyr::left_join(coord_grid, by = c('xanthos_id')) %>%
  dplyr::left_join(purpose %>% dplyr::mutate(xanthos_id = row_number()), 
                   by = 'xanthos_id') %>% 
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d'),
                value = if_else(is.na(value), 0, value)) %>% 
  ## Only considers inflows to irrigation and flood control reservoirs
  dplyr::filter(purpose %in% c(2, 3)) %>% 
  ## add all the cells within each basin for each month
  dplyr::group_by(basin_id, basin_name, time) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(month = lubridate::month(time),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = 2020 + (period - 1) * 5) %>% 
  ## get average monthly value (e.g., mean Jan, mean Feb, ...) 
  ## for each 5-year interval, unit: m3/s
  dplyr::group_by(basin_id, basin_name, period, month) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin_id, basin_name, period) %>% 
  dplyr::mutate(value_total = sum(value),
                value = value / value_total) %>% 
  dplyr::ungroup()

# Get average channel flow behavior in each basin
avgchflow_basin <- avgchflow %>%
  dplyr::select(c('id', all_of(ts_select))) %>% 
  dplyr::rename(xanthos_id = id) %>%
  dplyr::left_join(coord_grid, by = c('xanthos_id')) %>% 
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d'),
                value = if_else(is.na(value), 0, value)) %>%
  ## add all the cells within each basin for each month
  dplyr::group_by(basin_id, basin_name, time) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(month = lubridate::month(time),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = 2020 + (period - 1) * 5) %>% 
  ## get average monthly value (e.g., mean Jan, mean Feb, ...), unit: m3/s
  dplyr::group_by(basin_id, basin_name, period, month) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin_id, basin_name, period) %>% 
  dplyr::mutate(value_total = sum(value),
                value = value / value_total) %>% 
  dplyr::ungroup()

# VIC runoff profile
vic_basin <- vic %>% 
  dplyr::group_by(basin, month) %>% 
  dplyr::summarise(value = mean(q)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin) %>% 
  dplyr::mutate(value_total = sum(value),
                value = value / value_total) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(basin_id = basin) %>% 
  dplyr::left_join(basin, by = c('basin_id')) %>% 
  dplyr::select(basin_id, basin_name, month, value)

# finalize the inflow profile for each period
# Combine inflow and channel flow profiles
# if there is no value in inflow_reservoir, replace with avgchflow_basin
# if still no value, use vic_basin
inflow_profile <- avgchflow_basin %>% 
  dplyr::left_join(inflow_reservoir,
                   by = c('basin_id', 'basin_name', 'period', 'month'),
                   suffix = c('.basin', '.reservoir')) %>% 
  dplyr::left_join(vic_basin %>% dplyr::rename(value.vic = value),
                   by = c('basin_id', 'basin_name', 'month')) %>% 
  dplyr::mutate(note = if_else(is.na(value.reservoir), 'Basin Average*', 'Reservoir*'),
                value = if_else(is.na(value.reservoir), value.basin, value.reservoir),
                note = if_else(is.na(value), 'Basin VIC', note),
                value = if_else(is.na(value), value.vic, value)) %>% 
  dplyr::select(basin_id, basin_name, period, month, value, note)


# ------------------------------------------------------------------------------
# PET
# ------------------------------------------------------------------------------
pet_reservoir <- pet %>% 
  dplyr::select(c('id', all_of(ts_select))) %>% 
  dplyr::rename(xanthos_id = id) %>%
  dplyr::left_join(coord_grid, by = c('xanthos_id')) %>%
  dplyr::left_join(purpose %>% dplyr::mutate(xanthos_id = row_number()), by = 'xanthos_id') %>% 
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d')) %>% 
  dplyr::filter(purpose %in% c(1, 2, 3)) %>% 
  dplyr::group_by(basin_id, basin_name, time) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(month = lubridate::month(time),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = 2020 + (period - 1) * 5) %>% 
  dplyr::group_by(basin_id, basin_name, period, month) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin_id, basin_name, period) %>% 
  dplyr::mutate(value_total = sum(value),
                value = value / value_total) %>% 
  dplyr::ungroup()

pet_basin <- pet %>% 
  dplyr::select(c('id', all_of(ts_select))) %>% 
  dplyr::rename(xanthos_id = id) %>%
  dplyr::left_join(coord_grid, by = c('xanthos_id')) %>%
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts_select)) %>% 
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d')) %>%
  dplyr::group_by(basin_id, basin_name, time) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(month = lubridate::month(time),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = 2020 + (period - 1) * 5) %>% 
  dplyr::group_by(basin_id, basin_name, period, month) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin_id, basin_name, period) %>% 
  dplyr::mutate(value_total = sum(value),
                value = value / value_total) %>% 
  dplyr::ungroup()

pet_profile <- pet_basin %>% 
  dplyr::left_join(pet_reservoir,
                   by = c('basin_id', 'basin_name', 'period', 'month'),
                   suffix = c('.basin', '.reservoir')) %>% 
  dplyr::mutate(note = if_else(is.na(value.reservoir), 'Basin Average*', 'Reservoir*'),
                value = if_else(is.na(value.reservoir), value.basin, value.reservoir)) %>% 
  dplyr::select(basin_id, basin_name, period, month, value, note)


# ------------------------------------------------------------------------------
# Sectoral Water Demand
# ------------------------------------------------------------------------------
# using a saved data. since water demand is historical, it does not matter which
# gcm/rcp we are using
wd_file_name <- 'wd_tethys_2005_2010_gcam5p3-stash_GFDL-ESM2M_rcp2p6_235basins.RDS'

if(file.exists(wd_file_name)){
  
  # saveRDS(wd_tethys, file = wd_file_name)
  wd_tethys_235basins <- readRDS(file = wd_file_name)
  
  # Calculate sectoral demand profiles
  demand_sector_profile <- wd_tethys_235basins %>% 
    dplyr::select(subRegion, class, x, value) %>% 
    dplyr::rename(time = x) %>% 
    dplyr::mutate(time = as.Date(paste0(time, '01'), format = '%Y%m%d')) %>% 
    dplyr::mutate(month = lubridate::month(time)) %>% 
    dplyr::group_by(subRegion, class, month) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(subRegion, class) %>% 
    dplyr::mutate(value_total = sum(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(value = value / value_total)
  
  demand_sector_profile <- demand_sector_profile %>% 
    dplyr::rename(basin_name = subRegion) %>% 
    dplyr::left_join(basin_rmap, by = c('basin_name'))
  
  # calculate demand profiles with all sectors together
  demand_profile <- wd_tethys_235basins %>% 
    dplyr::select(subRegion, class, x, value) %>% 
    # dplyr::filter(!class %in% 'Electric') %>% 
    dplyr::rename(time = x) %>% 
    dplyr::mutate(time = as.Date(paste0(time, '01'), format = '%Y%m%d')) %>% 
    dplyr::group_by(subRegion, time) %>% 
    dplyr::summarise(value = sum(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(month = lubridate::month(time)) %>% 
    dplyr::group_by(subRegion, month) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(subRegion) %>% 
    dplyr::mutate(value_total = sum(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(value = value / value_total)
  
  demand_profile <- demand_profile %>% 
    dplyr::rename(basin_name = subRegion) %>% 
    dplyr::left_join(basin_rmap, by = c('basin_name'))
  
}

demand_sector_profile_spread <- demand_sector_profile %>% 
  dplyr::select(-value_total) %>% 
  dplyr::mutate(class = tolower(class)) %>% 
  tidyr::spread(key = class, value = value)

demand_sector_total <- demand_sector_profile %>% 
  dplyr::select(basin_id, basin_name, class, value_total) %>% 
  unique() %>% 
  dplyr::mutate(class = tolower(class)) %>% 
  dplyr::rename(demand_ann = value_total,
                sector = class)

write.csv(demand_sector_total, 
          file = file.path(output.dir, 'LP_demand_sector_mean_annual_hist.csv'), 
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Format output
# ------------------------------------------------------------------------------

output <- inflow_profile %>% 
  dplyr::rename(inflow = value) %>% 
  dplyr::select(basin_id, period, month, inflow) %>% 
  dplyr::left_join(pet_profile %>% 
                     dplyr::rename(pet = value) %>% 
                     dplyr::select(basin_id, period, month, pet),
                   by = c('basin_id', 'period', 'month')) %>% 
  dplyr::left_join(demand_sector_profile_spread,
                   by = c('basin_id', 'month')) %>% 
  dplyr::arrange(basin_id, basin_name, period, month) %>% 
  dplyr::select(basin_id, basin_name, period, month, inflow, pet, 
                domestic, electric, industry, irrigation, livestock, mining)

write.csv(output, 
          file = file.path(output.dir, 
                           paste0('LP_fraction_profile_', gcm, '_', rcp, '.csv')), 
          row.names = FALSE)
