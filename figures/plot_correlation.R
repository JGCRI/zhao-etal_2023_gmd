################################################################################
# Correlation between P, PET, COst, Yield, Slope, etc
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2023-02
################################################################################

library(data.table)
library(dplyr)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(rmap)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
setwd(work.dir)

# set paths
figure.dir <- file.path(work.dir, 'outputs', 'correlation')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}
output.dir <- file.path(work.dir, 'outputs')
output.glory.dir <- file.path(work.dir, 'outputs_glory')
output.workflow.dir <- file.path(dirname(work.dir), 'workflow', 'outputs')

# monthly xanthos streamflow output
xanthos.dir <- file.path(dirname(work.dir), 'workflow', 'inputs', 'hydrology_xanthos')
gcm <- 'MIROC-ESM-CHEM'
rcp <- 'rcp6p0'
scenario_name <- paste(gcm, rcp, sep = '_')
f_avgchflow <- paste0('avgchflow_m3persec_', scenario_name, '_1950_2099.csv')

# monthly maxSubResource runoff and runoff profile
f_runoff <- file.path(output.workflow.dir, 'LP_climate_MIROC-ESM-CHEM_rcp6p0.csv')
f_runoff_profile <- file.path(output.workflow.dir, 'LP_fraction_profile_MIROC-ESM-CHEM_rcp6p0.csv')

runoff <- data.table::fread(f_runoff) %>%
  dplyr::select(-pet_km)
runoff_profile <- data.table::fread(f_runoff_profile) %>%
  dplyr::select(basin_id, basin_name, period, month, inflow)


#  Experiment
exp <- data.frame(run = c('20230226_11h16m', # A + B: feedback on
                          '20230226_13h14m'), # A + C: feedback off
                  scenario = c('Climate_FB_ON',
                               'Climate_FB_OFF'))  



# ------------------------------------------------------------------------------
# Set scenarios names basin on run name, Mapping
# ------------------------------------------------------------------------------

gcam_basin_mapping <- rmap::mapping_gcambasins %>% 
  dplyr::left_join(rmap::mapping_tethys_grid_basin_region_country %>% 
                     dplyr::select(basinID, basinName) %>%
                     unique(),
                   by = c('subRegionMap' = 'basinName')) %>% 
  dplyr::mutate(resource = paste0(subRegion, '_water withdrawals')) %>%
  dplyr::rename(basin_name = subRegionMap,
                basin_name_gcam = subRegion,
                basin_id = basinID) %>%  
  dplyr::mutate_all(function(x) gsub('Madasgacar', 'Madagascar', x))
gcam_basin_mapping$basin_id <- as.integer(gcam_basin_mapping$basin_id)

coord_grid <- rmap::mapping_tethys_grid_basin_region_country %>% 
  dplyr::select(lat, lon, xanthos_id = gridID, 
                basin_id = basinID, basin_name = basinName)

# ------------------------------------------------------------------------------
# Define Functions
# ------------------------------------------------------------------------------

# Function to read before run
read_files <- function(file, year_loc, scen, header = TRUE, sep = 'auto', ...){
  data <- data.table::fread(file, header = header, sep = sep, ...)
  filename <- basename(file[1])
  year <- strsplit(filename, '\\_|\\.')[[1]][year_loc]
  data$year <- as.numeric(year)
  data$scenario <- scen
  return(data)
}

# ------------------------------------------------------------------------------
# Load data and Prepare Data
# ------------------------------------------------------------------------------

# demand -----------------------------------------------------------------------
# Read demand profile output from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_demand_profile',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

demand_profile <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::mutate(year = as.integer(year)) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('basin_id')) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, month, fraction)

# Read total demand
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'gcam_outputs', run),
                       pattern = 'gcam_demand_source',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

demand_source <- rbindlist(all_data, idcol = FALSE) %>% 
  tidyr::separate(col = subresource, sep = ' ', into = c('subresource', 'grade')) %>% 
  dplyr::group_by(scenario, resource, year) %>% 
  dplyr::summarise(value = sum(`physical-output`)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('resource')) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, value) %>% 
  tibble::as_tibble()

# get monthly demand
demand_monthly <- demand_profile %>% 
  dplyr::left_join(demand_source ,
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  dplyr::mutate(demand_monthly = fraction * value) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, month, demand_monthly)

# runoff -----------------------------------------------------------------------
# runoff data 1: to calculate deficit between runoff and demand
# get monthly runoff
runoff_monthly <- runoff_profile %>% 
  dplyr::left_join(runoff, by = c('basin_id', 'basin_name', 'period')) %>% 
  dplyr::mutate(runoff_monthly = inflow * smooth_runoff_km3) %>% 
  dplyr::select(basin_id, basin_name, year = period, month, runoff_monthly)


# combine runoff and demand
monthly <- demand_monthly %>% 
  dplyr::left_join(runoff_monthly,
                   by = c('basin_id', 'basin_name', 'year', 'month'))

# calculate annual deficit
deficit <- monthly %>% 
  dplyr::mutate(deficit = dplyr::if_else(runoff_monthly < demand_monthly,
                                         demand_monthly - runoff_monthly,
                                         0),
                deficit_count = dplyr::if_else(runoff_monthly < demand_monthly,
                                               1,
                                               0),
                surplus = dplyr::if_else(runoff_monthly < demand_monthly,
                                         0,
                                         runoff_monthly - demand_monthly)) %>% 
  dplyr::group_by(scenario, basin_id, basin_name, year) %>% 
  dplyr::summarise(deficit = sum(deficit),
                   surplus = sum(surplus),
                   demand = sum(demand_monthly),
                   runoff = sum(runoff_monthly),
                   duration = sum(deficit_count)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(deficit_ratio = dplyr::if_else(demand == 0, 0, deficit / demand),
                satisfy_ratio = 1 - deficit_ratio,
                drought_intensity = dplyr::if_else(duration > 0,
                                                   deficit / duration,
                                                   0))

# runoff data 2: to calculate Standardized Runoff Index (SRI)
# average channel flow (m3/s)
avgchflow <- data.table::fread(file.path(xanthos.dir, f_avgchflow))

# selected time series
start_yr <- 2016
ts_select <- format(seq.Date(from = as.Date(paste0(start_yr, '/1/1')),
                             to = as.Date('2050/12/1'),
                             by = 'month'),
                    format = '%Y%m')

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
                year = lubridate::year(time),
                period = cut(time, breaks = '5 years', labels = FALSE),
                period = (start_yr + 4) + (period - 1) * 5)


# Read demand and yield from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_gcam_capacity',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

demand_yield <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('basin_id')) %>% 
  dplyr::select(scenario, basin_id, basin_name, year,
                demand_gcam, yield_min, yield_max)

# Read demand output from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'gcam_outputs', run),
                       pattern = 'gcam_demand_source',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}


demand_source <- rbindlist(all_data, idcol = FALSE) %>% 
  tidyr::separate(col = subresource, sep = ' ', into = c('subresource', 'grade')) %>% 
  dplyr::group_by(scenario, resource, subresource, year) %>% 
  dplyr::summarise(value = sum(`physical-output`)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('resource')) %>% 
  dplyr::select(scenario, basin_id, basin_name, resource = subresource, year, value) %>%
  tidyr::pivot_wider(values_from = value, names_from = resource) %>% 
  dplyr::rename(demand_runoff = runoff,
                demand_gw = groundwater) %>% 
  tibble::as_tibble()

# capacity - yield info --------------------------------------------------------
# load LP_reservoir.csv to get reservoir expansion info
reservoir <- data.table::fread(
  file.path(dirname(output.workflow.dir), 'LP_reservoir.csv')
)
reservoir <- reservoir %>% 
  dplyr::select(basin_id, mean_cap_km3, expan_cap_km3, total_cap_km3)

# Read capacity yield curve from LP model
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_capacity_yield_curve',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 5, scenario)
  all_data <- c(all_data, data)
}


capacity_yield <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('basin_id')) %>% 
  dplyr::left_join(reservoir, by = 'basin_id') %>% 
  dplyr::mutate(expan_cap_km3 =
                  dplyr::if_else(expan_cap_km3 > total_cap_km3,
                                 expan_cap_km3,
                                 total_cap_km3)) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, 
                capacity, yield, expan_cap_km3, total_cap_km3) %>% 
  tibble::as_tibble()

curve_slope <- capacity_yield %>% 
  dplyr::left_join(demand_yield, 
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>%  
  dplyr::group_by(scenario, basin_id, basin_name, year) %>% 
  dplyr::mutate(delta_capacity = capacity - lag(capacity, default = capacity[1]),
                delta_yield = yield - lag(yield, default = yield[1]),
                slope = dplyr::if_else(delta_capacity == 0, 0, delta_yield/delta_capacity),
                slope = round(slope, digits = 4),
                yield_lag = dplyr::lag(yield, default = 0),
                target = dplyr::if_else(demand_gcam >= yield_lag & demand_gcam < yield, 1, 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(scenario, basin_id, basin_name, year) %>% 
  dplyr::summarise(slope_median = median(slope),
                   slope_target = sum(target * slope)) %>% 
  dplyr::ungroup()


# reservoir cost ---------------------------------------------------------------
# Read unit cost output from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_unit_cost',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

unit_cost <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::left_join(reservoir, by = c('basin_id')) %>% 
  # million USD per expansion
  dplyr::mutate(expan_cost = 10^9 * unit_cost * mean_cap_km3 / 10^6) %>% 
  dplyr::select(scenario, basin_id, year, expan_cost, unit_cost)


# Low-cost Renewable Water -----------------------------------------------------
# Read updated supply curves
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_supply_curve',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

supply_curve <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::mutate(year = as.integer(year)) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('resource')) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, 
                grade, available, extractioncost)

# Read maxsubresource
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'lp_outputs', run),
                       pattern = 'LP_maxsubresource',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 3, scenario)
  all_data <- c(all_data, data)
}

maxsubresource <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::mutate(year = as.integer(year)) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('resource')) %>% 
  dplyr::select(scenario, basin_id, basin_name, year, maxSubResource)

# for selected scenarios abd 2050
lowcost_water <- data.frame()
for(id in unique(gcam_basin_mapping$basin_id)) {
  for(scen in exp$scenario) {
    for(yr in seq(2020, 2050, 5)) {
      
      temp <- supply_curve %>% 
        dplyr::filter(basin_id == id, scenario == scen, year == yr)
      
      val_maxsubresource <- maxsubresource$maxSubResource[
        maxsubresource$basin_id == id & 
          maxsubresource$scenario == scen & 
          maxsubresource$year == yr
        ]
      
      if(length(val_maxsubresource) == 0){
        val_maxsubresource <- 0
      }
      
      cheap_fraction <- spline(x = temp$extractioncost, y = temp$available,
                               xout = 0.001, method = 'hyman')
      cheap_fraction_df <- data.frame(scenario = scen,
                                      basin_id = id,
                                      basin_name = unique(temp$basin_name),
                                      year = yr,
                                      cheap_fraction = min(cheap_fraction$y, 1)) %>% 
        dplyr::mutate(cheap_vol = cheap_fraction * val_maxsubresource)
      
      lowcost_water <- dplyr::bind_rows(lowcost_water, cheap_fraction_df)
      
    }
    
  }
  
}

# GCAM solved price ------------------------------------------------------------
# Read price output from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.glory.dir, 'gcam_outputs', run),
                       pattern = 'price_',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 2, scenario)
  all_data <- c(all_data, data)
}


price <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::mutate(year = as.integer(year)) %>% 
  dplyr::left_join(
    gcam_basin_mapping %>% 
      dplyr::mutate(
        market = paste0(basin_name_gcam, basin_name_gcam, '_water withdrawals')), 
    by = c('market')) %>% 
  dplyr::mutate(label = paste0(basin_id, ' - ', basin_name),
                label = forcats::fct_reorder(label, basin_id),
                price.new = dplyr::if_else(price.new < 0, 0, price.new),
                price.original = dplyr::if_else(price.original < 0, 0, price.original),
                price.diff = price.new - price.original,
                price.diff.pct = dplyr::if_else(price.original != 0,
                                                100 * (price.new - price.original) / ((price.new + price.original) / 2) ,
                                                0)) %>% 
  dplyr::select(scenario, basin_id, basin_name, label, year, 
                price.original, price.new, price.diff, price.diff.pct) %>% 
  tibble::as_tibble()


# join all variables -----------------------------------------------------------
df <- curve_slope %>% 
  dplyr::left_join(demand_yield, 
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  dplyr::left_join(demand_source, 
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  dplyr::left_join(unit_cost,
                   by = c('scenario', 'basin_id', 'year')) %>%
  dplyr::left_join(price %>% 
                     dplyr::select(scenario, basin_id, year, price = price.new),
                   by = c('scenario', 'basin_id', 'year')) %>% 
  dplyr::left_join(deficit,
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  dplyr::left_join(lowcost_water,
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  tidyr::drop_na()


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
library(GGally)
# library(caret)

# Change between Scenarios ----
# different methods of scaling is referred to at 
# https://www.digitalocean.com/community/tutorials/normalize-data-in-r

select_scenario <- c('Climate_FB_ON', 'Climate_FB_OFF')
cor_bin <- function(data, mapping, ...) {
  
  ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_point(...) +
    ggplot2::geom_smooth(method = 'lm', se = FALSE, 
                         method.args = list(family = "symmetric"))
}

for (select_year in seq(2020, 2050, 5)){
 
  df_plot <- df %>% 
    dplyr::filter(scenario %in% select_scenario, year %in% select_year,
                  # drought_intensity > 0,
                  # slope_target > 0,
                  # price < 0.5, expan_cost > 50
                  price > 0
    ) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(scenario = dplyr::case_when(scenario == 'Climate_FB_ON' ~ 'Feedback',
                                              scenario == 'Climate_FB_OFF' ~ 'No Feedback')) %>% 
    # Type: normalize with log
    dplyr::mutate_at(
    vars(price, demand, demand_runoff, demand_gw, 
         unit_cost, expan_cost, 
         slope_target, slope_median, drought_intensity,
         cheap_vol
         ), 
    ~log(as.data.frame(.))[[1]]) %>% 
    dplyr::select(scenario, basin_id, basin_name,
                  price, demand, demand_runoff, demand_gw, 
                  unit_cost, expan_cost, 
                  slope_target, slope_median, drought_intensity,
                  cheap_vol)
  
  
  pm <- 
    GGally::ggpairs(df_plot, aes(color = scenario, alpha = 0.8),
                    columns = c(
                      # 'scenario','demand_gw', 'demand',
                      'price', 'demand_runoff', 
                      'expan_cost',
                      'slope_median', 'drought_intensity',
                      'cheap_vol'
                      ),
                    columnLabels = c(
                      # 'Scenario','Groundwater Withdrawals', 'Total Withdrawals',
                      'Water Price', 'Surface Water Withdrawals', 
                      'Expansion Cost',
                      'Yield Gain/Unit Expansion', 'Socio-Eco Drought Intensity',
                      'Low-Cost Water Withdrawals'
                      ),
                    # lower = list(continuous = 'density'),
                    lower = list(continuous = cor_bin),
                    # upper = list(continuous = cor_bin),
                    upper = list(continuous = 'cor')
    ) +
    ggplot2::scale_colour_manual(values = c("darkorange","cyan4")) +
    ggplot2::scale_fill_manual(values = c("darkorange","cyan4")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = 'gray90'),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
    )

  
  ggplot2::ggsave(
    plot = pm,
    file.path(figure.dir, 
              paste0(
                'scenario_cor_axis_',
                # 'scenario_cor_',
                # 'scenario_point_',
                paste(select_scenario, collapse = '-'),
                '_', select_year, '_all.png')),
    height = 10, width = 12, unit = 'in', dpi = 600)
}


