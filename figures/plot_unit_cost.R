################################################################################
# Plot Unit Cost
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(sf)
library(ggnewscale)
library(patchwork)
library(paletteer)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
setwd(work.dir)

# set paths
figure.dir <- file.path(work.dir, 'outputs', 'unit_cost')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}

# Output paths
output.dir <- file.path(work.dir, 'outputs_glory')

#  Experiment
exp <- data.frame(run = c('climate_feedback_on', # A + B: feedback on
                          'climate_feedback_off'), # A + C: feedback off
                  scenario = c('Climate_FB_ON',
                               'Climate_FB_OFF'))  
select_scenario <- 'Climate_FB_ON'

# ------------------------------------------------------------------------------
# Spatial SF Data
# ------------------------------------------------------------------------------
# Load RDS data
world <- readRDS(file.path(dirname(figure.dir), 'world_basemap_robinson_projection.rds'))

# spatial reference code
# set epsg code 4326: https://spatialreference.org/ref/epsg/wgs-84/
epsg_code <- 4326

#  basin simple feature
basin_sf <- rmap::mapGCAMBasins %>% 
  dplyr::select(basin_id = subRegionAlt,
                basin_name = subRegion,
                geometry) %>% 
  sf::st_transform(epsg_code)

#  basin area
basin_area <- rmap::mapGCAMBasins %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(basin_name = subRegion, area) %>% 
  dplyr::mutate(basin_area_km2 = as.numeric(area)/10^6) %>% 
  dplyr::select(basin_name, basin_area_km2)

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

# Read demand output from GCAM
all_data <- list()
for(run in exp$run){
  f_demand_source <- list.files(file.path(output.dir, 'lp_outputs', run),
                                pattern = 'LP_unit_cost',
                                full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_demand_source, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

unit_cost <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::filter(scenario == select_scenario,
                year == 2020) 
 

# convert unit cost from 1975USD to 2020USD
# World bank GDP deflator https://data.worldbank.org/indicator/NY.GDP.DEFL.ZS?locations=US
# 1975 = 28.5, 2020 = 108.6
unit_cost <- unit_cost %>% 
  dplyr::mutate(unit_cost_2020USD = 108.6 * unit_cost / 28.5)

# make sf object
unit_cost_sf <- basin_sf %>% 
  dplyr::left_join(unit_cost,by = 'basin_id') 


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
# customize the breaks for original 1975USD
breaks <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 
            1.2, 1.6, 2.0, 3.0, 5.0, 7.0, Inf)
# labels <- paste(breaks[1:length(breaks) - 1], ' to ', breaks[2:length(breaks)])
unit_cost_sf <- unit_cost_sf %>% 
  dplyr::mutate(unit_cost_group = cut(unit_cost, breaks, right = FALSE))
title <- bquote('Unit Cost of Reservoir Construction (1975 USD/'~m^3~')')

# if using 2020 USD
# breaks <- c(0, 0.3, 0.6, 1.2, 1.8, 2.4, 3.0,
#             5.0, 7.0, 10, 15, 20, Inf)
breaks <- c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0,
            5.0, 7.0, 10, 15, 20, 40, 60.0)
unit_cost_sf <- unit_cost_sf %>% 
  dplyr::mutate(unit_cost_group = cut(unit_cost_2020USD, breaks, right = FALSE))
title <- bquote('Unit Cost of Reservoir Construction (2020 USD/'~m^3~')')

# Plot
world +
  ggnewscale::new_scale_fill() +
  # ggplot2::ggplot() +
  ggplot2::geom_sf(data = unit_cost_sf,
                   ggplot2::aes(fill = unit_cost_group),
                   color = 'grey30', linetype = 1, lwd = 0.25) +
  ggplot2::scale_fill_manual(
    values = paste0(rev(paletteer_c("grDevices::Purple-Yellow", length(breaks)))),
    guide = ggplot2::guide_colorsteps(
      barwidth = 20,
      barheight = 0.5,
      show.limits = T,
      frame.colour = 'gray50',
      ticks = TRUE,
      ticks.linewidth = 0.5,
      ticks.colour = 'gray50',
      title = title,
      title.position = "top"
    )) +
  ggplot2::coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = 'white', color = NA),
                 panel.background = ggplot2::element_rect(fill = 'white', color = NA),
                 legend.position = "bottom",
                 legend.margin = ggplot2::margin(t = 10, b = 2),
                 legend.title = ggplot2::element_text(size = 10),
                 legend.text = ggplot2::element_text(angle = 0,
                                                     margin = ggplot2::margin(t = 5)))

ggplot2::ggsave(file.path(figure.dir, paste('unit_cost_', select_scenario, '_2020.png')),
                height = 6, width = 10, unit = 'in', dpi = 600)
