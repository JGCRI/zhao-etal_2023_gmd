################################################################################
# LP capacity yield curves
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(snow)
library(rmap)
library(rchart)
library(paletteer)
library(ggnewscale)
library(ggstatsplot)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
setwd(work.dir)

# set paths
figure.dir <- file.path(work.dir, 'outputs', 'capacity_yield')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}

output.dir <- file.path(work.dir, 'outputs_glory')
output.workflow.dir <- file.path(dirname(work.dir), 'workflow', 'outputs')

#  Experiment
exp <- data.frame(run = c('climate_feedback_on', # A + B: feedback on
                          'climate_feedback_off'), # A + C: feedback off
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
  dplyr::mutate(market = paste0(subRegion, subRegion, '_water withdrawals')) %>% 
  dplyr::rename(basin_name = subRegionMap,
                basin_name_gcam = subRegion,
                basin_id = basinID) %>%  
  dplyr::mutate_all(function(x) gsub('Madasgacar', 'Madagascar', x))
gcam_basin_mapping$basin_id <- as.integer(gcam_basin_mapping$basin_id)


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

# load LP_reservoir.csv to get reservoir expansion info
reservoir <- data.table::fread(
  file.path(output.workflow.dir, 'LP_reservoir.csv')
)
reservoir <- reservoir %>% 
  dplyr::select(basin_id, expan_cap_km3, total_cap_km3)

# Read capacity yield curve from LP model
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.dir, 'lp_outputs', run),
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


# Read demand and yield from GCAM
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.dir, 'lp_outputs', run),
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
  

# -----------------------------------------------
# Plot: Capacity - Yield curve Slope
# -----------------------------------------------
# Set up the path and scenarios
save_path <- file.path(figure.dir, 'expansion_slope')
if(!dir.exists(save_path)){
  dir.create(save_path)
}

scenario_select <- 'Climate_FB_ON'
year_select <- c(2050)
basin_select <- c(79, 86, 87, 89, 131, 217)

df <- capacity_yield %>% 
  dplyr::left_join(demand_yield, 
                   by = c('scenario', 'basin_id', 'basin_name', 'year')) %>% 
  dplyr::filter(scenario == scenario_select, year == year_select) 

df_slope <- df %>% 
  dplyr::group_by(scenario, basin_id, basin_name, year) %>% 
  dplyr::mutate(delta_capacity = capacity - lag(capacity, default = capacity[1]),
                delta_yield = yield - lag(yield, default = yield[1]),
                slope = dplyr::if_else(delta_capacity == 0, 0, delta_yield/delta_capacity),
                slope = round(slope, digits = 4),
                yield_lag = dplyr::lag(yield, default = 0)) %>% 
  dplyr::ungroup()

# line chart: slope of the entire capacity - yield curve
chart_path <- file.path(save_path, paste0(scenario_select, '_', year_select))
if(!dir.exists(chart_path)){
  dir.create(chart_path)
}

for(id in basin_select){
  df_plot <- df_slope %>% 
    dplyr::filter(basin_id == id)
  
  stage_1 <- unique(df_plot$slope)[c(2,3, 4)]
  
  df_plot <- df_plot %>% 
    dplyr::mutate(slope = dplyr::if_else(slope == 0, stage_1[1], slope), 
                  slope_stage = as.factor(
                    dplyr::if_else(slope > quantile(df_plot$slope, 0.8)[[1]], 1, 2)))
  
  basin_name <- gcam_basin_mapping$basin_name[gcam_basin_mapping$basin_id == id]
  file_name <- paste(id, basin_name, sep = ' - ')
  
  y_axis_coeff <- min(df_plot$yield) / max(df_plot$slope)
  
  point_demand <- unique(df_plot$demand_gcam)
  if(point_demand <= min(df_plot$yield)){
    point_capacity <- data.frame(capacity_point = 0,
                                 gcam_demand = point_demand)
  } else {
    point <- spline(x = df_plot$yield, y = df_plot$capacity,
                    xout = c(point_demand), method = 'hyman') 
    point_capacity <- data.frame(capacity_point = point$y,
                                 gcam_demand = point_demand)
  }
  
  
  ggplot2::ggplot() +
    ggplot2::geom_area(data = df_plot, aes(x = capacity, y = slope * y_axis_coeff),
                       alpha = 0.8, lwd = 0.5, fill = '#91D1C2', color = '#666666') +
    ggplot2::geom_line(data = df_plot, ggplot2::aes(x = capacity, y = yield), 
                       color = 'black', size = 1) +
    ggplot2::labs(title = gsub('_', ' ', basin_name),
                  x = bquote('Storage Capacity ('~km^3~')')) +
    ggplot2::scale_y_continuous(
      name = bquote('Yield ('~km^3~')'),
      sec.axis = ggplot2::sec_axis(
        ~ ./y_axis_coeff,
        name = bquote('Yield Gain Per Unit Capacity Expansion'))
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = point_capacity,
                        ggplot2::aes(x = capacity_point, y = gcam_demand, 
                                     fill = paste0('GCAM Solved (', year_select, ')')), 
                        size = 3.5, shape = 21, color = 'black', 
                        stroke = 1.5, show.legend = T) +
    ggplot2::labs(fill = '') +
    ggplot2::scale_fill_manual(values = c('#FF7F00')) +
    ggplot2::theme(
      title = ggplot2::element_text(size = 18),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = 'black', fill = NA, size = 1.2),
      panel.grid.major.x = ggplot2::element_line(colour = 'grey'),
      axis.ticks = ggplot2::element_line(linewidth = 1.2),
      axis.text = ggplot2::element_text(colour = 'black', size = 14),
      axis.title = ggplot2::element_text(size = 14),
      axis.title.y = element_text(margin = ggplot2::margin(t=0, r=5, b=0, l=0)),
      axis.title.y.right = element_text(margin = ggplot2::margin(t=0, r=0, b=0, l=5)),
      strip.background = ggplot2::element_rect(fill = NA),
      strip.placement = 'outside',
      legend.position = c(0.04, 1),
      legend.justification = c('left', 'top'),
      legend.box = 'horizontal',
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1,'cm'),
      legend.box.margin = ggplot2::margin(t=0, r=0,b=0,l=0,unit='pt'),
      legend.margin = ggplot2::margin(0,0,0,0),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      plot.subtitle = ggplot2::element_text(size = 12)
    )
  
  ggplot2::ggsave(file.path(chart_path, paste0(file_name, '.png')),
                  height = 3.8, width = 6.5, unit = 'in', dpi = 300)
}


# Map: median slope
df_plot <- df_slope %>% 
  dplyr::group_by(scenario, basin_id, basin_name, year) %>% 
  dplyr::summarise(slope = median(slope)) %>% 
  dplyr::ungroup()

plot_name <- paste0('median slope_', scenario_select, '_', 
                    year_select, '.png')


# make sf object
breaks <- c(seq(0, 14, 1)) # change based on the max value
pal <- paletteer_c("grDevices::Blue-Yellow", length(breaks))
sf_plot <- basin_sf %>% 
  dplyr::left_join(df_plot %>% dplyr::select(basin_id, slope), 
                   by = c('basin_id')) %>% 
  dplyr::mutate(slope = dplyr::if_else(is.na(slope), 0, slope))


# Plot
world +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_sf(data = sf_plot,
                   ggplot2::aes(fill = slope),
                   color = 'grey30', linetype = 1, lwd = 0.25) +
  ggplot2::scale_fill_gradientn(
    colours = c('white', paste0(rev(paletteer_c("grDevices::Blue-Yellow", length(breaks) - 1)))),
    breaks = breaks,
    guide = ggplot2::guide_colorbar(
      barwidth = 20,
      barheight = 0.4,
      draw.llim = TRUE,
      frame.colour = 'gray50',
      direction = 'horizontal',
      label.position = 'bottom',
      ticks = TRUE,
      ticks.linewidth = 1.0,
      ticks.colour = 'gray50',
      title = bquote('Water Yield Gain Per Unit Expansion of Reservoir Storage Capacity ('~km^3~'/'~km^3~')'),
      title.position = "top"
    )) +
  ggplot2::coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = NA, color = NA),
                 plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                 panel.background = ggplot2::element_rect(fill = NA, color = NA),
                 legend.position = c(0.81, 0.15),
                 legend.justification = c('right', 'bottom'),
                 legend.box = 'horizontal',
                 legend.background = ggplot2::element_blank(),
                 legend.margin = ggplot2::margin(t = 10, b = 2),
                 legend.title = ggplot2::element_text(size = 8.5),
                 legend.text = ggplot2::element_text(angle = 0, size = 7.5,
                                                     margin = ggplot2::margin(t = 2)))

ggplot2::ggsave(file.path(save_path, plot_name),
                height = 6, width = 10, unit = 'in', dpi = 600)
