################################################################################
# Plot historical sub-annual patterns
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(rmap)
library(paletteer)
library(tidyverse)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly

figure.dir <- file.path(dirname(work.dir), 'outputs', 'historical_supply_demand')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}

# historical monthly demand
data.workflow.dir <- 'zhao-etal_2023_gmd/workflow/data'
demand_rds <- file.path(work.dir, 
                        'wd_tethys_2005_2010_gcam5p3-stash_GFDL-ESM2M_rcp2p6_235basins.RDS')
demand <- readRDS(demand_rds)

# historical monthly runoff
xanthos.output.dir <- file.path(data.workflow.dir, 'hydrology_xanthos')
f_runoff <- file.path(xanthos.output.dir, 'Basin_runoff_km3permonth_watch+wfdei_1970_2010.csv')
runoff <- data.table::fread(f_runoff, header = T)


# Historical (Reference) surface water withdrawals and total water demand
# output from postprocess_gcam_output_water_ag.R
output_path <- file.path(work.dir, 'outputs_glory', 'gcam_proc')
exp_name <- 'gcam_32regions_waterag_allscenarios'
rdata_name <- paste0(exp_name, '.rds')
dataGCAM <- readRDS(file.path(output_path, exp_name, rdata_name))

# ------------------------------------------------------------------------------
# Set scenarios names basin on run name, Mapping
# ------------------------------------------------------------------------------

gcam_basin_mapping <- rmap::mapping_gcambasins %>% 
  dplyr::left_join(rmap::mapping_tethys_grid_basin_region_country %>% 
                     dplyr::select(basinID, basinName) %>%
                     unique(),
                   by = c('subRegionMap' = 'basinName')) %>% 
  dplyr::rename(basin_name = subRegionMap,
                basin_name_gcam = subRegion,
                basin_id = basinID) %>%  
  dplyr::mutate_all(function(x) gsub('Madasgacar', 'Madagascar', x))
gcam_basin_mapping$basin_id <- as.integer(gcam_basin_mapping$basin_id)


# ------------------------------------------------------------------------------
# Spatial SF Data
# ------------------------------------------------------------------------------
# Load RDS data
world <- readRDS(file.path(data.workflow.dir, 'world_basemap_robinson_projection.rds'))

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
# Process Data
# ------------------------------------------------------------------------------

# get 2010 monthly total demand ---
demand_total <- demand %>% 
  dplyr::group_by(subRegion, x) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(time = as.Date(paste0(x, '01'), format = '%Y%m%d'),
                year = as.integer(year(time)),
                month = as.integer(month(time))) %>% 
  dplyr::filter(year == 2010) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('subRegion' = 'basin_name')) %>% 
  dplyr::arrange(basin_id, month) %>% 
  dplyr::select(basin_id, basin_name = subRegion, month, value)


# get average monthly runoff ---
ts <- format(seq.Date(from = as.Date('1970/1/1'),
                      to = as.Date('2010/12/1'),
                      by = 'month'),
             format = '%Y%m')
names(runoff) <- c('basin_id', 'name', ts)

runoff_avg <- runoff %>% 
  tidyr::gather(key = 'yearmon', value = 'value', all_of(ts)) %>%
  dplyr::mutate(time = as.Date(paste0(yearmon, '01'), format = '%Y%m%d'),
                month = month(time),
                value = if_else(is.na(value), 0, value)) %>% 
  dplyr::group_by(basin_id, name, month) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('basin_id')) %>% 
  dplyr::select(basin_id, basin_name, month, value)


df_all <- demand_total %>% 
  dplyr::left_join(runoff_avg, by = c('basin_id', 'basin_name', 'month'),
                   suffix = c('.demand', '.runoff')) 

df_demand_met <- df_all%>% 
  dplyr::mutate(demand_met = dplyr::if_else(value.runoff > value.demand,
                                            value.demand,
                                            value.runoff)) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::summarise(demand_met = sum(demand_met),
                   demand_ann = sum(value.demand)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(demand_met_frac = demand_met / demand_ann)

df_socioeco_drought <- df_all%>% 
  dplyr::mutate(deficit = dplyr::if_else(value.runoff > value.demand,
                                            0,
                                            value.demand - value.runoff)) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::mutate(n = dplyr::if_else(deficit > 0, 1, 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(basin_id, basin_name) %>%
  dplyr::summarise(deficit = sum(deficit),
                   duration = sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(drought_intensity = dplyr::if_else(duration > 0, 
                                                   deficit / duration,
                                                   0),
                drought_intensity_log = ifelse(drought_intensity == 0,
                                               NA,
                                               log(drought_intensity)))


# Process GCAM output data ---
df_hist <- dataGCAM %>% 
  dplyr::filter(x <= 2015, scenario == 'Reference')
source <- df_hist %>% 
  dplyr::filter(param %in% 'waterWithdrawROGW') %>% 
  dplyr::select(subRegion, year = x, class1, value) %>% 
  tidyr::separate(col = class1, sep = ' ', into = c('subresource', 'grade')) %>% 
  dplyr::mutate(subresource = paste0(subresource, '_withdrawal')) %>% 
  dplyr::group_by(subRegion, subresource, year) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(value = dplyr::if_else(is.na(value), 0, value))

runoff <- df_hist %>% 
  dplyr::filter(param %in% c('watSupRunoffBasin')) %>% 
  dplyr::mutate(subresource = 'total_runoff') %>% 
  dplyr::group_by(subRegion, subresource, year = x) %>% 
  dplyr::summarise(value = mean(value)) %>% 
  dplyr::ungroup()

df_hist_mean <- dplyr::bind_rows(source, runoff) %>% 
  tidyr::spread(key = 'subresource', value = 'value') %>% 
  dplyr::mutate(total_withdrawal = groundwater_withdrawal + runoff_withdrawal,
                historical_scarcity = dplyr::if_else(
                  total_runoff == 0,
                  0,
                  total_withdrawal / total_runoff),
                historical_supply_mix = dplyr::if_else(
                  total_withdrawal == 0,
                  0,
                  runoff_withdrawal / total_withdrawal)) %>% 
  dplyr::filter(year == 2015)

df_hist_mean <- gcam_basin_mapping %>% 
  dplyr::left_join(df_hist_mean, by = c('basin_name_gcam' = 'subRegion')) %>% 
  dplyr::select(basin_name, year, historical_supply_mix, historical_scarcity)

df_hist_mean_sf <- basin_sf %>% 
  dplyr::left_join(df_hist_mean, by = 'basin_name')

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

# Theme
fig_height <-  4.5 # in
fig_width <- 7 # in
theme <- ggplot2::theme(
  title = ggplot2::element_text(size = 16),
  panel.background = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(colour = 'black', fill = NA, size = 1.2),
  panel.grid.major.x = ggplot2::element_line(colour = 'grey'),
  axis.text = ggplot2::element_text(colour = 'black', size = 14),
  axis.title = ggplot2::element_text(size = 14),
  strip.text = ggplot2::element_text(size = 14),
  strip.placement = 'outside',
  legend.position = 'bottom',
  legend.box = 'horizontal',
  legend.background = ggplot2::element_blank(),
  legend.key = ggplot2::element_blank(),
  legend.key.width = ggplot2::unit(1.3,'cm'),
  legend.box.margin = ggplot2::margin(t=0, r=0,b=0,l=0,unit='pt'),
  legend.margin = ggplot2::margin(0,0,0,0),
  legend.title = ggplot2::element_blank(),
  legend.text = ggplot2::element_text(size = 12),
  axis.ticks = ggplot2::element_line(size = 1.2),
  plot.margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0),
  plot.subtitle = ggplot2::element_text(size = 12)
)



# Plot 1: Fraction of demand being met by runoff (Figure 1.2)

socioeco_drought_sf <- basin_sf %>% 
  dplyr::left_join(df_socioeco_drought, 
                   by = c('basin_id', 'basin_name'))

breaks <- seq(0, 1, 0.1)


world +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_sf(data = socioeco_drought_sf,
                   ggplot2::aes(fill = drought_intensity_log),
                   color = 'grey30', linetype = 1, lwd = 0.25) +
  ggplot2::scale_fill_gradientn(
    colours = rev(c(paletteer_c("grDevices::OrRd", length(breaks)-1), 'white')),
    na.value = 'white',
    # breaks = breaks,
    guide = ggplot2::guide_colorbar(
      barwidth = 0.7,
      barheight = 8,
      draw.llim = TRUE,
      frame.colour = 'gray50',
      direction = 'vertical',
      label.position = 'right',
      ticks = TRUE,
      ticks.linewidth = 1.0,
      ticks.colour = 'gray50',
      title = 'Log Transformed\nSocioeconomic\nDrought\nIntensity',
      title.position = "top"
    )) +
  ggplot2::coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = NA, color = NA),
                 plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                 panel.background = ggplot2::element_rect(fill = NA, color = NA),
                 legend.position = c(0.13, 0.7),
                 legend.justification = c('left', 'top'),
                 legend.box = 'horizontal',
                 legend.background = ggplot2::element_blank(),
                 legend.margin = ggplot2::margin(t = 10, b = 2),
                 legend.title = ggplot2::element_text(size = 8.5),
                 legend.text = ggplot2::element_text(angle = 0, size = 7.5,
                                                     margin = ggplot2::margin(t = 2)))
plot_name <- 'historical log transformed socioeconomic drought intensity.png'
ggplot2::ggsave(file.path(figure.dir, plot_name),
                height = 6, width = 10, unit = 'in', dpi = 600)


# Plot 2: deficit between inflow and demand
save_path <- file.path(figure.dir, 'runoff_demand_deficit')
if(!dir.exists(save_path)){
  dir.create(save_path)
}

select_basins <- c(217, 86, 79, 131, 87, 89)

for (id in select_basins){
  
  file_name <- paste0(id, ' - ',
                      gcam_basin_mapping$basin_name[gcam_basin_mapping$basin_id==id])
  
  df_basin <- df_all %>% 
    dplyr::filter(basin_id == id) %>% 
    dplyr::rename(a = value.runoff,
                  b = value.demand,
                  x = month) %>% 
    dplyr::mutate(basin_name = gsub('_', ' ', basin_name))
  
  title <- unique(df_basin$basin_name)
  
  df_basin_gather <- df_basin %>% 
    tidyr::gather(key = 'f', value = 'y', c('a', 'b'))
  bounds <- df_basin %>% 
    dplyr::mutate(ymax = pmax(a, b),
                  ymin = pmin(a, b),
                  fill = a > b)
  intervals <- bounds %>% 
    dplyr::filter(ymax > ymin) %>% 
    dplyr::select(-a, -b)
  intersections <- bounds %>% 
    dplyr::mutate(lag_fill = lag(fill), lead_fill = lead(fill)) %>%
    dplyr::filter(ymax == ymin) %>%
    dplyr::select(-a, -b, -fill) %>%
    tidyr::pivot_longer(lag_fill:lead_fill, names_to = NULL, values_to = "fill") %>%
    dplyr::filter(!is.na(fill)) %>%
    dplyr::distinct()
  
  other_intersections <- bounds %>%
    dplyr::transmute(
      x1 = x,       y1 = a,
      x2 = lead(x), y2 = lead(a),
      x3 = x,       y3 = b,
      x4 = lead(x), y4 = lead(b)
    ) %>%
    dplyr::filter(((y1 > y3) & (y2 < y4)) | ((y1 < y3) & (y2 > y4))) %>%  # only rows where an intersection occurs between two x
    dplyr::mutate(
      d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4),  # denominator
      u = x1 * y2 - y1 * x2,
      v = x3 * y4 - y3 * x4,
      x = (u * (x3 - x4) - v * (x1 - x2)) / d,
      y = (u * (y3 - y4) - v * (y1 - y2)) / d
    ) %>%
    dplyr::select(x, ymax = y, ymin = y)
  
  ribbons <- dplyr::bind_rows(
    intervals,
    intersections,
    dplyr::mutate(other_intersections, fill = TRUE),
    dplyr::mutate(other_intersections, fill = FALSE)
  ) %>%
    dplyr::arrange(x)
  
  ggplot2::ggplot(data = df_basin_gather) +
    ggplot2::geom_line(ggplot2::aes(x, y, linetype = f), size = 1.2) +
    ggplot2::geom_ribbon(data = ribbons,
                         ggplot2::aes(x, ymin = ymin, ymax = ymax, fill = fill)) +
    ggplot2::scale_fill_manual(values = c('FALSE' = 'firebrick4', 'TRUE' = 'cadetblue4'),
                               labels = c('Deficit', 'Surplus')) +
    ggplot2::scale_linetype_manual(values = c('solid', 'dashed'),
                                   labels = c('Natural Inflow', 'Water Demands')) +
    ggplot2::labs(x = NULL, 
                  y = bquote('Volume ('~km^3~'/month)'),
                  title = title) +
    ggplot2::scale_x_continuous(breaks = seq(1, 12, 1)) +
    theme

  
  ggplot2::ggsave(file.path(save_path, paste0(file_name, '.png')),
                  height = 4.5, width = 7, unit = 'in', dpi = 300)
}
