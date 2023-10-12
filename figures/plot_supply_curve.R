################################################################################
# Plot Supply Curves
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(ggnewscale)
library(patchwork)
library(paletteer)
library(snow)
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(forcats) # reorder 
library(rgcam)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
setwd(work.dir)

# set paths
figure.dir <- file.path(work.dir, 'outputs', 'supply_curve')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}
output.dir <- file.path(work.dir, 'outputs_glory')

# Read supply curves info
prj_ref <- rgcam::loadProject(
  file.path(output.dir, 'gcam_output_database', 'gcam_32regions_market_ref.dat')
)

#  Experiment
exp <- data.frame(run = c('20220807_21h29m', # A: climate impact reference MIROC-ESM-CHEM rcp6.0
                          '20230226_11h16m', # A + B: feedback on
                          '20230226_13h14m'), # A + C: feedback off
                  scenario = c('Climate',
                               'Climate_FB_ON',
                               'Climate_FB_OFF'))  

select_scenario <- 'Climate_FB_ON'

# ------------------------------------------------------------------------------
# Set scenarios names basin on run name, Mapping
# ------------------------------------------------------------------------------

gcam_basin_mapping <- rmap::mapping_gcambasins %>% 
  dplyr::left_join(rmap::mapping_tethys_grid_basin_region_country %>% 
                     dplyr::select(basinID, basinName) %>%
                     unique(),
                   by = c('subRegionMap' = 'basinName')) %>% 
  dplyr::mutate(market = paste0(subRegion, subRegion, '_water withdrawals'),
                resource = paste0(subRegion, '_water withdrawals')) %>% 
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

# Get Reference supply curves
avail <- prj_ref[['Reference']]$`resource supply curves` %>%
  dplyr::filter(subresource == 'runoff') %>%
  dplyr::mutate(scenario = 'Reference')
cost <- prj_ref[['Reference']]$`resource supply curve prices` %>%
  dplyr::filter(subresource == 'runoff') %>%
  dplyr::mutate(scenario = 'Reference')
supply_curve_orig <- avail %>% 
  dplyr::select(scenario, resource, year, grade, available = value) %>% 
  dplyr::left_join(
    cost %>%
      dplyr::select(scenario, resource, year, grade, extractioncost = value),
    by = c('scenario', 'resource', 'year', 'grade')) %>% 
  dplyr::left_join(gcam_basin_mapping, by = 'resource') %>% 
  dplyr::mutate(group = 'Original') %>% 
  dplyr::select(scenario, basin_id, basin_name, year, 
                grade, available, extractioncost, group) %>% 
  dplyr::filter(year >= 2020, year <= 2050)


# Read updated supply curves
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.dir, 'lp_outputs', run),
                       pattern = 'LP_supply_curve',
                       full.names = TRUE)
  scenario <- exp$scenario[exp$run == run]
  data <- lapply(f_list, read_files, 4, scenario)
  all_data <- c(all_data, data)
}

supply_curve <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::mutate(year = as.integer(year)) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('resource')) %>% 
  dplyr::mutate(group = 'Updated') %>% 
  dplyr::select(scenario, basin_id, basin_name, year, 
                grade, available, extractioncost, group) %>% 
  dplyr::bind_rows(supply_curve_orig)


# Read maxsubresource
all_data <- list()
for(run in exp$run){
  f_list <- list.files(file.path(output.dir, 'lp_outputs', run),
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

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
# Customize themes

theme <- ggplot2::theme(
  title = ggplot2::element_text(size = 16),
  panel.background = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(colour = 'black', fill = NA, linewidth = 1.2),
  panel.grid.major.x = ggplot2::element_line(colour = 'grey'),
  axis.text = ggplot2::element_text(colour = 'black', size = 14),
  axis.title = ggplot2::element_text(size = 14),
  strip.text = ggplot2::element_text(size = 14),
  strip.placement = 'outside',
  legend.position = 'right',
  legend.box = 'horizontal',
  legend.background = ggplot2::element_blank(),
  legend.key = ggplot2::element_blank(),
  legend.key.width = ggplot2::unit(1,'cm'),
  legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = 'pt'),
  legend.margin = ggplot2::margin(0, 0, 0, 0),
  legend.title = ggplot2::element_blank(),
  legend.text = ggplot2::element_text(size = 12),
  axis.ticks = ggplot2::element_line(linewidth = 1.2),
  plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0),
  plot.subtitle = ggplot2::element_text(size = 12)
)

pal_1 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'))
pal_2 <-  colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))
pal <- c(pal_1(6), pal_2(6), rep('black', 6))

# Lines ===

# Supply curves of all example basins c(79, 86, 87, 89, 131, 217)
if(T){
  
  library(ggh4x)
  df_plot <- supply_curve %>% 
    dplyr::filter(basin_id %in% c(79, 86, 87, 89, 131, 217), 
                  scenario %in% c('Climate_FB_ON', 'Climate_FB_OFF'),
                  year >= 2020 ) %>% 
    dplyr::bind_rows(supply_curve %>% 
                       dplyr::filter(basin_id %in% c(79, 86, 87, 89, 131, 217), 
                                     scenario %in% c('Reference'),
                                     year == 2050)) %>%
    dplyr::mutate(scenario = dplyr::case_when(scenario == 'Climate_FB_ON' ~ 'Feedback',
                                              scenario == 'Climate_FB_OFF' ~ 'No Feedback',
                                              scenario == 'Reference' ~ 'Reference/Climate Impact'),
                  basin_name = gsub('_', ' ', basin_name),
                  highlight = dplyr::case_when(year %in% c(2020, 2025, 2035, 2040, 2045) ~ 'other periods',
                                               year == 2030 ~ '2030',
                                               year == 2050 ~ '2050')) %>% 
    dplyr::filter(available <= 1)
  
  
  scales_y <- list(
    basin_name == 'Africa North West Coast' ~ scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02)),
    basin_name == 'Huang He' ~ scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05)),
    basin_name == 'Nile' ~ scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)),
    basin_name == 'California River' ~ scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)),
    basin_name == 'Indus' ~ scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02)),
    basin_name == 'Rio Lerma' ~ scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5))
  )
  
  file_name <- paste0('00_selected_basins_available lt 0p5.png')
  
  ggplot2::ggplot(df_plot, ggplot2::aes(
    x = available, y = extractioncost,
    group = interaction(scenario, year, highlight))) +
    ggplot2::geom_line(ggplot2::aes(color = scenario,
                                    size = highlight,
                                    linetype = highlight)) +
    ggplot2::facet_wrap(vars(basin_name), ncol = 2, scales = 'free_y') +
    ggh4x::facetted_pos_scales(y = scales_y) +
    ggplot2::coord_cartesian(xlim = c(0, 0.8)) +
    ggplot2::labs(x = 'Available Fraction of Total Runoff',
                  y = bquote('Water Price (1975 USD/'~m^3~')'),
                  title = NULL,
                  color = 'Scenario', size = 'Period', linetype = 'Period') +
    ggplot2::scale_color_manual(
      values = c('darkorange', 'darkseagreen', '#000000')
      ) +
    ggplot2::scale_size_manual(
      values = c(1, 1, 0.2)
    ) +
    ggplot2::scale_linetype_manual(
      values = c('longdash', 'solid', 'solid')
    ) +
    theme +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      strip.background = ggplot2::element_rect(fill = NA, size = 1.2),
      strip.text = ggplot2::element_text(size = 18),
      legend.position = 'bottom',
      legend.title = ggplot2::element_text(size = 12, face = 'bold'),
      legend.spacing = ggplot2::unit(2, 'cm')
    )
  
  ggplot2::ggsave(file.path(figure.dir, paste0(file_name, '_new.png')),
                  height = 10, width = 12, unit = 'in', dpi = 600)
}


# example of original supply curve from GCAM
df_plot <- supply_curve_orig %>% 
  dplyr::filter(scenario == 'Reference', basin_id == 217, year == 2020) %>% 
  dplyr::mutate(group = dplyr::case_when(grade == 'grade1' ~ 'origin',
                                         grade == 'grade10' ~ 'accessible point',
                                         grade == 'grade20' ~ 'max point',
                                         TRUE ~ 'interpolated point'))

ggplot2::ggplot(df_plot, ggplot2::aes(x = available, y = extractioncost)) +
  ggplot2::geom_line(linewidth = 1, color = 'black') +
  ggplot2::geom_point(ggplot2::aes(fill = group, size = group), shape = 21, color = 'black', 
                      stroke = 1) +
  ggplot2::scale_fill_manual(
    values = c('origin' = '#B6DBFF',
               'accessible point' = '#FFA319',
               'max point' = '#ABCD72',
               'interpolated point' = 'black'),
    guide = ggplot2::guide_legend(
      ncol = 4)) +
  ggplot2::scale_size_manual(
    values = c('origin' = 3.5,
               'accessible point' = 3.5,
               'max point' = 3.5,
               'interpolated point' = 2)
  ) +
  ggplot2::labs(x = 'Available Fraction',
                y = bquote('Price (1975 USD/'~m^3~')'),
                title = NULL) +
  theme +
  ggplot2::theme(
    legend.position = 'bottom'
  )

ggplot2::ggsave(file.path(figure.dir, '00_example_gcam_supply_curve.png'),
                height = 4.5, width = 7, unit = 'in', dpi = 300)


# Maps ===
# Fraction of total runoff that is cheaply available (Figure 1.3b)
# global map that classifies each basin's cost curve into bins
# showing how much surface water is cheaply accessible (<= $0.001/m3)
# Using Climate-FB-ON as example

# for single scenario and year
df <- supply_curve %>% 
  dplyr::filter(scenario == 'Reference', year == 2020) %>% 
  dplyr::select(basin_id, basin_name, available, extractioncost)

df <- supply_curve %>% 
  dplyr::filter(scenario == select_scenario, year == 2020) %>% 
  dplyr::select(basin_id, basin_name, available, extractioncost)
maxsubresource_df <- maxsubresource %>% 
  dplyr::filter(scenario == select_scenario, year == 2020)

df_plot <- data.frame()
for(id in unique(df$basin_id)) {
  df_basin <- df %>% 
    dplyr::filter(basin_id == id)
  
  val_maxsubresource <- maxsubresource_df$maxSubResource[
    maxsubresource_df$basin_id == id]
  if(length(val_maxsubresource) == 0){
    val_maxsubresource <- 0
  }
  
  cheap_fraction <- spline(x = df_basin$extractioncost, y = df_basin$available,
                           xout = 0.001, method = 'hyman')
  cheap_fraction_df <- data.frame(basin_id = id,
                                  basin_name = unique(df_basin$basin_name),
                                  cheap_fraction = min(cheap_fraction$y, 1)) %>% 
    dplyr::mutate(cheap_vol = cheap_fraction * val_maxsubresource)
  
  df_plot <- dplyr::bind_rows(df_plot, cheap_fraction_df)
  
}

# for year comparison
year_list <- c(2020, 2050)
cheap_price <- 0.001 # $/m3
i <- 1
df_plot <- data.frame()
for(yr in year_list){
  df <- supply_curve %>% 
    dplyr::filter(scenario == select_scenario, year == yr) %>% 
    dplyr::select(basin_id, basin_name, available, extractioncost)
  maxsubresource_df <- maxsubresource %>% 
    dplyr::filter(scenario == select_scenario, year == yr)
  
  cheap_df <- data.frame()
  for(id in unique(df$basin_id)) {
    df_basin <- df %>% 
      dplyr::filter(basin_id == id)
    
    val_maxsubresource <- maxsubresource_df$maxSubResource[
      maxsubresource_df$basin_id == id]
    if(length(val_maxsubresource) == 0){
      val_maxsubresource <- 0
    }
    
    cheap_fraction <- spline(x = df_basin$extractioncost, y = df_basin$available,
                             xout = cheap_price, method = 'hyman')
    cheap_fraction_df <- data.frame(basin_id = id,
                                    basin_name = unique(df_basin$basin_name),
                                    cheap_fraction = min(cheap_fraction$y, 1)) %>% 
      dplyr::mutate(cheap_vol = cheap_fraction * val_maxsubresource)
    
    cheap_df <- dplyr::bind_rows(cheap_df, cheap_fraction_df)
    
  }
  
  if(i == 1){
    df_plot <- cheap_df
  } else {
    df_plot <- df_plot %>% 
      dplyr::left_join(cheap_df, by = c('basin_id', 'basin_name'),
                       suffix = c(paste0('.', year_list[i-1]), 
                                  paste0('.', year_list[i])))
  }
  
  i <- i + 1
  
}

df_plot <- df_plot %>% 
  dplyr::mutate(cheap_fraction.diff = cheap_fraction.2050 - cheap_fraction.2020,
                cheap_vol.diff = cheap_vol.2050 - cheap_vol.2020)

save_path <- file.path(figure.dir, 'cheap_surface_water')
if(!dir.exists(save_path)){
  dir.create(save_path)
}


# plot as nicer maps
# prepare for sf data
cheap_fraction_sf <- basin_sf %>% 
  dplyr::left_join(df_plot, 
                   by = 'basin_id')

# customize the breaks
breaks_frac <- seq(0, 1, 0.1) # for cheap_fraction at a period
breaks_frac_diff <- seq(-0.4, 0.4, 0.05) # for abs diff of cheap_fraction
breaks_vol <- c(0, 5, 10, seq(20,100, 20), seq(150, 300, 50), seq(400, 600, 100)) # for cheap_vol
breaks_vol_diff <- c(-450,-300,-150,seq(-50,50,10),150,300,450) # for abs diff of cheap_vol between years

cheap_fraction_sf <- cheap_fraction_sf %>% 
  dplyr::mutate(cheap_fraction_group.2020 = cut(cheap_fraction.2020, breaks_frac, right = TRUE),
                cheap_fraction_group.2050 = cut(cheap_fraction.2050, breaks_frac, right = TRUE),
                cheap_fraction_group.diff = cut(cheap_fraction.diff, breaks_frac_diff, right = TRUE),
                cheap_vol_group.2020 = cut(cheap_vol.2020, breaks_vol, right = FALSE),
                cheap_vol_group.2050 = cut(cheap_vol.2050, breaks_vol, right = FALSE),
                cheap_vol_group.diff = cut(cheap_vol.diff, breaks_vol_diff, right = TRUE))

pal_n <- length(unique(cheap_fraction_sf$cheap_fraction_group.diff))

world +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_sf(data = cheap_fraction_sf,
                   ggplot2::aes(fill = cheap_fraction_group.diff),
                   color = 'grey30', linetype = 1, lwd = 0.25) +
  ggplot2::scale_fill_manual(
    # values = paste0(rev(paletteer_c("grDevices::Purple-Yellow", pal_n))),
    values = paletteer_c("ggthemes::Red-Blue Diverging", pal_n), # for abs diff
    guide = ggplot2::guide_colorsteps(
      barwidth = 20,
      barheight = 0.4,
      draw.llim = TRUE,
      frame.colour = 'gray50',
      direction = 'horizontal',
      label.position = 'bottom',
      ticks = TRUE,
      ticks.linewidth = 1.0,
      ticks.colour = 'gray50',
      show.limits = T,
      # title = 'Fraction of Low-cost ($0.001/m3) Renewable Water in 2050',
      title = 'Fraction Difference of Low-cost ($0.001/m3) Renewable Water from 2020 to 2050',
      # title = 'Low-cost ($0.001/m3) Renewable Water in 2050 (km3)',
      # title = 'Volume Difference of Low-cost ($0.001/m3) Renewable Water from 2020 to 2050 (km3)',
      title.position = "top"
    )) +
  ggplot2::coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = NA, color = NA),
                 plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                 panel.background = ggplot2::element_rect(fill = NA, color = NA),
                 legend.position = c(0.4, 0.15),
                 legend.justification = c('left', 'bottom'),
                 legend.box = 'horizontal',
                 legend.background = ggplot2::element_blank(),
                 legend.margin = ggplot2::margin(t = 10, b = 2),
                 legend.title = ggplot2::element_text(size = 8.5),
                 legend.text = ggplot2::element_text(angle = 0, size = 7.5,
                                                     margin = ggplot2::margin(t = 2)))

file_name <- paste0('cheap fraction_', select_scenario, '_2050_sf.png')
file_name <- paste0('cheap fraction_', select_scenario, '_2050-2020_sf.png')
file_name <- paste0('cheap volume_', select_scenario, '_2050_sf.png')
file_name <- paste0('cheap volume_', select_scenario, '_2050-2020_sf.png')
ggplot2::ggsave(file.path(save_path, file_name),
                height = 6, width = 10, unit = 'in', dpi = 600)
