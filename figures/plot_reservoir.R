################################################################################
# Plot Global Reservoir Purposes
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(rmap)
library(sf)
library(giscoR)
library(paletteer)
library(ggnewscale)
library(rgdal)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
# set paths
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
setwd(work.dir)

output.workflow.dir <- 'zhao-etal_2023_gmd/workflow/outputs'

data.workflow.dir <- 'zhao-etal_2023_gmd/workflow/data'

figure.dir <- file.path(work.dir, 'outputs', 'reservoirs')


# ------------------------------------------------------------------------------
# Load data and Prepare Data
# ------------------------------------------------------------------------------
# Load World Basemap
world <- readRDS(file.path(data.workflow.dir, 'world_basemap_robinson_projection.rds'))

data <- readRDS(file.path(output.workflow.dir, 'LP_inputs_reservoir_20221222.rds'))
dam_info <- data$dam_info

capacity <- dam_info %>% 
  dplyr::mutate(purpose = dplyr::if_else(purpose != 'hydropower',
                                         'non_hydropower', purpose)) %>% 
  dplyr::group_by(basin_id, basin_name, purpose) %>% 
  dplyr::summarise(value = sum(CAP_KM3)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = purpose, values_from = value) %>% 
  tidyr::replace_na(list(hydropower = 0, non_hydropower = 0)) %>% 
  dplyr::mutate(total = hydropower + non_hydropower)

# Get capacity layers
capacity <- capacity %>% 
  dplyr::mutate(hydropower_top = dplyr::if_else(hydropower < non_hydropower, hydropower, NULL),
                hydropower_bottom = dplyr::if_else(hydropower >= non_hydropower, hydropower, NULL),
                non_hydropower_top = dplyr::if_else(non_hydropower < hydropower, non_hydropower, NULL),
                non_hydropower_bottom = dplyr::if_else(non_hydropower >= hydropower, non_hydropower, NULL)) %>% 
  dplyr::mutate_at(c('hydropower_top',
                     'hydropower_bottom',
                     'non_hydropower_top',
                     'non_hydropower_bottom'),
                   ~ na_if(., 0) )


# ------------------------------------------------------------------------------
# Spatial Referencing
# ------------------------------------------------------------------------------
# set epsg code 4326: https://spatialreference.org/ref/epsg/wgs-84/
epsg_code <- 4326

# create basin simple feature
basin_sf <- rmap::mapGCAMBasins %>% 
  dplyr::select(basin_id = subRegionAlt,
                basin_name = subRegion,
                geometry) %>% 
  st_transform(epsg_code)

# basin dataframe
basin_df <- basin_sf %>% 
  sf::st_drop_geometry()

# get capacity sf
capacity_sf <- basin_sf %>% 
  dplyr::left_join(capacity %>% dplyr::select(basin_id, total), 
                   by = c('basin_id'))
  
# get bubble position - centroids
symbol_pos <- sf::st_centroid(basin_sf, of_largest_polygon = T)


# add basin point

capacity_point <- symbol_pos %>% 
  dplyr::left_join(capacity %>% 
                     dplyr::select(basin_id, 
                                   hydro_top = hydropower_top, 
                                   non_hydro_top = non_hydropower_top,
                                   hydro_bottom = hydropower_bottom, 
                                   non_hydro_bottom = non_hydropower_bottom,
                                   total), 
                   by = c('basin_id')) %>% 
  tidyr::gather(key = 'group', value = 'value', 
                c('hydro_top', 'non_hydro_top', 'hydro_bottom', 'non_hydro_bottom')) %>% 
  dplyr::arrange(match(group, c('hydro_bottom', 'non_hydro_bottom', 'hydro_top', 'non_hydro_top'))) %>%
  dplyr::mutate(group = factor(group, levels = c('hydro_bottom', 'non_hydro_bottom', 'hydro_top', 'non_hydro_top')))



# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

# plot bubble map
world +
  ggnewscale::new_scale_fill() +
# ggplot2::ggplot() +
  ggplot2::geom_sf(data = capacity_sf,
                   ggplot2::aes(fill = total),
                   color = 'grey30', linetype = 1, lwd = 0.25) +
  ggplot2::scale_fill_gradientn(colours = hcl.colors(3, "GnBu", rev = TRUE),
                                # labels = label_fun,
                                n.breaks = 8,
                                guide = ggplot2::guide_colorsteps(
                                  barwidth = 15,
                                  barheight = 0.5,
                                  show.limits = T,
                                  frame.colour = 'gray50',
                                  ticks = TRUE,
                                  ticks.linewidth = 1.0,
                                  ticks.colour = 'gray50',
                                  title = bquote('Basin Total Storage Capacity '(italic(km^3))),
                                  title.position = "top")) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_sf(data = capacity_point, 
                   ggplot2::aes(size = value, fill = group), show.legend = T,
                   pch = 21,  col = "grey20", alpha = 0.7) +

  ggplot2::scale_size(
    range = c(0, 13),
    guide = ggplot2::guide_legend(direction = "horizontal",
                                  nrow = 1,
                                  label.position = "bottom")) +
  ggplot2::scale_fill_manual(values = c('hydro_top' = 'yellow',
                                        'non_hydro_top' = 'firebrick3',
                                        'hydro_bottom' = 'yellow',
                                        'non_hydro_bottom' = 'firebrick3'),
                             breaks = c('hydro_top', 'non_hydro_top'),
                             labels = c('hydro_top' = 'Hydropower', 
                                        'non_hydro_top' = 'Non Hydropower')) +
  ggplot2::guides(
    size = ggplot2::guide_legend(
      title = bquote('Storage Capacity '(italic(km^3))),
      title.position = 'top',
      order = 2),
    fill = ggplot2::guide_legend(
      title = 'Reservoir Category',
      title.position = 'top',
      ncol = 1,
      override.aes = list(size = 3),
      order = 1)
    ) +
  ggplot2::coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = 'white', color = NA),
                 panel.background = ggplot2::element_rect(fill = 'white', color = NA),
                 legend.position = "bottom",
                 legend.spacing = ggplot2::unit(1.0, 'cm'),
                 legend.margin = ggplot2::margin(t = 10, b = 2),
                 legend.title = ggplot2::element_text(size = 10),
                 legend.text = ggplot2::element_text(angle = 0,
                                                     margin = ggplot2::margin(t = 5)))


ggplot2::ggsave(file.path(figure.dir, 'global_reservoirs.png'),
                height = 6, width = 10, unit = 'in', dpi = 300)

