################################################################################
# Plot Surface water Demand vs Groundwater Demand
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(gcamextractor)
library(rchart)
library(rmap)
library(paletteer)
library(ggplot2)
library(jgcricolors)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/figures' # update correspondingly
output.dir <- file.path(work.dir, 'outputs_glory')
db_path <- file.path(output.dir, 'gcam_output_database')

figure.dir <- file.path(work.dir, 'outputs', 'water_withdrawal')
if(!dir.exists(figure.dir)){
  dir.create(figure.dir)
}

# set experiment name
exp_name <- 'gcam_32regions_waterag_allscenarios'
rdata_name <- paste0(exp_name, '.rds')


exp <- data.frame(
  database = c('database_basexdb_default', 'database_basexdb_reference',
               'database_basexdb_feedback_on', 'database_basexdb_feedback_off'),
  scenario_orig = c('Reference', 
                    'Reference', 'feedback_on', 'feedback_off'),
  scenario_new = c('Reference',
                   'Climate', 'Climate_FB_ON', 'Climate_FB_OFF')
)

paramsSelect_water <- c('waterWithdrawROGW')
paramsSelect_i <- c(paramsSelect_water)



# ------------------------------------------------------------------------------
# Process GCAM output
# ------------------------------------------------------------------------------
dataGCAM <- tibble::tibble()
if(file.exists(file.path(output_path, exp_name, rdata_name))){
  
  dataGCAM <- readRDS(file.path(output_path, exp_name, rdata_name))
  
} else {
  
  for(db in 1:length(exp$database)){
    
    dataGCAM_temp <- gcamextractor::readgcam(
      gcamdatabase = file.path(db_path, exp$database[db]),
      regionsSelect = NULL,
      paramsSelect = paramsSelect_i,
      reReadData = T,
      scenOrigNames = exp$scenario_orig[db],
      scenNewNames = exp$scenario_new[db],
      nameAppend = paste0('_', exp$scenario_new[db]),
      folder = file.path(work.dir, 'gcam_proc', exp_name),
      save = F)
    
    dataGCAM <- dataGCAM %>%
      dplyr::bind_rows(dataGCAM_temp$dataAll)%>%
      dplyr::select(scenario, region, subRegion, param, class1, class2, x, value, aggregate)
  }
  
  saveRDS(dataGCAM, file.path(output_path, exp_name, paste0(exp_name, '.rds')))
}

df <- dataGCAM %>% 
  dplyr::filter(x %in% seq(2015, 2050, 5))

df_hist <- dataGCAM %>% 
  dplyr::filter(x <= 2015, scenario == 'Reference')

df_select <- dataGCAM %>% 
  dplyr::filter(x %in% seq(2020, 2050, 5), 
                scenario %in% c('Climate', 'Climate_FB_ON', 'Climate_FB_OFF'))

# ------------------------------------------------------------------------------
# Set scenarios names basin on run name, Mapping
# ------------------------------------------------------------------------------
# Plot comparison of parameters across scenarios on maps
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
scen_order <- c('Reference',
                'Climate',
                'Climate_FB_ON',
                'Climate_FB_OFF')


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

# Customize theme
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
  legend.key.width = ggplot2::unit(1,'cm'),
  legend.box.margin = ggplot2::margin(t=0, r=0,b=0,l=0,unit='pt'),
  legend.margin = ggplot2::margin(0,0,0,0),
  legend.title = ggplot2::element_blank(),
  legend.text = ggplot2::element_text(size = 12),
  axis.ticks = ggplot2::element_line(size = 1.2),
  plot.margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0),
  plot.subtitle = ggplot2::element_text(size = 12)
)


# runoff withdrawal with Groundwater ---------------------
# plot selected basins for the paper c(79, 86, 87, 89, 131, 217)
df_plot <- df %>% 
  dplyr::filter(scenario %in% c("Reference", "Climate", "Climate_FB_ON", "Climate_FB_OFF"), 
                x %in% seq(2020, 2050, 5),
                param %in% 'waterWithdrawROGW') %>%
  dplyr::select(scenario, subRegion, year = x, class1, value) %>% 
  tidyr::separate(col = class1, sep = ' ', into = c('subresource', 'grade')) %>% 
  dplyr::group_by(scenario, subRegion, subresource, year) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('subRegion' = 'basin_name_gcam')) %>% 
  dplyr::select(-subRegion) %>% 
  dplyr::rename(source = subresource) %>% 
  dplyr::mutate(scenario = dplyr::case_when(scenario == 'Reference' ~ 'Reference',
                                            scenario == 'Climate' ~ 'Climate Impact',
                                            scenario == 'Climate_FB_ON' ~ 'Feedback',
                                            scenario == 'Climate_FB_OFF' ~ 'No Feedback'),
                basin_name = gsub('_', ' ', basin_name),
                basin_name = gsub('North West', 'NW', basin_name),
                source = dplyr::case_when(source == 'runoff' ~ 'Suface Water',
                                          source == 'groundwater' ~ 'Groundwater')) %>% 
  dplyr::filter(basin_id %in% c(79, 86, 87, 89, 131, 217)) %>% 
  dplyr::mutate(across(scenario, factor, 
                       levels = c('Reference', 'Climate Impact', 'No Feedback', 'Feedback')))

file_name <- '00_selected basins'

ggplot2::ggplot(data = df_plot) +
  ggplot2::geom_bar(ggplot2::aes(x = year, y = value, fill = source), 
                    stat = 'identity', color = 'black') +
  ggplot2::facet_grid(rows = vars(basin_name), cols = vars(scenario), scales = 'free_y') +
  ggplot2::scale_fill_manual(values = c('Suface Water'='#6CA9C3', 
                                        'Groundwater' = '#CCFFFF')) +
  ggplot2::labs(title = NULL,
                x = NULL,
                y = bquote('Water Withdrawals ('~km^3~'/year)') )+
  theme +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = 'gray30', 
                                                          colour = 'gray30'),
                 strip.text = ggplot2::element_text(colour = 'white'))

ggplot2::ggsave(file.path(figure.dir, paste0(file_name, '.png')),
                height = 12, width = 12, unit = 'in', dpi = 600)
