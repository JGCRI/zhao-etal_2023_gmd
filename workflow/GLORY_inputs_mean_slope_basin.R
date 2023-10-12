################################################################################
# Calculate average slope for each basin
# Contact: Mengqi Zhao (mengqi.zhao@pnnl.gov)
# Last Update: 2022-08
################################################################################


# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(rmap)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/workflow' # update correspondingly
setwd(work.dir)
data.dir <- file.path(work.dir, 'data')
output.dir <- file.path(work.dir, 'outputs')

data_reservoir <- readRDS(file.path(work.dir, 'outputs', 'LP_inputs_reservoir_20221222.rds'))

# basin id and name
basin_rmap <- data.frame(basin_id = rmap::mapGCAMBasins$subRegionAlt,
                         basin_name = rmap::mapGCAMBasins$subRegion)

# use 50km slope (about 0.5 degree grid)
slope_remapped_file <- file.path(data.dir, 'slope_remapped.csv')
slope_remapped <- data.table::fread(slope_remapped_file)


# Get the mean slope for the expansion grids
expan_category_grid <- data_reservoir$expan_category_grid
slope_basin <- slope_remapped %>% 
  dplyr::left_join(expan_category_grid, by = c('lat', 'lon')) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::mutate(expan_grid = dplyr::if_else(expan_category > 0, 1, 0),
                count = sum(expan_grid),
                slope_basin = mean(value)) %>% 
  dplyr::summarise(slope = dplyr::if_else(count > 0,
                                          mean(value[expan_grid == 1]),
                                          slope_basin)) %>% 
  dplyr::ungroup() %>% 
  unique()


# write output
write.csv(x = slope_basin,
          file = file.path(output.dir, 'LP_mean_slope_basin.csv'),
          row.names = F)
