################################################################################
# For Linear Programming Model                                                 
# Maximum Reservoir Storage Capacity Potential                                 
# Author: Mengqi Zhao                                            
# Email: mengqi.zhao@pnnl.gov                                     
# Last Update: 2021-02-25                                        
################################################################################

# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(tidyr)
library(data.table)
library(sf) # geocomputation
library(rgis)
library(rmap)
library(rgdal)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(snow) # for parallel
library(minpack.lm) # non linear model

# display options of display of numbers and other data in a tibble
options(pillar.sigfig = 5)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/workflow' # update correspondingly
setwd(work.dir)

# set paths
data.dir <- file.path(work.dir, 'data')
output.dir <- file.path(work.dir, 'outputs')


# basin and grid mapping
basin_rmap <- data.frame(basin_id = rmap::mapGCAMBasins$subRegionAlt,
                         basin_name = rmap::mapGCAMBasins$subRegion)
gcam_basin_mapping <- rmap::mapping_gcambasins %>% 
  dplyr::rename(basin_name = subRegionMap,
                gcam_basin_name = subRegion) %>% 
  # correct Madagarscar
  dplyr::mutate(basin_name = gsub('Madasgacar', 'Madagascar', basin_name),
                gcam_basin_name = gsub('Madasgacar', 'Madagascar', gcam_basin_name))
grid_basin_mapping <- rmap::mapping_tethys_grid_basin_region_country %>% 
  dplyr::select(lat, lon, basinID, basinName) %>% 
  dplyr::rename(basin_id = basinID,
                basin_name = basinName) %>% 
  # correct Madagarscar
  dplyr::mutate(basin_name = gsub('Madasgacar', 'Madagascar', basin_name))

# HydroLAKES Dataset (including GranD and GLWD and it is newer)
hydrolakes <-  data.table::fread(
  file.path(data.dir, 'HydroLAKES_to_xanthos.csv'),
  header = TRUE)

# remapped hydroLakes
# default volume unit is million cubic meter
# lake type: 1 - lake, 2 - reservoir, 3 - lake control (natural lake with regulation structure)
hydrolakes_georef <- hydrolakes %>% 
  dplyr::mutate(LAKE_VOL_KM3 = Vol_total / 10^3,
                CAP_KM3 = Vol_res / 10^3) %>% 
  dplyr::select(xanthos_id, lat, lon, basin_id, basin_name,
                GRID_AREA_KM2 = GRID_AREA_SKM,
                GRAND_ID = Grand_id,
                LAKE_TYPE = Lake_type,
                LAKE_AREA_KM2 = Lake_area, 
                LAKE_VOL_KM3, 
                CAP_KM3)


# Read georeferenced and remapped-to-xanthos GranD reservoirs (remapped to xanthos grid)
# I created a new remapped GranD database based on the latest GranD_dams_v1_3.csv
# The new remapped file is called GranD_v1.3_remap_to_xanthos.csv
# File is created from remap_to_xanthos.py
# This file includes all the purposes so in total there are 7320 reservoirs included
# (Originally, we were using georeferencing output from Guta: GranD_Reservoirs_Georefrenced_Remapped.csv)
grand_georef <- data.table::fread(
  file.path(data.dir, 'GranD_v1.3_remap_to_xanthos.csv'), header = TRUE) %>% 
  dplyr::mutate(CAP_KM3 = dplyr::if_else(CAP_MCM == -99, 0, CAP_MCM / 1000),
                AREA_KM2 = dplyr::if_else(AREA_SKM == -99, 0, AREA_SKM),
                purpose = case_when(purpose == 0 ~ 'other',
                                    purpose == 1 ~ 'hydropower',
                                    purpose == 2 ~ 'irrigation',
                                    purpose == 3 ~ 'flood control')) 
unique(grand_georef$MAIN_USE)

# combine HydroLAKES and GranD
# By individual lakes and dams
dam_info <- hydrolakes_georef %>% 
  dplyr::left_join(grand_georef %>% 
                     dplyr::select(GRAND_ID, AREA_KM2, MAIN_USE, purpose),
                   by = c('GRAND_ID')) %>% 
  tidyr::replace_na(list(AREA_KM2 = 0,
                         MAIN_USE = 'None',
                         purpose = 'none'))
dam_purpose_grid <- dam_info %>% 
  dplyr::select(xanthos_id, purpose, CAP_KM3) %>% 
  dplyr::group_by(xanthos_id, purpose) %>%
  dplyr::summarise(CAP_KM3 = sum(CAP_KM3)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(xanthos_id) %>% 
  dplyr::mutate(MAX_CAP = max(CAP_KM3),
                purpose_main = unique(purpose[CAP_KM3 == MAX_CAP])) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(xanthos_id, purpose_main) %>% 
  unique()

# aggregate to xanthos grid
dam_grid_info <- dam_info %>% 
  dplyr::select(xanthos_id, lat, lon, basin_id, basin_name, GRID_AREA_KM2,
                AREA_KM2, CAP_KM3, LAKE_AREA_KM2, LAKE_VOL_KM3) %>% 
  dplyr::group_by(basin_id, basin_name, xanthos_id, lat, lon) %>% 
  dplyr::summarise(AREA_KM2 = sum(AREA_KM2),
                   CAP_KM3 = sum(CAP_KM3),
                   LAKE_AREA_KM2 = sum(LAKE_AREA_KM2),
                   LAKE_VOL_KM3 = sum(LAKE_VOL_KM3),
                   GRID_AREA_KM2 = mean(GRID_AREA_KM2)) %>% 
  dplyr::ungroup()


# ------------------------------------------------------------------------------
# Global Lakes and Wetlands
# ------------------------------------------------------------------------------
# GLWD Level 3 dataset
# Including type of each grid (0.008333)
glwd.file <- file.path(data.dir, 'GLWD_level3', 'glwd_3')
glwd_raster <- raster::raster(file.path(glwd.file, 'w001001.adf'))
glwd_raster <- raster::calc(glwd_raster, fun=function(x){if_else(x<=3 & x>1, 1, 0)})

# aggregate raster from 30arc second to 0.5 degree
glwd_raster_0p5 <- raster::aggregate(glwd_raster, fact = 3600*0.5/30, fun = mean)
plot(glwd_raster_0p5, col = rainbow(2))

raster_base <- raster::raster(resolution = 0.5,
                              xmn = -180, xmx = 180, ymn = -56, ymx = 84,
                              crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
glwd_raster_0p5resample <- raster::resample(glwd_raster_0p5, raster_base, method = 'ngb')
plot(glwd_raster_0p5resample, col = rainbow(2))

glwd_df <- as.data.frame(glwd_raster_0p5resample, xy = TRUE, na.rm = TRUE) %>% 
  dplyr::rename(water = layer,
                lon = x,
                lat = y)

# ------------------------------------------------------------------------------
# Population
# ------------------------------------------------------------------------------
# get population netcdf dir. Population is in capita
ssp <- 'ssp2'
year <- 2010
folder.name <- paste('popdynamics-1-8th-pop-base-year-projection-ssp-2000-2100-rev01-proj', ssp, 'netcdf', sep = '-')
pop.dir <- file.path(data.dir, 'Population', folder.name, toupper(ssp), 'Total', 'NetCDF')
pop.file <- file.path(pop.dir, paste0(ssp, '_', year, '.nc'))

pop_raster <- rgis::import_ncdf_to_raster(pop.file)
pop_raster <- subset(pop_raster, 1)
grid_area_raster <- raster::area(pop_raster) # km2
pop_density_raster <- pop_raster/grid_area_raster # capita/km2


# Aggregate raster from 0.125 degree to 0.5 degree
pop_raster_0p5 <- raster::aggregate(pop_density_raster, fact = 4, fun = mean)
raster::area(pop_raster_0p5)
raster_base <- raster::raster(resolution = 0.5,
                              xmn = -180, xmx = 180, ymn = -56, ymx = 84,
                              crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
pop_raster_0p5resample <- raster::resample(pop_raster_0p5, raster_base, method = 'bilinear')

# Calculate population density (capita/km2)
pop_density_df <- as.data.frame(pop_raster_0p5resample, xy = TRUE, na.rm = TRUE) %>% 
  dplyr::rename(pop_density = layer,
                lon = x,
                lat = y)
  
            
# ------------------------------------------------------------------------------
# Protected Areas
# ------------------------------------------------------------------------------
# WDPA shapefile path
version <- 'WDPA_Feb2022_Public_shp'
wdpa.dir <- file.path(data.dir, 'WDPA', version)
shp_list <- list.files(wdpa.dir, pattern = '\\polygons.shp$', full.names = TRUE, recursive = TRUE)

# Import WDPA shapfiles
# wdpa_raster_file <- 'wdpa_raster_0p5degree.rds' # The 0.5 degree process turns out to lose may Protected Areas
wdpa_raster_file <- 'wdpa_raster_0p125degree.rds'
if(file.exists(wdpa_raster_file)){
  wdpa_raster <- readRDS(wdpa_raster_file)
} else {
  shp_files <- lapply(shp_list, rgis::import_shapefile)
  wdpa_poly <- sf::st_as_sf(data.table::rbindlist(shp_files))
  wdpa_poly$VALUE <- 1
  
  if(grepl(wdpa_raster_file, pattern='0p5')){
    # 0.5 degree raster (same as rmap crs and general demeter gridcells lat lon)
    raster_base <- raster::raster(resolution = 0.5,
                                  xmn = -180, xmx = 180, ymn = -56, ymx = 84,
                                  crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  } else if(grepl(wdpa_raster_file, pattern='0p125')){
    # 0.125 degree raster
    raster_base <- raster::raster(ncols = 2880, nrows = 1117,
                                  xmn = -180, xmx = 180, ymn = -55.875, ymx = 83.75,
                                  crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  }

  # Rasterize WDPA polygon data. Take about 10 min
  wdpa_raster <- raster::rasterize(x = wdpa_poly,
                                   y = raster_base,
                                   field = 'VALUE',
                                   fun = 'mean')
  saveRDS(wdpa_raster, wdpa_raster_file)
}


# Aggregate raster from 0.125 degree to 0.5 degree.
wdpa_raster_0p5 <- raster::aggregate(wdpa_raster, fact = 4)
wdpa_raster_0p5resample <- raster::resample(wdpa_raster_0p5, raster_base, method = 'ngb')
raster::area(wdpa_raster_0p5resample)

# Get data frame for wdpa
wdpa_df <- as.data.frame(wdpa_raster_0p5resample, xy = TRUE, na.rm = TRUE) %>% 
  dplyr::rename(protected_area = layer,
                lon = x,
                lat = y) 

# ------------------------------------------------------------------------------
# Land Use
# ------------------------------------------------------------------------------
# Import landuse change from Demeter
demeter_output <- file.path(data.dir, 'landuse_demeter')
# GCAM v5.3
select_folders <- c('gcam5p3-stash_GFDL-ESM2M_rcp2p6_Reference_2021-05-10_16h55m37s')
select_output <- paste(demeter_output, select_folders, sep = '/')
output_files <- list.files(select_output, pattern = '0p5deg_2010', recursive = TRUE, full.names = TRUE)


# Function to read before run
read_files <- function(file, header = TRUE, sep = 'auto', ...){
  data <- data.table::fread(file, header = header, sep = sep, ...)
  filename <- basename(file[1])
  scenario <- tolower(strsplit(basename(filename), '\\_|\\.')[[1]][2])
  year <- strsplit(basename(filename), '\\_|\\.')[[1]][4]
  data$year <- year
  data$scenario <- scenario
  # data$file_type <- filename
  return(data)
}

# basename is a function to choose file by name
cl <- snow::makeSOCKcluster(4)
land_use <- snow::parLapply(cl, output_files, read_files)
stopCluster(cl)

col_gather <- c('water', 'forest', 'shrub', 'grass',
                'urban', 'snow', 'sparse', 'corn_irr',
                'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
                'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
                'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
                'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
                'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
                'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
                'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')

land_use <- data.table::rbindlist(land_use, idcol = FALSE) %>% 
  dplyr::rename(gridcode = pkey_0p5_deg, lon = longitude, lat = latitude) %>%
  tidyr::gather(key = 'crop', value = 'value', col_gather)

crop_land <- c('corn_irr',
               'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
               'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
               'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 
               'biomass_grass_irr','biomass_tree_irr',
               'corn_rfd',
               'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
               'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
               'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
               'biomass_grass_rfd', 'biomass_tree_rfd')

# Calculate exclusion zones. 0 - non-exclusion, 1 - exlusion zone
land_use_df <- land_use %>% 
  dplyr::mutate(zone = dplyr::if_else(crop %in% crop_land, 'crop', 'non_crop')) %>% 
  dplyr::group_by(scenario, year, gridcode, lon, lat, zone) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(values_from = value, names_from = 'zone') %>% 
  dplyr::mutate(crop_land = dplyr::if_else(crop > 0.1, 1, 0)) %>% 
  dplyr::select(lat, lon, crop_land)


# ------------------------------------------------------------------------------
# Calculate Exploitable Area and Volume
# ------------------------------------------------------------------------------
expan_zone <- land_use_df %>% 
  dplyr::left_join(wdpa_df, by = c('lat', 'lon')) %>% 
  dplyr::left_join(pop_density_df, by = c('lat', 'lon')) %>% 
  dplyr::left_join(glwd_df, by = c('lat', 'lon')) %>% 
  dplyr::mutate(protected_area = dplyr::if_else(is.na(protected_area), 0, protected_area),
                pop_density = dplyr::if_else(is.na(pop_density), 0, pop_density),
                water = dplyr::if_else(is.na(water), 0, water),
                expan_zone = dplyr::if_else((crop_land == 1 |
                                               protected_area == 1 |
                                               (pop_density == 0 |
                                                  pop_density > 1244) |
                                               water == 0), 0, 1))

expan_zone_df <- expan_zone %>% 
  dplyr::select(lat, lon, expan_zone) %>% 
  dplyr::rename(x = lon,
                y = lat,
                z = expan_zone) %>% 
  dplyr::select(x, y, z)
  
# Create raster from gridded exclusion zone data
expan_zone_raster <- rasterFromXYZ(expan_zone_df, res = 0.5,
                                   crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
plot(expan_zone_raster, col=c('lightgrey', 'black'))

# calculate the exploitable area in km2
expan_area <- as.data.frame(raster::area(expan_zone_raster), xy = TRUE, na.rm = TRUE) %>% 
  dplyr::rename(expan_area = layer,
                lon = x,
                lat = y)


# Exploitable areas by grid and by basin
# expan_category = 0: no expansion
# expan_category = 1: no expansion but has existing reservoirs in the grid
# expan_category = 2: expansion with existing reservoirs in the grid
# expan_category = 3: new expansion without existing reservoirs in the grid
expan_grid <- expan_zone %>% 
  dplyr::select(lat, lon, expan_zone) %>% 
  # grid cell area
  dplyr::left_join(expan_area, by = c('lat', 'lon')) %>% 
  dplyr::left_join(grid_basin_mapping, by = c('lat', 'lon')) %>% 
  dplyr::left_join(dam_grid_info %>% 
                     dplyr::select(lat, lon, AREA_KM2, CAP_KM3, LAKE_AREA_KM2, LAKE_VOL_KM3), 
                   by = c('lat', 'lon')) %>% 
  tidyr::replace_na(list(AREA_KM2 = 0,
                         CAP_KM3 = 0,
                         LAKE_AREA_KM2 = 0,
                         LAKE_VOL_KM3 = 0)) %>% 
  dplyr::mutate(expan_category =
                  dplyr::case_when(
                    expan_zone == 0 & CAP_KM3 == 0 ~ 0,
                    expan_zone == 0 & CAP_KM3 > 0 ~ 1,
                    expan_zone == 1 & CAP_KM3 >= LAKE_VOL_KM3 & LAKE_VOL_KM3 > 0 ~ 1,
                    expan_zone == 1 & CAP_KM3 > 0 & LAKE_VOL_KM3 > CAP_KM3  ~ 2,
                    expan_zone == 1 & CAP_KM3 == 0 ~ 3)) %>% 
  dplyr::mutate(expan_cap_km3 = 
                  dplyr::case_when(expan_category == 0 ~ 0,
                                   expan_category == 1 ~ CAP_KM3,
                                   expan_category == 2 ~ LAKE_VOL_KM3,
                                   expan_category == 3 ~ LAKE_VOL_KM3),
                expan_area_km2 = 
                  dplyr::case_when(expan_category == 0 ~ 0,
                                   expan_category == 1 ~ AREA_KM2,
                                   expan_category == 2 ~ LAKE_AREA_KM2,
                                   expan_category == 3 ~ LAKE_AREA_KM2)
) %>% 
  dplyr::select(lat, lon, basin_id, basin_name, expan_zone, expan_category,
                expan_cap_km3, expan_area_km2)

output_expan_grid <- expan_grid %>% 
  dplyr::select(lat, lon, basin_id, basin_name, expan_cap_km3, expan_area_km2)

# aggregate to basin
output_expan_basin <- output_expan_grid %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::summarise(expan_cap_km3 = sum(expan_cap_km3),
                   expan_area_km2 = sum(expan_area_km2)) %>% 
  dplyr::ungroup()


# ------------------------------------------------------------------------------
# Derive Nonlinear Area-Volume relationship V = cA^b (Liu et al., 2018)
# ------------------------------------------------------------------------------
# get reservoir capacity and reservoir surface area
grand_xanthos <- grand_georef %>% 
  dplyr::select(basin_id, basin_name, AREA_KM2, CAP_KM3)

res_basin <- grand_xanthos %>%
  dplyr::filter(!is.na(CAP_KM3) & AREA_KM2 > 0) %>% 
  dplyr::rename(AREA = AREA_KM2,
                CAP = CAP_KM3)

output_param <- data.frame(basin_id = numeric(),
                           basin_name = character(),
                           b = numeric(),
                           c = numeric())

# global reservoir param from GranD doc v1.3, converted to km2 (Area) and km3 (Capacity)
granD_c <- 0.030684
granD_b <- 0.9578
                            
for(id in 1:235){
  basin_df <- res_basin %>% 
    dplyr::filter(basin_id == id)
  basin_name <- basin_rmap$basin_name[basin_rmap$basin_id == id]
  
  if(nrow(basin_df) <= 5){
    # Use global reservoir param from GranD doc v1.3, converted to km2 (Area) and km3 (Capacity)
    c <- granD_c
    b <- granD_b
  } else {
    nls_model <- minpack.lm::nlsLM(formula = CAP ~ c * AREA ^ b,
                                   data = basin_df,
                                   start = list(c = 0.03, b = 1))
    
    summary(nls_model)
    c <- coef(nls_model)[[1]]
    b <- coef(nls_model)[[2]]

    # calculate linear regression by taking log of x and y.  y = m + n*x
    x <- log10(basin_df$AREA)
    y <- log10(basin_df$CAP)
    basin_log <- data.frame(x = x, y = y)
    lm_model <- lm(formula = y ~ x, data = basin_log)
    m <- format(unname(coef(lm_model)[1]), digits = 3)
    n <- format(unname(coef(lm_model)[2]), digits = 3)
    r2 <- format(summary(lm_model)$r.squared, digits = 3)
    eq_label_log <- substitute(italic(y) == m + n %.% italic(x)*','~~italic(r)^2~'='~r2,
                               list(m = m, n = n, r2 = r2))
    
    if(r2 < 0.2){
      c <- granD_c
      b <- granD_b
    }
    
  }
  # output parameter b and c: V = c*A^b. V is in km3, A is in km2
  param_temp <- data.frame(basin_id = id,
                           basin_name = basin_name,
                           b = b,
                           c = c)
  output_param <- output_param %>% 
    dplyr::bind_rows(param_temp)

  eq_label <- substitute(italic(y) == c%.%italic(x)^b, list(c = c, b = b))

  
  basin_log <- data.frame()
  
  
}


# ------------------------------------------------------------------------------
# Basin Historical Capacity and Future Mean Potential Capacity Expansion
# ------------------------------------------------------------------------------
# This is to estimate the unit size of reservoir expansion for each basin
# Calculate world mean storage capacity
expan_category_grid <- expan_grid %>% 
  dplyr::select(lat, lon, expan_category)
expan_dam <- dam_info %>% 
  dplyr::left_join(expan_category_grid,
                   by = c('lat', 'lon')) %>% 
  dplyr::left_join(dam_purpose_grid %>% 
                     dplyr::select(xanthos_id, purpose_main),
                   by = c('xanthos_id')) %>% 
  dplyr::mutate(expan_cap_km3 = 
                  dplyr::case_when(expan_category == 0 ~ 0,
                                   expan_category == 1 ~ CAP_KM3,
                                   expan_category == 2 ~ LAKE_VOL_KM3,
                                   expan_category == 3 ~ LAKE_VOL_KM3))

capacity_median_expan <- median(expan_dam$LAKE_VOL_KM3[expan_dam$GRAND_ID == 0 & expan_dam$expan_zone == 1])
capacity_mean_hist <- mean(expan_dam$CAP_KM3[expan_dam$GRAND_ID > 0])

# Get basin average storage capacity
# Use dam data point to estimate the mean expansion size of the reservoir
# Historical + future potential mean expansion capacity based on lake volumes in expandable grids
capacity_future_all <- expan_dam %>%
  dplyr::filter(expan_category > 0) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::summarise(mean_cap_km3 = mean(expan_cap_km3)) %>% 
  dplyr::ungroup()

capacity_future_no_constraint <- expan_dam %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::summarise(qt_cap_all_km3 = quantile(LAKE_VOL_KM3, c(0.9)),
                   mean_cap_all_km3 = mean(LAKE_VOL_KM3),
                   select_cap_km3 = min(qt_cap_all_km3, mean_cap_all_km3)) %>% 
  dplyr::ungroup()

capacity_future_all <- capacity_future_no_constraint %>% 
  dplyr::left_join(capacity_future_all, by = c('basin_id', 'basin_name')) %>% 
  dplyr::mutate(mean_cap_km3 = dplyr::if_else(is.na(mean_cap_km3), 
                                              select_cap_km3, 
                                              mean_cap_km3)) %>% 
  dplyr::select(basin_id, basin_name, mean_cap_km3)

# historical total dam capacity and surface area within each basin
capacity_hist_all <- expan_dam %>% 
  dplyr::filter(GRAND_ID > 0) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::summarise(total_cap_km3 = sum(CAP_KM3)) %>% 
  dplyr::ungroup()
capacity_hist_nonhydro <- expan_dam %>% 
  dplyr::filter(GRAND_ID > 0,
                !purpose_main %in% c('hydropower')) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::summarise(nonhydro_cap_km3 = sum(CAP_KM3),
                   nonhydro_area_km2 = sum(AREA_KM2),
                   count = mean(count)) %>% 
  dplyr::ungroup()
capacity_hist_all <- capacity_hist_all %>% 
  dplyr::full_join(capacity_hist_nonhydro, by = c('basin_id', 'basin_name'))

# manually reduce mean_cap_km3 for the following basins because 
# their supply curves' max price tend to be lower than GCAM solved price
# based on output from diagnostic_check_price.R
problem_basins_tier1 <- c(77, 87, 97, 107, 221)
problem_basins_tier2 <- c(89)
problem_basins_tier3 <- c(104)
output_capacity <- basin_rmap %>%
  dplyr::left_join(capacity_future_all, by = c('basin_id', 'basin_name')) %>%
  dplyr::left_join(capacity_hist_all, by = c('basin_id', 'basin_name')) %>%
  tidyr::replace_na(list(nonhydro_cap_km3 = 0,
                         nonhydro_area_km2 = 0,
                         total_cap_km3 = 0)) %>% 
  dplyr::group_by(basin_id, basin_name) %>% 
  dplyr::mutate(mean_cap_km3 = dplyr::case_when(
    basin_id %in% problem_basins_tier1 ~ mean_cap_km3 / 5,
    basin_id %in% problem_basins_tier2 ~ mean_cap_km3 / 20,
    basin_id %in% problem_basins_tier3 ~ mean_cap_km3 / 50,
    TRUE ~ mean_cap_km3)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-count)


# ------------------------------------------------------------------------------
# Aggregate Outputs for LP and Save
# ------------------------------------------------------------------------------
output <- output_expan_basin %>% 
  dplyr::left_join(output_param, by = c('basin_id', 'basin_name')) %>% 
  dplyr::left_join(output_capacity, by = c('basin_id', 'basin_name')) %>% 
  dplyr::left_join(gcam_basin_mapping, by = c('basin_name')) %>% 
  dplyr::rename(subsector = gcam_basin_name) %>% 
  dplyr::select(basin_id, basin_name, subsector, 
                # historical data
                mean_cap_km3, nonhydro_cap_km3, nonhydro_area_km2, total_cap_km3, 
                # future potential
                expan_cap_km3, expan_area_km2,
                # reservoir V = c*(A^b) from historical data
                b, c)

write.csv(x = output, 
          file = file.path(work.dir, 'outputs_LP', 'LP_reservoir.csv'),
          row.names = FALSE)


# Save certain params that will be used in other scripts
saveRDS(list(dam_info = dam_info, 
             dam_grid_info = dam_grid_info, 
             dam_purpose_grid = dam_purpose_grid,
             expan_dam = expan_dam,
             expan_grid = expan_grid, 
             expan_category_grid = expan_category_grid, 
             output_expan_basin = output_expan_basin, 
             output = output),
        file = file.path(output.dir, 'LP_inputs_reservoir_20221222.rds'))
