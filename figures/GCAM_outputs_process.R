################################################################################
# Post-process GCAM outputs and produce easy to access rds and dat files
# Author: Mengqi Zhao
# Email: mengqi.zhao@pnnl.gov
# Last Update: 2022-10
################################################################################

# Load libraries
library(data.table)
library(dplyr)
library(gcamextractor)
library(rgcam)

# ------------------------------------------------------------------------------
# Set working directory and data path
# ------------------------------------------------------------------------------
work.dir <- 'zhao-etal_2023_gmd/workflow' # update correspondingly
output.dir <- file.path(work.dir, 'outputs_glory')
db.dir <- file.path(output.dir, 'gcam_output_database')


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
# Process GCAM output and Extract selected Parameters
# ------------------------------------------------------------------------------
dataGCAM <- tibble::tibble()
if(file.exists(file.path(output_path, exp_name, rdata_name))){
  
  dataGCAM <- readRDS(file.path(output_path, exp_name, rdata_name))
  
} else {
  
  for(db in 1:length(exp$database)){
    
    dataGCAM_temp <- gcamextractor::readgcam(
      gcamdatabase = file.path(db.dir, exp$database[db]),
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


# ------------------------------------------------------------------------------
# Read GCAM output with rgcam and save to dat
# ------------------------------------------------------------------------------
f_query <- file.path(db.dir, 'queries_resource_supply_curves.xml')

# set experiment name
model <- 'gcam'
region <- '32regions'
sector <- 'market'
scenario <- 'allscenarios'
exp_name <- paste(model, region, sector, scenario, sep = '_')
rdata_name <- paste0(exp_name, '.rds')
prj_file <- file.path(db.dir, paste0(exp_name, '.dat'))

if(file.exists(prj_file)){
  prj <- rgcam::loadProject(prj_file)
} else {
  prj <- list()
  prj_ref <- list()
  for (i in 1:length(exp$database)){
    conn <- rgcam::localDBConn(dbPath = db.dir,
                               dbFile = exp$database[i])
    
    if(i == 1){
      prj_ref <- rgcam::addScenario(
        conn,
        file.path(db.dir,
                  paste0(paste(model, region, sector, 'ref', sep = '_'), '.dat')),
        scenario = exp$scenario_orig[i],
        queryFile = f_query)
    } else {
      prj <- rgcam::addScenario(
        conn,
        file.path(db.dir,
                  paste0(exp_name, '.dat')),
        scenario = exp$scenario_orig[i],
        queryFile = f_query)
    }
    
  }
}

# remove duplicated scenario and save project
prj <- prj[-2]
rgcam::saveProject(prj, file = prj_file)