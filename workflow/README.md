## Inputs and Outputs for Each Script

Each script processes input data and produces required files for running the GLORY model. 

1. Download the input data from [input data (DOI needed)](<Zenodo DOI>).
2. Download the publicly available data listed in [Table 1](../README.md#table1) into the `data` folder from Step 1. Please unzip the downloaded data if zipped.[Table 1](#table1) shows the data required for each script.


<a name="table1"></a>
**Table 1:** Input data.

| Script Name | Inputs | Outputs |
| --- | --- | --- |
| `GLORY_inputs_climate_future.R` | **Hydrologic Data** <br> * [Xanthos] monthly basin runoff `Basin_runoff_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` <br> * [Xanthos] smoothed annual basin runoff `runoff_impacts_MIROC-ESM-CHEM_rcp6p0.csv` <br> * [Xanthos] monthly gridded PET ` pet_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` <br><br> **Dam and Reservoir Data** <br> * groreferenced GranD `GranD_v1.3_remap_to_xanthos.csv` <br><br> **Reference Data** <br> * grid area `Grid_Areas_ID.csv` <br> * grid basin mapping `basin.csv` |  * GLORY climate inputs for annual runoff and ET depth `LP_climate_GCM_rcp.csv` |
| `GLORY_inputs_monthly_profiles_future.R` |  **Hydrologic Data** <br> * [Xanthos] monthly gridded streamflow `avgchflow_m3persec_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` <br> * [Xanthos] monthly gridded PET `pet_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` <br> * VIC historical monthly runoff `vic_watch_basin_km3_1971_2001_monthly.csv` <br><br> **Water Demand Data** <br> * tethys water demand `wd_tethys_2005_2010_gcam5p3-stash_GFDL-ESM2M_rcp2p6_235basins.RDS` <br><br>  **Dam and Reservoir Data** <br> * dam purpose `reservoir_purpose.csv` <br><br> **Reference Data** <br> * grid coordinates `coordinates.csv` <br> * grid basin mapping `basin.csv` <br> * basin id name mapping `basin_ID.csv` <br> * contributing grids to basin outlet `outlet.csv` | * monthly annual and sectoral demand profile `LP_fraction_profile_GCM_rcp.csv` <br> * historical mean annual sectoral demand `LP_demand_sector_mean_annual_hist.csv` |
| `GLORY_inputs_reservoir.R` | **Dam and Reservoir Data** <br> * GranD v1.3 <br> * HydroLAKES georeferenced to xanthos grid `HydroLAKES_to_xanthos.csv` <br><br> **Exclusion Area Data** <br> * GLWD level3 data ([Table 1](../README.md#table1)) <br> * population ([Table 1](../README.md#table1)) <br> * WDPA ([Table 1](../README.md#table1)) <br> * landuse from Demeter ([Table 1](../README.md#table1)) | * reservoir capacity and exploitable potential `LP_reservoir.csv` <br> * RData from the script `LP_inputs_reservoir_20221222.RData` |
| `GLORY_inputs_mean_slope_basin.R` | **Slope Data** <br> * 	EarthEnv remapped slope data `slope_remapped.csv` | * basin mean slope `LP_mean_slope_basin.csv` |
| `GLORY_gcamwrapper.py` | All the outputs from above | * GCAM output database <br> * selected intermediate outputs |