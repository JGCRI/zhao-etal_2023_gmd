## Inputs and Outputs for Each Script

Each script processes input data and produces required files for running the GLORY model. 

1. Clone the GitHub repository to your local device.
```
git clone https://github.com/JGCRI/zhao-etal_2023_gmd.git
```
2. Download the input data from [input data](https://doi.org/10.5281/zenodo.8436685), and put those data under cloned GitHub repository `zhao-etal_2023_gmd/workflow/data`. These input data does not include public available data from other studies (see step 3).
3. Download the publicly available data listed in [Table 1](../README.md#table1) into the `zhao-etal_2023_gmd/workflow/data` folder from Step 2. Please unzip the downloaded data if zipped. [Table 1](#table1) shows the data required for each script.


<a name="table1"></a>
**Table 1:** Input data.

| Script Name | Inputs | Outputs |
| --- | --- | --- |
| `GLORY_inputs_climate_future.R` | **Hydrologic Data** <br> <ul><li>[Xanthos] monthly basin runoff `Basin_runoff_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` </li><li>[Xanthos] smoothed annual basin runoff `runoff_impacts_MIROC-ESM-CHEM_rcp6p0.csv` </li><li>[Xanthos] monthly gridded PET ` pet_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` </li></ul><br> **Dam and Reservoir Data** <ul><li>groreferenced GranD `GranD_v1.3_remap_to_xanthos.csv` </li></ul><br> **Reference Data** <br><ul><li>grid area `Grid_Areas_ID.csv` </li><li>grid basin mapping `basin.csv`</li></ul> |  <ul><li>GLORY climate inputs for annual runoff and ET depth `LP_climate_GCM_rcp.csv`/li></ul> |
| `GLORY_inputs_monthly_profiles_future.R` |  **Hydrologic Data** <br> <ul><li>[Xanthos] monthly gridded streamflow `avgchflow_m3persec_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` </li><li>[Xanthos] monthly gridded PET `pet_km3permonth_MIROC-ESM-CHEM_rcp6p0_1950_2099.csv` </li><li>VIC historical monthly runoff `vic_watch_basin_km3_1971_2001_monthly.csv` </li></ul><br>**Water Demand Data** <br> <ul><li>tethys water demand `wd_tethys_2005_2010_gcam5p3-stash_GFDL-ESM2M_rcp2p6_235basins.RDS` </li></ul><br>**Dam and Reservoir Data** <br> <ul><li>dam purpose `reservoir_purpose.csv` </li></ul><br> **Reference Data** <br> <ul><li>grid coordinates `coordinates.csv` </li><li>grid basin mapping `basin.csv` </li><li>basin id name mapping `basin_ID.csv` </li><li>contributing grids to basin outlet `outlet.csv` </li></ul>| <ul><li>monthly annual and sectoral demand profile `LP_fraction_profile_GCM_rcp.csv`</li><li> historical mean annual sectoral demand `LP_demand_sector_mean_annual_hist.csv`</li></ul> |
| `GLORY_inputs_reservoir.R` | **Dam and Reservoir Data** <br> <ul><li>GranD v1.3 </li>HydroLAKES georeferenced to xanthos grid `HydroLAKES_to_xanthos.csv` </li></ul><br> **Exclusion Area Data** <ul><li>GLWD level3 data ([Table 1](../README.md#table1)) </li><li>population ([Table 1](../README.md#table1)) </li><li>WDPA ([Table 1](../README.md#table1))</li><li>landuse from Demeter ([Table 1](../README.md#table1)) </li></ul>| <ul><li>reservoir capacity and exploitable potential `LP_reservoir.csv` </li><li>RData from the script `LP_inputs_reservoir_20221222.RData`</li></ul> |
| `GLORY_inputs_mean_slope_basin.R` | **Slope Data** <br><ul><li>EarthEnv remapped slope data `slope_remapped.csv`</li></ul> | <ul><li>basin mean slope `LP_mean_slope_basin.csv`</li></ul> |
| `GLORY_gcamwrapper.py` | All the outputs from above | <ul><li>GCAM output database</li><li>selected intermediate outputs</li></ul> |