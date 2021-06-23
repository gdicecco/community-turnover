# community-turnover

## Anthropogenic drivers of avian community turnover across scales

This study examines how avian community turnover is predicted by climate and land cover change across North America using data from the Breeding Bird Survey. We examine the relative importance of anthropogenic change drivers in explaining turnover from local to regional scales. We also explore whether turnover can be attritbuted to species groups (trophic group, foraging guild, migratory status, and breeding biome) or is driven by individual, influential species.

### `/code` contains scripts to run analyses and produce figures

#### Main analysis scripts:
- `data_cleaning.R`: cleans Breeding Bird Survey data for target survey locations
- `landcover_change.R`: calculates land cover change metrics at Breeding Bird Survey routes
- `trajectory_analysis.R`: contains code for main analyses and manuscript figures

#### Other scripts:
- `scale_model_input_vars_spatialautocorr.R`: calculates spatial autocorrelation of input anthropogenic change variables
- `trajectory_sensitivity_analyses.R`: examines how sensitive turnover metric is to species pool, time windows, and other analysis choices
- `directionality_sim_datasets.R`: explores how turnover metric behaves in different cases with simulated data

### `/figures` contains manuscript figures

### `/data` contains raw data used in analyses

#### `/derived_data` contains data produced by analysis scripts
