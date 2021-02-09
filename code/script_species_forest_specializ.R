### Calculate forest specialization of species using mean forest cover of BBS routes
### Create traits table used for trait models in model/script_main_analyses.R

library(tidyverse)

## Make master data table

abund_trend <- read.csv("data/derived_data//BBS_abundance_trends.csv", stringsAsFactors = F)

env_change <- read.csv("data/derived_data/bbs_route_env_change.csv", stringsAsFactors = F)

clim_hab_poptrend <- abund_trend %>%
  left_join(env_change, by = "stateroute") %>%
  filter(!is.na(propForest))

## Read in forest cover data for BBS routes

frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)
frags_ca <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
  group_by(year, stateroute) %>% 
  mutate(n_zones = as.factor(n_distinct(file))) %>%
  filter(n_zones == 1)

# Landcover legend US
newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Landcover legend Canada
canada_code <- read.csv("data/canada_landcover_classification.csv", stringsAsFactors = F)
colnames(canada_code)[1:2] <- c("legend", "code")

#### Species traits ####

# Use all routes for cover_breadth

## Breadth of forest cover for all routes
forest_ed_allroutes <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2016)

forest_ed_canada <- frags_ca %>%
  left_join(canada_code, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2010)

forest_allroutes <- forest_ed_allroutes %>%
  bind_rows(forest_ed_canada)

clim_hab_pop_allroutes <- abund_trend %>%
  left_join(forest_allroutes, by = "stateroute") %>%
  left_join(dplyr::select(env_change, stateroute, year, tmax, tmin), by = "stateroute") %>%
  filter(!is.na(propForest))

env_breadth_allroutes <- clim_hab_pop_allroutes %>%
  dplyr::group_by(aou) %>%
  summarize(propFor = mean(propForest),
            min_for = quantile(propForest, c(0.05))[[1]],
            max_for = quantile(propForest, c(0.95))[[1]],
            patchArea = mean(meanPatchArea),
            std_patch = sd(meanPatchArea),
            for_range = max_for - min_for)

write.csv(env_breadth_allroutes, "data/derived_data/spp_forest_traits.csv", row.names = F)

