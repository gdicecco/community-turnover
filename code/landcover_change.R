### Land cover change metrics

library(tidyverse)
library(purrr)
library(tmap)

## read in data

bbs_halfroutes <- read.csv("/Volumes/hurlbertlab/DiCecco/data/bbs_half_routes_fragstats.csv", stringsAsFactors = F)

## Measure land cover change as max delta at each route

# Function to calculate the class with the maximum change in proportion landscape from earliest to late time period
max_delta <- function(df) {
  delta <- df %>%
    filter(year == min_year | year == max_year) %>%
    mutate(year_name = case_when(year == min_year ~ "year1",
                                 year == max_year ~ "year2")) %>%
    select(year_name, class, prop.landscape) %>%
    distinct() %>%
    pivot_wider(names_from = year_name, values_from = prop.landscape) %>%
    replace_na(list(year1 = 0, year2 = 0)) %>%
    mutate(deltaCover = year2 - year1) %>%
    mutate(absCover = abs(deltaCover)) %>%
    filter(absCover == max(absCover)) %>%
    select(class, deltaCover)
  
  return(delta)
}

possibly_max_delta <- possibly(max_delta, data.frame(class = NA, deltaCover = NA))

bbs_maxdeltas <- bbs_halfroutes %>%
  mutate(min_year = case_when(country == "US" ~ 1992,
                              TRUE ~ 1990),
         max_year = case_when(country == "US" ~ 2016,
                              TRUE ~ 2010)) %>%
  group_by(country, stateroute, stops) %>%
  nest() %>%
  mutate(maxDelta = map(data, ~possibly_max_delta(.))) %>%
  select(-data) %>%
  unnest()
# write.csv(bbs_maxdeltas, "data/bbs_half_route_max_land_change.csv", row.names = F)

# Correlation between max land cover values for half routes

bbs_maxdeltas_complete <- bbs_maxdeltas %>%
  group_by(stateroute) %>%
  filter(n() == 2)

cor(filter(bbs_maxdeltas_complete, stops == "1-25")$class, filter(bbs_maxdeltas_complete, stops == "26-50")$class, use = "pairwise.complete.obs")
# 0.9

cor(filter(bbs_maxdeltas_complete, stops == "1-25")$deltaCover, filter(bbs_maxdeltas_complete, stops == "26-50")$deltaCover, use = "pairwise.complete.obs")
# 0.5

# For single routes - max land cover change

landcover_us <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F) %>%
  mutate(country = "US")
landcover_ca <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
  mutate(country = "Canada")

landcover_na <- bind_rows(landcover_us, landcover_ca)

oneroute_maxdeltas <- landcover_na %>%
  mutate(min_year = case_when(country == "US" ~ 1992,
                              TRUE ~ 1990),
         max_year = case_when(country == "US" ~ 2016,
                              TRUE ~ 2010)) %>%
  group_by(country, stateroute) %>%
  nest() %>%
  mutate(maxDelta = map(data, ~possibly_max_delta(.))) %>%
  select(-data) %>%
  unnest()
# write.csv(oneroute_maxdeltas, "data/bbs_route_max_landcover_change.csv", row.names = F)

# Category of maximum land cover change

# Landcover legend US
newcode <- data.frame(Code = seq(1,9), 
                      Label = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Landcover legend Canada
ca_codes <- read.csv("data/canada_landcover_classification.csv", stringsAsFactors = F) %>%
  select(-Definition)

landcover_legend <- bind_rows(newcode, ca_codes)

oneroute_plot <- oneroute_maxdeltas %>%
  na.omit() %>%
  left_join(landcover_legend, by = c("class" = "Code")) %>%
  left_join(routes) %>%
  st_as_sf(coords = c("longitude", "latitude"))

na_map <- read_sf("data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(iso_a2 == "CA" | iso_a2 == "US") %>%
  filter(name != "Hawaii")

max_delta_map <- tm_shape(na_map) + tm_borders() +
  tm_shape(oneroute_plot) + tm_dots(col = "Label", size = 0.1, palette = "Set2")
tmap_save(max_delta_map, "figures/max_landcover_change_bbs_routes.pdf")

## Measure land cover change as PCA values of land cover class deltas

# This doesn't work - need to remove routes with duplicates
bbs_landcover_wide <- bbs_halfroutes %>%
  filter(country == "US") %>%
  mutate(min_year = 1992,
         max_year = 2016) %>%
  filter(year == min_year | year == max_year) %>%
  mutate(year_name = case_when(year == min_year ~ "year1",
                               year == max_year ~ "year2")) %>%
  select(year_name, class, prop.landscape) %>%
  distinct() %>%
  pivot_wider(names_from = year_name, values_from = prop.landscape) %>%
  replace_na(list(year1 = 0, year2 = 0)) %>%
  mutate(deltaCover = year2 - year1)

bbs_landcover_wide_ca <- bbs_maxdeltas %>%
  filter(country == "Canada") %>%
  pivot_wider(names_from = class, values_from = deltaCover)
