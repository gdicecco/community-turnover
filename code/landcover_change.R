### Land cover change metrics

library(tidyverse)
library(purrr)

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

## Measure land cover change as PCA values (also use directionality?)