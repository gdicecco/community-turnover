### Subset BBS data for community trajectory analysis

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)

### Read in data #####

## NA map

na <- world %>%
  filter(continent == "North America")

## BBS 2017 Version
routes <- read.csv("\\\\BioArk\\HurlbertLab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# species list
species_list <- species %>%  
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010) %>%
  filter(sporder != "Accipitriformes", 
         sporder != "Falconiformes", 
         sporder != "Anseriformes",
         sporder != "Cathartiformes")

counts.subs <- counts %>%
  filter(rpid == 101) %>%
  filter(aou %in% species_list$aou) %>%
  right_join(RT1.routes, by = c("countrynum", "statenum", "stateroute", "year"))

### sample sizes for route time series

rtes_per_year <- counts.subs %>%
  group_by(year) %>%
  summarize(n_rtes = n_distinct(stateroute))

### routes sampled continuously 1970-present

cont_routes <- counts.subs %>%
  filter(year >= 1970) %>%
  mutate(year_bin = 5*floor(year/5)) %>%
  group_by(stateroute) %>%
  mutate(n_bins = n_distinct(year_bin)) %>%
  filter(n_bins == 10)

route_sf <- cont_routes %>%
  ungroup() %>%
  distinct(stateroute, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))

bbs_map <- tm_shape(na) + tm_polygons() + tm_shape(route_sf) + tm_dots()
# tmap_save(bbs_map, "figures/bbs_route_map_1970-2016.pdf")

## Annual rank abundance distributions
# min_rank gives ties the same value - e.g. 1,2,2,4

rank_abund <- cont_routes %>%
  group_by(stateroute, year) %>%
  mutate(rank = min_rank(desc(speciestotal))) %>%
  dplyr::select(stateroute, year, aou, rank)

# write.csv(rank_abund, "data/bbs_subset_1970-2016_ranks.csv", row.names = F)
