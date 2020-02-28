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
# write.csv(species_list, "data/species_list.csv", row.names = F)

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
  group_by(countrynum, stateroute) %>%
  mutate(n_bins = n_distinct(year_bin)) %>%
  filter(n_bins == 10)

route_sf <- cont_routes %>%
  ungroup() %>%
  distinct(stateroute, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))

bbs_map <- tm_shape(na) + tm_polygons() + tm_shape(route_sf) + tm_dots()
# tmap_save(bbs_map, "figures/bbs_route_map_1970-2016.pdf")

### Route density 1990-2016 - 1-5 years in every five year time window

counts_landcover_years <- counts.subs %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016),
         max_bins = case_when(countrynum == 124 ~ 5,
                            countrynum == 840 ~ 6))

counts_per_window <- function(surveys_per_window) {
  cont_routes <- counts_landcover_years %>%
    filter(year >= y1, year <= y2) %>%
    mutate(year_bin = 5*floor(year/5)) %>%
    group_by(max_bins, stateroute, year_bin) %>%
    summarize(n_years = n_distinct(year)) %>%
    filter(n_years >= surveys_per_window) %>%
    group_by(stateroute) %>%
    mutate(n_bins = n_distinct(year_bin)) %>%
    filter(n_bins == max_bins)
    
  return(list(df = cont_routes, n_routes = length(unique(cont_routes$stateroute))))
}

bbs_density <- data.frame(years_per_window = 1:5) %>%
  mutate(routes = map_dbl(years_per_window, ~{
      output <- counts_per_window(.)
      output$n_routes
    }))

ggplot(bbs_density, aes(x = years_per_window, y = routes)) + geom_point() + geom_line(cex = 1) + theme_classic(base_size = 15)

pdf("figures/bbs_route_map_1990-2016_surveys_per_window.pdf")
for(i in c(1:2)) {
  output <- counts_per_window(i)
  cont_routes <- output$df
  routes_short <- route_sf %>%
    filter(stateroute %in% cont_routes$stateroute)
  map <- tm_shape(na) + tm_polygons() + tm_shape(routes_short) + tm_dots() + 
    tm_layout(title = paste0("Surveys per window = ", i))
  print(map)
}
dev.off()

counts_subs <- counts_per_window(1)
routes_subs <- counts_subs$df

# write.csv(routes_subs, "data/bbs_route_subset_1990-2016.csv", row.names = F)

## Log abundance

log_abund <- cont_routes %>%
  group_by(countrynum, stateroute, year) %>%
  mutate(log_abund = log10(speciestotal)) %>%
  dplyr::select(stateroute, year, aou, log_abund)

# write.csv(log_abund, "data/bbs_subset_1970-2016_logabund.csv", row.names = F)
