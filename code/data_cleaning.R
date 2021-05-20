### Subset BBS data for community trajectory analysis

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)
library(ecospat)

### Read in data #####

## NA map

na <- world %>%
  filter(continent == "North America")

## BBS 2017 Version

# Append correct BioArk path
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\ad.unc.edu\\bio")

routes <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_routes_20170712.csv"))
counts <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_counts_20170712.csv"))
species <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_species_20170712.csv"))
weather <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_weather_20170712.csv"))

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

# four letter species codes 

fourletter_codes <- read.csv("data/four_letter_codes_birdspp.csv", stringsAsFactors = F) %>%
  left_join(species_list, by = c("COMMONNAME" = "english_common_name")) %>%
  filter(!is.na(aou))
# write.csv(fourletter_codes, "data/four_letter_codes_aous.csv", row.names = F)

# Filter BBS to rpid = 101, runtype = 1, land birds != birds of prey
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

### Route density 1990-2016 - 1-4 years in every four year time window
# countrynum 124 = Canada, 840 = US
counts_landcover_years <- counts.subs %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016),
         max_bins = case_when(countrynum == 124 ~ 5,
                            countrynum == 840 ~ 6))

## Count number of routes that meet a threshhold (surveys_per_window) of surveys per five year time window across study period
counts_per_window <- function(surveys_per_window) {
  cont_routes <- counts_landcover_years %>%
    filter(year >= y1, year <= y2) %>%
    group_by(countrynum) %>%
    nest() %>%
    mutate(year_bins = map2(countrynum, data, ~{
      country <- .x
      df <- .y
      
      if(country == 124) {
        df %>%
          mutate(year_bin = case_when(year >= 1990 & year <= 1993 ~ 1990,
                                      year >= 1994 & year <= 1997 ~ 1994,
                                      year >= 1998 & year <= 2001 ~ 1998,
                                      year >= 2002 & year <= 2005 ~ 2002,
                                      TRUE ~ 2006))
      } else {
        df %>%
          mutate(year_bin = case_when(year >= 1992 & year <= 1995 ~ 1992,
                                      year >= 1996 & year <= 1999 ~ 1996,
                                      year >= 2000 & year <= 2003 ~ 2000,
                                      year >= 2004 & year <= 2007 ~ 2004,
                                      year >= 2008 & year <= 2011 ~ 2008,
                                      year >= 2012 & year <= 2016 ~ 2012))
      }
    })) %>%
    select(-data) %>%
    unnest(cols = c(year_bins)) %>%
    group_by(max_bins, stateroute, year_bin) %>%
    summarize(n_years = n_distinct(year)) %>%
    filter(n_years >= surveys_per_window) %>%
    group_by(stateroute) %>%
    mutate(n_bins = n_distinct(year_bin)) %>%
    filter(n_bins == max_bins)
    
  return(list(df = cont_routes, n_routes = length(unique(cont_routes$stateroute))))
}

bbs_density <- data.frame(years_per_window = 1:4) %>%
  mutate(routes = map_dbl(years_per_window, ~{
      output <- counts_per_window(.)
      output$n_routes
    }))

ggplot(bbs_density, aes(x = years_per_window, y = routes)) + geom_point() + geom_line(cex = 1) + theme_classic(base_size = 15)
ggsave("figures/bbs_route_density.pdf")

pdf("figures/bbs_route_map_1990-2016_surveys_per_window.pdf")
for(i in c(1:4)) {
  output <- counts_per_window(i)
  cont_routes <- output$df
  routes_short <- route_sf %>%
    filter(stateroute %in% cont_routes$stateroute)
  map <- tm_shape(na) + tm_polygons() + tm_shape(routes_short) + tm_dots() + 
    tm_layout(title = paste0("Surveys per window = ", i))
  print(map)
}
dev.off()

counts_subs <- counts_per_window(3)
routes_subs <- counts_subs$df

# write.csv(routes_subs, "data/bbs_route_subset_1990-2016.csv", row.names = F)

bbs_subset <- counts.subs %>%
  filter(stateroute %in% routes_subs$stateroute) %>%
  filter(year >= 1990)

# write.csv(bbs_subset, "data/bbs_counts_subset_1990-2016.csv", row.names = F)

## Sample size per BCR

bcr_n <- bbs_subset %>%
  group_by(bcr) %>%
  summarize(n_routes = n_distinct(stateroute))

## Figure: map of BCRs with BBS routes, number of routes per BCR that meet sampling density thresholds

bcrs <- read_sf(paste0(bioark, "/HurlbertLab/DiCecco/bcr_terrestrial_shape/BCR_terrestrial_master.shp")) %>%
  filter(COUNTRY == "USA" | COUNTRY == "CANADA") %>%
  st_crop(c(xmin = -178, ymin = 18.9, xmax = -53, ymax = 60)) %>%
  filter(PROVINCE_S != "ALASKA" & PROVINCE_S != "HAWAIIAN ISLANDS" & PROVINCE_S != "NUNAVUT" & PROVINCE_S != "NORTHWEST TERRITORIES" & PROVINCE_S != "YUKON") %>%
  filter(WATER == 3) %>%
  left_join(bcr_n, by = c("BCR" = "bcr"))

bcr_route_density <- tm_shape(bcrs) + tm_polygons(col = "n_routes", palette = "YlGnBu", title = "BBS routes") +
  tm_layout(scale = 2)
tmap_save(bcr_route_density, "figures/bcr_route_density.pdf")

## Log abundance

log_abund <- cont_routes %>%
  group_by(countrynum, stateroute, year) %>%
  mutate(log_abund = log10(speciestotal)) %>%
  dplyr::select(stateroute, year, aou, log_abund)

# write.csv(log_abund, "data/derived_data/bbs_subset_1970-2016_logabund.csv", row.names = F)
