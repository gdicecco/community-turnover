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
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

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

# Filter BBS to rpid = 101, runtype = 1, land birds != birds of prey
counts.subs <- counts %>%
  filter(rpid == 101) %>%
  filter(aou %in% species_list$aou) %>%
  right_join(RT1.routes, by = c("countrynum", "statenum", "stateroute", "year"))

# BBS 50 stop data for 1/2 route estimates
# Data through 2018
wd <- paste0(bioark, "/HurlbertLab/Databases/BBS/FiftyStopDataThru2018/")

files <- list.files(path = wd)
data_files <- files[grepl("fifty[-0-9]", files)]

fifty_stop <- data.frame(filename = data_files) %>%
  group_by(filename) %>%
  nest() %>%
  mutate(data = map(filename, ~read.csv(paste0(wd, .)))) %>%
  unnest(data) %>%
  ungroup()

fifty_stop$stops1_25 <- rowSums(select(fifty_stop, Stop1:Stop25))
fifty_stop$stops26_50 <- rowSums(select(fifty_stop, Stop26:Stop50))

half_route_counts <- fifty_stop %>%
  select(RouteDataID, CountryNum, StateNum, Route, RPID, Year, AOU, stops1_25, stops26_50) %>%
  mutate(stateroute = StateNum*1000 + Route) %>%
  filter(RPID == 101, AOU %in% species_list$aou) %>%
  filter(Year >= 1990) %>%
  inner_join(RT1.routes, by = c("CountryNum" = "countrynum", "StateNum" = "statenum", "stateroute", "Year" = "year")) %>%
  pivot_longer(stops1_25:stops26_50, names_to = "stops", values_to = "count")

# All species
log_abund_wider <- half_route_counts %>%
  mutate(y1 = case_when(CountryNum == 124 ~ 1990,
                        CountryNum == 840 ~ 1992),
         y2 = case_when(CountryNum == 124 ~ 2010,
                        CountryNum == 840 ~ 2016),
         max_bins = case_when(CountryNum == 124 ~ 5,
                              CountryNum == 840 ~ 6)) %>%
  filter(Year >= y1, Year <= y2) %>%
  group_by(CountryNum) %>%
  nest() %>%
  mutate(year_bins = map2(CountryNum, data, ~{
    country <- .x
    df <- .y
    
    if(country == 124) {
      df %>%
        mutate(year_bin = case_when(Year >= 1990 & Year <= 1993 ~ 1990,
                                    Year >= 1994 & Year <= 1997 ~ 1994,
                                    Year >= 1998 & Year <= 2001 ~ 1998,
                                    Year >= 2002 & Year <= 2005 ~ 2002,
                                    TRUE ~ 2006))
    } else {
      df %>%
        mutate(year_bin = case_when(Year >= 1992 & Year <= 1995 ~ 1992,
                                    Year >= 1996 & Year <= 1999 ~ 1996,
                                    Year >= 2000 & Year <= 2003 ~ 2000,
                                    Year >= 2004 & Year <= 2007 ~ 2004,
                                    Year >= 2008 & Year <= 2011 ~ 2008,
                                    Year >= 2012 & Year <= 2016 ~ 2012))
    }
  })) %>%
  select(-data) %>%
  unnest(cols = c(year_bins)) %>%
  group_by(stateroute, stops, AOU, year_bin) %>%
  summarize(mean_abund = mean(count) + 1,
            log_abund = log10(mean_abund)) %>%
  dplyr::select(stateroute, stops, AOU, year_bin, log_abund) %>%
  pivot_wider(names_from = AOU, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

# Excluding transient species
log_abund_core <- half_route_counts %>%
  mutate(y1 = case_when(CountryNum == 124 ~ 1990,
                        CountryNum == 840 ~ 1992),
         y2 = case_when(CountryNum == 124 ~ 2010,
                        CountryNum == 840 ~ 2016),
         max_bins = case_when(CountryNum == 124 ~ 5,
                              CountryNum == 840 ~ 6)) %>%
  filter(Year >= y1, Year <= y2) %>%
  group_by(CountryNum) %>%
  nest() %>%
  mutate(year_bins = map2(CountryNum, data, ~{
    country <- .x
    df <- .y
    
    if(country == 124) {
      df %>%
        mutate(year_bin = case_when(Year >= 1990 & Year <= 1993 ~ 1990,
                                    Year >= 1994 & Year <= 1997 ~ 1994,
                                    Year >= 1998 & Year <= 2001 ~ 1998,
                                    Year >= 2002 & Year <= 2005 ~ 2002,
                                    TRUE ~ 2006))
    } else {
      df %>%
        mutate(year_bin = case_when(Year >= 1992 & Year <= 1995 ~ 1992,
                                    Year >= 1996 & Year <= 1999 ~ 1996,
                                    Year >= 2000 & Year <= 2003 ~ 2000,
                                    Year >= 2004 & Year <= 2007 ~ 2004,
                                    Year >= 2008 & Year <= 2011 ~ 2008,
                                    Year >= 2012 & Year <= 2016 ~ 2012))
    }
  })) %>%
  select(-data) %>%
  unnest(cols = c(year_bins)) %>%
  group_by(stateroute, stops, AOU, year_bin) %>%
  summarize(n_years = n_distinct(Year),
            mean_abund = mean(count) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, stops, AOU, year_bin, log_abund) %>%
  pivot_wider(names_from = AOU, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

half_route_dir <- half_route_counts %>%
  ungroup() %>%
  distinct(stateroute, stops) %>%
  mutate(scale = 0.5,
         dir_vals = map2(stateroute, stops, ~{
          strte <- .x
          stps <- .y
        
           log_abund <- log_abund_wider %>%
             filter(stateroute == strte, stops == stps) %>%
             ungroup() %>%
             select(-stateroute, -stops) %>%
             group_by(year_bin) %>%
             summarize_all(mean, na.rm = T)
             
           if(nrow(log_abund) > 3) {
             abund_dist <- dist(log_abund[, -1])
             dir_all <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(log_abund)), surveys = log_abund$year_bin)
           } else {data.frame(dir_all = NA, dir_core = NA)}
             
             log_core <- log_abund_core %>%
               filter(stateroute == strte, stops == stps) %>%
               ungroup() %>%
               select(-stateroute, -stops) %>%
               group_by(year_bin) %>%
               summarize_all(mean, na.rm = T)
             
             if(nrow(log_core) > 3) {
             abund_dist_core <- dist(log_core[, -1])
             dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
             
             data.frame(dir_all = dir_all, dir_core = dir_core)
             } else {data.frame(dir_all = NA, dir_core = NA)}

         }))

half_route_dir_write <- half_route_dir %>%
  unnest(dir_vals) %>%
  filter(!is.na(dir_all))
write.csv(half_route_dir_write, "data/half_route_directionality.csv", row.names = F)

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

# write.csv(log_abund, "data/bbs_subset_1970-2016_logabund.csv", row.names = F)
