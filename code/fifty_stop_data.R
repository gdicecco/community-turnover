## 50 stop BBS data - aggregate species counts to 1/2 routes

library(tidyverse)

species_list <- read.csv("data/species_list.csv")

wd <- "/Volumes/hurlbertlab/Databases/BBS/FiftyStopData/"

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
  select(RouteDataID, countrynum, statenum, Route, RPID, year, AOU, stops1_25, stops26_50) %>%
  mutate(stateroute = statenum*1000 + Route) %>%
  filter(RPID == 101, AOU %in% species_list$aou)
