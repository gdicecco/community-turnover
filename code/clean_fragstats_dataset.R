### Clean longleaf output: half-BBS route FRAGSTATS

library(tidyverse)
library(purrr)

# US routes

us_path <- "/proj/hurlbertlab/gdicecco/nlcd_bbs_half_route_fragstats/"

us_routes <- list.files(us_path) %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename(filename = ".")

us_fragstats <- us_routes %>%
  mutate(country = "US",
         year = word(filename, 4, 4, sep = "_"),
         stateroute = word(filename, 6, 6, sep = "_"),
         stops = word(filename, 7, 7, sep = "_")) %>%
  mutate_at(c("stops"), ~word(., 1, 1, sep = "\\.")) %>%
  group_by(filename, country, year, stateroute, stops) %>%
  nest() %>%
  mutate(data = map(filename, ~read_csv(paste0(us_path, .)))) %>%
  unnest(cols = c(data))

# Canada routes

ca_path <- "/proj/hurlbertlab/gdicecco/canada_bbs_half_route_fragstats/"

ca_routes <- list.files(ca_path) %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename(filename = ".")

ca_fragstats <- ca_routes %>%
  mutate(country = "Canada",
         year = word(filename, 4, 4, sep = "_"),
         stateroute = word(filename, 11, 11, sep = "_"),
         stops = word(filename, 12, 12, sep = "_")) %>%
  mutate_at(c("stops"), ~word(., 1, 1, sep = "\\.")) %>%
  group_by(filename, country, year, stateroute, stops) %>%
  nest() %>%
  mutate(data = map(filename, ~read_csv(paste0(ca_path, .)))) %>%
  unnest(cols = c(data))

# Join for North America dataset

fragstats <- bind_rows(us_fragstats, ca_fragstats)
write.csv(fragstats, "/proj/hurlbertlab/gdicecco/bbs_half_routes_fragstats.csv", row.names = F)

