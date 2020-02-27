### Species occurrence analysis: colonizations and extinctions

library(tidyverse)
library(purrr)

## Plotting theme

theme_set(theme_classic(base_size = 15))

## Read in data

log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)
routes_subs <- read.csv("data/bbs_route_subset_1990-2016.csv", stringsAsFactors = F)

log_abund_subs <- log_abund %>%
  filter(stateroute %in% routes_subs$stateroute) %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016)) %>%
  filter(year >= y1, year <= y2)

## Species occurrence over time - split time series in half to measure colonizations and extinctions

spp_occ <- log_abund_subs %>%
  mutate(decade = case_when(countrynum == 124 & year >= 1990 & year <= 1999 ~ "d1",
                        countrynum == 840 & year >= 1992 & year <= 2001 ~ "d1",
                        countrynum == 124 & year >= 2001 & year <= 2010 ~ "d2",
                        countrynum == 840 & year >= 2007 & year <= 2016 ~ "d2",
                        TRUE ~ "NA")) %>%
  filter(decade == "d1" | decade == "d2") %>%
  group_by(aou, stateroute, decade) %>%
  summarize(n_years = n_distinct(year),
            occ = n_years/10)
