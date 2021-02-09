### Route-level abundance trends for full species list

library(tidyverse)
library(purrr)
library(broom)

# Append correct BioArk path
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

# Bird population data
## BBS 2017 Version
routes <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_routes_20170712.csv"))
counts <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_counts_20170712.csv"))
species <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_species_20170712.csv"))
weather <- read.csv(paste0(bioark, "/HurlbertLab/Databases/BBS/2017/bbs_weather_20170712.csv"))

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# Diurnal land birds
species_list <- read.csv("data/species_list.csv", stringsAsFactors = F)

# Routes surveyed 3/4 years in every 4 year window
bbs_subset <- read.csv("data/bbs_route_subset_1990-2016.csv", stringsAsFactors = F)

counts.subs <- counts %>%
  filter(aou %in% species_list$aou) %>%
  merge(filter(RT1.routes, stateroute %in% bbs_subset$stateroute), by = c("stateroute", "year")) %>%
  filter(rpid == 101) %>%
  filter(year >= 1990, year < 2017)

# first year observers

first_year <- function(obsns) {
  first_yrs <- c()
  for(i in 1:length(obsns)) {
    ob <- obsns[i]
    obs <- obsns[1:i]
    first_yrs <- c(first_yrs, ifelse(length(obs[obs == ob]) == 1, 1, 0))
  }
  return(first_yrs)
}

obs_years <- weather %>%
  group_by(stateroute) %>%
  arrange(year) %>%
  nest() %>%
  mutate(first_yr = purrr::map(data, ~{
    df <- .
    first_year(df$obsn)
  })) %>%
  unnest() %>%
  dplyr::select(stateroute, year, first_yr)

# abundance trends read in
# 1990-2010 for Canada, 1992-2016 for US
abund_trend <- counts.subs %>%
  left_join(obs_years, by = c('stateroute', 'year')) %>%
  group_by(aou, stateroute) %>%
  nest() %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    country <- unique(df$countrynum.x)
    if(country == 840){
      df <- df %>%
        filter(year >= 1992 & year <= 2016) %>%
        unique()
      nrow(df)
    } else {
      df <- df %>%
        filter(year >= 1990 & year <= 2010) %>%
        unique()
      nrow(df)
    }
  })) %>%
  filter(nObs > 9) %>%
  mutate(lmFit = purrr::map(data, ~{
    df <- .
    country <- unique(df$countrynum.x)
    if(country == 840){
      df.short <- df %>%
        dplyr::select(year, first_yr, speciestotal) %>%
        filter(year >= 1992 & year <= 2016) %>%
        unique()
      glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
    } else {
      df.short <- df %>%
        dplyr::select(year, first_yr, speciestotal) %>%
        filter(year >= 1990 & year <= 2010) %>%
        unique()
      glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
    }
  }))  %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(abundTrend = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[2]
  })) %>%
  mutate(trendInt = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[1]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[2]
  })) %>%
  mutate(obsTrend = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[3]
  })) %>%
  mutate(obsPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[3]
  }))

abund_trend_write <- abund_trend %>%
  select(-data, -lmFit, -lm_broom)
write.csv(abund_trend_write, "data/derived_data/BBS_abundance_trends.csv", row.names = F)
