### Species mean temperature range (5th percentile to 95th percentile)

library(tidyverse)

## List of species observed in BCRs of interest during time window (1990-present)
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, select = c("statenum","stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

landbirds <- species %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

# Observations with RT=1, 1990-present, diurnal land birds
counts.subs <- counts %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  filter(rpid == 101) %>%
  filter(year >= 1990) %>%
  filter(aou %in% landbirds$aou)

length(unique(counts.subs$aou))

require(raster)
require(maps)

# WorldClim
climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
climatelayers_ss = climatelayers[[10]]

climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))

## Species observed on 10 or more routes
birds.subs <- counts.subs %>%
  group_by(aou) %>%
  summarize(nRoutes = length(unique(stateroute))) %>%
  filter(nRoutes >= 10) %>%
  left_join(counts.subs)

bird_temps <- birds.subs %>%
  group_by(aou) %>%
  dplyr::select(stateroute) %>%
  unique() %>%
  left_join(routes) %>%
  dplyr::select(stateroute, longitude, latitude) %>%
  nest() %>%
  mutate(temp_range = purrr::map_dbl(data, ~{
    df <- .
    climate <- raster::extract(climatelayers_ss_cropped, dplyr::select(df, longitude, latitude))
    quants <- quantile(climate, c(0.05, 0.95), na.rm = T)
    print(quants)
    (quants[[2]] - quants[[1]])/10
    }),
    temp_mean = purrr::map_dbl(data, ~{
      df <- .
      climate <- raster::extract(climatelayers_ss_cropped, dplyr::select(df, longitude, latitude))
      mean(climate/10, na.rm = T)
    })) %>%
  dplyr::select(-data)

write.csv(bird_temps, "data/derived_data/bbs_aou_temp_range.csv", row.names = F)
