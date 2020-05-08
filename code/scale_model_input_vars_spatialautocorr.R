## Scale model input spatial autocorrelation

library(tidyverse)
library(ncf)
library(purrr)

# Read in data

scale_model_input <- read.csv("/proj/hurlbertlab/gdicecco/community-turnover/data/scale_model_input.csv", stringsAsFactors = F)

routes <- read.csv("/proj/hurlbertlab/bbs/bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# Moran's I for env variables, each scale

spatial_corr <- scale_model_input %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  group_by(scale) %>%
  nest() %>%
  mutate(moransI = map(data, ~{
    df <- .
    cor_lc <- correlog(df$longitude, df$latitude, df$max_lc, increment = 250, latlon = T, na.rm = T)
    cor_tmax <- correlog(df$longitude, df$latitude, df$trend_tmax, increment = 250, latlon = T, na.rm = T)
    cor_tmin <- correlog(df$longitude, df$latitude, df$trend_tmin, increment = 250, latlon = T, na.rm = T)

    data.frame(mean.of.class = cor_lc$mean.of.class,
               lc = cor_lc$correlation,
               tmax = cor_tmax$correlation,
               tmin = cor_tmin$correlation)

  })) %>%
  select(-data) %>%
  unnest(moransI)

write.csv(spatial_corr, "/proj/hurlbertlab/gdicecco/community-turnover/data/scale_model_vars_moransI.csv", row.names = F)
