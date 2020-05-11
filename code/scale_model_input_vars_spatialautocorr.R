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

## Plot spatial autocorrelation for env vars, lines for different scales, facet by variable

spatial_corr <- read.csv("data/scale_model_vars_moransI.csv", stringsAsFactors = F)

spatial_corr_long <- spatial_corr %>%
  pivot_longer(lc:tmin, names_to = "env", values_to = "moransI") %>%
  filter(mean.of.class < 4000)

ggplot(spatial_corr_long, aes(x = mean.of.class, y = moransI, col = scale, group = scale)) + 
  geom_line() + 
  facet_wrap(~env) + theme_bw(base_size = 15) +
  geom_hline(yintercept = 0) +
  labs(x = "Distance (km)", y = "Moran's I", col = "Scale (routes)")
ggsave("figures/env_vars_moransI.pdf", units = "in", height = 6, width = 18)
