### Community trajectory analysis
## Compare 1970-2016 to 1990-2016

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)

## Plotting theme

theme_set(theme_classic(base_size = 15))

## Read in data

log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)

routes <- read.csv("/Volumes/hurlbertlab/databases/BBS/2017/bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# Land cover and climate data

bbs_landcover <- read.csv("data/bbs_route_max_landcover_change.csv", stringsAsFactors = F)
bbs_climate <- read.csv("data/bbs_routes_climate_trends.csv", stringsAsFactors = F)

## Community trajectories
# https://cran.r-project.org/web/packages/vegclust/vignettes/CTA.html

# long to wide

logabund_wide <- log_abund %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir50 = purrr::map_dbl(data, ~{
    df <- .
    abund_dist <- dist(df)
    trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year)
  })) %>%
  mutate(dir25 = purrr::map_dbl(data, ~{
    df <- . 
    df <- df %>%
      filter(year >= 1990)
    if(max(df$year) - min(df$year) > 15) {
      abund_dist <- dist(df)
      trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year)
    } else(NA)
  }))

r <- round(cor(logabund_wide$dir50, logabund_wide$dir25), 2)

ggplot(logabund_wide, aes(x = dir50, y = dir25)) + geom_point() + 
  labs(x = "Directionality 1970-2016", y = "Directionality 1990-2016") + 
  geom_abline(intercept = 0, slope = 1, cex = 1) +
  annotate(geom= "text", x = 0.78, y = 0.7, label = paste0("r = ", r), size = 8)
ggsave("figures/directionality_time_series_comparison.pdf")

## Model of directionality ~ land cover change + climate change

dir_model <- logabund_wide %>%
  select(-data, -dir50) %>%
  left_join(bbs_climate) %>%
  left_join(bbs_landcover) %>%
  filter(env == "tmin") %>%
  na.omit() # why are there NAs - should be able to get all routes for climate
# 629 routes
# write.csv(dir_model, "data/directionality_mod_input.csv", row.names = F)

# Spatial CAR model -- Talk to Allen/James - this needs to be GLM (directionality is bounded 0-1)
## Possibly CARBayes?
# get species weights matrix

# Coordinates of BBS routes

coords_bbs <- dir_model %>%
  ungroup() %>%
  dplyr::select(stateroute) %>%
  distinct() %>%
  left_join(routes) %>%
  dplyr::select(stateroute, longitude, latitude)

coords_mat <- as.matrix(coords_bbs[, -1])

# Calculate nearest neighbors and make spatial neighborhood
k0 <- knearneigh(coords_mat, longlat=T, k=1)
k1 <- knn2nb(k0)

# Find maximum neighbor distance and use this distance to define neighborhoods
max.d <- max(unlist(nbdists(k1, coords_mat, longlat=T)))
nb0.birds <- dnearneigh(coords_mat, 0, max.d, longlat=T)
plot(nb0.birds, coords_mat)

# Using a distance threshold of 100km
nb1.birds <- dnearneigh(coords_mat,1,100,longlat=T)
plot(nb1.birds, coords_mat)

# Create spatial weights based on linear distance decay
glist.birds <- nbdists(nb0.birds, coords_mat, longlat=T)
glist.birds <- lapply(glist.birds, function(x) 1-(x/ceiling(max.d))) # Round up max.d so that all point get weights
wt.birds <- nb2listw(nb0.birds, style='B', glist=glist.birds)

dir_mod <- spautolm(dir25 ~ climateTrend + deltaCover, data = dir_model, wt.birds, na.action = na.fail, family = "CAR")
