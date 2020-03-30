### Community trajectory analysis
## Compare 1970-2016 to 1990-2016

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)
library(ecospat)

## Plotting theme

theme_set(theme_classic(base_size = 15))

## Read in data

log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)

routes <- read.csv("/Volumes/hurlbertlab/databases/BBS/2017/bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# Land cover and climate data - 1 route scale

bbs_landcover <- read.csv("data/bbs_route_max_landcover_change.csv", stringsAsFactors = F)
bbs_climate <- read.csv("data/bbs_routes_climate_trends.csv", stringsAsFactors = F)

## Land cover and climate data for calculating multiple scales

landcover_us <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F) %>%
  mutate(country = "US")
landcover_ca <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
  mutate(country = "Canada")

landcover_na <- bind_rows(landcover_us, landcover_ca) %>%
  mutate(min_year = case_when(country == "US" ~ 1992,
                              TRUE ~ 1990),
         max_year = case_when(country == "US" ~ 2016,
                              TRUE ~ 2010))

# Function to calculate the class with the maximum change in proportion landscape from earliest to late time period
## If multiple stateroutes are present (for scale model) - sums area across classes pooling stateroutes, calculates prop landscape for summed area
max_delta <- function(df) {
  delta <- df %>%
    filter(year == min_year | year == max_year) %>%
    mutate(year_name = case_when(year == min_year ~ "year1",
                                 year == max_year ~ "year2")) %>%
    group_by(year_name, class) %>%
    summarize(sum.total.area = sum(total.area)) %>%
    group_by(year_name) %>%
    mutate(prop.landscape = sum.total.area/sum(sum.total.area)) %>%
    select(year_name, class, prop.landscape) %>%
    distinct() %>%
    pivot_wider(names_from = year_name, values_from = prop.landscape) %>%
    replace_na(list(year1 = 0, year2 = 0)) %>%
    mutate(deltaCover = year2 - year1) %>%
    mutate(absCover = abs(deltaCover)) %>%
    filter(absCover == max(absCover)) %>%
    select(class, deltaCover)
  
  return(delta)
}

possibly_max_delta <- possibly(max_delta, data.frame(class = NA, deltaCover = NA))

## Climate trends

bbs_climate_avgs <- read.csv("data/bbs_routes_breeding_season_climate.csv", stringsAsFactors = F) %>%
  left_join(routes) %>%
  mutate(min_year = case_when(countrynum == 840 ~ 1992,
                              countrynum == 124 ~ 1990),
         max_year = case_when(countrynum == 840 ~ 2016,
                              countrynum == 124 ~ 2010)) %>%
  filter(countrynum == 840 | countrynum == 124) %>%
  filter(year >= min_year & year <= max_year)

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
  na.omit() ## Need to ID why there are NAs in this - all routes should be present in climate/land cover
# 629 routes

# GLM model -- Talk to Allen/James - this needs to be spatial CAR + GLM (directionality is bounded 0-1)
## Possibly CARBayes?

# Coordinates of BBS routes

coords_bbs <- dir_model %>%
  ungroup() %>%
  dplyr::select(stateroute) %>%
  distinct() %>%
  left_join(routes) %>%
  dplyr::select(stateroute, bcr, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))

## Find nearest grouped BBS routes within BCR

dist_bcr <- function(strte) {
  focal_rte <- coords_bbs %>%
    filter(stateroute == strte)
  
  bcr_rtes <- coords_bbs %>%
    filter(bcr == focal_rte$bcr)
  
  dist <- st_distance(bcr_rtes, focal_rte)
  
  dist_df <- bcr_rtes %>%
    st_set_geometry(NULL) %>%
    mutate(distance = dist[, 1], 
           focal_rte = strte)
  
  return(dist_df)
}

bbs_route_distances <- map_dfr(unique(dir_model$stateroute), ~dist_bcr(.))

## Model for 1/2 routes

## Model from 1 route up to nearest 30 routes

log_abund_wider <- log_abund %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean))

stateroutes <- data.frame(stateroute = unique(dir_model$stateroute))

## At each scale (1 route, up to nearest 30 routes within BCR) 
## Get max land cover delta from raw land cover data and get trend in Tmin and trend in Tmax
## Calculate directionality for each grouping of routes: average log(abundance) across routes when pooled
## Fit models across BBS routes for each scale (1:30)
## Variance partitioning of land cover and climate variables explaining trajectory 
## Directionality ~ tmin + tmax + max(deltaLandCover)
## Variance partitioning with ecospat.varpat(model.1, model.2, model.12) in ecospat

scale_model_input <- stateroutes %>%
  mutate(scale = map(stateroute, ~data.frame(scale = rep(1:30)))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = map(model_input, ~{
    df <- .
    
    land_cover <- landcover_na %>%
      filter(stateroute %in% df$stateroute)
      
    max_delta <- possibly_max_delta(land_cover)
    
    climate <- bbs_climate_avgs %>%
      filter(stateroute %in% df$stateroute)
    
    trend_tmax <- coef(lm(mean_tmax ~ year, data = climate))[[2]]
    trend_tmin <- coef(lm(mean_tmin ~ year, data = climate))[[2]]
    
    log_abund <- log_abund_wider %>%
      filter(stateroute %in% df$stateroute) %>%
      dplyr::select(-countrynum) %>%
      group_by(year) %>%
      summarize_all(mean, na.rm = T)
    abund_dist <- dist(log_abund)
    dir25 <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(log_abund)), surveys = log_abund$year)
    
    data.frame(focal_rte = unique(df$focal_rte), max_lc = max_delta, trend_tmax = trend_tmax, trend_tmin = trend_tmin, dir25 = dir25)
  
  }))

scale_model_variables <- scale_model_input %>%
  dplyr::select(scale, input_vars) %>%
  unnest(cols = c(input_vars)) %>%
  group_by(scale) %>%
  nest()

# scale_model_variables_unnest <- scale_model_variables %>%
#   unnest(cols = c(data))
# write.csv(scale_model_variables_unnest, "data/scale_model_input.csv", row.names = F)

scale_model_output <- data.frame(scale = c(), part1 = c(), part2 = c(), joined = c(), unexpl = c())
for(i in 1:30) {
  model_input <- scale_model_variables$data[[i]]
  
  mod1 <- glm(dir25 ~ trend_tmax + trend_tmin, family = "binomial", data = model_input)
  mod2 <- glm(dir25 ~ max_lc.deltaCover, family = "binomial", data = model_input)
  mod12 <- glm(dir25 ~ trend_tmax + trend_tmin + max_lc.deltaCover, family = "binomial", data = model_input)
  
  varpart <- ecospat.varpart(model.1 = mod1, model.2 = mod2, model.12 = mod12)
  
  scale_model_output <- rbind(scale_model_output, 
        data.frame(scale = i, part1 = varpart[[1]], part2 = varpart[[2]],
                   joined = varpart[[3]], unexpl = varpart[[4]]))
}
# write.csv(scale_model_output, "data/scale_model_output_deviance.csv", row.names = F)

scale_model_plot <- scale_model_output %>%
  pivot_longer(part1:unexpl, names_to = "deviance")

ggplot(filter(scale_model_plot, deviance == "part1" | deviance == "part2"), aes(x = scale, y = value, fill = deviance)) +
  geom_col(position = "stack") +
  scale_fill_discrete(name = "Deviance explained", labels = c("Climate", "Land cover")) +
  labs(x = "Aggregated routes", y = "Deviance")
ggsave("figures/scale_model_deviance.pdf")
