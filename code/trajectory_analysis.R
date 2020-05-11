### Community trajectory analysis
## Compare 1970-2016 to 1990-2016
## Scale model

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)
library(ecospat)

#### Set up ####

## Plotting theme

theme_set(theme_classic(base_size = 15))

## Read in data

# BBS sampled every 5 years from 1970 to 2016
log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)

# BBS sampled 3/4 years every 4 year time window
bbs_subset <- read.csv("data/bbs_counts_subset_1990-2016.csv", stringsAsFactors = F)

routes <- read.csv("/Volumes/hurlbertlab/databases/BBS/2017/bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# Land cover and climate data - 1/2 route scale

bbs_half_landcover <- read.csv('data/bbs_half_route_max_land_change.csv', stringsAsFactors = F)
bbs_half_climate <- read.csv("data/bbs_half_route_breeding_season_climate.csv", stringsAsFactors = F)

# Land cover and climate data - 1 route scale

bbs_landcover <- read.csv("data/bbs_route_max_landcover_change.csv", stringsAsFactors = F)
bbs_climate <- read.csv("data/bbs_routes_climate_trends.csv", stringsAsFactors = F)

## Raw land cover and climate data for calculating multiple scales

landcover_us <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F) %>%
  mutate(country = "US")
landcover_ca <- read.csv("/Volumes/hurlbertlab/dicecco/data/fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
  mutate(country = "Canada")

landcover_na <- bind_rows(landcover_us, landcover_ca) %>%
  mutate(min_year = case_when(country == "US" ~ 1992,
                              TRUE ~ 1990),
         max_year = case_when(country == "US" ~ 2016,
                              TRUE ~ 2010))

## Land cover legend - classes in common between US and Canadian classifications

landcover_common_legend <- read.csv("data/landcover_code_common_legend_US_Canada.csv", stringsAsFactors = F)

# Function to calculate the class with the maximum change in proportion landscape from earliest to late time period
## If multiple stateroutes are present (for scale model) - sums area across classes pooling stateroutes, calculates prop landscape for summed area
max_delta <- function(df) {
  delta <- df %>%
    filter(year == min_year | year == max_year) %>%
    mutate(year_name = case_when(year == min_year ~ "year1",
                                 year == max_year ~ "year2")) %>%
    left_join(landcover_common_legend, by = c("class" = "Code")) %>%
    filter(!(is.na(class))) %>%
    group_by(year_name, common_label) %>%
    summarize(sum.total.area = sum(total.area)/(mean(max_year) - mean(min_year))) %>%
    group_by(year_name) %>%
    mutate(prop.landscape = sum.total.area/sum(sum.total.area)) %>%
    select(year_name, common_label, prop.landscape) %>%
    distinct() %>%
    pivot_wider(names_from = year_name, values_from = prop.landscape) %>%
    replace_na(list(year1 = 0, year2 = 0)) %>%
    mutate(deltaCover = year2 - year1) %>%
    mutate(absCover = abs(deltaCover)) %>%
    filter(absCover == max(absCover)) %>%
    select(common_label, deltaCover)
  
  return(delta)
}

possibly_max_delta <- possibly(max_delta, data.frame(common_label = NA, deltaCover = NA))

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
    abund_dist <- dist(df[, -c(1:2)])
    trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year)
  })) %>%
  mutate(dir25 = purrr::map_dbl(data, ~{
    df <- . 
    df <- df %>%
      filter(year >= 1990)
    if(max(df$year) - min(df$year) > 15) {
      abund_dist <- dist(df[, -c(1:2)])
      trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year)
    } else(NA)
  }))

r <- round(cor(logabund_wide$dir50, logabund_wide$dir25), 2)

ggplot(logabund_wide, aes(x = dir50, y = dir25)) + geom_point() + 
  labs(x = "Directionality 1970-2016", y = "Directionality 1990-2016") + 
  geom_abline(intercept = 0, slope = 1, cex = 1) +
  annotate(geom= "text", x = 0.35, y = 0.3, label = paste0("r = ", r), size = 8)
# ggsave("figures/directionality_time_series_comparison.pdf")

## Model of directionality ~ land cover change + climate change
# 935 routes

# Coordinates of BBS routes

coords_bbs <- bbs_subset %>%
  ungroup() %>%
  dplyr::select(stateroute) %>%
  distinct() %>%
  left_join(routes) %>%
  dplyr::select(stateroute, bcr, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4236)

bcr_sample_size <- coords_bbs %>%
  group_by(bcr) %>%
  summarize(n_routes = n_distinct(stateroute))

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

# units in m
bbs_route_distances <- map_dfr(unique(bbs_subset$stateroute), ~dist_bcr(.))

stateroutes <- data.frame(stateroute = unique(bbs_subset$stateroute))
  
bbs_distances_km <- bbs_route_distances %>%
  mutate(dist_km = distance/1000)

mean_bbs_distances <- stateroutes %>%
  mutate(scale = map(stateroute, ~data.frame(scale = rep(1:25)))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_distances_km %>%
      filter(focal_rte == .x) %>%
      arrange(dist_km) %>%
      slice(1:.y)
  }),
  mean_dist = map_dbl(model_input, ~mean(.$dist_km)))
  
mean_distances <- mean_bbs_distances %>%
  group_by(scale) %>%
  summarize(mean_dist_rtes = mean(mean_dist))

mean_bcr_distances <- mean_bbs_distances %>%
  left_join(coords_bbs) %>%
  group_by(bcr, scale) %>%
  summarize(mean_dist_bcr = mean(mean_dist))

bcr_subset <- bcr_sample_size %>%
  filter(n_routes > 25)

ggplot(filter(mean_bcr_distances, bcr %in% bcr_subset$bcr), 
       aes(x = scale, y = mean_dist_bcr, col = bcr, group = bcr)) + 
  geom_line() +
  scale_color_viridis_c() +
  labs(x = "Scale (routes)", y = "Mean distance between routes (km)", col = "BCR")
# ggsave("figures/scale_aggregated_routes.pdf")

# % overlap of aggregated stateroutes at each spatial scale x bcr

pct_overlap <- mean_bbs_distances %>%
  filter(scale > 1) %>%
  select(-stateroute) %>%
  unnest(model_input) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  group_by(scale, bcr) %>%
  count(stateroute) %>%
  left_join(bcr_subset) %>%
  mutate(pct_reps = n/n_routes)

mean_pct_overlap <- pct_overlap %>%
  group_by(bcr, scale) %>%
  summarize(mean_reps = mean(pct_reps))

ggplot(mean_pct_overlap, aes(x = scale, y = mean_reps, col = as.factor(bcr), group = as.factor(bcr))) + 
  geom_line() +
  labs(x = "Scale (routes)", y = "Avg percent of aggregates stateroute occurs in", col = "BCR")
ggsave("figures/percent_overlap_aggregates.pdf")

## Model for 1/2 routes

#### Need half BBS data at stop level!!

## Model from 1 route up to nearest 25 routes
## One data point per 4 year time window
## ID core and transient species

# All species
log_abund_wider <- bbs_subset %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016),
         max_bins = case_when(countrynum == 124 ~ 5,
                              countrynum == 840 ~ 6)) %>%
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
  group_by(stateroute, aou, year_bin) %>%
  summarize(mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  dplyr::select(stateroute, aou, year_bin, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

# Excluding transient species
log_abund_core <- bbs_subset %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016),
         max_bins = case_when(countrynum == 124 ~ 5,
                              countrynum == 840 ~ 6)) %>%
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
  group_by(stateroute, aou, year_bin) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, year_bin, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

#### Scale model ####

## At each scale (1 route, up to nearest 25 routes within BCR) 
## Get max land cover delta from raw land cover data and get trend in Tmin and trend in Tmax
## Calculate directionality for each grouping of routes: average log(abundance) across routes when pooled
## Fit models across BBS routes for each scale (1:25)
## Variance partitioning of land cover and climate variables explaining trajectory 
## Directionality ~ tmin + tmax + max(deltaLandCover)
## Variance partitioning with ecospat.varpat(model.1, model.2, model.12) in ecospat

na_climate <- filter(bbs_climate_avgs, is.na(mean_tmax) | is.na(mean_tmin))

climate_trend <- function(climate_df) {
  trend_tmax <- coef(lm(mean_tmax ~ year, data = climate_df))[[2]]
  trend_tmin <- coef(lm(mean_tmin ~ year, data = climate_df))[[2]]
  
  return(list(trend_tmax = trend_tmax, trend_tmin = trend_tmin))
}
possibly_climate_trend <- possibly(climate_trend, list(trend_tmax = NA, trend_tmin = NA))

scale_model_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = map(stateroute, ~data.frame(scale = rep(1:25)))) %>%
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

    max_lc_delta <- possibly_max_delta(land_cover)

    climate <- bbs_climate_avgs %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year) %>%
      summarize(mean_tmax = mean(mean_tmax, na.rm = T),
                mean_tmin = mean(mean_tmin, na.rm = T))

    climate_trend <- possibly_climate_trend(climate)
    
    log_abund <- log_abund_wider %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin) %>%
      summarize_all(mean, na.rm = T) %>%
      dplyr::select(-stateroute)
    abund_dist <- dist(log_abund[, -1])
    dir_all <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(log_abund)), surveys = log_abund$year_bin)
    
    log_core <- log_abund_core %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin) %>%
      summarize_all(mean, na.rm = T) %>%
      dplyr::select(-stateroute)
    abund_dist_core <- dist(log_core[, -1])
    dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
    
    data.frame(focal_rte = unique(df$focal_rte), max_lc_class = max_lc_delta$common_label, max_lc = max_lc_delta$deltaCover,
               trend_tmax = climate_trend$trend_tmax, trend_tmin = climate_trend$trend_tmin,
               dir_all = dir_all, dir_core = dir_core)
  }))

scale_model_variables <- scale_model_input %>%
  dplyr::select(scale, input_vars) %>%
  unnest(cols = c(input_vars)) %>%
  group_by(scale) %>%
  nest()

# Join 1/2 route data

scale_model_variables_unnest <- scale_model_variables %>%
  unnest(cols = c(data))
# write.csv(scale_model_variables_unnest, "data/scale_model_input.csv", row.names = F)

ggplot(scale_model_variables_unnest, aes(x = scale, y = dir_all, col = focal_rte, group = focal_rte)) + 
  geom_line(alpha = 0.1) + labs(x = "Scale (routes)", y = "Directionality (all spp.)", col = "Stateroute")
ggsave("figures/directionality_scale.pdf")

ggplot(filter(scale_model_variables_unnest, scale == 25), aes(x = dir_all, y = dir_core)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Directionality (all spp.)", y = "Directionality (excl. transients)", title = "25 routes")
ggsave("figures/directionality_all_vs_core_25.pdf")

scale_model_output <- data.frame(scale = c(), part1 = c(), part2 = c(), joined = c(), unexpl = c())
for(i in 1:25) {
  model_input <- scale_model_variables$data[[i]]
  
  mod1 <- lm(dir_core ~ trend_tmax + trend_tmin, data = model_input)
  mod2 <- lm(dir_core ~ max_lc, data = model_input)
  mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  mod1.r2 <- summary(mod1)$r.squared
  mod2.r2 <- summary(mod2)$r.squared
  mod12.r2 <- summary(mod12)$r.squared
  
  part1 <- mod12.r2 - mod2.r2 # climate alone
  part2 <- mod12.r2 - mod1.r2 # land cover alone
  joined <- mod1.r2 - part1 # shared variance
  unexpl <- 1 - mod12.r2 # unexplained variance
  
  scale_model_output <- rbind(scale_model_output, 
        data.frame(scale = i, part1 = part1, part2 = part2,
                   joined = joined, unexpl = unexpl))
}
# write.csv(scale_model_output, "data/scale_model_output_deviance.csv", row.names = F)

scale_model_plot <- scale_model_output %>%
  pivot_longer(part1:unexpl, names_to = "deviance")

ggplot(filter(scale_model_plot, deviance == "part1" | deviance == "part2"), aes(x = scale, y = value, fill = deviance)) +
  geom_col(position = "stack") +
  scale_fill_discrete(name = "Variance explained", labels = c("Climate", "Land cover")) +
  labs(x = "Aggregated routes", y = "Variance")
ggsave("figures/scale_model_variance.pdf")
