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
ggsave("figures/directionality_time_series_comparison.pdf")

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
ggsave("figures/scale_aggregated_routes.pdf")

## Model for 1/2 routes

#### Need half route climate, land cover, BBS data at stop level!!

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

# Directionality at 1 route scale, core vs. transients

dir_all_spp <- log_abund_wider %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_all = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$year_bin)
  }))

dir_core_spp <- log_abund_core %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_core = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$year_bin)
  }))

dir_compare <- dir_all_spp %>%
  dplyr::select(-data) %>%
  left_join(dir_core_spp) %>%
  dplyr::select(-data)

ggplot(dir_compare, aes(x = dir_all, y = dir_core)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) + labs(x = "Directionality (all spp)", y = "Directionality (excl. transients)")
ggsave("figures/directionality_all_vs_core.pdf")

### How much variance is there in directionality based on delineation of time windows?
## Pull stateroutes with full time series, try multiple four year time windows

all_years <- bbs_subset %>% 
  group_by(stateroute) %>% 
  summarize(n_years = n_distinct(year)) %>%
  filter(n_years == 26)

yearbins_test <- bbs_subset %>%
  filter(stateroute %in% all_years$stateroute) %>%
  mutate(yearbin1 = case_when(year >= 1992 & year <= 1995 ~ 1992,
                              year >= 1996 & year <= 1999 ~ 1996,
                              year >= 2000 & year <= 2003 ~ 2000,
                              year >= 2004 & year <= 2007 ~ 2004,
                              year >= 2008 & year <= 2011 ~ 2008,
                              year >= 2012 & year <= 2015 ~ 2012),
  yearbin2 = case_when(year >= 1991 & year <= 1994 ~ 1992,
                       year >= 1995 & year <= 1998 ~ 1996,
                       year >= 1999 & year <= 2002 ~ 2000,
                       year >= 2003 & year <= 2006 ~ 2004,
                       year >= 2007 & year <= 2010 ~ 2008,
                       year >= 2011 & year <= 2014 ~ 2012),
  yearbin3 = case_when(year >= 1990 & year <= 1993 ~ 1992,
                       year >= 1994 & year <= 1997 ~ 1996,
                       year >= 1999 & year <= 2001 ~ 2000,
                       year >= 2002 & year <= 2005 ~ 2004,
                       year >= 2006 & year <= 2009 ~ 2008,
                       year >= 2010 & year <= 2013 ~ 2012),
  yearbin4 = case_when(year >= 1993 & year <= 1996 ~ 1992,
                       year >= 1997 & year <= 2000 ~ 1996,
                       year >= 2001 & year <= 2004 ~ 2000,
                       year >= 2005 & year <= 2008 ~ 2004,
                       year >= 2009 & year <= 2012 ~ 2008,
                       year >= 2013 & year <= 2016 ~ 2012))

dir_yearbin1 <- yearbins_test %>%
  filter(!is.na(yearbin1)) %>%
  group_by(stateroute, aou, yearbin1) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, yearbin1, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_yb1 = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$yearbin1)
  }))

dir_yearbin2 <- yearbins_test %>%
  filter(!is.na(yearbin2)) %>%
  group_by(stateroute, aou, yearbin2) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, yearbin2, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_yb2 = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$yearbin2)
  }))

dir_yearbin3 <- yearbins_test %>%
  filter(!is.na(yearbin3)) %>%
  group_by(stateroute, aou, yearbin3) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, yearbin3, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_yb3 = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$yearbin3)
  }))

dir_yearbin4 <- yearbins_test %>%
  filter(!is.na(yearbin4)) %>%
  group_by(stateroute, aou, yearbin4) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, yearbin4, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir_yb4 = map_dbl(data, ~{
    df <- .
    dist <- dist(df[, -1])
    trajectoryDirectionality(dist, sites = rep(1, nrow(df)), surveys = df$yearbin4)
  }))

compare_yearbins <- dir_yearbin1 %>%
  left_join(select(dir_yearbin2, -data)) %>%
  left_join(select(dir_yearbin3, -data)) %>%
  left_join(select(dir_yearbin4, -data)) %>%
  select(-data) %>%
  pivot_longer(2:5, names_to = "yearbin", values_to = "dir") %>%
  group_by(stateroute) %>%
  summarize(dir_mean = mean(dir, na.rm = T),
            dir_min = min(dir, na.rm = T), 
            dir_max = max(dir, na.rm = T))

ggplot(compare_yearbins, aes(x = fct_reorder(as.factor(stateroute), dir_mean), y = dir_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = dir_min, ymax = dir_max)) +
  theme(axis.text.x = element_blank()) +
  labs(x = "BBS routes", y = "Mean directionality")
ggsave("figures/directionality_range_by_yearbin.pdf", units = "in", width = 12, height = 8)

### How does species richness impact directionality value?

dir_spp_rich <- dir_core_spp %>%
  mutate(spp_rich = map_dbl(data, ~{
    df <- .
    
    df_long <- pivot_longer(df, 2:311, names_to = "aou", values_to = "logabund") %>%
      filter(logabund > 0)
    
    length(unique(df_long$aou))
    
  }))

dir_all_spp_rich <- dir_all_spp %>%
  mutate(spp_rich = map_dbl(data, ~{
    df <- .
    
    df_long <- pivot_longer(df, 2:311, names_to = "aou", values_to = "logabund") %>%
      filter(logabund > 0)
    
    length(unique(df_long$aou))
    
  }))

excl_trans <- ggplot(dir_spp_rich, aes(x = spp_rich, y = dir_core)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Species richness", y = "Directionality (excl. transients)")

all_spp <- ggplot(dir_all_spp_rich, aes(x = spp_rich, y = dir_all)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Species richness", y = "Directionality (all spp.)")

plot_grid(excl_trans, all_spp, nrow = 1)
ggsave("figures/directionality_vs_spprich.pdf", units = "in", height = 6, width = 12)

### Pop trends at top and bottom 3% of routes by directionality

top3 <- quantile(dir_core_spp$dir_core, 0.97)
bottom3 <- quantile(dir_core_spp$dir_core, 0.03)

min_max_dir <- dir_core_spp %>%
  mutate(pctl = case_when(dir_core >= top3 ~ "Top 3%",
                          dir_core <= bottom3 ~ "Bottom 3%",
                          TRUE ~ "middle")) %>%
  filter(pctl != "middle")

abund_trends <- read.csv("/Users/gracedicecco/Desktop/git/poptrends_envchange/model/BBS_abundance_trends.csv", stringsAsFactors = F)

min_max_dir_abund <- min_max_dir %>%
  left_join(abund_trends) %>%
  rename(Directionality = "pctl")

min_max_abundtrend_means <- min_max_dir_abund %>%
  group_by(pctl) %>%
  summarize(meantrend = mean(abundTrend, na.rm = T))

ggplot(min_max_dir_abund, aes(x = abundTrend, group = stateroute, col = Directionality, fill = Directionality)) + geom_density(alpha = 0.2) +
  geom_vline(col = "#F8766D", xintercept = min_max_abundtrend_means$meantrend[min_max_abundtrend_means$Directionality == "bottom 3%"], cex = 2, lty = 2) +
  geom_vline(col = "#00BFC4", xintercept = min_max_abundtrend_means$meantrend[min_max_abundtrend_means$Directionality == "top 3%"], cex = 2, lty = 2) +
  labs(x = "Abundance trend", y = "Count")
ggsave("figures/abundtrend_by_directionality.pdf")

### Compare 5 time points to all years directionality for one route

dir_compare_timepts <- dir_all_spp %>%
  dplyr::select(-data) %>%
  left_join(logabund_wide) %>%
  na.omit() %>%
  mutate(n_years = map_dbl(data, ~{
    df <- .
    length(unique(df$year[df$year >= 1990]))
  }))

ggplot(dir_compare_timepts, aes(x = n_years, y = dir25)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Years sampled", y = "Directionality")
ggsave("figures/directionality_by_years_sampled.pdf")

ggplot(dir_compare_timepts, aes(x = dir_all, y = dir25)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Directionality (4 year time windows)", y = "Directionality (all years)")
ggsave("figures/directionality_25yr_vs_5timepts.pdf")

### Compare trajectoryDir to manually calculated trajectory distances from PCoA

directionality_manual <- function(points) {
  if(nrow(points) == 6) {
    totaldist <- sqrt((points[6,1] - points[1,1])^2 + 
                        (points[6,2] - points[1,2])^2 +
                        (points[6,3] - points[1,3])^2 +
                        (points[6,4] - points[1,4])^2 +
                        (points[6,5] - points[1,5])^2)
    
    cumdist <- sum(sqrt((points[2,1] - points[1,1])^2 + 
                          (points[2,2] - points[1,2])^2 +
                          (points[2,3] - points[1,3])^2 +
                          (points[2,4] - points[1,4])^2 +
                          (points[2,5] - points[1,5])^2) +
                     sqrt((points[3,1] - points[2,1])^2 + 
                            (points[3,2] - points[2,2])^2 +
                            (points[3,3] - points[2,3])^2 +
                            (points[3,4] - points[2,4])^2 +
                            (points[3,5] - points[2,5])^2) +
                     sqrt((points[4,1] - points[3,1])^2 + 
                            (points[4,2] - points[3,2])^2 +
                            (points[4,3] - points[3,3])^2 +
                            (points[4,4] - points[3,4])^2 +
                            (points[4,5] - points[3,5])^2) +
                     sqrt((points[5,1] - points[4,1])^2 + 
                            (points[5,2] - points[4,2])^2 +
                            (points[5,3] - points[4,3])^2 +
                            (points[5,4] - points[4,4])^2 +
                            (points[5,5] - points[4,5])^2) +
                     sqrt((points[6,1] - points[5,1])^2 + 
                            (points[6,2] - points[5,2])^2 +
                            (points[6,3] - points[5,3])^2 +
                            (points[6,4] - points[5,4])^2 +
                            (points[6,5] - points[5,5])^2))
    dir <- totaldist/cumdist
    return(dir) }
  else {
    totaldist <- sqrt((points[5,1] - points[1,1])^2 + 
                        (points[5,2] - points[1,2])^2 +
                        (points[5,3] - points[1,3])^2 +
                        (points[5,4] - points[1,4])^2)
    
    cumdist <- sum(sqrt((points[2,1] - points[1,1])^2 + 
                          (points[2,2] - points[1,2])^2 +
                          (points[2,3] - points[1,3])^2 +
                          (points[2,4] - points[1,4])^2) +
                     sqrt((points[3,1] - points[2,1])^2 + 
                            (points[3,2] - points[2,2])^2 +
                            (points[3,3] - points[2,3])^2 +
                            (points[3,4] - points[2,4])^2) +
                     sqrt((points[4,1] - points[3,1])^2 + 
                            (points[4,2] - points[3,2])^2 +
                            (points[4,3] - points[3,3])^2 +
                            (points[4,4] - points[3,4])^2) +
                     sqrt((points[5,1] - points[4,1])^2 + 
                            (points[5,2] - points[4,2])^2 +
                            (points[5,3] - points[4,3])^2 +
                            (points[5,4] - points[4,4])^2))
    dir <- totaldist/cumdist
    return(dir) }
  }

dir_manual_all_spp <- log_abund_wider %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(pcoa_all = map(data, ~{
    df <- .
    dist <- dist(df[, -1])
    pcoa <- trajectoryPCoA(dist, sites = rep(1, nrow(df)), surveys = df$year_bin)
  }),
  dir_manual = map_dbl(pcoa_all, ~{
    pcoa <- .
    points <- pcoa$points
    directionality_manual(points)
  }))

dir_compare_manual <- dir_manual_all_spp %>%
  dplyr::select(-data) %>%
  left_join(dir_all_spp)

ggplot(dir_compare_manual, aes(x = dir_manual, y = dir_all)) + 
  geom_point() + 
  labs(x = "Distance between first and last points/Cumulative distance", y = "Directionality (all species)")
ggsave("figures/manual_directionality_vs_trajectory.pdf")
# cor(dir_compare_manual$dir_manual, dir_compare_manual$dir_all) = -0.37

ggplot(dir_compare_manual, aes(x = dir_manual)) + 
  geom_histogram(col = "white") + 
  labs(x = "Distance between first and last points/Cumulative distance", y = "Count")
ggsave("figures/manual_directionality_histogram.pdf")

### Pop trends at top and bottom 3% of routes for manual directionality

top3_manual <- quantile(dir_manual_all_spp$dir_manual, 0.97)
bottom3_manual <- quantile(dir_manual_all_spp$dir_manual, 0.03)

manual_dir_min_max <- dir_manual_all_spp %>%
  mutate(pctl = case_when(dir_manual >= top3_manual ~ "Top 3%",
                          dir_manual <= bottom3_manual ~ "Bottom 3%",
                          TRUE ~ "middle")) %>%
  filter(pctl != "middle")

min_max_dir_manual_abund <- manual_dir_min_max %>%
  left_join(abund_trends) %>%
  rename(Directionality = "pctl")

min_max_manual_means <- min_max_dir_manual_abund %>%
  group_by(Directionality) %>%
  summarize(meantrend = mean(abundTrend, na.rm = T))

ggplot(min_max_dir_manual_abund, aes(x = abundTrend, group = stateroute, col = Directionality, fill = Directionality)) + geom_density(alpha = 0.2) +
  geom_vline(col = "#F8766D", xintercept = min_max_manual_means$meantrend[min_max_manual_means$Directionality == "Bottom 3%"], cex = 2, lty = 2) +
  geom_vline(col = "#00BFC4", xintercept = min_max_manual_means$meantrend[min_max_manual_means$Directionality == "Top 3%"], cex = 2, lty = 2) +
  labs(x = "Abundance trend", y = "Count")
ggsave("figures/abundtrend_by_manual_directionality.pdf")

  
## Determine trajectoryDirectionality for time series sensitivity
## Subsample time points: 5 years up to 15-20

dir_sample_sens <- log_abund %>%
  select(-countrynum) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(samplesize = map(data, ~{
    df <- .
    res <- data.frame(n = c(), dir = c())
    for(i in 5:nrow(df)) {
      df_sample <- df %>%
        sample_n(i)
      
      abund_dist <- dist(df_sample[, -1])
      dir <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df_sample)), surveys = df_sample$year)
      
      res <- rbind(res, data.frame(n = i, dir = dir))
      
    }
    
    res
  }))

dir_sample_plot <- dir_sample_sens %>%
  select(-data) %>%
  unnest(cols = c(samplesize))

ggplot(dir_sample_plot, aes(x = n, y = dir, group = stateroute, col = stateroute)) + geom_smooth(se = F, alpha = 0.1) +
  labs(x = "Years", y = "Directionality")
ggsave("figures/directionality_sample_size_sensitivity.pdf")

# All species
log_abund_rsample <- bbs_subset %>%
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
  nest() %>%
  mutate(log_abund = map(data, ~{
    df <- .
    sample <- sample_n(df, 1)
    log10(sample$speciestotal + 1)
  }))

log_abund_rsample_wide <- log_abund_rsample %>%
  dplyr::select(stateroute, aou, year_bin, log_abund) %>%
  unnest(cols = c(log_abund)) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

dir_rsample_wide <- log_abund_rsample_wide %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir = map_dbl(data, ~{
    df <- .
      
    abund_dist <- dist(df[, -1])
    dir <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year_bin)

  }))

dir_rsample_plot <- dir_rsample_wide %>%
  select(-data) %>%
  left_join(dir_all_spp)

ggplot(dir_rsample_plot, aes(x = dir, y = dir_all)) + geom_point() + 
  geom_abline(intercept = 0, slope= 1) +
  labs(x = "Directionality (random sample count)", y = "Directionality (average count)")
ggsave("figures/directionality_averaging_sensitivity.pdf")

## How does taking time window averages impact abundances
## Directionality with presence-absence for comparison?

pres_wide <- log_abund_wider %>%
  mutate_all(~ifelse(. > 0 & . < 1000, 1, .))

pres_dir <- pres_wide %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir = map_dbl(data, ~{
    df <- .
    
    abund_dist <- dist(df[, -1])
    dir <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(df)), surveys = df$year_bin)
    
  }))

abund_pres_dir <- pres_dir %>%
  select(-data) %>%
  left_join(dir_all_spp)

ggplot(abund_pres_dir, aes(x = dir, y = dir_all)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Directionality (presence-absence)", y = "Directionality (average count)")
ggsave("figures/directionality_abund_vs_presabs.pdf")

## At each scale (1 route, up to nearest 25 routes within BCR) 
## Get max land cover delta from raw land cover data and get trend in Tmin and trend in Tmax
## Calculate directionality for each grouping of routes: average log(abundance) across routes when pooled
## Fit models across BBS routes for each scale (1:25)
## Variance partitioning of land cover and climate variables explaining trajectory 
## Directionality ~ tmin + tmax + max(deltaLandCover)
## Variance partitioning with ecospat.varpat(model.1, model.2, model.12) in ecospat

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
    
    # land_cover <- landcover_na %>%
    #   filter(stateroute %in% df$stateroute)
    #   
    # max_delta <- possibly_max_delta(land_cover)
    # 
    # climate <- bbs_climate_avgs %>%
    #   filter(stateroute %in% df$stateroute)
    # 
    # trend_tmax <- coef(lm(mean_tmax ~ year, data = climate))[[2]]
    # trend_tmin <- coef(lm(mean_tmin ~ year, data = climate))[[2]]
    
    log_abund <- log_abund_wider %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin) %>%
      summarize_all(mean, na.rm = T) %>%
      dplyr::select(-stateroute)
    abund_dist <- dist(log_abund[, -1])
    dir_all <- trajectoryDirectionality(abund_dist, sites = rep(1, nrow(log_abund)), surveys = log_abund$year_bin)
    
    log_abund_core <- log_abund_core %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin) %>%
      summarize_all(mean, na.rm = T) %>%
      dplyr::select(-stateroute)
    abund_dist_core <- dist(log_abund_core[, -1])
    dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_abund_core)), surveys = log_abund_core$year_bin)
    
    # data.frame(focal_rte = unique(df$focal_rte), max_lc = max_delta, 
    #            trend_tmax = trend_tmax, trend_tmin = trend_tmin, 
    #            dir_all = dir_all, dir_core = dir_core)
    
    data.frame(focal_rte = unique(df$focal_rte), dir_all = dir_all, dir_core = dir_core)
  
  }))

scale_model_variables <- scale_model_input %>%
  dplyr::select(scale, input_vars) %>%
  unnest(cols = c(input_vars)) %>%
  group_by(scale) %>%
  nest()

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
