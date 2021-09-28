### Community trajectory analysis
## 1990-2016
## Scale model

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)
library(ecospat)
library(cowplot)
library(broom)
library(grid)

#### Set up ####

## Plotting theme

theme_set(theme_classic(base_size = 15))

# Append correct BioArk path
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\ad.unc.edu\\bio")

## North America map

na <- read_sf("data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN") %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

## Read in data

# Species list 

species_list <- read.csv("data/species_list.csv", stringsAsFactors = F)
fourletter_codes <- read.csv("data/four_letter_codes_aous.csv", stringsAsFactors = F)

# Species traits/guilds

bird_traits <- read.csv("data/Master_RO_Correlates_20110610.csv", stringsAsFactors = F) %>%
  select(AOU, CommonName, Foraging, Trophic.Group, migclass)

habitat_guilds <- read.csv("data/aaw1313_Data_S1.csv", stringsAsFactors = F)

# BBS sampled every 5 years from 1970 to 2016
log_abund <- read.csv("data/derived_data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)

# BBS sampled 3/4 years every 4 year time window
bbs_subset <- read.csv("data/bbs_counts_subset_1990-2016.csv", stringsAsFactors = F)

routes <- read.csv(paste0(bioark, "/hurlbertlab/databases/BBS/2017/bbs_routes_20170712.csv"), stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# Land cover and climate data - 1 route scale

bbs_landcover <- read.csv("data/derived_data/bbs_route_max_landcover_change.csv", stringsAsFactors = F)
bbs_climate <- read.csv("data/derived_data/bbs_routes_climate_trends.csv", stringsAsFactors = F)

## Raw land cover and climate data for calculating multiple scales

landcover_us <- read.csv(paste0(bioark, "/hurlbertlab/dicecco/data/fragmentation_indices_nlcd_simplified.csv"), stringsAsFactors = F) %>%
  mutate(country = "US")
landcover_ca <- read.csv(paste0(bioark, "/hurlbertlab/dicecco/data/fragmentation_indices_canada.csv"), stringsAsFactors = F) %>%
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

bbs_climate_avgs <- read.csv("data/derived_data/bbs_routes_breeding_season_climate.csv", stringsAsFactors = F) %>%
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
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(1:25)))) %>%
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

#### BCR map and aggregation suppl fig ####

## BCR map

bcr <- read_sf(paste0(bioark,"/HurlbertLab/DiCecco/bcr_terrestrial_shape/BCR_Terrestrial_master.shp")) %>%
  filter(BCR %in% bcr_subset$bcr) %>%
  mutate_at(c("BCR"), ~as.factor(.)) %>%
  st_transform(st_crs(na)) %>%
  st_crop(st_bbox(na))

study_routes <- bbs_subset %>%
  ungroup() %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  distinct(stateroute, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

bcr_map <- tm_shape(na) + tm_polygons(col = "gray50") +
  tm_shape(bcr) + tm_polygons(col = "BCR") +
  tm_shape(study_routes) + tm_dots(col = "black", size = 0.05) + 
  tm_layout(legend.text.size = 1.25, legend.title.size = 1.5, outer.margins = c(0.01,0,0.01,0),
            inner.margins = c(0.0, 0.08, 0.0, 0.0), legend.position = c("left", "bottom"),
            main.title = "A", title.size = 4) +
  tm_compass(type = "8star", position = c("right", "top")) +
  tm_scale_bar(breaks = c(0, 500, 1000), text.size = 1.5)

bcr_22 <- bcr %>%
  filter(BCR == 22)

na_22 <- na %>%
  st_transform(st_crs(bcr_22)) %>%
  st_crop(st_bbox(bcr_22))

rtes_22 <- study_routes %>%
  st_set_crs(st_crs(bcr)) %>%
  st_intersection(bcr_22) 

focal_rte <- study_routes %>%
  filter(stateroute == 34020)

rtes_focal <- mean_bbs_distances %>%
  filter(stateroute == 34020, scale == 25)

nearest_rte <- study_routes %>%
  filter(stateroute %in% rtes_focal$model_input[[1]]$stateroute) %>%
  filter(stateroute != 34020)

circle_buffer <- st_buffer(focal_rte, dist = 2.5) %>%
  st_set_crs(st_crs(focal_rte))

agg_panel <- tm_shape(na_22) + tm_polygons(col = "gray50") +
  tm_shape(bcr_22) + tm_polygons(col = "BCR", legend.show = F) +
  tm_shape(rtes_22) + tm_dots(col = "black", size = 0.1) + 
  tm_shape(focal_rte) + tm_symbols(col = "blue", shape = 15, size = 0.2) +
  tm_shape(nearest_rte) + tm_symbols(col = "purple", shape = 17, size = 0.2) +
  tm_shape(circle_buffer) + tm_borders(col = "black") +
  tm_layout(legend.text.size = 1.25, legend.title.size = 1.5, outer.margins = c(0.01,0,0.01,0),
            inner.margins = c(0.02, 0.02, 0.02, 0.02), legend.position = c("left", "bottom"),
            main.title = "B", title.size = 4) +
  tm_compass(type = "8star", position = c("right", "bottom")) +
  tm_scale_bar(breaks = c(0, 100, 200), text.size = 0.9)

bcr_panels <- tmap_arrange(bcr_map, agg_panel, nrow = 2)
tmap_save(bcr_panels, "figures/bcr_aggregation.pdf", units = "in", height = 10, width = 8)

bcr_palette <- tmaptools::get_brewer_pal("Set3", n = 15)
names(bcr_palette) <- as.factor(unique(bcr$BCR))

mean_dist_plot <- ggplot(filter(mean_bcr_distances, bcr %in% bcr_subset$bcr), 
                         aes(x = scale, y = mean_dist_bcr, col = as.factor(bcr), group = as.factor(bcr))) + 
  geom_line(cex = 1) +
  scale_color_manual(values = bcr_palette) +
  labs(x = "Scale (routes)", y = "Mean distance between routes (km)", col = "BCR") 

legend <- get_legend(mean_dist_plot)

pct_plot <- ggplot(mean_pct_overlap, aes(x = scale, y = mean_reps*100, col = as.factor(bcr), group = as.factor(bcr))) + 
  geom_line(cex = 1) +
  scale_color_manual(values = bcr_palette) +
  labs(x = "Scale (routes)", y = "Avg percent of aggregates route occurs in", col = "BCR") +
  theme(legend.position = "none")

plot_grid(mean_dist_plot+theme(legend.position = "none"), pct_plot, legend,
                         labels = c("A", "B", ""), rel_widths = c(0.45, 0.45, 0.1),
                         ncol = 3)
ggsave(paste0(getwd(),"/figures/bcr_aggregation_multipanel.pdf"), units = "in", height = 5, width = 10)

## Log abund from 1 route up to nearest 25 routes
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
  dplyr::select(-data) %>%
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
  dplyr::select(-data) %>%
  unnest(cols = c(year_bins)) %>%
  group_by(stateroute, aou, year_bin) %>%
  summarize(n_years = n_distinct(year),
            mean_abund = mean(speciestotal) + 1,
            log_abund = log10(mean_abund)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, year_bin, log_abund) %>%
  pivot_wider(names_from = aou, values_from = log_abund, values_fn = list(log_abund = mean), values_fill = list(log_abund = 0))

#### Scale model ####

## At each scale (0.5 route, 1 route, up to nearest 25 routes within BCR) 
## Get max land cover delta from raw land cover data and get trend in Tmin and trend in Tmax
## Calculate directionality for each grouping of routes: average log(abundance) across routes when pooled
## Fit models across BBS routes for each scale (1:25)
## Variance partitioning of land cover and climate variables explaining trajectory 
## Directionality ~ tmin + tmax + max(deltaLandCover)
## Variance partitioning 

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

scale_model_variables_unnest <- scale_model_variables %>%
  unnest(cols = c(data))
# write.csv(scale_model_variables_unnest, "data/derived_data/scale_model_input.csv", row.names = F)

scale_model_variables_unnest <- read.csv("data/derived_data/scale_model_input.csv", stringsAsFactors = F)

scale_model_variables <- read.csv("data/derived_data/scale_model_input.csv", stringsAsFactors = F) %>%
  group_by(scale) %>%
  nest()

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

# write.csv(scale_model_output, "data/derived_data/scale_model_output_variance.csv", row.names = F)
scale_model_output <- read.csv("data/derived_data/scale_model_output_variance.csv")

# Scale model figures 

scale_model_plot <- scale_model_output %>%
  pivot_longer(part1:unexpl, names_to = "variance") %>%
  mutate_at(c("value"), .funs = ~ifelse(. < 0, 0, .))

all_routes <- ggplot(filter(scale_model_plot, variance != "unexpl"), aes(x = scale, y = value, fill = variance)) +
  geom_col(position = "stack", col = "white") +
  scale_fill_manual(values = c("gray", "#92C5DE", "#A6D854"), name = "Variance explained", labels = c("Shared", "Climate", "Land cover")) +
  labs(x = "Spatial scale (no. of routes)", y = "Variance")

# Scale model predictor effects

scale_model_variables <-  read.csv("data/derived_data/scale_model_input.csv", stringsAsFactors = F) %>%
  group_by(scale) %>%
  nest()

scale_model_ests <- data.frame(scale = c(), term = c(), estimate = c(), std.error = c(), 
                                 statistic = c(), p.value = c(), conf_lo = c(), conf_hi = c())
for(i in 1:25) {
  model_input <- scale_model_variables$data[[i]]
  
  mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  tidy_mod <- tidy(mod12) %>%
    mutate(scale = i,
           conf_lo = confint(mod12)[, 1],
           conf_hi = confint(mod12)[, 2])
  
  scale_model_ests <- rbind(scale_model_ests, 
                              tidy_mod)
}

# Compare effect estimates w/ beta regression

scale_model_output_beta <- data.frame(scale = c(), term = c(), estimate = c(), std.error = c(), 
                                      statistic = c(), p.value = c(), 
                                      component = c())
for(i in 1:25) {
  model_input <- scale_model_variables$data[[i]]
  
  mod12 <- betareg(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  tidy_mod <- tidy(mod12) %>%
    mutate(scale = i)
  
  scale_model_output_beta <- rbind(scale_model_output_beta, 
                              tidy_mod)
}

mod_compare <- scale_model_ests %>%
  left_join(scale_model_output_beta, by = c("scale", "term"), suffix = c("_lm", "_beta")) %>%
  filter(term != "(Intercept)")

ggplot(mod_compare, aes(x = estimate_lm, y = estimate_beta)) + geom_point() + geom_abline(slope = 1, intercept= 0) +
  facet_wrap(~term) + labs(x = "LM estimate", y = "Betareg estimate")
ggsave("figures/lm_betareg_ests.pdf", units = "in", height = 4, width = 10)

betareg_r2 <- data.frame(scale = c(), r2 = c())
for(i in 1:25) {
  model_input <- scale_model_variables$data[[i]]
  
  mod12 <- betareg(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  tidy_mod <- data.frame(r2 = summary(mod12)$pseudo.r.squared) %>%
    mutate(scale = i)
  
  betareg_r2 <- rbind(betareg_r2, 
                                   tidy_mod)
}

r2_compare <- scale_model_output %>%
  left_join(betareg_r2) %>%
  mutate(expl = 1 - unexpl)

ggplot(r2_compare, aes(x = expl, y = r2)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  labs(x = "Linear model total R2", y = "Betareg total R2")
ggsave("figures/lm_betareg_r2.pdf")

# Compare distance between routes with directionality

dist_dir_cor <- scale_model_variables_unnest %>%
  filter(scale == 25) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  left_join(mean_bcr_distances, by = c("bcr", "scale"))

ggplot(dist_dir_cor, aes(x = mean_dist_bcr, y = dir_core)) + geom_point() + 
  labs(x = "Mean distance between routes", y = "Turnover")
ggsave("figures/distance_btw_routes_dir.pdf")

mod <- betareg(dir_core ~ mean_dist_bcr, data = dist_dir_cor)
summary(mod)

#### Scale model with low/no overlap ####

scale_25 <- mean_bbs_distances %>%
  filter(scale == 25) %>%
  select(-stateroute) %>%
  unnest(model_input) %>%
  group_by(focal_rte, bcr, scale) %>%
  nest()

# What is the maximum number of focal routes per BCR assuming 25 aggregated routes, no more than 40% overlap (10 routes)?
max_focal_rtes <- bcr_subset %>%
  st_set_geometry(NULL) %>%
  mutate(max_focal_rtes = floor(1 + (n_routes - 25)/15))
# 33 data points total

large_bcrs <- max_focal_rtes %>%
  filter(max_focal_rtes > 1)

# For BCRs with more than one focal route, sample x routes and calculate pairwise overlap of aggregated 25 nearest routes 100x
# Select x focal routes with min pairwise overlap

two_rtes <- large_bcrs %>%
  filter(max_focal_rtes == 2)

min_overlap_rtes <- data.frame(bcr = rep(NA, 7), focal_rte1 = rep(NA, 7), focal_rte2 = rep(NA, 7), overlap = rep(NA, 7))

for(i in 1:7) {
  b <- two_rtes[[i, 1]]
  
  rtes <- filter(scale_25, bcr == b)
  
  res <- data.frame(stateroute1 = rep(NA, 100), stateroute2 = rep(NA, 100), overlap = rep(NA, 100))
  
  for(j in 1:100) {
    focal_rtes <- sample(rtes$focal_rte, 2)
    rte1 <- focal_rtes[[1]]
    rte2 <- focal_rtes[[2]]
    
    rte_list1 <- scale_25 %>%
      filter(focal_rte == rte1) %>%
      unnest(data)
    
    rte_list2 <- scale_25 %>%
      filter(focal_rte == rte2) %>%
      unnest(data)
    
    ovlp <- sum(rte_list1$stateroute %in% rte_list2$stateroute)/25
    
    res[j, 1] <- rte1
    res[j, 2] <- rte2
    res[j, 3] <- ovlp
  }
  
  best_pair <- filter(res, overlap == min(overlap))
  
  min_overlap_rtes[i, 1] <- b
  min_overlap_rtes[i, 2] <- best_pair$stateroute1[[1]]
  min_overlap_rtes[i, 3] <- best_pair$stateroute2[[1]]
  min_overlap_rtes[i, 4] <- best_pair$overlap[[1]]
  
}


three_rtes <- large_bcrs %>%
  filter(max_focal_rtes == 3)

min_overlap_3rtes <- data.frame(bcr = rep(NA, 2), 
                               focal_rte1 = rep(NA, 2), focal_rte2 = rep(NA, 2), focal_rte3 = rep(NA, 2),
                               overlap = rep(NA, 2))

for(i in 1:2) {
  b <- three_rtes[[i, 1]]
  
  rtes <- filter(scale_25, bcr == b)
  
  res <- data.frame(stateroute1 = rep(NA, 1000), stateroute2 = rep(NA, 1000), stateroute3 = rep(NA, 1000), 
                    max.overlap = rep(NA, 1000), mean.overlap = rep(NA, 1000))
  
  for(j in 1:1000) {
    focal_rtes <- sample(rtes$focal_rte, 3)
    rte1 <- focal_rtes[[1]]
    rte2 <- focal_rtes[[2]]
    rte3 <- focal_rtes[[3]]
    
    rte_list1 <- scale_25 %>%
      filter(focal_rte == rte1) %>%
      unnest(data)
    
    rte_list2 <- scale_25 %>%
      filter(focal_rte == rte2) %>%
      unnest(data)
    
    rte_list3 <- scale_25 %>%
      filter(focal_rte == rte3) %>%
      unnest(data)
    
    max.ovlp <- max(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                       sum(rte_list1$stateroute %in% rte_list3$stateroute), 
                       sum(rte_list2$stateroute %in% rte_list3$stateroute)))
    mean.ovlp <- mean(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                        sum(rte_list1$stateroute %in% rte_list3$stateroute), 
                        sum(rte_list2$stateroute %in% rte_list3$stateroute)))/25
    
    res[j, 1] <- rte1
    res[j, 2] <- rte2
    res[j, 3] <- rte3
    res[j, 4] <- max.ovlp
    res[j, 5] <- mean.ovlp
  }
  
  best_pair <- filter(res, max.overlap <= 10, mean.overlap == min(mean.overlap))
  
  min_overlap_3rtes[i, 1] <- b
  min_overlap_3rtes[i, 2] <- best_pair$stateroute1[[1]]
  min_overlap_3rtes[i, 3] <- best_pair$stateroute2[[1]]
  min_overlap_3rtes[i, 4] <- best_pair$stateroute3[[1]]
  min_overlap_3rtes[i, 5] <- best_pair$mean.overlap[[1]]
  
}

# Four routes - BCR 22

rtes <- filter(scale_25, bcr == 22)

res <- data.frame(stateroute1 = rep(NA, 1000), stateroute2 = rep(NA, 1000), stateroute3 = rep(NA, 1000), stateroute4 = rep(NA, 1000), 
                  max.overlap = rep(NA, 1000), mean.overlap = rep(NA, 1000))

for(j in 1:1000) {
  focal_rtes <- sample(rtes$focal_rte, 4)
  rte1 <- focal_rtes[[1]]
  rte2 <- focal_rtes[[2]]
  rte3 <- focal_rtes[[3]]
  rte4 <- focal_rtes[[4]]
  
  rte_list1 <- scale_25 %>%
    filter(focal_rte == rte1) %>%
    unnest(data)
  
  rte_list2 <- scale_25 %>%
    filter(focal_rte == rte2) %>%
    unnest(data)
  
  rte_list3 <- scale_25 %>%
    filter(focal_rte == rte3) %>%
    unnest(data)
  
  rte_list4 <- scale_25 %>%
    filter(focal_rte == rte4) %>%
    unnest(data)
  
  max.ovlp <- max(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                    sum(rte_list1$stateroute %in% rte_list3$stateroute),
                    sum(rte_list1$stateroute %in% rte_list4$stateroute),
                    sum(rte_list2$stateroute %in% rte_list3$stateroute),
                    sum(rte_list2$stateroute %in% rte_list4$stateroute),
                    sum(rte_list3$stateroute %in% rte_list4$stateroute)))
  mean.ovlp <- mean(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                      sum(rte_list1$stateroute %in% rte_list3$stateroute),
                      sum(rte_list1$stateroute %in% rte_list4$stateroute),
                      sum(rte_list2$stateroute %in% rte_list3$stateroute),
                      sum(rte_list2$stateroute %in% rte_list4$stateroute),
                      sum(rte_list3$stateroute %in% rte_list4$stateroute)))/25
  
  res[j, 1] <- rte1
  res[j, 2] <- rte2
  res[j, 3] <- rte3
  res[j, 4] <- rte4
  res[j, 5] <- max.ovlp
  res[j, 6] <- mean.ovlp
}

best_pair_4 <- filter(res, max.overlap <= 10, mean.overlap == min(mean.overlap)) %>%
  slice(1:1)

# Five routes - BCR 28

rtes <- filter(scale_25, bcr == 28)

res <- data.frame(stateroute1 = rep(NA, 10000), stateroute2 = rep(NA, 10000), 
                  stateroute3 = rep(NA, 10000), stateroute4 = rep(NA, 10000), stateroute5 = rep(NA, 10000),
                  max.overlap = rep(NA, 10000), mean.overlap = rep(NA, 10000))

for(j in 1:10000) {
  focal_rtes <- sample(rtes$focal_rte, 5)
  rte1 <- focal_rtes[[1]]
  rte2 <- focal_rtes[[2]]
  rte3 <- focal_rtes[[3]]
  rte4 <- focal_rtes[[4]]
  rte5 <- focal_rtes[[5]]
  
  rte_list1 <- scale_25 %>%
    filter(focal_rte == rte1) %>%
    unnest(data)
  
  rte_list2 <- scale_25 %>%
    filter(focal_rte == rte2) %>%
    unnest(data)
  
  rte_list3 <- scale_25 %>%
    filter(focal_rte == rte3) %>%
    unnest(data)
  
  rte_list4 <- scale_25 %>%
    filter(focal_rte == rte4) %>%
    unnest(data)
  
  rte_list5 <- scale_25 %>%
    filter(focal_rte == rte5) %>%
    unnest(data)
  
  max.ovlp <- max(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                    sum(rte_list1$stateroute %in% rte_list3$stateroute),
                    sum(rte_list1$stateroute %in% rte_list4$stateroute),
                    sum(rte_list1$stateroute %in% rte_list5$stateroute),
                    sum(rte_list2$stateroute %in% rte_list3$stateroute),
                    sum(rte_list2$stateroute %in% rte_list4$stateroute),
                    sum(rte_list2$stateroute %in% rte_list5$stateroute),
                    sum(rte_list3$stateroute %in% rte_list4$stateroute),
                    sum(rte_list3$stateroute %in% rte_list5$stateroute),
                    sum(rte_list4$stateroute %in% rte_list5$stateroute)))
  mean.ovlp <- mean(c(sum(rte_list1$stateroute %in% rte_list2$stateroute), 
                      sum(rte_list1$stateroute %in% rte_list3$stateroute),
                      sum(rte_list1$stateroute %in% rte_list4$stateroute),
                      sum(rte_list1$stateroute %in% rte_list5$stateroute),
                      sum(rte_list2$stateroute %in% rte_list3$stateroute),
                      sum(rte_list2$stateroute %in% rte_list4$stateroute),
                      sum(rte_list2$stateroute %in% rte_list5$stateroute),
                      sum(rte_list3$stateroute %in% rte_list4$stateroute),
                      sum(rte_list3$stateroute %in% rte_list5$stateroute),
                      sum(rte_list4$stateroute %in% rte_list5$stateroute)))/25
  
  res[j, 1] <- rte1
  res[j, 2] <- rte2
  res[j, 3] <- rte3
  res[j, 4] <- rte4
  res[j, 5] <- rte5
  res[j, 6] <- max.ovlp
  res[j, 7] <- mean.ovlp
}

best_pair_5 <- filter(res, max.overlap <= 10, mean.overlap == min(mean.overlap)) %>%
  slice(1:1)

## List of focal routes with low overlap

one_focal_rte <- max_focal_rtes %>%
  filter(max_focal_rtes == 1)

focal_routes_1 <- scale_25 %>%
  filter(bcr %in% one_focal_rte$bcr) %>%
  group_by(bcr) %>%
  sample_n(1) %>%
  select(-scale, -data) %>%
  mutate(rte = "stateroute1")

best_routes_2 <- min_overlap_rtes %>%
  pivot_longer(2:3, names_to = "rte", values_to = "focal_rte")

best_routes_3 <- min_overlap_3rtes %>%
  pivot_longer(2:4, names_to = "rte", values_to = "focal_rte")

best_routes_4 <- best_pair_4 %>%
  mutate(bcr = 22) %>%
  select(-max.overlap) %>%
  rename(overlap = "mean.overlap") %>%
  pivot_longer(1:4, names_to = "rte", values_to = "focal_rte")

best_routes_5 <- best_pair_5  %>%
  mutate(bcr = 28) %>%
  select(-max.overlap) %>%
  rename(overlap = "mean.overlap") %>%
  pivot_longer(1:5, names_to = "rte", values_to = "focal_rte")

low_overlap_focal_routes <- bind_rows(focal_routes_1, best_routes_2, best_routes_3, best_routes_4, best_routes_5)
# write.csv(low_overlap_focal_routes, "data/derived_data/low_overlap_focal_routes.csv", row.names = F)

low_overlap_focal_routes <- read.csv("data/derived_data/low_overlap_focal_routes.csv")

## Run models for just focal routes 

scale_model_variables_unnest <- read.csv("data/derived_data/scale_model_input.csv", stringsAsFactors = F)

low_overlap_model_variables <- scale_model_variables_unnest %>%
  filter(focal_rte %in% low_overlap_focal_routes$focal_rte) %>%
  group_by(scale) %>%
  nest()

LO_scale_model_output <- data.frame(scale = c(), part1 = c(), part2 = c(), joined = c(), unexpl = c())
for(i in 1:25) {
  model_input <- low_overlap_model_variables$data[[i]]
  
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
  
  LO_scale_model_output <- rbind(LO_scale_model_output, 
                              data.frame(scale = i, part1 = part1, part2 = part2,
                                         joined = joined, unexpl = unexpl))
}
# write.csv(LO_scale_model_output, "data/derived_data/low_overlap_scale_model_output_deviance.csv", row.names = F)
LO_scale_model_output <- read.csv("data/derived_data/low_overlap_scale_model_output_deviance.csv", stringsAsFactors = F)

LO_scale_model_plot <- LO_scale_model_output %>%
  pivot_longer(part1:unexpl, names_to = "variance")  %>%
  mutate_at(c("value"), .funs = ~ifelse(. < 0, 0, .))

lo_overlap <- ggplot(filter(LO_scale_model_plot, variance != "unexpl"), aes(x = scale, y = value, fill = variance)) +
  geom_col(position = "stack", width = 0.9, col = "white") +
  scale_fill_manual(values = c("gray", "#92C5DE", "#A6D854"), name = "Variance explained", labels = c("Shared", "Climate", "Land cover")) +
  labs(x = "Spatial scale (no. of routes)", y = "Variance")

## Variance partitioning multipanel figure
legend <- get_legend(all_routes)

plot_grid(all_routes + theme(legend.position = "none"), lo_overlap + theme(legend.position = "none"), legend,
          ncol = 3, rel_widths = c(0.4, 0.4, 0.2), labels = c("A", "B", ""))
ggsave("figures/variance_partitioning_multipanel.pdf", units = "in", height = 5, width = 12)

# Scale model lo overlap predictor effects

low_overlap_model_variables <- scale_model_variables_unnest %>%
  filter(focal_rte %in% low_overlap_focal_routes$focal_rte) %>%
  group_by(scale) %>%
  nest()

scale_model_output <- data.frame(scale = c(), part1 = c(), part2 = c(), joined = c(), unexpl = c())
for(i in 1:25) {
  model_input <- low_overlap_model_variables$data[[i]]
  
  mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  tidy_mod <- tidy(mod12) %>%
    mutate(scale = i,
           conf_lo = confint(mod12)[, 1],
           conf_hi = confint(mod12)[, 2])
  
  scale_model_output <- rbind(scale_model_output, 
                              tidy_mod)
}

#### Map of directionality values ####

dir_sf <- scale_model_variables_unnest %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  filter(statenum != 3) %>%
  st_as_sf(coords = c("longitude", "latitude"))

one_route <- tm_shape(na) + tm_polygons(col = "gray50") +
  tm_shape(filter(dir_sf, scale == 1)) + 
  tm_dots(col = "dir_core", title = "Turnover", palette = "YlGnBu", size = 0.3, legend.show = F) +
  tm_compass(type = "8star", position = c("right", "top")) +
  tm_scale_bar(breaks = c(0, 500, 1000), text.size = 1.5, position = c("left", "top")) +
  tm_layout(main.title = "A", main.title.size = 1.5, inner.margins = c(0.15, 0.02, 0.02, 0.02),
            title = "Spatial scale:\n1 route", title.position = c(0.8, 0.1))
# col breaks: 0.25-0.30, ... 0.5-0.55 by 0.5
# 6 colors

# histogram legend
one_route_cols <- RColorBrewer::brewer.pal("YlGnBu", n = 6)
names(one_route_cols) <- seq(1:6)

one_route_vals <- dir_sf %>%
  filter(scale == 1) %>%
  mutate(color = case_when(dir_core > 0.25 & dir_core <= 0.30 ~ 1,
                           dir_core > 0.30 & dir_core <= 0.35 ~ 2,
                           dir_core > 0.35 & dir_core <= 0.40 ~ 3,
                           dir_core > 0.40 & dir_core <= 0.45 ~ 4,
                           dir_core > 0.45 & dir_core <= 0.50 ~ 5,
                           dir_core > 0.50 ~ 6))

one_legend <- ggplot(one_route_vals, aes(x = dir_core, fill = as.factor(color))) +
  geom_histogram(breaks = seq(0.25, 0.55, by = 0.01)) + 
  scale_fill_manual(values = one_route_cols) +
  labs(x = "Turnover", y = "Count") + theme(legend.position = "none", 
                                                  panel.background = element_rect(fill = "transparent"), # bg of the panel
                                                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot


tf_route <-  tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(dir_sf, scale == 25)) + 
  tm_dots(col = "dir_core", title = "Turnover", palette = "YlGnBu", size = 0.3, legend.show = F) +
  tm_layout(main.title = "B", main.title.size = 1.5, inner.margins = c(0.15, 0.02, 0.02, 0.02),
            title = "Spatial scale:\n25 routes", title.position = c(0.8, 0.1))
# col breaks: 0.25 to 0.65 by 0.5
# 8 colors

tf_route_cols <- RColorBrewer::brewer.pal("YlGnBu", n = 8)

names(tf_route_cols) <- seq(1:8)

tf_route_vals <- dir_sf %>%
  filter(scale == 25) %>%
  mutate(color = case_when(dir_core > 0.25 & dir_core <= 0.30 ~ 1,
                           dir_core > 0.30 & dir_core <= 0.35 ~ 2,
                           dir_core > 0.35 & dir_core <= 0.40 ~ 3,
                           dir_core > 0.40 & dir_core <= 0.45 ~ 4,
                           dir_core > 0.45 & dir_core <= 0.50 ~ 5,
                           dir_core > 0.50 & dir_core <= 0.55 ~ 6,
                           dir_core > 0.55 & dir_core <= 0.60 ~ 7,
                           dir_core > 0.6 ~ 8))

tf_legend <- ggplot(tf_route_vals, aes(x = dir_core, fill = as.factor(color))) +
  geom_histogram(breaks = seq(0.25, 0.65, by = 0.01)) + 
  scale_fill_manual(values = tf_route_cols) +
  labs(x = "Turnover", y = "Count") + theme(legend.position = "none",
                                                  panel.background = element_rect(fill = "transparent"), # bg of the panel
                                                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot

vp_one <- viewport(0.1, 0.25, width = 0.175, height = 0.25)
vp_tf <- viewport(0.6, 0.25, width = 0.175, height = 0.25)

wd <- getwd()
pdf(paste0(wd, "/figures/directionality_map.pdf"), height = 8, width = 16)
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
print(one_route, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tf_route, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(one_legend, vp = vp_one)
print(tf_legend, vp = vp_tf)
dev.off()

#### High leverage species ####

# Directionality values at 25 route scales
# Regional high leverage species

scale_model_variables_unnest <- read.csv("data/derived_data/scale_model_input.csv", stringsAsFactors = F)
low_overlap_focal_routes <- read.csv("data/derived_data/low_overlap_focal_routes.csv")

regional_rtes <- bbs_route_distances %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(grouped_rtes = purrr::map(data, ~{
    df <- .
    
      df %>%
      arrange(distance) %>%
      slice(1:25) %>%
      dplyr::select(stateroute)
  })) %>%
  dplyr::select(-data) %>%
  unnest(grouped_rtes)

# High leverage species with leave-one-out calculations of directionality (@25 route scales)

core_spp <- bbs_subset %>%
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
  group_by(stateroute) %>%
  distinct(aou) %>%
  filter(!is.na(aou))

spp_dir_deltas <- regional_rtes %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(spp_list = map(data, ~{
    df <- .
    core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou)
  }),
  n_spp = map_dbl(spp_list, ~nrow(.))) %>%
  filter(n_spp > 0) %>%
  mutate(spp_dir = map2(data, spp_list, ~{
    routes <- .x
    spp <- .y
    spp$excl_dir <- NA
    
    for(sp in spp$aou) {
      log_core <- log_abund_core %>%
        filter(stateroute %in% routes$stateroute) %>%
        select(-contains(as.character(sp))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)

      if(nrow(log_core) > 3) {
      abund_dist_core <- dist(log_core[, -1])
      dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
      
      spp$excl_dir[spp$aou == sp] <- dir_core
      }
      
      else {
        spp$excl_dir[spp$aou == sp] <- NA
      }
    }
    
    spp

  }))

spp_dir_unnest <- spp_dir_deltas %>%
  select(-data, -spp_list, -n_spp) %>%
  unnest(spp_dir)
# write.csv(spp_dir_unnest, "data/derived_data/species-leave-one-out-directionality.csv", row.names = F)

spp_dir_unnest <- read.csv("data/derived_data/species-leave-one-out-directionality.csv", stringsAsFactors = F)

regional_dir_core <- scale_model_variables_unnest %>%
  filter(scale == 25) %>%
  select(focal_rte, dir_core)

spp_dir_diffs <- spp_dir_unnest %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(dir_diff > 0) %>%
  group_by(focal_rte) %>%
  filter(dir_diff == max(dir_diff, na.rm = T)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))

high_impact_spp <- spp_dir_diffs %>% 
  group_by(aou) %>% 
  summarize(n_regions = n()) %>%
  arrange(desc(n_regions)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))
# write.csv(high_impact_spp, "data/derived_data/spec-LOO-directionality-hi-impact.csv", row.names = F)

# Local high leverage species - how well do these match with regional

local_rtes <- bbs_route_distances %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(grouped_rtes = purrr::map(data, ~{
    df <- .
    
    df %>%
      arrange(distance) %>%
      slice(1:3) %>%
      dplyr::select(stateroute)
  })) %>%
  dplyr::select(-data) %>%
  unnest(grouped_rtes)

spp_dir_deltas_local <- local_rtes %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(spp_list = map(data, ~{
    df <- .
    core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou)
  }),
  n_spp = map_dbl(spp_list, ~nrow(.))) %>%
  filter(n_spp > 0) %>%
  mutate(spp_dir = map2(data, spp_list, ~{
    routes <- .x
    spp <- .y
    spp$excl_dir <- NA
    
    for(sp in spp$aou) {
      log_core <- log_abund_core %>%
        filter(stateroute %in% routes$stateroute) %>%
        select(-contains(as.character(sp))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        spp$excl_dir[spp$aou == sp] <- dir_core
      }
      
      else {
        spp$excl_dir[spp$aou == sp] <- NA
      }
    }
    
    spp
    
  }))

spp_dir_unnest_local <- spp_dir_deltas_local %>%
  select(-data, -spp_list, -n_spp) %>%
  unnest(spp_dir)
# write.csv(spp_dir_unnest_local, "data/derived_data/spec-leave-one-out-directionality-local.csv", row.names = F)

spp_dir_unnest_local <- read.csv("data/derived_data/spec-leave-one-out-directionality-local.csv", stringsAsFactors = F)

local_dir_core <- scale_model_variables_unnest %>%
  filter(scale == 3) %>%
  select(focal_rte, dir_core)

spp_dir_diffs_local <- spp_dir_unnest_local %>%
  left_join(local_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(dir_diff > 0) %>%
  group_by(focal_rte) %>%
  filter(dir_diff == max(dir_diff, na.rm = T)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))

high_impact_spp_local <- spp_dir_diffs_local %>% 
  group_by(aou) %>% 
  summarize(n_regions = n()) %>%
  arrange(desc(n_regions)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))
# write.csv(high_impact_spp_local, "data/derived_data/spec-LOO-directionality-hi-impact-local.csv", row.names = F)

local_v_regional <- spp_dir_diffs_local %>%
  left_join(spp_dir_diffs, by = c("focal_rte"), suffix = c('_local', "_regional"))
ggplot(local_v_regional, aes(x = dir_diff_local, y = dir_diff_regional)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + labs(x = "Directionality impact local", y = "Directionality impact regional")
ggsave("figures/spp_impacts_local_regional.pdf")

sum(local_v_regional$SPEC_local == local_v_regional$SPEC_regional, na.rm = T)/nrow(local_v_regional)
# 30% match between local/regional scale

# make local & regional maps

# get list of combined top tens from local/regional

high_impact_spp <- read.csv("data/derived_data/spec-LOO-directionality-hi-impact.csv", stringsAsFactors = F)

top_ten <- high_impact_spp %>%
  slice(1:10)

top_ten_local <- high_impact_spp_local %>%
  slice(1:10)

top_spp <- top_ten %>%
  full_join(top_ten_local, 
            by = c("aou", "SPEC", "COMMONNAME", "SCINAME", "SPEC6", "species_id", "french_common_name", 
                   "spanish_common_name", "sporder", "family", "genus", "species"), 
            suffix = c("_regional", "_local")) %>%
  slice(1:11)

# color palette - gray is "Other"

cols <- RColorBrewer::brewer.pal(12, "Set3")

cols_graylast <- c(cols[cols != "#D9D9D9"], "#D9D9D9")

# regional map

spp_dir_diffs_map <- spp_dir_diffs %>%
  mutate(plot_spp = case_when(aou %in% top_spp$aou ~ COMMONNAME,
                              TRUE ~ "Other")) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  mutate(sign = case_when(SPEC %in% c("CLSW", "EUCD", "CORA") ~ "POS",
                          SPEC == "HOFI" & longitude < -75 ~ "POS",
                          SPEC == "TRES" & latitude < 41 ~ "POS",
                          TRUE ~ "NEG")) %>%
  filter(longitude > -130) %>%
  st_as_sf(coords = c("longitude", "latitude"))

spp_dir_diffs_map$spp_fct <- fct_relevel(as.factor(spp_dir_diffs_map$plot_spp),
                                         "Other", after = Inf)
spp_signs <- spp_dir_diffs_map %>%
  ungroup() %>%
  dplyr::select(sign) %>%
  filter(sign == "POS")

hi_lev_spp <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(spp_dir_diffs_map) + tm_bubbles(col = "spp_fct", size = "dir_diff", 
                                           title.col = "Species", title.size = "Turnover impact",
                                           palette = cols_graylast) +
  tm_shape(spp_signs) + tm_symbols(shape = 3, size = 0.05, col = "black", alpha = 0.5) +
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 500, 1000), text.size = 1.5, position = c("left", "bottom")) +
  tm_layout(legend.show = F, main.title = "B. Scale: regional")

# local map

spp_dir_diffs_map_local <- spp_dir_diffs_local %>%
  mutate(plot_spp = case_when(aou %in% top_spp$aou ~ COMMONNAME,
                              TRUE ~ "Other")) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  mutate(sign = case_when(SPEC %in% c("CLSW", "EUCD", "CORA") ~ "POS",
                          SPEC == "HOFI" & longitude < -75 ~ "POS",
                          SPEC == "TRES" & latitude < 41 ~ "POS",
                          TRUE ~ "NEG")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

spp_dir_diffs_map_local$spp_fct <- fct_relevel(as.factor(spp_dir_diffs_map_local$plot_spp),
                                               "Other", after = Inf)

spp_signs_local <- spp_dir_diffs_map_local %>%
  ungroup() %>%
  dplyr::select(sign) %>%
  filter(sign == "POS")

hi_lev_spp_local <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(spp_dir_diffs_map_local) + tm_bubbles(col = "spp_fct", size = "dir_diff", 
                                                 title.col = "Species", title.size = "Turnover impact",
                                                 palette = cols_graylast) +
  tm_shape(spp_signs_local) + tm_symbols(shape = 3, size = 0.05, col = "black", alpha = 0.5) +
  tm_layout(legend.position=c(0.82, 0.02), main.title = "A. Scale: local")

# multi-panel

spp_loo_maps <- tmap_arrange(hi_lev_spp_local, hi_lev_spp, ncol = 1)
tmap_save(spp_loo_maps, "figures/spec-LOO-dir_map.pdf", units = "in", height = 11.5, width = 8.25)

#### Suppl: hi lev species violins by spp ####

reg_spp <- spp_dir_unnest %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC)) %>%
  filter(aou %in% top_spp$aou) %>%
  dplyr::select(focal_rte, aou, SPEC, excl_dir)

local_spp <- spp_dir_unnest_local %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC)) %>%
  filter(aou %in% top_spp$aou) %>%
  dplyr::select(focal_rte, aou, SPEC, excl_dir)

loo_viol <- reg_spp %>%
  left_join(local_spp, by = c("focal_rte", "aou", "SPEC"), suffix = c("_regional", "_local")) %>%
  pivot_longer(names_to = "scale", values_to = "excl_dir", excl_dir_regional:excl_dir_local)

ggplot(loo_viol, aes(x = fct_reorder(SPEC, excl_dir), y = excl_dir, fill = scale)) + geom_violin(draw_quantiles = c(0.5), trim = F)+ 
  coord_flip() +
  labs(x = "Species", y = "Turnover impact", fill = "Scale") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("excl_dir_regional" = "Regional", "excl_dir_local" = "Local")) +
ggsave("figures/spp_loo_violins.pdf", units = "in", height = 6, width = 9)

### Leave-one-out directionality species guilds ####
# local and regional scales

guild_core_spp <- core_spp %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  left_join(dplyr::select(habitat_guilds, species, sci_name, Breeding.Biome), by = c("CommonName" = "species"))

# write trait table for MS supplement
trait_list <- guild_core_spp %>%
  ungroup() %>%
  dplyr::select(aou, CommonName, sci_name, Breeding.Biome, migclass, Foraging, Trophic.Group) %>%
  distinct() %>%
  filter(!(is.na(CommonName)))
# write.csv(trait_list, "data/derived_data/species_trait_list_ms.csv", row.names = F)

# add missing species from Rosenberg (Bicknell's Thrush, Gunnison Sage-Grouse)
# update common names

guild_spp_list <- read.csv("data/derived_data/species_trait_list_ms_updated_names.csv", stringsAsFactors = F) %>%
  right_join(core_spp)

guild_list <- guild_spp_list %>%
  ungroup() %>%
  dplyr::select(aou, Breeding.Biome, migclass, Foraging, Trophic.Group) %>%
  distinct()

# leave one out directionality calculations
# local and regional scales

local_grped <- local_rtes %>%
  group_by(focal_rte) %>%
  nest(data_local = c("stateroute"))

guild_excl_dirs <- regional_rtes %>%
  group_by(focal_rte) %>%
  nest() %>%
  left_join(local_grped) %>%
  mutate(spp_list = map(data, ~{
    df <- .
    guild_core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou, Foraging, Trophic.Group, migclass, Breeding.Biome)
  }),
  spp_list_local = map(data_local, ~{
    df <- .
    guild_core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou, Foraging, Trophic.Group, migclass, Breeding.Biome)
  }),
  n_spp = map_dbl(spp_list, ~nrow(.)),
  n_spp_local = map_dbl(spp_list_local, ~nrow(.))) %>%
  filter(n_spp > 0, n_spp_local > 0) %>%
  mutate(forage_dir = map2(data, spp_list, ~{
    r <- .x
    spp <- .y
    
    foraging <- unique(spp$Foraging)
    
    res <- data.frame(Foraging = foraging, excl_dir = NA)
    
    for(sp in foraging) {
      spec <- spp %>%
        filter(Foraging == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Foraging == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Foraging == sp] <- NA
      }
    }
    
    res
    
  }),
  trophic_dir = map2(data, spp_list, ~{
    r <- .x
    spp <- .y
    
    trophic <- unique(spp$Trophic.Group)
    
    res <- data.frame(Trophic.Group = trophic, excl_dir = NA)
    
    for(sp in trophic) {
      spec <- spp %>%
        filter(Trophic.Group == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Trophic.Group == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Trophic.Group == sp] <- NA
      }
    }
    
    res
    
  }),
  mig_dir = map2(data, spp_list, ~{
    r <- .x
    spp <- .y
    
    migclass <- unique(spp$migclass)
    
    res <- data.frame(migclass = migclass, excl_dir = NA)
    
    for(sp in migclass) {
      spec <- spp %>%
        filter(migclass == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$migclass == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$migclass == sp] <- NA
      }
    }
    
    res
    
  }),
  nesting_dir = map2(data, spp_list, ~{
    r <- .x
    spp <- .y
    
    nesting <- unique(spp$Breeding.Biome)
    
    res <- data.frame(Breeding.Biome = nesting, excl_dir = NA)
    
    for(sp in nesting) {
      spec <- spp %>%
        filter(Breeding.Biome == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Breeding.Biome == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Breeding.Biome == sp] <- NA
      }
    }
    
    res
    
  }),
  forage_dir_local = map2(data_local, spp_list_local, ~{
    r <- .x
    spp <- .y
    
    foraging <- unique(spp$Foraging)
    
    res <- data.frame(Foraging = foraging, excl_dir = NA)
    
    for(sp in foraging) {
      spec <- spp %>%
        filter(Foraging == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Foraging == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Foraging == sp] <- NA
      }
    }
    
    res
    
  }),
  trophic_dir_local = map2(data_local, spp_list_local, ~{
    r <- .x
    spp <- .y
    
    trophic <- unique(spp$Trophic.Group)
    
    res <- data.frame(Trophic.Group = trophic, excl_dir = NA)
    
    for(sp in trophic) {
      spec <- spp %>%
        filter(Trophic.Group == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Trophic.Group == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Trophic.Group == sp] <- NA
      }
    }
    
    res
    
  }),
  mig_dir_local = map2(data_local, spp_list_local, ~{
    r <- .x
    spp <- .y
    
    migclass <- unique(spp$migclass)
    
    res <- data.frame(migclass = migclass, excl_dir = NA)
    
    for(sp in migclass) {
      spec <- spp %>%
        filter(migclass == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$migclass == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$migclass == sp] <- NA
      }
    }
    
    res
    
  }),
  nesting_dir_local = map2(data_local, spp_list_local, ~{
    r <- .x
    spp <- .y
    
    nesting <- unique(spp$Breeding.Biome)
    
    res <- data.frame(Breeding.Biome = nesting, excl_dir = NA)
    
    for(sp in nesting) {
      spec <- spp %>%
        filter(Breeding.Biome == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$Breeding.Biome == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$Breeding.Biome == sp] <- NA
      }
    }
    
    res
    
  }))

all_dir_core <- regional_dir_core %>%
  left_join(local_dir_core, by = c("focal_rte"), suffix = c("_regional", "_local"))

forage_over5 <- guild_list %>%
  group_by(Foraging) %>%
  summarize(n_spp = n_distinct(aou)) %>%
  filter(n_spp > 5)

forage_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, forage_dir, forage_dir_local) %>%
  pivot_longer(forage_dir:forage_dir_local, names_to = "scale", values_to = "data") %>%
  unnest(data) %>%
  left_join(all_dir_core) %>%
  mutate(dir_diff = ifelse(scale == "forage_dir", dir_core_regional - excl_dir, dir_core_local - excl_dir)) %>%
  filter(!is.na(Foraging), Foraging %in% forage_over5$Foraging) %>%
  left_join(forage_over5) %>%
  mutate(forage_plot = paste0(str_to_title(Foraging), " (", n_spp, ")"))
# write.csv(forage_dir_diffs, "data/derived_data/guild_LOO_dir_impact_foraging.csv", row.names = F)
forage_dir_diffs <- read.csv("data/derived_data/guild_LOO_dir_impact_foraging.csv", stringsAsFactors = F)

trophic_over5 <- guild_list %>%
  group_by(Trophic.Group) %>%
  summarize(n_spp = n_distinct(aou)) %>%
  filter(n_spp > 5)

trophic_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, trophic_dir, trophic_dir_local) %>%
  pivot_longer(trophic_dir:trophic_dir_local, names_to = "scale", values_to = "data") %>%
  unnest(data) %>%
  left_join(all_dir_core) %>%
  mutate(dir_diff = ifelse(scale == "trophic_dir", dir_core_regional - excl_dir, dir_core_local - excl_dir)) %>%
  filter(!is.na(Trophic.Group), Trophic.Group %in% trophic_over5$Trophic.Group) %>%
  left_join(trophic_over5) %>%
  mutate_at(c("Trophic.Group"), ~ifelse(. == "insct/om", "insect/omnivore", .)) %>%
  mutate(trophic_plot = paste0(str_to_title(Trophic.Group), " (", n_spp, ")"))
# write.csv(trophic_dir_diffs, "data/derived_data/guild_LOO_dir_impact_trophic.csv", row.names = F)
trophic_dir_diffs <- read.csv("data/derived_data/guild_LOO_dir_impact_trophic.csv", stringsAsFactors = F)

mig_over5 <- guild_list %>%
  group_by(migclass) %>%
  summarize(n_spp = n_distinct(aou)) %>%
  filter(n_spp > 5)

mig_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, mig_dir, mig_dir_local) %>%
  pivot_longer(mig_dir:mig_dir_local, names_to = "scale", values_to = "data") %>%
  unnest(data) %>%
  left_join(all_dir_core) %>%
  mutate(dir_diff = ifelse(scale == "mig_dir", dir_core_regional - excl_dir, dir_core_local - excl_dir)) %>%
  filter(!is.na(migclass)) %>%
  left_join(mig_over5) %>%
  mutate_at(c("migclass"), ~case_when(. == "short" ~ "short distance",
                                      . == "resid" ~ "resident",
                                      . == "neotrop" ~ "long distance")) %>%
  mutate(mig_plot = paste0(str_to_title(migclass), " (", n_spp, ")"))
# write.csv(mig_dir_diffs, "data/derived_data/guild_LOO_dir_impact_migclass.csv", row.names = F)
mig_dir_diffs <- read.csv("data/derived_data/guild_LOO_dir_impact_migclass.csv", stringsAsFactors = F)

hab_over5 <- guild_list %>%
  group_by(Breeding.Biome) %>%
  summarize(n_spp = n_distinct(aou)) %>%
  filter(n_spp > 5)

hab_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, nesting_dir, nesting_dir_local) %>%
  pivot_longer(nesting_dir:nesting_dir_local, names_to = "scale", values_to = "data") %>%
  unnest(data) %>%
  left_join(all_dir_core) %>%
  mutate(dir_diff = ifelse(scale == "nesting_dir", dir_core_regional - excl_dir, dir_core_local - excl_dir)) %>%
  filter(!is.na(Breeding.Biome), Breeding.Biome %in% hab_over5$Breeding.Biome) %>%
  left_join(hab_over5) %>%
  mutate(hab_plot = paste0(str_to_title(Breeding.Biome), " (", n_spp, ")"))
# write.csv(hab_dir_diffs, "data/derived_data/guild_LOO_dir_impact_habitat.csv", row.names = F)
hab_dir_diffs <- read.csv("data/derived_data/guild_LOO_dir_impact_habitat.csv", stringsAsFactors = F)

foraging_plot <- ggplot(forage_dir_diffs, aes(x = forage_plot, y = dir_diff, fill = scale)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, col = NA) +
  stat_summary(aes(group = scale), fun.y = "median", fun.ymin = "median", fun.ymax = "median",
               geom = "crossbar", width = 0.8, position = position_dodge(width = 0.9), col = "gray50", show.legend = F) +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  annotate("segment", x = 3, xend = 3,  y = 0.06, yend = 0.08, arrow = arrow(), size = 2, col = "firebrick") +
  annotate("text", 2, 0.075, label = "Higher turnover\nincluding group", size = 5, col = "firebrick") +
  labs(x = "Foraging guild", y = "", fill = "") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("forage_dir" = "Regional", "forage_dir_local" = "Local")) +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 16)) + coord_flip() 

trophic_plot <- ggplot(trophic_dir_diffs, aes(x = trophic_plot, y = dir_diff, fill = scale)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, col = NA) +
  stat_summary(aes(group = scale), fun.y = "median", fun.ymin = "median", fun.ymax = "median",
               geom = "crossbar", width = 0.8, position = position_dodge(width = 0.9), col = "gray50", show.legend = F) +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  annotate("segment", x = 2.25, xend = 2.25,  y = -0.04, yend = -0.06, arrow = arrow(), size = 2, col = "dodgerblue") +
  annotate("text", 1.25, -0.05, label = "Higher turnover\nexcluding group", size = 5, col = "dodgerblue") +
  labs(x = "Trophic group", y = "", fill = "") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("trophic_dir" = "Regional", "trophic_dir_local" = "Local")) +
  theme(legend.position = c(0.2, 0.9), axis.text = element_text(size = 15), axis.title = element_text(size = 16), legend.text = element_text(size = 14)) + coord_flip() 

mig_plot <- ggplot(mig_dir_diffs, aes(x = mig_plot, y = dir_diff, fill= scale)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, col = NA) +
  stat_summary(aes(group = scale), fun.y = "median", fun.ymin = "median", fun.ymax = "median",
               geom = "crossbar", width = 0.8, position = position_dodge(width = 0.9), col = "gray50") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Migration distance", y = "Turnover impact") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("mig_dir" = "Regional", "mig_dir_local" = "Local")) +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 16)) + coord_flip() 

# merge forest, remove introduced/wetland groups to simplify hab_plot

hab_dir_diffs_few <- hab_dir_diffs %>%
  filter(Breeding.Biome != "Wetland", Breeding.Biome != "Introduced") %>%
  mutate(biome_few = case_when(grepl("Forest", Breeding.Biome) ~ "Forest",
                               TRUE ~ Breeding.Biome)) %>%
  group_by(biome_few) %>%
  mutate(n_few = sum(unique(n_spp))) %>%
  mutate(hab_plot_few = paste0(biome_few, " (", n_few, ")"))

hab_plot <- ggplot(hab_dir_diffs_few, aes(x = hab_plot_few, y = dir_diff, fill = scale)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, col = NA) +
  stat_summary(aes(group = scale), fun.y = "median", fun.ymin = "median", fun.ymax = "median",
               geom = "crossbar", width = 0.8, position = position_dodge(width = 0.9), col = "gray50") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Breeding biome", y = "Turnover impact") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("nesting_dir" = "Regional", "nesting_dir_local" = "Local")) +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 16)) + coord_flip() 

plot_grid(trophic_plot, foraging_plot, mig_plot, hab_plot, ncol = 2, 
          labels = c("A", "B", "C", "D"), label_size = 16)
ggsave("figures/guild_LOO_directionality.pdf", units = "in", height = 8, width = 13)

### Supplemental habitat groups plot with all groups

ggplot(hab_dir_diffs, aes(x = hab_plot, y = dir_diff, fill = scale)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, col = NA) +
  stat_summary(aes(group = scale), fun.y = "median", fun.ymin = "median", fun.ymax = "median",
               geom = "crossbar", width = 0.8, position = position_dodge(width = 0.9), col = "gray50", show.legend = F) +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Breeding biome", y = "Turnover impact", fill = "") +
  scale_fill_manual(values = c("skyblue3", "gray"), labels = c("nesting_dir" = "Regional", "nesting_dir_local" = "Local")) +
  theme(legend.position = c(0.8, 0.15)) + coord_flip() 
ggsave("figures/breeding_LOO_directionality.pdf")

#### Route env change figures ####

## Climate trends at 25 route scale

climate_25route <- scale_model_variables_unnest %>%
  filter(scale == 25) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  select(focal_rte, longitude, latitude, trend_tmax, trend_tmin) %>%
  filter(!is.na(trend_tmax), !is.na(trend_tmin)) %>%
  pivot_longer(trend_tmax:trend_tmin, names_to = "climate", values_to = "trend") %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)


climate_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(climate_25route) + 
  tm_dots(col = "trend", palette = "-RdBu", size = 0.5, title = "Trend (deg/year)") + 
  tm_layout(legend.text.size = 1, legend.title.size =2) +
  tm_facets(by = "climate")
tmap_save(climate_map, "figures/regional_climate_trends_map.pdf", units = "in", height = 10, width = 8)

## Land cover values at 1 route scale

lc_1route <- scale_model_variables_unnest %>%
  filter(scale == 1) %>%
  filter(!is.na(max_lc_class), max_lc_class != "Water") %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

theme_set(theme_classic(base_size = 20))
class_change <- ggplot(lc_1route, aes(x = fct_rev(max_lc_class), y = max_lc, fill = max_lc_class)) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "", y = "Change in proportion cover", fill = "", title = "B") + 
  coord_flip() +
  theme(plot.title = element_text(size = 18, hjust = -0.25), 
        plot.margin = unit(c(1.25, 0, 0, 0), "cm"))

lc_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "max_lc_class", size = 0.5, title = "Land cover") + 
  tm_compass(type = "8star", position = c("left", "top")) +
  tm_scale_bar(breaks = c(0, 500, 1000), text.size = 1.5, position = c("left", "bottom")) +
  tm_layout(main.title = "A", title.size = 4, legend.text.size = 1.25, legend.title.size =2, legend.position = c("right", "bottom"), outer.margins = c(0.01,0.01,0.01,0.01))

color_scale <- data.frame(color = c(1:4), 
                          temp_hex = c("#92C5DE", "#FDDBC7", "#EF8A62", "#B2182B"), stringsAsFactors = F)

tmin_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "trend_tmin", size = 0.5, title = "Trend in Tmin", breaks = quantile(lc_1route$trend_tmin, na.rm = T), palette = color_scale$temp_hex) + 
  tm_layout(main.title = "C", title.size = 4, legend.show = F, 
            inner.margins = c(0.12, 0.02, 0.02, 0.02), 
            outer.margins = c(0.01,0.01,0.01,0.01))


tmax_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "trend_tmax", size = 0.5, title = "Trend in Tmax", 
          breaks = quantile(lc_1route$trend_tmax, na.rm = T), palette = color_scale$temp_hex) + 
  tm_layout(main.title = "D", title.size = 4, legend.show = F, 
            inner.margins = c(0.12, 0.02, 0.02, 0.02), 
            outer.margins = c(0.01,0.01,0.01,0.01))

tmin_quant <- quantile(lc_1route$trend_tmin, na.rm = T)
tmax_quant <- quantile(lc_1route$trend_tmax, na.rm = T)

hist_breaks <- function(quantiles, bins) {
  breaks <- seq(min(quantiles), max(quantiles), length.out = bins)
  
  breaks_total <- c(breaks, quantiles) %>% sort()
  
  return(unique(breaks_total))
}

tmin_breaks <- hist_breaks(tmin_quant, 30)
tmax_breaks <- hist_breaks(tmax_quant, 30)

quantile_group <- function(quantiles, value) {
  case_when(
    value <= quantiles[2] ~ 1,
    value > quantiles[2] & value <= quantiles[3] ~ 2,
    value > quantiles[3] & value <= quantiles[4] ~ 3,
    value > quantiles[4] & value <= quantiles[5] ~ 4)
}

lc_1route_plots <- lc_1route %>%
  mutate(tmin_color = quantile_group(tmin_quant, trend_tmin),
         tmax_color = quantile_group(tmax_quant, trend_tmax))

tmin_hist <- ggplot(lc_1route_plots, aes(x = trend_tmin, fill = as.factor(tmin_color))) + 
  geom_histogram(breaks = tmin_breaks) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  labs(x = "Trend in Tmin", y = "Count") +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_fill_manual(values = color_scale$temp_hex) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(r = 0)),
        legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  )

tmax_hist <- ggplot(lc_1route_plots, aes(x = trend_tmax, fill = as.factor(tmax_color))) + 
  geom_histogram(breaks = tmax_breaks) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  labs(x = "Trend in Tmax", y = "Count") +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1)) +
  scale_fill_manual(values = color_scale$temp_hex) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(r = 0)),
        legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  )

vp_tmin <- viewport(0.095, 0.09, width = 0.19, height = 0.15)
vp_tmax <- viewport(0.595, 0.09, width = 0.19, height = 0.15)

grid.newpage()
pdf(paste0(getwd(), "/figures/max_landcover_multipanel.pdf"), height = 12, width = 15)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
print(lc_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(class_change, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(tmin_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(tmax_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(tmin_hist, vp = vp_tmin)
print(tmax_hist, vp = vp_tmax)
dev.off()
