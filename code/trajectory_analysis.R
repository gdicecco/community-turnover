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
library(cowplot)
library(broom)
library(grid)

#### Set up ####

## Plotting theme

theme_set(theme_classic(base_size = 15))

## North America map

na <- read_sf("data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN") %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

## Read in data

# Append correct BioArk path
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

# Species list 

species_list <- read.csv("data/species_list.csv", stringsAsFactors = F)
fourletter_codes <- read.csv("data/four_letter_codes_aous.csv", stringsAsFactors = F)

# Species traits/guilds

bird_traits <- read.csv("data/Master_RO_Correlates_20110610.csv", stringsAsFactors = F) %>%
  select(AOU, CommonName, Foraging, Trophic.Group, migclass)

habitat_guilds <- read.csv("data/habitat_guilds_new_aous.csv", stringsAsFactors = F)

# BBS sampled every 5 years from 1970 to 2016
log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)

# BBS sampled 3/4 years every 4 year time window
bbs_subset <- read.csv("data/bbs_counts_subset_1990-2016.csv", stringsAsFactors = F)

routes <- read.csv(paste0(bioark, "/hurlbertlab/databases/BBS/2017/bbs_routes_20170712.csv"), stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

# BBS half-route directionality

half_route_dir <- read.csv("data/half_route_directionality.csv", stringsAsFactors = F)

# Land cover and climate data - 1/2 route scale

bbs_half_landcover <- read.csv('data/bbs_half_route_max_land_change.csv', stringsAsFactors = F)
bbs_half_climate <- read.csv("data/bbs_half_route_breeding_season_climate.csv", stringsAsFactors = F)

# Land cover and climate data - 1 route scale

bbs_landcover <- read.csv("data/bbs_route_max_landcover_change.csv", stringsAsFactors = F)
bbs_climate <- read.csv("data/bbs_routes_climate_trends.csv", stringsAsFactors = F)

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

mean_dist_plot <- ggplot(filter(mean_bcr_distances, bcr %in% bcr_subset$bcr), 
       aes(x = scale, y = mean_dist_bcr, col = bcr, group = bcr)) + 
  geom_line(cex = 1) +
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

pct_plot <- ggplot(mean_pct_overlap, aes(x = scale, y = mean_reps, col = bcr, group = bcr)) + 
  geom_line(cex = 1) +
  scale_color_viridis_c() +
  labs(x = "Scale (routes)", y = "Avg percent of aggregates route occurs in", col = "BCR")
# ggsave("figures/percent_overlap_aggregates.pdf")

## Aggregation methods fig

bcr <- read_sf("\\\\BioArk//HurlbertLab//DiCecco//bcr_terrestrial_shape//BCR_Terrestrial_master.shp") %>%
  filter(BCR %in% bcr_subset$bcr) %>%
  mutate_at(c("BCR"), ~as.factor(.))

study_routes <- bbs_subset %>%
  ungroup() %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  distinct(stateroute, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

bcr_map <- tm_shape(na) + tm_polygons(col = "gray50") +
  tm_shape(bcr) + tm_polygons(col = "BCR") +
  tm_shape(study_routes) + tm_dots(col = "black", size = 0.05) + 
  tm_layout(legend.text.size = 1.25, legend.title.size = 1.5, outer.margins = c(0.01,0,0.01,0))

bcr_palette <- tmaptools::get_brewer_pal("Set3", n = 15)
names(bcr_palette) <- as.factor(unique(bcr$BCR))

mean_dist_plot <- ggplot(filter(mean_bcr_distances, bcr %in% bcr_subset$bcr), 
                         aes(x = scale, y = mean_dist_bcr, col = as.factor(bcr), group = as.factor(bcr))) + 
  geom_line(cex = 1) +
  scale_color_manual(values = bcr_palette) +
  labs(x = "Scale (routes)", y = "Mean distance between routes (km)", col = "BCR") +
  theme(legend.position = "none")

pct_plot <- ggplot(mean_pct_overlap, aes(x = scale, y = mean_reps, col = as.factor(bcr), group = as.factor(bcr))) + 
  geom_line(cex = 1) +
  scale_color_manual(values = bcr_palette) +
  labs(x = "Scale (routes)", y = "Avg percent of aggregates route occurs in", col = "BCR") +
  theme(legend.position = "none")

line_panels <- plot_grid(mean_dist_plot, pct_plot, ncol = 2)

grid.newpage()
pdf(paste0(getwd(),"/figures/bcr_aggregation_multipanel.pdf"), height = 12, width = 12)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
print(bcr_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(line_panels, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

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
# write.csv(scale_model_variables_unnest, "data/scale_model_input.csv", row.names = F)

scale_model_variables_unnest <- read.csv("data/scale_model_input.csv", stringsAsFactors = F)

scale_model_variables <- read.csv("data/scale_model_input.csv", stringsAsFactors = F) %>%
  group_by(scale) %>%
  nest()

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

# 1/2 route model

half_climate_trends <- bbs_half_climate %>%
  left_join(routes) %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016),
         max_bins = case_when(countrynum == 124 ~ 5,
                              countrynum == 840 ~ 6)) %>%
  filter(year >= y1, year <= y2) %>%
  group_by(stateroute, stops) %>%
  nest() %>%
  mutate(trends = map(data, ~{
    df <- .
    climate_trend <- possibly_climate_trend(df)
    data.frame(trend_tmax = climate_trend$trend_tmax, trend_tmin = climate_trend$trend_tmin)
  })) %>%
  dplyr::select(-data) %>%
  unnest(trends)

half_route_mod_input <- half_route_dir %>%
  mutate_at("stops", ~case_when(. == "stops1_25" ~ "1-25",
                                . == "stops26_50" ~ "26-50")) %>%
  left_join(half_climate_trends) %>%
  left_join(bbs_half_landcover) %>%
  left_join(routes) %>%
  filter(stateroute %in% bbs_subset$stateroute, bcr %in% bcr_subset$bcr) %>%
  rename(max_lc = "deltaCover")
# write.csv(half_route_mod_input, "data/half_route_scale_model_input.csv", row.names = F)

mod1 <- lm(dir_core ~ trend_tmax + trend_tmin, data = half_route_mod_input)
mod2 <- lm(dir_core ~ max_lc, data = half_route_mod_input)
mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = half_route_mod_input)

mod1.r2 <- summary(mod1)$r.squared
mod2.r2 <- summary(mod2)$r.squared
mod12.r2 <- summary(mod12)$r.squared

part1 <- mod12.r2 - mod2.r2 # climate alone
part2 <- mod12.r2 - mod1.r2 # land cover alone
joined <- mod1.r2 - part1 # shared variance
unexpl <- 1 - mod12.r2 # unexplained variance

half_route_output <- data.frame(scale = 0.5, part1 = part1, part2 = part2,
                                       joined = joined, unexpl = unexpl)  

scale_model_output_all <- bind_rows(scale_model_output, half_route_output)
# write.csv(scale_model_output, "data/scale_model_output_variance.csv", row.names = F)
scale_model_output <- read.csv("data/scale_model_output_variance.csv")

  
# Scale model figures 

scale_model_plot <- scale_model_output %>%
  pivot_longer(part1:unexpl, names_to = "variance")

all_routes <- ggplot(filter(scale_model_plot, variance == "part1" | variance == "part2"), aes(x = scale, y = value, fill = variance)) +
  geom_col(position = "stack", col = "white") +
  scale_fill_manual(values = c("#92C5DE", "#A6D854"), name = "Variance explained", labels = c("Climate", "Land cover")) +
  labs(x = "Aggregated routes", y = "Variance")
ggsave("figures/scale_model_variance.pdf")

# Scale model predictor effects

scale_model_variables <-  read.csv("data/scale_model_input.csv", stringsAsFactors = F) %>%
  group_by(scale) %>%
  nest()

scale_model_output <- data.frame(scale = c(), part1 = c(), part2 = c(), joined = c(), unexpl = c())
for(i in 1:25) {
  model_input <- scale_model_variables$data[[i]]
  
  mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = model_input)
  
  tidy_mod <- tidy(mod12) %>%
    mutate(scale = i,
           conf_lo = confint(mod12)[, 1],
           conf_hi = confint(mod12)[, 2])
  
  scale_model_output <- rbind(scale_model_output, 
                              tidy_mod)
}

LC <- ggplot(filter(scale_model_output, term == "max_lc"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of change in max land cover")

tmax <- ggplot(filter(scale_model_output, term == "trend_tmax"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of trend in Tmax")

tmin <- ggplot(filter(scale_model_output, term == "trend_tmin"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of trend in Tmin")

plot_grid(tmax, tmin, LC, nrow = 2)
ggsave("figures/scale_model_effect_ests.pdf", units = "in", height = 9, width = 12)

## Scale model partial effects plots

pdf("figures/scale_model_tmax_effects.pdf")
for(scale in scale_model_variables$scale) {
  df <- scale_model_variables$data[[scale]]
  print(ggplot(df, aes(x = trend_tmax, y = dir_core)) + geom_point() + geom_smooth(method = "lm", se = F) +
          labs(title = paste0("Routes = ", scale)))
}
dev.off()

pdf("figures/scale_model_tmin_effects.pdf")
for(scale in scale_model_variables$scale) {
  df <- scale_model_variables$data[[scale]]
  print(ggplot(df, aes(x = trend_tmin, y = dir_core)) + geom_point() + geom_smooth(method = "lm", se = F) +
          labs(title = paste0("Routes = ", scale)))
}
dev.off()


  
# Map of directionality values

dir_sf <- scale_model_variables_unnest %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  filter(statenum != 3) %>%
  st_as_sf(coords = c("longitude", "latitude"))

one_route <- tm_shape(na) + tm_polygons(col = "gray50") +
  tm_shape(filter(dir_sf, scale == 1)) + 
  tm_dots(col = "dir_core", title = "Directionality", palette = "YlGnBu", size = 0.3, legend.show = F) +
  tm_layout(main.title = "Route-level")
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
  labs(x = "Directionality", y = "Count") + theme(legend.position = "none")


tf_route <-  tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(dir_sf, scale == 25)) + 
  tm_dots(col = "dir_core", title = "Directionality", palette = "YlGnBu", size = 0.3, legend.show = F) +
  tm_layout(main.title = "25 Routes")
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
  labs(x = "Directionality", y = "Count") + theme(legend.position = "none")

dir_map <- tmap_arrange(one_route, tf_route, ncol = 2)

vp_one <- viewport(0.1, 0.15, width = 0.175, height = 0.25)
vp_tf <- viewport(0.6, 0.15, width = 0.175, height = 0.25)

wd <- getwd()
pdf(paste0(wd, "/figures/directionality_map.pdf"), height = 8, width = 16)
print(dir_map)
print(one_legend, vp = vp_one)
print(tf_legend, vp = vp_tf)
dev.off()

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
# write.csv(low_overlap_focal_routes, "data/low_overlap_focal_routes.csv", row.names = F)

low_overlap_focal_routes <- read.csv("data/low_overlap_focal_routes.csv")

## Run models for just focal routes 

scale_model_variables_unnest <- read.csv("data/scale_model_input.csv", stringsAsFactors = F)
half_route_mod_input <- read.csv("data/half_route_scale_model_input.csv", stringsAsFactors = F)

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
# write.csv(LO_scale_model_output, "data/low_overlap_scale_model_output_deviance.csv", row.names = F)
LO_scale_model_output <- read.csv("data/low_overlap_scale_model_output_deviance.csv", stringsAsFactors = F)

half_route_variables <- half_route_mod_input %>%
  filter(stateroute %in% low_overlap_focal_routes$focal_rte)

mod1 <- lm(dir_core ~ trend_tmax + trend_tmin, data = half_route_mod_input)
mod2 <- lm(dir_core ~ max_lc, data = half_route_mod_input)
mod12 <- lm(dir_core ~ trend_tmax + trend_tmin + max_lc, data = half_route_mod_input)

mod1.r2 <- summary(mod1)$r.squared
mod2.r2 <- summary(mod2)$r.squared
mod12.r2 <- summary(mod12)$r.squared

part1 <- mod12.r2 - mod2.r2 # climate alone
part2 <- mod12.r2 - mod1.r2 # land cover alone
joined <- mod1.r2 - part1 # shared variance
unexpl <- 1 - mod12.r2 # unexplained variance

half_route_output <- data.frame(scale = 0.5, part1 = part1, part2 = part2,
                                joined = joined, unexpl = unexpl) 

LO_scale_model_plot <- LO_scale_model_output %>%
#  bind_rows(half_route_output) %>%
  pivot_longer(part1:unexpl, names_to = "variance")

lo_overlap <- ggplot(filter(LO_scale_model_plot, variance == "part1" | variance == "part2"), aes(x = scale, y = value, fill = variance)) +
  geom_col(position = "stack", width = 0.9, col = "white") +
  scale_fill_manual(values = c("#92C5DE", "#A6D854"), name = "Variance explained", labels = c("Climate", "Land cover")) +
  labs(x = "Aggregated routes", y = "Variance")
ggsave("figures/low_overlap_scale_model_variance.pdf")

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

LC <- ggplot(filter(scale_model_output, term == "max_lc"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of change in max land cover")

tmax <- ggplot(filter(scale_model_output, term == "trend_tmax"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of trend in Tmax")

tmin <- ggplot(filter(scale_model_output, term == "trend_tmin"), aes(x = scale, y = estimate)) + 
  geom_hline(yintercept = 0, col = "darkgray", lty = 2, cex = 1) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = conf_lo, ymax = conf_hi), width = 0) +
  labs(x = "Scale", y = "Effect of trend in Tmin")

plot_grid(tmax, tmin, LC, nrow = 2)
ggsave("figures/scale_model_lo_overlap_effect_ests.pdf", units = "in", height = 9, width = 12)


## Scale model partial effects plots

pdf("figures/scale_model_lo_overlap_tmax_effects.pdf")
for(scale in low_overlap_model_variables$scale) {
  df <- low_overlap_model_variables$data[[scale]]
  print(ggplot(df, aes(x = trend_tmax, y = dir_core)) + geom_point() + geom_smooth(method = "lm", se = F) +
          labs(title = paste0("Routes = ", scale)))
}
dev.off()

pdf("figures/scale_model_lo_overlap_tmin_effects.pdf")
for(scale in low_overlap_model_variables$scale) {
  df <- low_overlap_model_variables$data[[scale]]
  print(ggplot(df, aes(x = trend_tmin, y = dir_core)) + geom_point() + geom_smooth(method = "lm", se = F) +
          labs(title = paste0("Routes = ", scale)))
}
dev.off()


#### Explaining high directionality ####

# Directionality values at 25 route scales
# Regional high leverage species

scale_model_variables_unnest <- read.csv("data/scale_model_input.csv", stringsAsFactors = F)
abund_trends <- read.csv("data/BBS_abundance_trends.csv", stringsAsFactors = F)
low_overlap_focal_routes <- read.csv("data/low_overlap_focal_routes.csv")

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

regional_abund_trends <- scale_model_variables_unnest %>%
  filter(focal_rte %in% low_overlap_focal_routes$focal_rte) %>%
  filter(scale == 25) %>%
  left_join(regional_rtes) %>%
  left_join(abund_trends) %>%
  group_by(focal_rte, dir_all, dir_core, aou) %>%
  summarize(abund_trend = mean(abundTrend, na.rm = T))

spp_cor <- regional_abund_trends %>%
  group_by(aou) %>%
  nest() %>%
  mutate(n_regions = map_dbl(data, ~nrow(.))) %>%
  filter(n_regions > 15) %>%
  mutate(abund_dir_r = map_dbl(data, ~cor(.$dir_core, .$abund_trend))) %>%
  select(-data) %>%
  ungroup() %>%
  left_join(species_list) %>%
  filter(!grepl("unid.", english_common_name)) %>%
  arrange(desc(abs(abund_dir_r))) %>%
  slice(1:10) %>%
  select(aou, n_regions, abund_dir_r, english_common_name)
# write.csv(spp_cor, "data/high_leverage_spp.csv", row.names = F)

# Compare correlations of directionality with habitat, trophic guild abund trends

foraging_cor <- regional_abund_trends %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  filter(!is.na(Foraging)) %>%
  group_by(Foraging) %>%
  nest() %>%
  mutate(cor_test = map(data, ~cor.test(.$dir_core, .$abund_trend)),
         abund_dir_r = map_dbl(cor_test, ~.$estimate),
         conf_low = map_dbl(cor_test, ~.$conf.int[[1]]),
         conf_hi = map_dbl(cor_test, ~.$conf.int[[2]])) %>%
  select(-data)

trophic_cor <- regional_abund_trends %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  filter(!is.na(Trophic.Group)) %>%
  group_by(Trophic.Group) %>%
  nest() %>%
  mutate(cor_test = map(data, ~cor.test(.$dir_core, .$abund_trend)),
         abund_dir_r = map_dbl(cor_test, ~.$estimate),
         conf_low = map_dbl(cor_test, ~.$conf.int[[1]]),
         conf_hi = map_dbl(cor_test, ~.$conf.int[[2]])) %>%
  select(-data)  

mig_cor <- regional_abund_trends %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  filter(!is.na(migclass)) %>%
  group_by(migclass) %>%
  nest() %>%
  mutate(cor_test = map(data, ~cor.test(.$dir_core, .$abund_trend)),
         abund_dir_r = map_dbl(cor_test, ~.$estimate),
         conf_low = map_dbl(cor_test, ~.$conf.int[[1]]),
         conf_hi = map_dbl(cor_test, ~.$conf.int[[2]])) %>%
  select(-data)

habitat_cor <- regional_abund_trends %>%
  left_join(habitat_guilds) %>%
  filter(!is.na(nesting_group)) %>%
  group_by(nesting_group) %>%
  nest() %>%
  mutate(cor_test = map(data, ~cor.test(.$dir_core, .$abund_trend)),
         abund_dir_r = map_dbl(cor_test, ~.$estimate),
         conf_low = map_dbl(cor_test, ~.$conf.int[[1]]),
         conf_hi = map_dbl(cor_test, ~.$conf.int[[2]])) %>%
  select(-data)

foraging_plot <- ggplot(foraging_cor, aes(x = Foraging, y = abund_dir_r)) +
  geom_hline(yintercept = 0, cex = 1, col = "red", lty = 2) +
  labs(x = "Foraging guild", y = " ") +
  geom_point(cex = 2) + 
  geom_errorbar(aes(ymin = conf_low, ymax = conf_hi), width = 0) +
  coord_flip() 

trophic_plot <- ggplot(trophic_cor, aes(x = Trophic.Group, y = abund_dir_r)) +
  geom_hline(yintercept = 0, cex = 1, col = "red", lty = 2) +
  labs(x = "Trophic Group", y = " ") +
  geom_point(cex = 2) + 
  geom_errorbar(aes(ymin = conf_low, ymax = conf_hi), width = 0) + 
  coord_flip() 

mig_plot <- ggplot(mig_cor, aes(x = migclass, y = abund_dir_r)) +
  geom_hline(yintercept = 0, cex = 1, col = "red", lty = 2) +
  labs(x = "Migration distance", y = "") +
  geom_point(cex = 2) + 
  geom_errorbar(aes(ymin = conf_low, ymax = conf_hi), width = 0) +
  coord_flip() 

hab_plot <- ggplot(habitat_cor, aes(x = nesting_group, y = abund_dir_r)) +
  geom_hline(yintercept = 0, cex = 1, col = "red", lty = 2) +
  labs(x = "Nesting habitat", y = "Correlation with directionality") +
  geom_point(cex = 2) + 
  geom_errorbar(aes(ymin = conf_low, ymax = conf_hi), width = 0) +
  coord_flip() 

plot_grid(foraging_plot, trophic_plot, mig_plot, hab_plot, ncol = 1)
ggsave("figures/guild_cor_directionality.pdf", units = "in", height = 10, width = 6)

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
# write.csv(spp_dir_unnest, "data/species-leave-one-out-directionality.csv", row.names = F)

spp_dir_unnest <- read.csv("data/species-leave-one-out-directionality.csv", stringsAsFactors = F)

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
# write.csv(high_impact_spp, "data/spec-LOO-directionality-hi-impact.csv", row.names = F)

high_impact_spp <- read.csv("data/spec-LOO-directionality-hi-impact.csv", stringsAsFactors = F)

top_ten <- high_impact_spp %>%
  slice(1:10)

spp_dir_diffs_map <- spp_dir_diffs %>%
  mutate(SPEC_10 = case_when(aou %in% top_ten$aou ~ SPEC,
                             TRUE ~ "Other")) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  mutate(sign = case_when(SPEC %in% c("CLSW", "EUCD") ~ "POS",
                          SPEC == "HOFI" & longitude < -75 ~ "POS",
                          SPEC == "TRES" & latitude < 41 ~ "POS",
                          TRUE ~ "NEG")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

spp_signs <- spp_dir_diffs_map %>%
  ungroup() %>%
  dplyr::select(sign) %>%
  filter(sign == "POS")

hi_lev_spp <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(spp_dir_diffs_map) + tm_bubbles(col = "SPEC_10", size = "dir_diff", title.col = "SPEC", title.size = "Change in dir.") +
  tm_shape(spp_signs) + tm_symbols(shape = 3, size = 0.05, col = "black", alpha = 0.5) +
  tm_layout(legend.position=c("right", "bottom"))
tmap_save(hi_lev_spp, "figures/spec-LOO-dir_map.pdf")

# To what degree is one species dominant in driving directionality - how many species w/in 15% of max dir_diff value

dominant_spp <- spp_dir_unnest %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(dir_diff > 0) %>%
  group_by(focal_rte) %>%
  filter(dir_diff >= max(dir_diff, na.rm = T)-0.15*max(dir_diff, na.rm = T)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))

dominant_spp_list <- dominant_spp %>%
  group_by(aou) %>% 
  summarize(n_regions = n()) %>%
  arrange(desc(n_regions)) %>%
  left_join(fourletter_codes) %>%
  filter(!is.na(SPEC))
# write.csv(dominant_spp_list, "data/spec-LOO-hi-directionality-ties.csv", row.names = F)

n_regions_tied <- dominant_spp %>%
  group_by(focal_rte) %>%
  summarize(n_spp = n_distinct(aou))

nrow(filter(n_regions_tied, n_spp > 1))/nrow(n_regions_tied)

max(n_regions_tied$n_spp)

# Leave-one-out directionality for species guilds

guild_core_spp <- core_spp %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  left_join(habitat_guilds)

guild_excl_dirs <- regional_rtes %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(spp_list = map(data, ~{
    df <- .
    guild_core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou, Foraging, Trophic.Group, migclass, nesting_group)
  }),
  n_spp = map_dbl(spp_list, ~nrow(.))) %>%
  filter(n_spp > 0) %>%
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
    
    nesting <- unique(spp$nesting_group)
    
    res <- data.frame(nesting_group = nesting, excl_dir = NA)
    
    for(sp in nesting) {
      spec <- spp %>%
        filter(nesting_group == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$nesting_group == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$nesting_group == sp] <- NA
      }
    }
    
    res
    
  }))

forage_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, forage_dir) %>%
  unnest(forage_dir) %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(!is.na(Foraging)) 
# write.csv(forage_dir_diffs, "data/guild_LOO_dir_impact_foraging.csv", row.names = F)
forage_dir_diffs <- read.csv("data/guild_LOO_dir_impact_foraging.csv", stringsAsFactors = F)

trophic_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, trophic_dir) %>%
  unnest(trophic_dir) %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(!is.na(Trophic.Group)) 
# write.csv(trophic_dir_diffs, "data/guild_LOO_dir_impact_trophic.csv", row.names = F)
trophic_dir_diffs <- read.csv("data/guild_LOO_dir_impact_trophic.csv", stringsAsFactors = F)

mig_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, mig_dir) %>%
  unnest(mig_dir) %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(!is.na(migclass)) 
# write.csv(mig_dir_diffs, "data/guild_LOO_dir_impact_migclass.csv", row.names = F)
mig_dir_diffs <- read.csv("data/guild_LOO_dir_impact_migclass.csv", stringsAsFactors = F)

hab_dir_diffs <- guild_excl_dirs %>%
  select(focal_rte, nesting_dir) %>%
  unnest(nesting_dir) %>%
  left_join(regional_dir_core) %>%
  mutate(dir_diff = dir_core - excl_dir) %>%
  filter(!is.na(nesting_group)) 
# write.csv(hab_dir_diffs, "data/guild_LOO_dir_impact_habitat.csv", row.names = F)

foraging_plot <- ggplot(forage_dir_diffs, aes(x = Foraging, y = dir_diff)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, fill = "gray") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Foraging guild", y = "Directionality impact") +
  coord_flip() 

trophic_plot <- ggplot(trophic_dir_diffs, aes(x = Trophic.Group, y = dir_diff)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, fill = "gray") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Trophic Group", y = "Directionality impact") +
  coord_flip() 

mig_plot <- ggplot(mig_dir_diffs, aes(x = migclass, y = dir_diff)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, fill = "gray") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Migration distance", y = "Directionality impact") +
  coord_flip() 

hab_plot <- ggplot(hab_dir_diffs, aes(x = nesting_group, y = dir_diff)) +
  geom_violin(draw_quantiles = c(0.5), trim = T, scale = "width", cex = 1, fill = "gray") +
  geom_hline(yintercept = 0, cex = 1, col = "black", lty = 2) +
  labs(x = "Nesting habitat", y = "Directionality impact") +
  coord_flip() 

plot_grid(foraging_plot, trophic_plot, mig_plot, hab_plot, ncol = 1)
ggsave("figures/guild_LOO_directionality.pdf", units = "in", height = 10, width = 6)

### Maps: max delta guild for each route

max_diff_guild_map <- function(df) {
  legend_title <- colnames(df)[[2]]
  max_sf <- df %>%
    group_by(focal_rte) %>%
    filter(dir_diff == max(dir_diff)) %>%
    left_join(routes, by = c("focal_rte" = "stateroute")) %>%
    st_as_sf(coords = c("longitude", "latitude"))
  
  map <- tm_shape(na) + tm_polygons(col = "gray50") + 
    tm_shape(max_sf) + tm_bubbles(col = legend_title, size = "dir_diff", title.size = "Impact on directionality") +
    tm_layout(legend.position=c("right", "bottom"))
  tmap_save(map, paste0("figures/guild_LOO_map_", legend_title, ".pdf"))
}

max_diff_guild_map(mig_dir_diffs)
max_diff_guild_map(forage_dir_diffs)
max_diff_guild_map(hab_dir_diffs)
max_diff_guild_map(trophic_dir_diffs)

### Can max impact guild be attributed to high leverage species?

# Highest impact species per focal_rte
spp_dir_diffs
hi_impact_spp <- unique(spp_dir_diffs$aou)

# Remove high impact species from communities
guild_core_spp <- core_spp %>%
  left_join(bird_traits, by = c("aou" = "AOU")) %>%
  left_join(habitat_guilds) %>%
  filter(!(aou %in% spp_dir_diffs$aou))

# Calculate guild LOO with high impact species excluded
guild_excl_dirs <- regional_rtes %>%
  group_by(focal_rte) %>%
  nest() %>%
  left_join(spp_dir_diffs) %>%
  mutate(spp_list = map(data, ~{
    df <- .
    guild_core_spp %>%
      filter(stateroute %in% df$stateroute) %>%
      ungroup() %>%
      distinct(aou, Foraging, Trophic.Group, migclass, nesting_group)
  }),
  n_spp = map_dbl(spp_list, ~nrow(.))) %>%
  filter(n_spp > 0) %>%
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
        select(-c(contains(as.character(hi_impact_spp)))) %>%
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
        select(-c(contains(as.character(hi_impact_spp)))) %>%
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
        select(-c(contains(as.character(hi_impact_spp)))) %>%
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
    
    nesting <- unique(spp$nesting_group)
    
    res <- data.frame(nesting_group = nesting, excl_dir = NA)
    
    for(sp in nesting) {
      spec <- spp %>%
        filter(nesting_group == sp)
      
      log_core <- log_abund_core %>%
        filter(stateroute %in% r$stateroute) %>%
        select(-c(contains(as.character(hi_impact_spp)))) %>%
        select(-c(contains(as.character(spec$aou)))) %>%
        group_by(year_bin) %>%
        summarize_all(mean, na.rm = T) %>%
        dplyr::select(-stateroute)
      
      if(nrow(log_core) > 3) {
        abund_dist_core <- dist(log_core[, -1])
        dir_core <- trajectoryDirectionality(abund_dist_core, sites = rep(1, nrow(log_core)), surveys = log_core$year_bin)
        
        res$excl_dir[res$nesting_group == sp] <- dir_core
      }
      
      else {
        res$excl_dir[res$nesting_group == sp] <- NA
      }
    }
    
    res
    
  }))


forage_dir_diffs <- guild_excl_dirs %>%
  rename(spp_excl_dir = "excl_dir") %>%
  select(focal_rte, spp_excl_dir, forage_dir) %>%
  unnest(forage_dir) %>%
  mutate(dir_diff = spp_excl_dir - excl_dir) %>%
  filter(!is.na(Foraging)) 
write.csv(forage_dir_diffs, "data/guild_LOO_noHIspp_dir_impact_foraging.csv", row.names = F)

trophic_dir_diffs <- guild_excl_dirs %>%
  rename(spp_excl_dir = "excl_dir") %>%
  select(focal_rte, spp_excl_dir, trophic_dir) %>%
  unnest(trophic_dir) %>%
  mutate(dir_diff = spp_excl_dir - excl_dir) %>%
  filter(!is.na(Trophic.Group)) 
write.csv(trophic_dir_diffs, "data/guild_LOO_noHIspp_dir_impact_trophic.csv", row.names = F)

mig_dir_diffs <- guild_excl_dirs %>%
  rename(spp_excl_dir = "excl_dir") %>%
  select(focal_rte, spp_excl_dir, mig_dir) %>%
  unnest(mig_dir) %>%
  mutate(dir_diff = spp_excl_dir - excl_dir) %>%
  filter(!is.na(migclass)) 
write.csv(mig_dir_diffs, "data/guild_LOO_noHIspp_dir_impact_migclass.csv", row.names = F)

hab_dir_diffs <- guild_excl_dirs %>%
  rename(spp_excl_dir = "excl_dir") %>%
  select(focal_rte, spp_excl_dir, nesting_dir) %>%
  unnest(nesting_dir) %>%
  mutate(dir_diff = spp_excl_dir - excl_dir) %>%
  filter(!is.na(nesting_group)) 
write.csv(hab_dir_diffs, "data/guild_LOO_noHIspp_dir_impact_habitat.csv", row.names = F)

max_diff_guild_map_exclHIspp <- function(df) {
  legend_title <- colnames(df)[[3]]
  max_sf <- df %>%
    group_by(focal_rte) %>%
    filter(dir_diff == max(dir_diff)) %>%
    left_join(routes, by = c("focal_rte" = "stateroute")) %>%
    st_as_sf(coords = c("longitude", "latitude"))
  
  map <- tm_shape(na) + tm_polygons(col = "gray50") + 
    tm_shape(max_sf) + tm_bubbles(col = legend_title, size = "dir_diff", title.size = "Impact on directionality") +
    tm_layout(legend.position=c("right", "bottom"))
  tmap_save(map, paste0("figures/guild_LOO_noHIspp_map_", legend_title, ".pdf"))
}

max_diff_guild_map_exclHIspp(mig_dir_diffs)
max_diff_guild_map_exclHIspp(forage_dir_diffs)
max_diff_guild_map_exclHIspp(hab_dir_diffs)
max_diff_guild_map_exclHIspp(trophic_dir_diffs)

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
  filter(!is.na(max_lc_class)) %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_crop(xmin = -130, ymin = 18, xmax = -57, ymax = 60)

lc_posneg <- scale_model_variables_unnest %>%
  filter(scale == 1) %>%
  filter(!is.na(max_lc_class)) %>%
  mutate(sign = case_when(max_lc > 0 ~ "Increase",
                          TRUE ~ "Decrease"))

theme_set(theme_classic(base_size = 17))
class_change <- ggplot(lc_posneg, aes(x = max_lc_class, fill = sign)) + geom_bar() + labs(x = "", y = "BBS Routes", fill = "") + 
  scale_fill_manual(values = c("Increase" = "#EF8A62", "Decrease" ="#92C5DE")) + theme(legend.position = c(0.9, 0.9)) +
  coord_flip()

lc_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "max_lc_class", size = 0.5, title = "Land cover") + 
  tm_layout(legend.text.size = 1.25, legend.title.size =2, legend.position = c("right", "bottom"), outer.margins = c(0.01,0.01,0.01,0.01))
# tmap_save(lc_map, "figures/max_landcover_map.pdf", units = "in", height = 6, width = 8)

color_scale <- data.frame(color = c(1:4), 
                          temp_hex = c("#92C5DE", "#FDDBC7", "#EF8A62", "#B2182B"), stringsAsFactors = F)

tmin_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "trend_tmin", size = 0.5, title = "Trend in Tmin", breaks = quantile(lc_1route$trend_tmin, na.rm = T), palette = color_scale$temp_hex) + 
  tm_layout(legend.show = F, outer.margins = c(0.01,0.01,0.01,0.01))
  #  tm_layout(legend.text.size = 1.25, legend.title.size =2, legend.position = c("right", "bottom"), outer.margins = c(0.01,0.01,0.01,0.01))

tmax_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(lc_1route) + 
  tm_dots(col = "trend_tmax", size = 0.5, title = "Trend in Tmax", 
          breaks = quantile(lc_1route$trend_tmax, na.rm = T), palette = color_scale$temp_hex) + 
  tm_layout(legend.show = F, outer.margins = c(0.01,0.01,0.01,0.01))
#  tm_layout(legend.text.size = 1.25, legend.title.size =2, legend.position = c("right", "bottom"), outer.margins = c(0.01,0.01,0.01,0.01))

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
  theme(text = element_text(size = 12), axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
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
  scale_fill_manual(values = color_scale$temp_hex) +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  )

vp_tmin <- viewport(0.42, 0.12, width = 0.16, height = 0.15)
vp_tmax <- viewport(0.92, 0.12, width = 0.16, height = 0.15)

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


##### Community weighted means ######

# Long format log abund data

# All species
log_abund_long <- bbs_subset %>%
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
  summarize(mean_abund = mean(speciestotal)) %>%
  dplyr::select(stateroute, aou, year_bin, mean_abund)

# Excluding transient species
log_abund_core_long <- bbs_subset %>%
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
            mean_abund = mean(speciestotal)) %>%
  filter(n_years > 1) %>%
  dplyr::select(stateroute, aou, year_bin, mean_abund)

# Trait data
bbs_aou_temp_range <- read.csv("data/bbs_aou_temp_range.csv", stringsAsFactors = F)
forest_range <- read.csv("data/spp_forest_traits.csv", stringsAsFactors = F)
body_size <- read.csv("data/Master_RO_Correlates_20110610.csv", stringsAsFactors = F) %>%
  dplyr::select(AOU, Mass.g., logMass)

cwm_traits <- bbs_aou_temp_range %>%
  left_join(forest_range) %>%
  left_join(body_size, by = c("aou" = "AOU")) %>%
  na.omit()

# CWM trait correlations

cwm_traits_plot <- cwm_traits %>%
  left_join(fourletter_codes)

cor(cwm_traits_plot[, c(2:4, 9)])

ggplot(cwm_traits_plot, aes(x = temp_mean, y = temp_range)) + geom_text(aes(label = SPEC)) +
  labs(x = "Mean temperature", y = "Temperature range")
ggsave("figures/temp_traits_spp.pdf", units = "in", height = 6, width = 8)

ggplot(cwm_traits_plot, aes(x = propFor, y = for_range)) + geom_text(aes(label = SPEC)) +
  labs(x = "Mean forest", y = "Forest range")
ggsave("figures/hab_traits_spp.pdf", units = "in", height = 6, width = 8)

# Climate trends
env_vars <- read.csv("data/scale_model_input.csv", stringsAsFactors = F)

# Takes a long time ~40 minutes
cwm_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(c(2,3,21))))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = purrr::map(model_input, ~{
    df <- .
    
    cwm_all <- log_abund_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      mutate(wtd_temp_range = mean_abund*temp_range,
             wtd_temp_mean = mean_abund*temp_mean,
             wtd_for_range = mean_abund*for_range,
             wtd_for_mean = mean_abund*propFor,
             wtd_logmass = mean_abund*logMass) %>%
      group_by(year_bin) %>%
      summarize(cwm_temp_range = mean(wtd_temp_range, na.rm = T),
                cwm_temp_mean = mean(wtd_temp_mean, na.rm = T),
                cwm_for_range = mean(wtd_for_range, na.rm = T),
                cwm_for_mean = mean(wtd_for_mean, na.rm = T),
                cwm_logmass = mean(wtd_logmass, na.rm =T)) %>%
      mutate(spp = "all",
             focal_rte = unique(df$focal_rte))
    
    cwm_core <- log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      mutate(wtd_temp_range = mean_abund*temp_range,
             wtd_temp_mean = mean_abund*temp_mean,
             wtd_for_range = mean_abund*for_range,
             wtd_for_mean = mean_abund*propFor,
             wtd_logmass = mean_abund*logMass) %>%
      group_by(year_bin) %>%
      summarize(cwm_temp_range = mean(wtd_temp_range, na.rm = T),
                cwm_temp_mean = mean(wtd_temp_mean, na.rm = T),
                cwm_for_range = mean(wtd_for_range, na.rm = T),
                cwm_for_mean = mean(wtd_for_mean, na.rm = T),
                cwm_logmass = mean(wtd_logmass, na.rm =T)) %>%
      mutate(spp = "no transients",
             focal_rte = unique(df$focal_rte))
    
    rbind(cwm_all, cwm_core)
    
  }))

cwm_unnest <- cwm_input %>%
  dplyr::select(-model_input) %>%
  unnest(cols = c(input_vars))

# write.csv(cwm_unnest, "data/cwm_traits_all_scales.csv", row.names = F)
cwm_unnest <- read.csv("data/cwm_traits_all_scales.csv", stringsAsFactors = F)

## CWM linear models

cwm_temp_plots <- cwm_unnest %>%
  filter(scale == 21, spp == "no transients") %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_mean ~ year_bin, data = df))
  }),
  body_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_logmass ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:body_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

temp_plots <- ggplot(filter(cwm_temp_plots, model %in% c("mean_mod", "range_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Regional") +
  scale_x_discrete(labels = c("mean_mod" = "Temperature mean", "range_mod" = "Temperature range")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)


body_plot <- ggplot(filter(cwm_temp_plots, model %in% c("body_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Regional") +
  scale_x_discrete(labels = c("body_mod" = "Log(body size (g))")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)

cwm_hab_plots <- cwm_unnest %>%
  filter(scale == 2, spp == "no transients") %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_mean ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:mean_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

for_plots <- ggplot(filter(cwm_hab_plots, model %in% c("mean_mod", "range_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Local") +
  scale_x_discrete(labels = c("mean_mod" = "% Forest mean", "range_mod" = "% Forest range")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)

plot_grid(temp_plots, body_plot, for_plots, ncol = 2)
# ggsave("figures/cwms_breadth_position.pdf", units = "in", height = 10, width = 10)

plot_grid(foraging_plot, trophic_plot, mig_plot, temp_plots, body_plot, for_plots, ncol = 3, nrow = 2,
          labels = c("A", "B", "C", "D", "E", "F"), label_size = 17)
ggsave("figures/guild_niche_impacts.pdf", units = "in", height = 10, width = 15)

# Maps of CWMs

cwm_temp_shape <- cwm_temp_plots %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

cwm_hab_shape <- cwm_hab_plots %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

temp_range_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cwm_temp_shape, model == "range_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Temp range")

temp_mean_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cwm_temp_shape, model == "mean_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Temp mean")

body_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cwm_temp_shape, model == "body_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Body size")

regional_maps <- tmap_arrange(temp_mean_map, temp_range_map, body_map, nrow = 2)
tmap_save(regional_maps, "figures/cwms_regional_scale.pdf", units = "in", height = 10, width = 10)

hab_range_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cwm_hab_shape, model == "range_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Forest range")

hab_mean_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cwm_hab_shape, model == "mean_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Forest mean")

local_maps <- tmap_arrange(hab_mean_map, hab_range_map, nrow = 1)
tmap_save(local_maps, "figures/cwms_local_scale.pdf", units = "in", height = 5, width = 10)


## CMs - no abundance weighting

cm_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(c(2,21))))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = purrr::map(model_input, ~{
    df <- .
    
    cwm_core <- log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      group_by(year_bin) %>%
      summarize(cwm_temp_range = mean(temp_range, na.rm = T),
                cwm_temp_mean = mean(temp_mean, na.rm = T),
                cwm_for_range = mean(for_range, na.rm = T),
                cwm_for_mean = mean(propFor, na.rm = T),
                cwm_logmass = mean(logMass, na.rm =T)) %>%
      mutate(spp = "no transients",
             focal_rte = unique(df$focal_rte))
    
    cwm_core
    
  }))

cm_unnest <- cm_input %>%
  dplyr::select(-model_input) %>%
  unnest(cols = c(input_vars))
# write.csv(cm_unnest, "data/community_means_noAbund.csv", row.names = F)


cm_temp_plots <- cm_unnest %>%
  filter(scale == 21, spp == "no transients") %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_mean ~ year_bin, data = df))
  }),
  body_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_logmass ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:body_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

temp_plots <- ggplot(filter(cm_temp_plots, model %in% c("mean_mod", "range_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Regional") +
  scale_x_discrete(labels = c("mean_mod" = "Temperature mean", "range_mod" = "Temperature range")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)


body_plot <- ggplot(filter(cm_temp_plots, model %in% c("body_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Regional") +
  scale_x_discrete(labels = c("body_mod" = "Log(body size (g))")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)

cm_hab_plots <- cm_unnest %>%
  filter(scale == 2, spp == "no transients") %>%
  group_by(focal_rte) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_mean ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:mean_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

for_plots <- ggplot(filter(cm_hab_plots, model %in% c("mean_mod", "range_mod")), aes(x = model, y = estimate)) + 
  geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
  labs(x = "", y = "Change over time", title = "Local") +
  scale_x_discrete(labels = c("mean_mod" = "% Forest mean", "range_mod" = "% Forest range")) +
  geom_hline(yintercept = 0, lty = 2, cex = 1)

plot_grid(temp_plots, body_plot, for_plots, ncol = 2)
# ggsave("figures/cms_noAbund_breadth_position.pdf", units = "in", height = 10, width = 10)

# Maps of CMs - no abundance

cm_temp_shape <- cm_temp_plots %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

cm_hab_shape <- cm_hab_plots %>%
  left_join(routes, by = c("focal_rte" = "stateroute")) %>%
  st_as_sf(coords = c("longitude", "latitude"))

temp_range_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cm_temp_shape, model == "range_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Temp range")

temp_mean_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cm_temp_shape, model == "mean_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Temp mean")

body_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cm_temp_shape, model == "body_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Body size")

regional_maps <- tmap_arrange(temp_mean_map, temp_range_map, body_map, nrow = 2)
tmap_save(regional_maps, "figures/cms_noAbund_regional_scale.pdf", units = "in", height = 10, width = 10)

hab_range_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cm_hab_shape, model == "range_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Forest range")

hab_mean_map <- tm_shape(na) + tm_polygons(col = "gray50") + 
  tm_shape(filter(cm_hab_shape, model == "mean_mod")) + tm_dots(col = "estimate", size = 0.25, title = "Forest mean")

local_maps <- tmap_arrange(hab_mean_map, hab_range_map, nrow = 1)
tmap_save(local_maps, "figures/cms_noAbund_local_scale.pdf", units = "in", height = 5, width = 10)


## Guild CWMs

cwm_foraging_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(c(2,21))))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = purrr::map(model_input, ~{
    df <- .
    
    cwm_core <- log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      left_join(bird_traits, by = c("aou" = "AOU")) %>%
      mutate(wtd_temp_range = mean_abund*temp_range,
             wtd_temp_mean = mean_abund*temp_mean,
             wtd_for_range = mean_abund*for_range,
             wtd_for_mean = mean_abund*propFor,
             wtd_logmass = mean_abund*logMass) %>%
      group_by(year_bin, Foraging) %>%
      summarize(cwm_temp_range = mean(wtd_temp_range, na.rm = T),
                cwm_temp_mean = mean(wtd_temp_mean, na.rm = T),
                cwm_for_range = mean(wtd_for_range, na.rm = T),
                cwm_for_mean = mean(wtd_for_mean, na.rm = T),
                cwm_logmass = mean(wtd_logmass, na.rm =T),
                abund = sum(mean_abund, na.rm = T)) %>%
      mutate(spp = "no transients",
             focal_rte = unique(df$focal_rte))
    
    cwm_core
    
  }))

cwm_foraging_unnest <- cwm_foraging_input %>%
  dplyr::select(-model_input) %>%
  unnest(cols = c(input_vars)) %>%
  mutate(trait_grp = "Foraging") %>%
  rename("trait_val" = "Foraging")

cwm_mig_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(c(2,21))))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = purrr::map(model_input, ~{
    df <- .
    
    cwm_core <- log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      left_join(bird_traits, by = c("aou" = "AOU")) %>%
      mutate(wtd_temp_range = mean_abund*temp_range,
             wtd_temp_mean = mean_abund*temp_mean,
             wtd_for_range = mean_abund*for_range,
             wtd_for_mean = mean_abund*propFor,
             wtd_logmass = mean_abund*logMass) %>%
      group_by(year_bin, migclass) %>%
      summarize(cwm_temp_range = mean(wtd_temp_range, na.rm = T),
                cwm_temp_mean = mean(wtd_temp_mean, na.rm = T),
                cwm_for_range = mean(wtd_for_range, na.rm = T),
                cwm_for_mean = mean(wtd_for_mean, na.rm = T),
                cwm_logmass = mean(wtd_logmass, na.rm =T),
                abund = sum(mean_abund, na.rm = T)) %>%
      mutate(spp = "no transients",
             focal_rte = unique(df$focal_rte))
    
    cwm_core
    
  }))

cwm_mig_unnest <- cwm_mig_input %>%
  dplyr::select(-model_input) %>%
  unnest(cols = c(input_vars)) %>%
  mutate(trait_grp = "migclass") %>%
  rename("trait_val" = "migclass")

cwm_trophic_input <- bbs_subset %>%
  ungroup() %>%
  distinct(stateroute, bcr) %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  dplyr::select(stateroute) %>%
  mutate(scale = purrr::map(stateroute, ~data.frame(scale = rep(c(2,21))))) %>%
  unnest(cols = c(scale)) %>%
  mutate(model_input = map2(stateroute, scale, ~{
    bbs_route_distances %>%
      filter(focal_rte == .x) %>%
      arrange(distance) %>%
      slice(1:.y)
  }),
  input_vars = purrr::map(model_input, ~{
    df <- .
    
    cwm_core <- log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      left_join(bird_traits, by = c("aou" = "AOU")) %>%
      mutate(wtd_temp_range = mean_abund*temp_range,
             wtd_temp_mean = mean_abund*temp_mean,
             wtd_for_range = mean_abund*for_range,
             wtd_for_mean = mean_abund*propFor,
             wtd_logmass = mean_abund*logMass) %>%
      group_by(year_bin, Trophic.Group) %>%
      summarize(cwm_temp_range = mean(wtd_temp_range, na.rm = T),
                cwm_temp_mean = mean(wtd_temp_mean, na.rm = T),
                cwm_for_range = mean(wtd_for_range, na.rm = T),
                cwm_for_mean = mean(wtd_for_mean, na.rm = T),
                cwm_logmass = mean(wtd_logmass, na.rm =T),
                abund = sum(mean_abund, na.rm = T)) %>%
      mutate(spp = "no transients",
             focal_rte = unique(df$focal_rte))
    
    cwm_core
    
  }))

cwm_trophic_unnest <- cwm_trophic_input %>%
  dplyr::select(-model_input) %>%
  unnest(cols = c(input_vars)) %>%
  mutate(trait_grp = "Trophic.Group") %>%
  rename("trait_val" = "Trophic.Group")

cwm_guild_unnest <- bind_rows(cwm_foraging_unnest, cwm_mig_unnest, cwm_trophic_unnest)
write.csv(cwm_guild_unnest, "data/cwm_guild_data_unnest.csv", row.names = F)

cwm_guild_temp_plots <- cwm_guild_unnest %>%
  filter(scale == 21, spp == "no transients") %>%
  group_by(focal_rte, trait_grp, trait_val) %>%
  filter(!is.na(trait_val)) %>%
  nest() %>%
  mutate(n_bins = purrr::map_dbl(data, ~{nrow(.)})) %>%
  filter(n_bins == 5 | n_bins == 6) %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_temp_mean ~ year_bin, data = df))
  }),
  body_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_logmass ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:body_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

guild_plot_labels <- data.frame(mod = c("mean_mod", "range_mod", "body_mod"), label = c("Mean temp", "Temp range", "Body size"))

forage_n <- bird_traits %>%
  group_by(Foraging) %>%
  count()

trophic_n <- bird_traits %>%
  group_by(Trophic.Group) %>%
  count()

mig_n <- bird_traits %>%
  group_by(migclass) %>%
  count()

for(i in 1:3) {
  mod <- guild_plot_labels$mod[i]
  lab <- guild_plot_labels$label[i]
  
  forage_plot <- ggplot(filter(cwm_guild_temp_plots, model == mod, trait_grp == "Foraging"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = forage_n$Foraging, y = max(filter(cwm_guild_temp_plots, model == mod, trait_grp == "Foraging")$estimate) - 1, label = forage_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  trophic_plot <- ggplot(filter(cwm_guild_temp_plots, model == mod, trait_grp == "Trophic.Group"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = trophic_n$Trophic.Group, y = max(filter(cwm_guild_temp_plots, model == mod, trait_grp == "Trophic.Group")$estimate) - 1, label = trophic_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  mig_plot <- ggplot(filter(cwm_guild_temp_plots, model == mod, trait_grp == "migclass"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = mig_n$migclass, y = max(filter(cwm_guild_temp_plots, model == mod, trait_grp == "migclass")$estimate) - 1, label = mig_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  plot_grid(forage_plot, trophic_plot, mig_plot, nrow = 1)
  ggsave(paste0('figures/cwm_guild_plot', mod, ".pdf"), units = "in", height = 5, width = 15)
}

cwm_guild_hab_plots <- cwm_guild_unnest %>%
  filter(scale == 2, spp == "no transients") %>%
  group_by(focal_rte, trait_grp, trait_val) %>%
  filter(!is.na(trait_val)) %>%
  nest() %>%
  mutate(n_bins = purrr::map_dbl(data, ~{nrow(.)})) %>%
  filter(n_bins == 5 | n_bins == 6) %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cwm_for_mean ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:mean_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

guild_hab_plot_labels <- data.frame(mod = c("mean_mod", "range_mod"), label = c("Mean forest", "Forest range"))

for(i in 1:2) {
  mod <- guild_hab_plot_labels$mod[i]
  lab <- guild_hab_plot_labels$label[i]
  
  forage_plot <- ggplot(filter(cwm_guild_hab_plots, model == mod, trait_grp == "Foraging"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = forage_n$Foraging, y = max(filter(cwm_guild_hab_plots, model == mod, trait_grp == "Foraging")$estimate) - 1, label = forage_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  trophic_plot <- ggplot(filter(cwm_guild_hab_plots, model == mod, trait_grp == "Trophic.Group"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = trophic_n$Trophic.Group, y = max(filter(cwm_guild_hab_plots, model == mod, trait_grp == "Trophic.Group")$estimate) - 1, label = trophic_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  mig_plot <- ggplot(filter(cwm_guild_hab_plots, model == mod, trait_grp == "migclass"), aes(x = trait_val, y = estimate)) +
    geom_violin(trim = T, draw_quantiles = c(0.5), cex = 1, fill = "gray") +
    annotate(geom = "text", x = mig_n$migclass, y = max(filter(cwm_guild_hab_plots, model == mod, trait_grp == "migclass")$estimate) - 1, label = mig_n$n) +
    labs(y = "Change in CWM", x = "", title = lab) + coord_flip() +
    geom_hline(yintercept = 0, lty = 2, cex = 1)
  
  plot_grid(forage_plot, trophic_plot, mig_plot, nrow = 1)
  ggsave(paste0('figures/cwm_guild_plot_hab_', mod, ".pdf"), units = "in", height = 5, width = 15)
}