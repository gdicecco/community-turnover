### Species occurrence analysis: colonizations and extinctions

library(tidyverse)
library(purrr)
library(broom)

## Plotting theme

theme_set(theme_classic(base_size = 15))

## Read in data

log_abund <- read.csv("data/bbs_subset_1970-2016_logabund.csv", stringsAsFactors = F)
routes_subs <- read.csv("data/bbs_route_subset_1990-2016.csv", stringsAsFactors = F) # only routes with at least 3 surveys per 5 year window
species_list <- read.csv("data/species_list.csv", stringsAsFactors = F)

log_abund_subs <- log_abund %>%
  left_join(species_list) %>%
  filter(!grepl("unid", english_common_name)) %>%
  filter(stateroute %in% routes_subs$stateroute) %>%
  mutate(y1 = case_when(countrynum == 124 ~ 1990,
                        countrynum == 840 ~ 1992),
         y2 = case_when(countrynum == 124 ~ 2010,
                        countrynum == 840 ~ 2016)) %>%
  filter(year >= y1, year <= y2)

## Species occurrence over time - split time series in half to measure colonizations and extinctions

spp_occ <- log_abund_subs %>%
  mutate(decade = case_when(countrynum == 124 & year >= 1990 & year <= 1999 ~ "d1",
                        countrynum == 840 & year >= 1992 & year <= 2001 ~ "d1",
                        countrynum == 124 & year >= 2001 & year <= 2010 ~ "d2",
                        countrynum == 840 & year >= 2007 & year <= 2016 ~ "d2",
                        TRUE ~ "NA")) %>%
  filter(decade == "d1" | decade == "d2") %>%
  group_by(aou, english_common_name, stateroute, decade) %>%
  summarize(n_years = n_distinct(year),
            occ = n_years/10)

spp_occ_deltas <- spp_occ %>%
  select(aou, english_common_name, stateroute, decade, occ) %>%
  pivot_wider(names_from = decade, values_from = occ) %>%
  replace_na(list(d1 = 0, d2 = 0)) %>%
  mutate(deltaOcc = d2-d1)

ggplot(spp_occ_deltas, aes(x = deltaOcc)) + geom_histogram(bins = 20, col = "white") + 
  geom_vline(aes(xintercept = mean(deltaOcc)), col = "red", cex = 1) +
  labs(y = "Count", x = "Change in occupancy")
ggsave("figures/hist_change_in_occ.pdf")

spp_colext <- spp_occ_deltas %>%
  group_by(aou, english_common_name) %>%
  nest() %>%
  mutate(rte_col = map_dbl(data, ~{
    df <- .
    n_distinct(df$stateroute[df$d1 == 0 & df$d2 > 0])/n_distinct(df$stateroute[df$d2 > 0])
  }),
  rte_ext = map_dbl(data, ~{
    df <- .
    n_distinct(df$stateroute[df$d1 > 0 & df$d2 == 0])/n_distinct(df$stateroute[df$d1 > 0])
  }),
  rte_dTrans = map_dbl(data, ~{
    df <- .
    trans_d1 <- n_distinct(df$stateroute[df$d1 <= 0.3 & df$d1 > 0])/n_distinct(df$stateroute[df$d1 > 0])
    trans_d2 <- n_distinct(df$stateroute[df$d2 <= 0.3 & df$d2 > 0])/n_distinct(df$stateroute[df$d2 > 0])
    trans_d2-trans_d1
  }),
  n_routes = map_dbl(data, nrow),
  colext = rte_col/rte_ext) %>%
  filter(n_routes > 20)

ggplot(spp_colext, aes(x = rte_col)) + geom_histogram(bins = 20, col = "white") +
  geom_vline(aes(xintercept = mean(rte_col)), col = "red", cex = 1) +
  labs(y = "Count", x = "Percent routes colonized")
ggsave("figures/hist_rtes_col.pdf")

ggplot(spp_colext, aes(x = rte_ext)) + geom_histogram(bins = 20, col = "white") +
  geom_vline(aes(xintercept = mean(rte_ext)), col = "red", cex = 1) +
  labs(y = "Count", x = "Percent routes extinct")
ggsave("figures/hist_rtes_ext.pdf")

ggplot(spp_colext, aes(x = rte_dTrans)) + geom_histogram(bins = 20, col = "white") +
  geom_vline(aes(xintercept = mean(rte_dTrans)), col = "red", cex = 1) +
  labs(y = "Count", x = "Change in percent routes transient")
ggsave("figures/hist_dTransient.pdf")

ggplot(spp_colext, aes(x = colext)) + geom_histogram(bins = 20, col = "white") + 
  geom_vline(aes(xintercept = mean(colext)), col = "red", cex = 1) +
  labs(x = "Colonization rate/extinction rate", y = "Count")
ggsave("figures/hist_ratio_col_ext.pdf")

## Species occupancy by year bin

spp_occ_year_bins <- log_abund_subs %>%
  group_by(countrynum) %>%
  nest() %>%
  mutate(year_bins = map2(countrynum, data, ~{
    country <- .x
    df <- .y
    
    if(country == 124) {
      df %>%
        mutate(year_bin = case_when(year >= 1990 & year <= 1994 ~ 1990,
                                    year >= 1995 & year <= 1999 ~ 1995,
                                    year >= 2000 & year <= 2004 ~ 2000,
                                    TRUE ~ 2005))
    } else {
      df %>%
        mutate(year_bin = case_when(year >= 1992 & year <= 1996 ~ 1992,
                                    year >= 1997 & year <= 2001 ~ 1997,
                                    year >= 2002 & year <= 2006 ~ 2002,
                                    year >= 2007 & year <= 2011 ~ 2007,
                                    year >= 2012 & year <= 2016 ~ 2012))
    }
  })) %>%
  select(-data) %>%
  unnest(cols = c(year_bins)) %>%
  group_by(aou, english_common_name, stateroute, year_bin) %>%
  summarize(years_obs = n_distinct(year)) %>%
  left_join(routes_subs, by = c("stateroute", "year_bin")) %>%
  mutate(occ = years_obs/n_years) %>%
  group_by(max_bins, aou, english_common_name, stateroute) %>%
  nest() %>%
  mutate(n_points = map_dbl(data, ~nrow(.))) %>%
  filter(n_points == max_bins) %>%
  mutate(occ_mod = map(data, ~glm(occ~year_bin, data = ., family = "binomial", weights = n_years)),
         tidy_mod = map(occ_mod, ~tidy(.))) %>%
  select(-data, -occ_mod) %>%
  unnest(cols = c("tidy_mod"))

# write.csv(spp_occ_year_bins, "data/occupancy_trend_mods.csv", row.names = F)

# Vis model output

spp_occ_year_bins <- read.csv("data/occupancy_trend_mods.csv", stringsAsFactors = F)

spp_occ_trends <- spp_occ_year_bins %>%
  group_by(aou, english_common_name) %>%
  nest() %>%
  mutate(mean_occ_trend = map_dbl(data, ~{
    df <- .
    trend <- df %>%
      filter(term == "year_bin")
    
    mean(trend$estimate, na.rm = T)
  }),
  sd_occ_trend = map_dbl(data, ~{
    df <- .
    trend <- df %>%
      filter(term == "year_bin")
    
    sd(trend$estimate, na.rm = T)
  }),
  routes = map_dbl(data, ~nrow(.))) %>%
  filter(routes > 50) %>%
  select(-data)
# write.csv(spp_occ_trends, "data/spp_occ_trends.csv", row.names = F)

ggplot(spp_occ_trends, aes(x = mean_occ_trend)) + geom_histogram(col = "white") +
  geom_vline(aes(xintercept = mean(mean_occ_trend)), col = "blue", cex = 1) +
  labs(x = "Mean trend in occupancy", y = "Count")
ggsave("figures/mean_occupancy_trend_hist.pdf")
