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

rank_abund <- read.csv("data/bbs_subset_1970-2016_ranks.csv", stringsAsFactors = F)

## Community trajectories
# https://cran.r-project.org/web/packages/vegclust/vignettes/CTA.html

# long to wide

rank_wide <- rank_abund %>%
  pivot_wider(names_from = aou, values_from = rank, values_fn = list(rank = mean)) %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(dir50 = purrr::map_dbl(data, ~{
    df <- .
    rank_dist <- dist(df)
    trajectoryDirectionality(rank_dist, sites = rep(1, nrow(df)), surveys = df$year)
  })) %>%
  mutate(dir25 = purrr::map_dbl(data, ~{
    df <- . 
    df <- df %>%
      filter(year >= 1990)
    if(max(df$year) - min(df$year) > 15) {
      rank_dist <- dist(df)
      trajectoryDirectionality(rank_dist, sites = rep(1, nrow(df)), surveys = df$year)
    } else(NA)
  }))

r <- round(cor(rank_wide$dir50, rank_wide$dir25), 2)

ggplot(rank_wide, aes(x = dir50, y = dir25)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Directionality 1970-2016", y = "Directionality 1990-2016") + 
  annotate(geom= "text", x = 0.7, y = 0.35, label = paste0("r = ", r), size = 8)
ggsave("figures/directionality_time_series_comparison.pdf")

