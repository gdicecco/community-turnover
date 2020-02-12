### Community trajectory analysis
## Compare 1970-2016 to 1990-2016

library(tidyverse)
library(purrr)
library(spdep)
library(tmap)
library(sf)
library(vegclust)

## Read in data

rank_abund <- read.csv("data/bbs_subset_1970-2016_ranks.csv", stringsAsFactors = F)

## Community trajectories

# https://cran.r-project.org/web/packages/vegclust/vignettes/CTA.html