### Trajectory directionality simulations
### How sensitive is this metric to changes

library(tidyverse)
library(purrr)
library(vegclust)

theme_set(theme_classic(base_size = 15))

### Sample community of 30 species

set.seed(10)
community <- data.frame(species = paste0("s", c(1:30)), t1 = round(exp(rlnorm(30, meanlog = 0.75, sd=0.5))))

### Large vs. small population trends
### Sample community at 5 time points
### Slopes 0-10 by 0.5 indiv/time point
### Test from 2-10 species with trends (half positive, half negative)

slopes <- seq(0,10, by= 0.5)
n_spec <- c(1:5)

poptrend_results <- data.frame(n_spec = c(), slope = c(), dir = c())
for(i in n_spec) {
  
  pos_trends <- sample_n(community, i)
  neg_trends <- sample_n(filter(community, !(species %in% pos_trends$species)), i)
  
  community_change <- community %>%
    mutate(trend = case_when(species %in% pos_trends$species ~ "pos",
                             species %in% neg_trends$species ~ "neg",
                             TRUE ~ "none"))
  
  for(j in slopes) {
    
    community_timeseries <- community_change %>%
      mutate(t2 = case_when(trend == "pos" ~ t1 + j*t1,
                            trend == "neg" ~ ifelse(t1 - j*t1 >= 0, t1 - j*t1, 0),
                            trend == "none" ~ t1),
             t3 = case_when(trend == "pos" ~ t2 + j*t2,
                            trend == "neg" ~ ifelse(t2 - j*t2 >= 0, t2 - j*t2, 0),
                            trend == "none" ~ t2),
             t4 = case_when(trend == "pos" ~ t3 + j*t3,
                            trend == "neg" ~ ifelse(t3 - j*t3 >= 0, t3 - j*t3, 0),
                            trend == "none" ~ t3),
             t5 = case_when(trend == "pos" ~ t4 + j*t4,
                            trend == "neg" ~ ifelse(t4 - j*t4 >= 0, t4 - j*t4, 0),
                            trend == "none" ~ t4)) %>%
      mutate_at(c('t1', "t2", "t3", "t4", "t5"), ~ifelse(.<1, 1, .)) %>%
      mutate_at(c("t1", "t2", "t3", "t4", "t5"), ~log(round(.)))
    
    community_names <- community_timeseries$species
    
    community_long <- as.data.frame(t(community_timeseries[,-c(1,3)]))
    colnames(community_long) <- community_names
    community_long$time <- row.names(community_long)

    dist_mat <- dist(community_long)  
    dir <- trajectoryDirectionality(dist_mat, sites = rep(1, nrow(community_long)), surveys = community_long$time)
    
    if(j > 0) {
      trajectoryPCoA(dist_mat, sites = rep(1, nrow(community_long)), surveys = community_long$time)
    }
    
    poptrend_results <- bind_rows(poptrend_results, data.frame(n_spec = i, slope = j, dir = dir))
  
  }
}

poptrend_plot <- poptrend_results %>%
  replace_na(list(dir = 0))

ggplot(poptrend_plot, aes(x = (2*n_spec)/30, y = slope, fill = dir)) + geom_tile() + scale_fill_viridis_c() +
  labs(x = "Proportion species with directional trends", y = "Population trend", fill = "Directionality")
ggsave("figures/directionality_sim_poptrends.pdf")


### Sample community at 10 time points
### Extinctions in t2, colonizations in t4
### Test from 4-20 species C&E (half colonize, half extinct)

n_spec <- c(1:5)

colext_results <- data.frame(n_spec = c(), dir = c())

for(i in n_spec) {
  col1_spec <- data.frame(species = paste0("s", c(31:35)), t1 = 0) %>%
    sample_n(i)
  col2_spec <- data.frame(species = paste0("s", c(36:40)), t1 = 0) %>%
    sample_n(i)
  
  ext1_spec <- community %>%
    sample_n(i)
  ext2_spec <- community %>%
    filter(!(species %in% ext1_spec$species)) %>%
    sample_n(i)
  
  community_timeseries <- community %>%
    bind_rows(col1_spec, col2_spec) %>%
    mutate(t2 = case_when(species %in% ext1_spec$species ~ 0,
                          species %in% col1_spec$species ~ 0,
                          TRUE ~ t1),
           t3 = case_when(species %in% ext2_spec$species ~ 0,
                          species %in% col2_spec$species ~ 0,
                          TRUE ~ t2),
           t4 = case_when(species %in% col1_spec$species ~ round(exp(rlnorm(1, meanlog = 0.75, sd=0.5))),
                          TRUE ~ t3),
           t5 = case_when(species %in% col2_spec$species ~ round(exp(rlnorm(1, meanlog = 0.75, sd=0.5))),
                          TRUE ~ t4)) %>%
    mutate_at(c('t1', "t2", "t3", "t4", "t5"), ~ifelse(.<1, 1, .)) %>%
    mutate_at(c("t1", "t2", "t3", "t4", "t5"), ~log(.))
  
  community_names <- community_timeseries$species
  
  community_long <- as.data.frame(t(community_timeseries[,-c(1)]))
  colnames(community_long) <- community_names
  community_long$time <- row.names(community_long)
  
  dist_mat <- dist(community_long)  
  dir <- trajectoryDirectionality(dist_mat, sites = rep(1, nrow(community_long)), surveys = community_long$time)
  trajectoryPCoA(dist_mat, sites = rep(1, nrow(community_long)), surveys = community_long$time)
  
  colext_results <- bind_rows(colext_results, data.frame(n_spec = i, dir = dir))
  
}

colext_results
## Always 0.5 if col/ext are balanced

