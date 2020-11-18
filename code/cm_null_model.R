### CM null models

### For each BBS route, how many species have changed between first and last time period
### What is the species pool to pick from? BCR list?
### Should we just consider net species or total loss and total gain

cm_null_input <- bbs_subset %>%
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
    
    log_abund_core_long %>%
      filter(stateroute %in% df$stateroute) %>%
      group_by(year_bin, aou) %>%
      summarize_all(mean, na.rm = T) %>%
      left_join(cwm_traits) %>%
      group_by(year_bin)
    
  }))

cm_null_input_unnest <- cm_null_input %>%
  dplyr::select(-model_input) %>%
  unnest() 

spp_changes <- cm_null_input %>%
  dplyr::select(-model_input) %>%
  unnest() %>%
  group_by(stateroute, scale) %>%
  mutate(n_bins = n_distinct(year_bin)) %>%
  filter(n_bins == 6) %>%
  nest() %>%
  mutate(turnover = map(data, ~{
    df <- .
    
    years <- unique(df$year_bin)
    
    start <- length(unique(df$aou[df$year_bin == min(df$year_bin, na.rm = T)]))
    
    transitions <- data.frame(transition = c(1:5), yr_start = c(rep(NA,5)),
                              yr_end = c(rep(NA, 5)), 
                              gain = c(rep(NA, 5)), loss = c(rep(NA, 5)))
    
    for(i in 1:5) {
      y1 <- years[i]
      y2 <- years[i + 1]
      
      list1 <- unique(df$aou[df$year_bin == y1])
      list2 <- unique(df$aou[df$year_bin == y2])
      
      g <- length(list2[!(list2 %in% list1)])
      l <- length(list1[!(list1 %in% list2)])
      
      transitions[i, 2] <- y1
      transitions[i, 3] <- y2
      transitions[i, 4] <- g
      transitions[i, 5] <- l
    }
    
    transitions$start_n <- start
    
    transitions
    
  }))


cm_transitions <- spp_changes %>%
  dplyr::select(-data) %>%
  unnest(cols = c("turnover"))

### Species list for each BCR
bcr_spp_list <- bbs_subset %>%
  ungroup() %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  group_by(bcr) %>%
  dplyr::select(bcr, aou) %>%
  distinct(aou) %>%
  left_join(cwm_traits)

cm_null_means <- cm_transitions %>%
  left_join(dplyr::select(routes, stateroute, bcr)) %>%
  group_by(stateroute, scale, bcr) %>%
  nest() %>%
  mutate(null_cm = map2(bcr, data, ~{
    b <- .x
    df <- .y
    
    start <- unique(df$start_n)
    
    start_y <- min(df$yr_start)
    
    res <- data.frame(sim = c(), year_bin = c(), cm_temp_mean = c(), 
                      cm_temp_range = c(),
                      cm_for_mean = c(),
                      cm_for_range = c(),
                      cm_bodysize = c())
    
    for(i in 1:999) {
      samp <- bcr_spp_list %>%
        filter(bcr == b) %>%
        sample_n(start)
      
      res <- rbind(res, data.frame(sim = i, year_bin = start_y,
                              cm_temp_mean = mean(samp$temp_mean, na.rm = T),
                              cm_temp_range = mean(samp$temp_range, na.rm = T),
                              cm_for_mean = mean(samp$propFor, na.rm = T),
                              cm_for_range = mean(samp$for_range, na.rm = T),
                              cm_bodysize = mean(samp$logMass, na.rm = T)))
      
      for(j in 1:5) {
        g <- df$gain[j]
        l <- df$loss[j]
        
        yr <- df$yr_end[j]
        
        g_spp <- bcr_spp_list %>%
          filter(bcr == b) %>%
          filter(!(aou %in% samp$aou)) %>%
          sample_n(g)
        
        l_spp <- samp %>%
          sample_n(l)
        
        samp <- samp %>%
          bind_rows(g_spp) %>%
          filter(!(aou %in% l_spp$aou))
        
        cm <- data.frame(sim = i, year_bin = yr,
                         cm_temp_mean = mean(samp$temp_mean, na.rm = T),
                         cm_temp_range = mean(samp$temp_range, na.rm = T),
                         cm_for_mean = mean(samp$propFor, na.rm = T),
                         cm_for_range = mean(samp$for_range, na.rm = T),
                         cm_bodysize = mean(samp$logMass, na.rm = T))
        
        res <- rbind(res, cm)
        
      }
      
    }
    
    res
    
  }))

local_null <- cm_null_means %>%
  filter(scale == 2) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("null_cm"))
write.csv(local_null, "bigdata/community_mean_null_mod_res_local_scale.csv", row.names = F)

regional_null <- cm_null_means %>%
  filter(scale == 21) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("null_cm"))
write.csv(regional_null, "bigdata/community_mean_null_mod_res_regional_scale.csv", row.names = F)


#### Figures ####

library(tidyverse)
library(purrr)
library(broom)

theme_set(theme_classic(base_size = 15))

# Null model output
local_null <- read.csv("bigdata/community_mean_null_mod_res_local_scale.csv", stringsAsFactors = F)
regional_null <- read.csv("bigdata/community_mean_null_mod_res_regional_scale.csv", stringsAsFactors = F)

# CM empirical
cm_unnest <- read.csv("data/community_means_noAbund.csv", stringsAsFactors = F)

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

# CM slope z scores

cm_local_z <- local_null %>%
  group_by(stateroute, sim) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cm_for_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cm_for_mean ~ year_bin, data = df))
  })) 

cm_local_z_mod <- cm_local_z %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:mean_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

cm_local_z_scores <- cm_local_z_mod %>%
  group_by(stateroute, model) %>%
  summarize(mean_est = mean(estimate),
            sd_est = sd(estimate)) %>%
  left_join(cm_hab_plots, by = c("stateroute" = "focal_rte", "model")) %>%
  mutate(z_est = (estimate - mean_est)/sd_est)

ggplot(cm_local_z_scores, aes(x = model, y = z_est)) + 
  geom_violin(draw_quantiles = c(0.5), fill = "gray") + geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Proportion forest cover", y = "Z-Community mean slope")
ggsave("figures/community_mean_habitat_z-score.pdf")

Sys.time()
cm_regional_z <- regional_null %>%
  group_by(stateroute, sim) %>%
  nest() %>%
  mutate(range_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cm_temp_range ~ year_bin, data = df))
  }),
  mean_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cm_temp_mean ~ year_bin, data = df))
  }),
  body_mod = purrr::map(data, ~{
    df <- .
    tidy(lm(cm_bodysize ~ year_bin, data = df))
  })) %>%
  dplyr::select(-data) %>%
  pivot_longer(range_mod:body_mod, names_to = "model", values_to = "table") %>%
  unnest(cols = c("table")) %>%
  filter(term == "year_bin")

cm_regional_z_scores <- cm_regional_z %>%
  group_by(stateroute, model) %>%
  summarize(mean_est = mean(estimate),
            sd_est = sd(estimate)) %>%
  left_join(cm_temp_plots, by = c("stateroute" = "focal_rte", "model")) %>%
  mutate(z_est = (estimate - mean_est)/sd_est)

ggplot(cm_regional_z_scores, aes(x = model, y = z_est)) + 
  geom_violin(draw_quantiles = c(0.5), fill = "gray") + geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Temperature niche/body size traits", y = "Z-Community mean slope")
ggsave("figures/community_mean_temperature_z-score.pdf")


