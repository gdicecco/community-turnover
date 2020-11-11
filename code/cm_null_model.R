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
      group_by(year_bin) %>%
      summarize(n_spp = n_distinct(aou))
    
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
  summarize(start_spp = n_spp[year_bin == min(year_bin, na.rm = T)],
            end_spp = n_spp[year_bin == max(year_bin, na.rm = T)],
    delta_spp = n_spp[year_bin == max(year_bin, na.rm = T)] - n_spp[year_bin == min(year_bin, na.rm = T)])

local_change <- ggplot(filter(spp_changes, scale == 2), aes(x = delta_spp)) + geom_histogram() +
  labs(x = "Net change in species", main = "Scale = 2 routes")

regional_change <- ggplot(filter(spp_changes, scale == 21), aes(x = delta_spp)) + geom_histogram() +
  labs(x = "Net change in species", main = "Scale = 21 routes")

plot_grid(local_change, regional_change, nrow = 1)
ggsave("figures/cm_change_in_spp.pdf", units = "in", width = 10, height = 5)

### Species list for each BCR
bcr_spp_list <- bbs_subset %>%
  ungroup() %>%
  filter(bcr %in% bcr_subset$bcr) %>%
  group_by(bcr) %>%
  dplyr::select(bcr, aou) %>%
  distinct(aou) %>%
  left_join(cwm_traits)

cm_null_means <- cm_null_input_unnest %>%
  left_join(routes) %>%
  mutate(null_cm = map2(bcr, n_spp, ~{
    b <- .x
    n <- .y
    
    res <- data.frame(sim = c(), cm_temp_mean = c(), 
                      cm_temp_range = c(),
                      cm_for_mean = c(),
                      cm_for_range = c(),
                      cm_bodysize = c())
    
    for(i in 1:999) {
      samp <- bcr_spp_list %>%
        filter(bcr == b) %>%
        sample_n(n)
      
      cm <- data.frame(sim = i, cm_temp_mean = mean(samp$temp_mean),
                       cm_temp_range = mean(samp$temp_range),
                       cm_for_mean = mean(samp$propFor),
                       cm_for_range = mean(samp$for_range),
                       cm_bodysize = mean(samp$logMass))
      
      res <- rbind(res, cm)
      
    }
    
    res
    
  }))
