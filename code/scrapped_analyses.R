### Scrapped analyses

#### Change in niche specialization over time ####
### For each species, what is temp range and mean, forest range and mean of occurrences in each time window
### Slope of these four variables with time

temp_range <- function(clims) {
  quants <- quantile(clims, c(0.05, 0.95), na.rm = T)
  quants[[2]] - quants[[1]]
}

bbs_mean_temp <- read.csv("data/bbs_routes_breeding_season_climate.csv", stringsAsFactors = F)
bbs_mean_temp$mean_temp <- rowMeans(bbs_mean_temp[,2:3])

bbs_yearbin_temp <- bbs_mean_temp %>%
  left_join(routes) %>%
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
  unnest(cols = c("year_bins")) %>%
  group_by(countrynum, stateroute, year_bin) %>%
  summarize(meanTemp = mean(mean_temp, na.rm = T))

spp_temp_niche <- log_abund_core_long %>%
  left_join(bbs_yearbin_temp, by = c("stateroute", "year_bin")) %>%
  group_by(year_bin, aou) %>%
  summarize(temp_mean = mean(meanTemp, na.rm = T),
            temp_range = temp_range(meanTemp),
            nRoutes = n_distinct(stateroute)) %>%
  filter(nRoutes >= 10) %>% # 228 species
  group_by(aou) %>%
  nest() %>%
  mutate(mean_mod = purrr::map(data, ~{
    df <- .
    lm(temp_mean ~ year_bin, data = df)
  }),
  range_mod = purrr::map(data, ~{
    df <- .
    lm(temp_range ~ year_bin, data = df)
  }),
  mean_tidy = purrr::map(mean_mod, ~tidy(.)),
  range_tidy = purrr::map(range_mod, ~tidy(.)),
  n_years = purrr::map_dbl(data, ~nrow(.)),
  delta_mean = purrr::map_dbl(data, ~{
    df <- .
    df$temp_mean[df$year_bin == max(df$year_bin)] - df$temp_mean[df$year_bin == min(df$year_bin)]
  }),
  delta_range = purrr::map_dbl(data, ~{
    df <- .
    df$temp_range[df$year_bin == max(df$year_bin)] - df$temp_range[df$year_bin == min(df$year_bin)]
  })) %>%
  dplyr::select(-data, -mean_mod, -range_mod) %>%
  pivot_longer(mean_tidy:range_tidy, names_to = "model", values_to = "tidy") %>%
  unnest(cols = c("tidy"))

# Fig: violin plot, delta ranges and delta mean temps

delta_plot <- spp_temp_niche %>%
  dplyr::select(aou, delta_range, delta_mean) %>%
  ungroup() %>%
  distinct() %>%
  pivot_longer(delta_range:delta_mean, names_to = "niche_prop", values_to = "delta")

ggplot(delta_plot, aes(x = niche_prop, y = delta)) + geom_violin(draw_quantiles = c(0.5), fill = "gray") +
  geom_hline(yintercept = 0, lty = 2, cex = 1) +
  labs(x = "", y = "Change betw first and last time period (deg. C)") +
  scale_x_discrete(labels = c("delta_mean" = "Mean temperature", "delta_range" = "Temperature range"))
ggsave("figures/temperature_niche_changes.pdf")

# Fig: violin plot, temp range and mean slopes

slope_plot <- spp_temp_niche %>%
  filter(term == "year_bin")

ggplot(slope_plot, aes(x = model, y = estimate)) + geom_violin(draw_quantiles = c(0.5), fill = "gray") +
  geom_hline(yintercept = 0, lty = 2, cex = 1) +
  labs(x = "", y = "Slope of temperature niche change (deg. C/year)") +
  scale_x_discrete(labels = c("mean_tidy" = "Mean temperature", "range_tidy" = "Temperature range"))
ggsave("figures/temperature_niche_slopes.pdf")

trait_vs_change <- spp_temp_niche %>%
  dplyr::select(aou, delta_range, delta_mean) %>%
  ungroup() %>%
  distinct() %>%
  left_join(cwm_traits)

temp_range <- ggplot(trait_vs_change, aes(x = temp_range, y = delta_range)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Temperature range", y = "Change in temperature range")
temp_mean <- ggplot(trait_vs_change, aes(x = temp_mean, y = delta_mean)) + geom_point() + geom_smooth(method = "lm", se = F) +
  labs(x = "Mean temperature", y = "Change in mean temperature")

plot_grid(temp_mean, temp_range, ncol = 2)
ggsave("figures/niche_vs_niche_shifts.pdf", units = "in", height = 6, width = 10)
