### Create 1/2 route path shapefiles for land cover analysis

library(tidyverse)
library(sf)
library(lwgeom)
library(rgdal)

### Canada and US BBS route paths with 5-km buffers

us_bbs <- read_sf("\\\\BioArk.bio.unc.edu/HurlbertLab/Databases/BBS/GPS_stoplocations/bbsrte_2012_alb/bbsrte_2012_alb.shp")

# Convert multilinestring to linestring

bbs_route_lines <- us_bbs %>%
  group_by(rteno) %>%
  summarize(sum_length = sum(ARC_LENGTH),
            rte_length = mean(rte_length))

us_bbs_contig <- bbs_route_lines %>%
  st_cast("LINESTRING") %>%
  group_by(rteno) %>%
  count() %>%
  filter(n == 1)

us_1to25 <- st_linesubstring(us_bbs_contig, from = 0, to = 0.5) %>%
  mutate(stops = "1-25")
us_26to50 <- st_linesubstring(us_bbs_contig, from = 0.5, to = 1) %>%
  mutate(stops = "26-50")

us_bbs_lines <- rbind(us_1to25, us_26to50)

us_bbs_buffers <- st_buffer(us_bbs_lines, 5000)
# write_sf(us_bbs_buffers, "C:/Users/gdicecco/Desktop/data/us_bbs_half_route_paths_5km.shp")


# Use Canada route paths from 2014 and 2019
# CA lines can't be converted to LINESTRINGS like US can - start and stop GPS where there are gaps in the route

ca_bbs <- readOGR(dsn = "\\\\BioArk.bio.unc.edu/HurlbertLab/Databases/BBS/GPS_stoplocations/BBS_routes_Canada_dec2019.gdb",
                  layer = "CURRENT_STOPS")

ca_bbs_2014 <- readOGR(dsn = "\\\\BioArk.bio.unc.edu/HurlbertLab/Databases/BBS/GPS_stoplocations/BBS_Routes_NAD83_15042014.gdb",
                       layer = "CURRENT_STOPS")

ca_bbs_sf <- st_as_sf(ca_bbs) %>%
  dplyr::select(-ProvRoute_, -Nbr_FullNa)
ca_bbs_2014_sf <- st_as_sf(ca_bbs_2014) %>%
  rename("Province_R" = Province_Route) %>%
  mutate(Year = 0)

ca_bbs_lines <- rbind(ca_bbs_sf, filter(ca_bbs_2014_sf, !(Province_R %in% ca_bbs_sf$Province_R))) %>%
  group_by(Province_R) %>%
  filter(n_distinct(Stop) == 50) %>%
  mutate(stops = case_when(Stop < 26 ~ "1-25",
                           Stop > 25 ~ "26-50")) %>%
  group_by(Province_R, stops) %>%
  summarize(nStops = n_distinct(Stop)) %>%
  st_cast("LINESTRING")

ca_bbs_lines$rte_length <- st_length(ca_bbs_lines)

len <- units::set_units(300000, "m")
ca_bbs_lines <- ca_bbs_lines %>%
  filter(rte_length < len) 

ca_bbs_buffers <- ca_bbs_lines %>%
  st_transform(st_crs(us_bbs_buffers)) %>%
  st_buffer(5000)
# write_sf(ca_bbs_buffers, "C:/Users/gdicecco/Desktop/data/canada_bbs_half_route_paths_5km.shp")


