### Create 1/2 route path shapefiles for land cover analysis

library(tidyverse)
library(sf)
library(lwgeom)

### Canada and US BBS route paths with 5-km buffers

bbs_route_paths <- read_sf("\\\\BioArk.bio.unc.edu/HurlbertLab/Databases/BBS/GPS_stoplocations/BBS_routes_USandCanada/bbs_routes_usandcanada.shp")

# Use new Canada route paths - 2019 data

ca_bbs <- read_sf("\\\\BioArk.bio.unc.edu/HurlbertLab/Databases/BBS/GPS_stoplocations/BBS_routes_Canada_dec2019.gdb")

# Convert multilinestring to linestring
# st_linesubstring(line, from = 0, to = 0.5)
# st_linesubstring(line, from = 0.5, to = 1)

bbs_route_lines <- st_cast(bbs_route_paths, "LINESTRING")

# Many BBS route paths can't be converted into linestrings (bbs_route_paths has 6000 rows, bbs_route_lines has 10000)

st_length(bbs_route_paths) # lengths of route paths
st_centroid(bbs_route_paths) # centroids of route paths

