### Climate change at 1/2 route scale - crop daymet and get breeding season averages for 1/2 routes
## Run on longleaf

## Libraries

library(sf)
library(rgdal)
library(rgeos)
library(daymetr)
library(tidyverse)
library(raster)
library(ncdf4)
library(lubridate)
library(units)
library(purrr)

## Read in data

# Half-route paths

ca_routes <- read_sf("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/canada_bbs_half_route_paths_5km.shp") %>%
  dplyr::select(rteno, stops, RTENAME, STATUS, geometry) %>%
  mutate(country = "Canada")
us_routes <- read_sf("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/us_bbs_half_route_paths_5km.shp") %>%
  dplyr::select(rteno, stops, RTENAME, STATUS, geometry) %>%
  mutate(country = "US")

na_routes <- rbind(ca_routes, us_routes)

# Climate data - DaymetR

# Function: get average temperature for one BBS route polygon

daymetMean <- function(stateroute) {
  rte <- filter(na_routes_transf, rteno == stateroute)
  
  daymet_crop <- raster::crop(daymet_mean, rte)
  daymet_mask <- raster::mask(daymet_crop, rte)
  
  local_extract <- raster::extract(daymet_mask, rte, fun = mean, na.rm = T, df = T)
  
  return(mean(local_extract$layer))
}

# Function errors if some routes are outside the extent of DAYMET raster - use possibly function to ignore those routes

possibly_daymetMean <- possibly(daymetMean, otherwise = NA)

# Calculate mean annual breeding season temp at BBS routes

setwd("/proj/hurlbertlab/gdicecco/")

years <- c(1990:2017)

for(y in years) {
  download_daymet_ncss(location = c(60, -145, 15, -52),
                       start = y,
                       end = y,
                       param = c("tmin", "tmax"), 
                       frequency = "monthly",
                       path = "daymet/")
  
  # Read in data
  files <- list.files("daymet/")
  
  for(f in files) {
    daymet_nc <- nc_open(paste0("daymet/", f))
    daymet_raster <- brick(paste0("daymet/", f))
    crs(daymet_raster) <- "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0"
    
    daymet_breeding <- daymet_raster[[5:7]]
    
    daymet_mean <- mean(daymet_breeding, na.rm = T)
    
    crs_daymet <- crs(daymet_breeding)
    
    na_routes_transf <- st_transform(na_routes, "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0")

    routeclim <- data.frame(stateroute = na_routes_transf$rteno) %>%
      mutate(mean_temp = purrr::map_dbl(stateroute, possibly_daymetMean))
    
    routeclim <- bind_rows(routeclim_us, routeclim_ca)
    
    write.csv(routeclim, paste0("/proj/hurlbertlab/gdicecco/half_bbs_route_daymet_out/", f, ".csv"), row.names = F)
    print(f)
    nc_close(daymet_nc)
  }
  
  print(y)
  sapply(paste0("daymet/", files), unlink)
  
}

# Read in indiv year files
# Join together

dir <- "/proj/hurlbertlab/gdicecco/half_bbs_route_daymet_out/"

routeDAYMET <- data.frame(stateroute = c(), year = c(), mean_tmax = c(), mean_tmin = c())

for(y in years) {
  files <- list.files(dir)
  files_y <- files[grepl(y, files)]
  
  tmax <- read.csv(paste0(dir, files_y[grepl("tmax", files_y)])) %>%
    group_by(stateroute) %>%
    summarize(mean_tmax = mean(mean_temp, na.rm = T))
  tmin <- read.csv(paste0(dir, files_y[grepl("tmin", files_y)])) %>%
    group_by(stateroute) %>%
    summarize(mean_tmin = mean(mean_temp, na.rm = T))
  
  tmp <- tmax %>%
    left_join(tmin) %>%
    mutate(year = y)
  print(nrow(tmp))
  
  routeDAYMET <- rbind(routeDAYMET, tmp)
}

setwd("/proj/hurlbertlab/gdicecco/")
write.csv(routeDAYMET, "bbs_half_routes_breeding_season_climate.csv", row.names = F)
