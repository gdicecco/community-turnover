# Calculate fragmentation measures for half BBS routes, US
# 1992

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(dplyr)
library(stringr)
library(SDMTools)

# Read in BBS routes shapefile

# 1992
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_1992_landcover_2018_08_31/")
nlcd <- raster("nlcd_1992_whole_simplified.tif")
routes <- readOGR("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/us_bbs_half_route_paths_5km.shp")
routes_tr <- spTransform(routes, crs(nlcd))

routenos <- routes_tr@data[ , 1]

setwd("/proj/hurlbertlab/gdicecco/nlcd_bbs_half_route_fragstats/")
for(i in 1:nrow(routes_tr@data)) {
  rte <- subset(routes_tr, rteno == routenos[i])
  rtenum <- routenos[i]
  
  for(i in 1:nrow(rte)){
    half <- rte[i, ]
    stp <- half$stops
    nlcd_crop <- crop(nlcd, half)
    nlcd_mask <- mask(nlcd_crop, half)
    class <- ClassStat(nlcd_mask)
    filename <- paste0("classStat_nlcd_30x30_1992_route_", rtenum,"_", stp, ".csv")
    write.csv(class, filename, row.names = F)
  }
  
}
