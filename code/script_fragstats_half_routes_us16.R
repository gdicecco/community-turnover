# Calculate fragmentation measures for half BBS routes, US
# 2016

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(SDMTools)

# Read in BBS routes shapefile

# 2016
setwd("/proj/hurlbertlab/nlcd_landcover/NLCD_2016_Land_Cover_L48_20190424/")
nlcd <- raster("nlcd_2016_whole_simplified.tif")
routes <- readOGR("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/us_bbs_half_route_paths_5km.shp")
routes_transf <- spTransform(routes, crs(nlcd))
routes_tr <- crop(routes_transf, extent(nlcd))

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
    filename <- paste0("classStat_nlcd_30x30_2016_route_", rtenum,"_", stp, ".csv")
    write.csv(class, filename, row.names = F)
  }

}
