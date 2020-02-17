# Canada frag stats

library(raster)
library(SDMTools)
library(rgdal)
library(sp)
library(rgeos)

setwd("/proj/hurlbertlab/gdicecco/canada_landuse")
dir <- getwd()

# read in Canada BBS routes

ca_routes <- readOGR("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/canada_bbs_half_route_paths_5km.shp")

# Crop out BBS route & compute fragstats - go zone by zone

folder <- "CanadaLandUse2000"

files <- list.files(paste0(dir, "/", folder))
files <- files[9:12]

for(f in files) {
  zone_raster <- raster(paste0(dir, "/", folder, "/", f, "/", f, ".tif"))
  
  zone_crs <- crs(zone_raster)
  
  ca_transf <- spTransform(ca_routes, zone_crs)
  
  routes_tr <- crop(ca_transf, extent(zone_raster))
  
  if(is.null(routes_tr)) {
    print(f)
  } else {
    
    routenos <- routes_tr@data[ , 1]
    
    for(i in 1:nrow(routes_tr@data)) {
      rte <- subset(routes_tr, Province_R == routenos[i])
      rtenum <- routenos[i]
      
      for(i in 1:nrow(rte)){
        half <- rte[i, ]
        stp <- half$stops
      
        nlcd_crop <- crop(zone_raster, half)
        nlcd_mask <- mask(nlcd_crop, half)
      
        mask_poly <- rasterToPolygons(nlcd_mask)
      
      if(is.null(mask_poly)) {
        print(f)
        print(rtenum)
      } else {
        
        class <- ClassStat(nlcd_mask)
        mask_poly <- rasterToPolygons(nlcd_mask)
        area_polygon <- area(mask_poly)
        class$area <- sum(area_polygon)
        filename <- paste0("/proj/hurlbertlab/gdicecco/canada_bbs_half_route_fragstats/classStat_ca_30x30_2000_zone_", 
                           f, "route_", rtenum, "_", stp, ".csv")
        write.csv(class, filename, row.names = F) 
      }
      }
    }
  }
  
}
