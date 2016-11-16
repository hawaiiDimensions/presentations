library(sp)
library(maptools)
library(rgdal)
library(rgeos)

setwd('~/Dropbox/hawaiiDimensions/presentations/calAcad_2016-11')

## ===================
## hawaii map in white
## ===================

## load flow layers
oldwd <- setwd('~/Dropbox/hawaiiDimensions/geoData/site_selection/Haw_St_shapefiles/Haw_St_geo_20070426_region')
hi.geo.poly <- readOGR('.', 'Haw_St_geo_20070426_region')

## get island outlines
islands <- gUnionCascaded(hi.geo.poly)
islands <- SpatialPolygons(list(
    Polygons(islands@polygons[[1]]@Polygons[sapply(islands@polygons[[1]]@Polygons, 
                                                   function(p) !p@hole & p@area > 1e+07)
                                            ], ID=1)), 
    proj4string = CRS(proj4string(hi.geo.poly)))

## plot it
setwd(oldwd)
pdf('fig_hawaiiMap.pdf', width = 5, height = 3)
par(mar = rep(0, 4))
plot(islands, col = 'white', border = 'white')
dev.off()


