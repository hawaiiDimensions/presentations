library(sp)
library(maptools)
library(rgdal)
library(maps)
library(mapdata)
library(RColorBrewer)

setwd('~/Dropbox/hawaiiDimensions/presentations/evol2015')

## load flow layers
oldwd <- setwd('~/Dropbox/hawaiiDimensions/site_selection/Haw_St_shapefiles/Haw_St_geo_20070426_region')
hi.geo.poly <- readOGR('.', 'Haw_St_geo_20070426_region')
setwd(oldwd)


## colors for flow ages
geo.col <- colorRampPalette(brewer.pal(9,"YlGnBu"))(max(hi.geo.poly@data$AGE_GROUP) + 1)

ka.geo <- hi.geo.poly[hi.geo.poly$ISLAND=='Kauai', ]
ma.geo <- hi.geo.poly[hi.geo.poly$ISLAND %in% c('Molokai', 'Maui'), ]
hi.geo <- hi.geo.poly[hi.geo.poly$ISLAND=='Hawaii', ]

pdf(file='fig_chrono.pdf', width=5, height=5)

map('worldHires', xlim=c(-160, -154), ylim=c(18, 23), resolution=0, mar=rep(0, 4))

plot(spTransform(hi.geo, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')), 
     col=geo.col[hi.geo$AGE_GROUP], border=geo.col[hi.geo$AGE_GROUP], add=TRUE)
plot(spTransform(ma.geo, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')), 
     col=geo.col[ma.geo$AGE_GROUP], border=geo.col[ma.geo$AGE_GROUP], add=TRUE)
plot(spTransform(ka.geo, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')), 
     col=geo.col[ka.geo$AGE_GROUP], border=geo.col[ka.geo$AGE_GROUP], add=TRUE)

map('worldHires', xlim=c(-160, -154), ylim=c(18, 23), resolution=0, add=TRUE)

dev.off()


## phylogeny
library(ape)

set.seed(123)
tre <- rphylo(30, 0.2, 0.15, fossils=TRUE)

pdf(file='fig_phylo.pdf', width=8, height=4)
par(mar=rep(0, 4))
plot(tre, show.tip.label=FALSE, edge.width=2)
dev.off()

## smaller phylo
set.seed(12)
tre <- rphylo(5, 0.2, 0.15, fossil=TRUE)
pdf(file='fig_phyloSml.pdf', width=8, height=4)
plot(tre, show.tip.label=FALSE, edge.width=2, edge.col='gray')
dev.off()

## bigger phylo
set.seed(123)
tre <- rphylo(50, 0.2, 0.15, fossil=FALSE)
pdf(file='fig_phyloBig.pdf', width=4, height=4)
plot(tre, show.tip.label=FALSE, edge.width=2)
dev.off()

## phylo for mete
set.seed(123)
tre <- rphylo(12, 0.2, 0.15, fossil=FALSE)
pdf(file='fig_phyloMETE.pdf', width=4, height=4)
plot(tre, show.tip.label=FALSE, edge.width=2)
dev.off()
