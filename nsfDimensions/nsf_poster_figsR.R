## =====================================
## script to make figures for NSF poster
## =====================================

setwd('~/Dropbox/hawaiiDimensions/presentations/nsfDimensions')
corepar <- list(fg='black', bg='white', cex.lab=5, cex.axis=5.5)

## =================================
## map of hawaii
## =================================

library(rgdal)
library(RColorBrewer)

## load flow layers
oldwd <- setwd('~/Dropbox/hawaiiDimensions/geoData/site_selection/Haw_St_shapefiles/Haw_St_geo_20070426_region')
hi.geo <- readOGR('.', 'Haw_St_geo_20070426_region')

chrono.age <- read.csv("../Haw_St_ageCode.csv",stringsAsFactors=FALSE)
chrono.age$age.low <- chrono.age$age.low*c(yr=10^-6,ka=10^-3,Ma=10^0)[chrono.age$unit]
chrono.age$age.hi <- chrono.age$age.hi*c(yr=10^-6,ka=10^-3,Ma=10^0)[chrono.age$unit]
chrono.age <- chrono.age[,-4]

chrono.age <- rbind(chrono.age, cbind(code=c(13:14), age.low=c(2, 4), age.hi=c(4, 6)))

setwd(oldwd)

## colors for flow ages
geo.col <- c('gray45', colorRampPalette(brewer.pal(9,"YlGnBu")[-9])(max(hi.geo@data$AGE_GROUP)))

jpeg(filename='fig_hawaii_flowAge.jpg', bg='black', width=4800, height=4800, quality=100)
par(corepar)
par(mar=rep(0, 4))
plot(hi.geo, col=geo.col[hi.geo$AGE_GROUP+1], border=geo.col[hi.geo$AGE_GROUP+1])
dev.off()

## legend
source('~/R_functions/logAxis.R')

jpeg(filename='fig_hawaii_flowLegend.jpg', width=450*2, height=150*2, quality=100, pointsize=24)
par(corepar)
par(mar=c(2.5, 1, 0, 1) + 1, mgp=c(2,0.5,0), xpd=NA, cex.lab=1.5, cex.axis=1)

plot(1,xlim=range(chrono.age[,-1])+0.0001,ylim=0:1,log="x",type="n",
     axes=FALSE,xlab="Substrate age (My)",ylab="",xaxs="i",yaxs="i")

apply(chrono.age,1,function(x) rect(x[2]+0.0001,0,x[3]+0.0001,1,col=geo.col[x[1]+1],border=NA))
box()
logAxis(1, labels=FALSE)
axis(side=1, at=c(10^(-3:0), 6))

dev.off()