## =====================================
## script to make figures for NSF poster
## =====================================

setwd('~/Dropbox/hawaiiDimensions/presentations/nsfDimensions')
corepar <- list(fg='black', bg='white', cex.lab=2, cex.axis=1.5)

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

jpeg(filename='fig_hawaii_flowAge.jpg', width=4800, height=4800, quality=100)
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


## =================================
## simulate LiDAR data
## =================================

library(sp)
library(raster)
library(rgdal)
library(maptools)

old.wd <- setwd('~/Dropbox/hawaiiDimensions/geodata')
site.poly <- readOGR('./sites', 'dimensions_plots')
setwd(old.wd)

laup <- site.poly[grep('laupLSAG', site.poly@data$name), ]
laup <- spTransform(laup, CRS('+proj=utm +zone=5 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
laupbbox <- bbox(laup)

x <- seq(laupbbox[1, 1]-1, laupbbox[1, 2]+1, by=1)
y <- seq(laupbbox[2, 1]-1, laupbbox[2, 2]+1, by=1)
z <- matrix(rnorm(length(x)*length(y), 13.5, 3), nrow=length(x), ncol=length(y))
z[z < 0] <- runif(sum(z < 0), 0, 25)
r <- raster(z, xmn=laupbbox[1, 1], xmx=laupbbox[1, 2], ymn=laupbbox[2, 1], ymx=laupbbox[2, 2], 
            crs=projection(laup))

writeRaster(r, filename='laup_sim.tif', format='GTiff', overwrite=TRUE)


## =================================
## metabarcoding trial data
## =================================

x <- read.csv('Spider_DNA_vs_Reads.csv', as.is=TRUE)

pdf(file='fig_metaB_res.pdf', width=5, height=5)
par(corepar)
par(mar=c(4, 4, 1, 1) + 0.1, mgp=c(2.75, 1, 0))
plot(x, xlab='Fraction DNA', ylab='Fraction reads', xlim=range(x), ylim=range(x), cex=2)
abline(0, 1)
r2 <- round(summary(lm(x$percentReads ~ x$percentDNA))$r.squared, 2)
legend('topleft', legend=bquote(R^2==.(r2)), bty='n')
dev.off()


## =================================
## METE concept fig
## =================================

Ebreak <- cumsum(sort(0.03 + (1-0.03*12)*rmultinom(1,1000,rexp(11))[,1]/1000,decreasing=TRUE))
Ebreak <- c(0,Ebreak,1)

indAt <- seq(0.1,0.9,length=12)
sppAt <- seq(0.3,0.7,by=0.2)

pdf(file='fig_meteExplc.pdf', width=3,height=4)
par(mar=c(4,0,0,0)+0.1,bg="transparent",fg="black",col.axis="black")
plot(1,type="n",xlim=c(0.4, 1), ylim=c(0,1),xlab="",axes=FALSE)

rect(0.9,0,1,1,col=hsv(0.42,0.6,0.6),border=NA)
segments(x0=0.9,x1=1,y0=Ebreak[2:12],col=hsv(0.42,0.6,0.8),lwd=2)

rect(0.65,indAt,0.75,indAt+0.03,col=hsv(0.6,0.3,0.7),border=NA)

for(i in 1:12) {
    polygon(x=c(0.75,0.75,0.9,0.9),
            y=c(indAt[i],indAt[i]+0.03,Ebreak[i+1],Ebreak[i]),
            col=hsv(0.42,0.6,0.8),border=hsv(0.42,0.6,0.8),lwd=0.2)
    
    if(i < 3) {
        polygon(x=c(0.65,0.65,0.45),
                y=c(indAt[i],indAt[i]+0.03,sppAt[1]),
                col=hsv(0.6,0.3,0.8),border=hsv(0.6,0.3,0.8))
    } else if(i < 7) {
        polygon(x=c(0.65,0.65,0.45),
                y=c(indAt[i],indAt[i]+0.03,sppAt[2]),
                col=hsv(0.6,0.3,0.8),border=hsv(0.6,0.3,0.8))
    } else {
        polygon(x=c(0.65,0.65,0.45),
                y=c(indAt[i],indAt[i]+0.03,sppAt[3]),
                col=hsv(0.6,0.3,0.8),border=hsv(0.6,0.3,0.8))
    }
}

points(rep(0.45,3),sppAt,cex=4,pch=21,bg=hsv(0.05,0.7,0.8))

axis(1,at=c(0.45,0.7,0.95),labels=c('Species','Individuals','Energy'))

dev.off()


## =================================
## METE concept fig
## =================================