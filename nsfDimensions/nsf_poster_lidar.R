## =====================================
## script to make LiDAR figure for NSF poster
## =====================================

corepar <- list(fg='black', bg='white', cex.lab=1.5, cex.axis=1.5)

## you'll need to change this wd to the directory where you have the CAO data
setwd('~/Dropbox/hawaiiDimensions/presentations/nsfDimensions')

## and of course change this to the correct raster file. I'm thinking the best one would be 
## `LaupahoehoeLSAG_14-65K'

r <- raster('laup_sim.tif')

## this should most likely work on your computer if you accept the `hawaiiDimensions' dropbox inviation
old.wd <- setwd('~/Dropbox/hawaiiDimensions/geodata')
site.poly <- readOGR('./sites', 'dimensions_plots')
setwd(old.wd)

laup <- site.poly[grep('laupLSAG', site.poly@data$name), ]
laup <- spTransform(laup, CRS(projection(r)))
laup.vals <- extract(r, laup)
laup.vals[is.na(laup.vals)] <- runif(sum(is.na(laup.vals)), min(laup.vals, na.rm=TRUE), max(laup.vals, na.rm=TRUE))


## this figure file will be in what ever directory this is:
getwd()

pdf(file='fig_canopyHeight.pdf', width=7, height=4)

par(mfrow=c(1, 2), oma=c(4, 1, 0, 1)+0.1, mar=rep(0.1, 4), mgp=c(2.25, 1, 0))
par(corepar)
plot(r, axes=FALSE, col=colorRampPalette(c(hsv(0.2, 0.6, 1), hsv(0.5, 0.8, 0.4)))(20), legend=FALSE)
points(laup, pch=21, bg='black', col='white', cex=2, lwd=2)

par(xpd=NA)
xyhist <- hist(values(r), plot=FALSE)
hist(values(r), ylim=c(-0.03*max(xyhist$counts), max(xyhist$counts)), 
     yaxt='n', main='', ylab='', xlab='Canopy heights (m)',
     xaxs='i')
rect(xleft=seq(xyhist$breaks[1], max(xyhist$breaks), length=20)[-20], 
     xright=seq(xyhist$breaks[1], max(xyhist$breaks), length=20)[-1], 
     ybottom=par('usr')[3], ytop=0, 
     col=colorRampPalette(c(hsv(0.2, 0.6, 1), hsv(0.5, 0.8, 0.4)))(20),
     border=NA)
rect(xleft=xyhist$breaks[1], xright=max(xyhist$breaks), ybottom=par('usr')[3], ytop=0)
segments(x0=laup.vals, y0=par('usr')[3], y1=0, lwd=3)

dev.off()
