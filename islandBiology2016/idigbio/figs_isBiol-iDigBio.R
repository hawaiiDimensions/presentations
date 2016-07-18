library(ape)
library(devtools)
library(sp)
library(maptools)
library(rgdal)
library(maps)
library(mapdata)
load_all('~/Dropbox/Research/meteR')
load_all('~/Dropbox/Research/socorro')

## place to save figures
setwd('~/Dropbox/hawaiiDimensions/presentations/islandBiology2016/idigbio')

## gbif and hdim

thiswd <- getwd()
eval(parse(text=c('{', readLines('~/Dropbox/hawaiiDimensions/geodata/maps/sites_map.R', n=29), '}')))
setwd(thiswd)

gbif <- occ_data(scientificName = c('Insecta', 'Arachnida'), 
                 decimalLatitude = '18,24', decimalLongitude = '-161,-154', 
                 limit = 200000)

rec <- rbind(as.data.frame(gbif[[1]][[2]][, c('decimalLongitude', 'decimalLatitude')]),
             as.data.frame(gbif[[2]][[2]][, c('decimalLongitude', 'decimalLatitude')]))

rec <- SpatialPoints(rec, proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
rec <- spTransform(rec, CRS(proj4string(islands)))
rec <- rec[!is.na(over(rec, islands))]

plots <- readOGR('/Users/ajr/Dropbox/hawaiiDimensions/geodata/sites', 'dimensions_plots')
plots <- spTransform(plots, CRS(proj4string(islands)))

pdf('fig_gbif1.pdf', width=5, height=3)
par(mar=rep(0, 4))
plot(islands)
points(plots, pch=16)
dev.off()

pdf('fig_gbif2.pdf', width=5, height=3)
par(mar=rep(0, 4))
plot(rec, pch=16, cex=0.5, col='gray')
plot(islands, add=TRUE)
points(plots, pch=16)
dev.off()

## SAD

this.sad <- sad(meteESF(S0=300, N0=3000))

pdf('fig_sad.pdf', width=3, height=3)
par(lwd=2, mgp=c(1, 0, 0), mar=c(2, 2, 0, 0) + 0.1)
plot(this.sad, ptype = 'rad', add.legend = FALSE, th.col = 'black', axes=FALSE, log='y')
par(lwd=1)
axis(1, labels=NA)
logAxis(2, labels=NA)
box()
dev.off()


## phylo
set.seed(0)
tre <- rphylo(30, 1, 0.8)

pdf('fig_phylo.pdf', width = 3, height = 3)
par(mar=rep(0.1, 4))
plot(tre, show.tip.label = FALSE, edge.width = 2)
dev.off()

## sequence
nnuc <- 20
nseq <- 8
mmat <- c('A', 'T', 'G', 'C')
names(mmat) <- c('T', 'A', 'C', 'G')

con <- sample(mmat, nnuc, rep=TRUE)

pdf('fig_seq.pdf', width = 6.5, height = 3.5)

par(mar=rep(0.2, 4), xpd=NA)
plot(1, xlim=c(1, nnuc), ylim=c(2, 1), type='n', axes=FALSE)
cxy <- abs(par('cxy'))
yy <- cumsum(c(1, rep(cxy[2]+cxy[2]/1.2, nseq)))

seqi <- con

for(i in 1:length(yy)) {
    if(i > 1) {
        mut <- sample(nnuc, sample(6, 1))
        seqi[mut] <- mmat[seqi[mut]]
        
        realMut <- which(seqi != con)
        rect(xleft = realMut - cxy[1]/1.2, ybottom = yy[i] - cxy[2]/1.2, 
             xright = realMut + cxy[1]/1.2, ytop = yy[i] + cxy[2],
             border = NA, col='gray35')
    }
    
    text(1:nnuc, rep(yy[i], nnuc), labels = seqi, col=c('white', 'black')[as.numeric(seqi == con) + 1])
    rect(xleft = 1 - cxy[1]/1.2, ybottom = yy[i] - cxy[2]/1.2, 
         xright= 1 - cxy[1]/4 + nnuc*cxy[1]*2 + cxy[1]/1.5, ytop = yy[i] + cxy[2])
}

dev.off()


## sampling

n <- 100
xy <- matrix(runif(2*n), ncol=2)
hs <- c(0.15, 0.3, 0.5, 0.6)
cols <- sample(hs, n, rep=TRUE)
detect <- sample(c(0, 0.9), n, rep=TRUE)

pdf('fig_sampComplete.pdf', width=4, height=4)
par(mar=rep(0.1, 4))
plot(xy, bg=hsv(cols), pch=21, cex=2, axes=FALSE)
box()
dev.off()

pdf('fig_sampIncomplete.pdf', width=4, height=4)
par(mar=rep(0.1, 4))
plot(xy, bg=hsv(cols, 1-detect, 1), col=hsv(0, 0, detect), pch=21, cex=2, axes=FALSE)
box()
dev.off()


## latent process

randLogistic <- function(x, b, k, s) {
    ceiling(k / (1 + exp(-b*x)) * exp(rnorm(length(x), 0, s)))
}

y1 <- randLogistic(seq(-10, 12, by=0.25), 0.5, 100, 0.05)
y2 <- randLogistic(seq(-10, 12, by=0.25), 0.5, 80, 0.05)
y3 <- randLogistic(seq(-10, 12, by=0.25), 0.75, 150, 0.075)
y4 <- randLogistic(seq(-10, 12, by=0.25), 0.25, 50, 0.025)

pdf('fig_latentMod.pdf', width=4, height=4)
par(mar=c(3, 3, 0, 0) + 0.1, mgp=c(1.5, 1, 0))
plot(y1, type='l', col=hsv(hs[4]), ylim=range(y1, y2, y3, y4), 
     lwd=2, axes=FALSE, xlab='Time', ylab='Population size')
box()
axisArrows(length=0.1)
lines(y2, col=hsv(hs[2]), lwd=2)
lines(y3, col=hsv(hs[3]), lwd=2)
lines(y4, col=hsv(hs[1]), lwd=2)
dev.off()


## network

circweb <- function(x, rowc, colc, ...) {
    x <- x[rowSums(x) > 0, ]
    x <- x[, colSums(x) > 0]
    xy <- cbind(cos(seq(0, 2*pi, length=1+sum(dim(x)))), sin(seq(0, 2*pi, length=1+sum(dim(x)))))[-1, ]
    rownode <- rep(1:nrow(x), ncol(x))
    colnode <- rep(1:ncol(x) + nrow(x), each=nrow(x))
    links <- as.vector(x)
    
    xnode1 <- xy[rownode[links > 0], 1]
    xnode2 <- xy[colnode[links > 0], 1]
    ynode1 <- xy[rownode[links > 0], 2]
    ynode2 <- xy[colnode[links > 0], 2]
    
    plot(xy, col=c(rep(rowc, nrow(x)), rep(colc, ncol(x))), 
         panel.first=segments(xnode1, ynode1, xnode2, ynode2),
         asp=1, axes=FALSE, xlab='', ylab='', ...)
    
}

nr <- 30
nc <- 20
bimat <- matrix(0, nrow=nr, ncol=nc)
bimat[upper.tri(bimat)] <- 1
bimat[upper.tri(bimat) & sample(c(TRUE, FALSE), nr*nc, prob=c(6, 1), rep=TRUE)] <- 0

pdf('fig_network.pdf', width=4, height=4)

close.screen(all.screens = TRUE)
par(mar=rep(0.1, 4), fg='black')
circweb(bimat, rowc=hsv(0.45, 0.5, 0.9), colc=hsv(0.1, 0.9, 0.7), pch=16, cex=2)

split.screen(c(1, 1), erase=FALSE)
par(fg='transparent')
circweb(bimat, rowc='black', colc='black', cex=2)

dev.off()
close.screen(all.screens = TRUE)
