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
