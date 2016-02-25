setwd('~/Dropbox/hawaiiDimensions/geb_paper/islandBiology')

## substrait age and assembly figure
pdf(width=4.5, height=3.8, file='fig_age-assembly.pdf')
par(mar=c(0.05, 2, 0.05, 0), mgp=c(1, 0, 0), bg='white')
curve(1*x + 0.08, lwd=3, axes=FALSE, frame.plot=TRUE, 
      ylab='Relative importance', xlab='', col.lab='gray',
      ylim=c(0, 1.1))
curve(0.13*exp(2*x) - 0.1, add=TRUE, lwd=3)
curve(0.1*x, add=TRUE, lwd=3)
curve(2.02/(1 + exp(-10*(x)))-0.9, add=TRUE, lwd=3)
mtext('Ecology', side=2, line=1, at=0.1)
mtext('Evolution', side=2, line=1, at=1)
dev.off()

## phylo figure
library(ape)
x <- rcoal(10)
pdf(width=3, height=3, file='fig_phylo_eg2.pdf')
par(mar=rep(0, 4))
plot(x, show.tip.label=FALSE, edge.width=3)
dev.off()


## mac-Wilson
x <- exp(-seq(0, 1, length=40))

pdf(width=4.5, height=3.8, file='fig_macWilson.pdf')
par(mar=c(2, 2.25, 0.05, 0), mgp=c(0.75, 0, 0))
plot(x, type='l', col='blue', lwd=3, 
     xlab='Species richness', ylab='Rates',
     axes=FALSE, frame.plot=TRUE, cex.lab=1.5)
lines(rev(x), col='red', lwd=3)
text(c(12, 30), c(0.9, 0.9), labels=c('Immigration', 'Extinction'),
     col=c('blue', 'red'), cex=1.2)
dev.off()

## fitness landscape
library(MASS)
x <- kde2d(c(rnorm(800, mean=0), rnorm(2000, mean=2), rnorm(1000, mean=1)),
           c(rnorm(800, mean=2), rnorm(2000, mean=4), rnorm(1000, mean=0)),
           n=64)

surf.colors <- function(x, col = terrain.colors(20)) {
  x <- x$z
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
             x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4

  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]

  return(colors)
}

pdf(width=4.5, height=3.8, file='fig_fitnessLandscape.pdf')
par(mar=rep(0, 4))
persp(x, col=surf.colors(x, col=topo.colors(32)), phi=30, theta=60,
      border='gray75', lwd=0.01, box=FALSE)
dev.off()

## function to plot bipartite networks

plot.bimat <- function(x, col=hsv(c(0.08, 0.45), c(1, 0.9), c(0.7, 0.7)), file=NULL) {
	plot(rep(1:2, dim(x)), c(1:nrow(x), 1:ncol(x)), 
	     ylim=c(max(dim(x)), 1), axes=FALSE, cex=4, xpd=NA,
	     panel.first={
	     	for(i in 1:nrow(x)) {
				segments(x0=1, x1=2, y0=i, y1=which(x[i, ] > 0), lwd=5, col='gray25')
			}
	     }, xlab='', ylab='', main='', 
	     pch=21, bg=rep(col, dim(x)))
}

## nestedness and modularity explained
n <- 5

net.nest <- matrix(0, nrow=n, ncol=n)
net.nest[upper.tri(net.nest, diag=TRUE)] <- 1
net.nest <- net.nest[, ncol(net.nest):1]

net.mod <- matrix(0, nrow=n, ncol=n)
net.mod[1:2, 1:2] <- 1
net.mod[2, 2] <- 0
net.mod[3:4, 3:4] <- 1
net.mod[4, 4] <- 0
net.mod[5:n, 5:n] <- 1

pdf(width=4, height=4, file='fig_nested.pdf')
par(mar=rep(0.1, 4))
plot.bimat(net.nest)
dev.off()

pdf(width=4, height=4, file='fig_nestedImm.pdf')
par(mar=rep(0.1, 4))
plot.bimat(net.nest[, -5])
points(2, 5, bg=hsv(0.45, 0.9, 0.7), pch=21, cex=4)
dev.off()

pdf(width=4, height=4, file='fig_mod.pdf')
par(mar=rep(0.1, 4))
plot.bimat(net.mod)
dev.off()

pdf(width=4, height=4, file='fig_modEvol.pdf')
par(mar=rep(0.1, 4))
plot.bimat(net.mod, col=rep('transparent', 2))
points(rep(1, 5), 1:5, bg=hsv(0.08, c(0.6, 0.6, 1, 1, 0.6), 
                              c(1, 1, 0.7, 0.7, 0.7)), 
       pch=21, cex=4)
points(rep(2, 5), 1:5, bg=hsv(0.45, c(0.6, 0.6, 1, 1, 0.4), 
                              c(1, 1, 0.5, 0.5, 0.6)), 
       pch=21, cex=4)
dev.off()


net.rand <- matrix(sample(0:1, size=10*10, rep=TRUE), nrow=10)
pdf(width=4, height=4, file='fig_exampleBipartite.pdf')
par(mar=rep(0.1, 4))
plot.bimat(net.rand)
dev.off()


pdf(width=4, height=4, file='fig_degreeDist.pdf')
par(mar=c(2, 2, 0, 0)+0.1, mgp=c(0.75, 1, 0))
curve(exp(-x), from=0, to=2, axes=FALSE, frame.plot=TRUE, xlab='Degree', ylab='Probability')
polygon(c(0, seq(0, 2, length=20)), c(exp(-2), exp(-seq(0, 2, length=20))), col='gray')
curve(exp(-x), lwd=3, add=TRUE)
dev.off()


pdf(width=4, height=4, file='fig_spaceTime.pdf')
par(mar=c(2, 2, 0, 0)+0.1, mgp=c(0.75, 1, 0))
plot(1, type='n', axes=FALSE, frame.plot=TRUE, xlab='Space', ylab='Time')
dev.off()

