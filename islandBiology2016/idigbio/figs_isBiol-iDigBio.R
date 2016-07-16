library(ape)
library(devtools)
load_all('~/Dropbox/Research/meteR')
load_all('~/Dropbox/Research/socorro')

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
