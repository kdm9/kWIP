#!/usr/bin/env Rscript
library(fields, quietly = T)
argv = commandArgs(trailingOnly=T)
dist.file = paste0(argv[1], ".dist")
kern.file = paste0(argv[1], ".kern")

base = basename(argv[1])
dist = as.matrix(read.delim(dist.file)[, -1])
kern = as.matrix(read.delim(kern.file)[, -1])
n = ncol(dist)

clust = hclust(as.dist(dist))
sam_ord = clust$order

pdf(paste0(base, ".pdf"))
image.plot(1:n, 1:n, dist[sam_ord, sam_ord],
           axes=FALSE,
           xlab="",
           ylab="",
           col=rainbow(80, start=0.2),
           main=paste("Distance of", base))
axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)

image.plot(1:n, 1:n, kern[sam_ord, sam_ord],
           axes=FALSE,
           xlab="",
           ylab="",
           col=rainbow(80, start=0.2),
           main=paste("Kernel of", base))
axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)

d = diag(kern)
geom.mean.d = sqrt(d %*% t(d))
normk = kern / geom.mean.d

image.plot(1:n, 1:n, normk[sam_ord, sam_ord],
           axes=FALSE,
           xlab="",
           ylab="",
           col=rainbow(80, start=0.2),
           main=paste("Normalised Kernel of", base))
axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)

plot(clust, cex=0.5, main=base)

fit <- cmdscale(dist, eig=TRUE, k=2) # k is the number of dim

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS", type='p')
#text(x, y, labels = colnames(dist), cex=.4)

dev.off()
