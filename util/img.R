#!/usr/bin/env Rscript
library(fields)
argv = commandArgs(trailingOnly=T)
dist.file = paste0(argv[1], ".dist")
kern.file = paste0(argv[1], ".kern")

base = basename(argv[1])
dist = as.matrix(read.delim(dist.file)[, -1])
kern = as.matrix(read.delim(kern.file)[, -1])
n = ncol(dist)

clust = hclust(as.dist(dist))
sam_ord = clust$order
sam_ord = 1:n

pdf(paste0(base, ".pdf"))
image.plot(1:n, 1:n, dist[sam_ord, sam_ord],
           axes=FALSE,
           xlab="",
           ylab="",
           col=rainbow(80, start=0.2),
           main=paste("Distance of", base))
axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.2, las=2, tick=F)
axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.2, las=2, tick=F)

image.plot(1:n, 1:n, kern[sam_ord, sam_ord],
           axes=FALSE,
           xlab="",
           ylab="",
           col=rainbow(80, start=0.2),
           main=paste("Kernel of", base))
axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.2, las=2, tick=F)
axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.2, las=2, tick=F)

plot(clust, cex=0.2, main=base)
dev.off()
