#!/usr/bin/env Rscript
library(RColorBrewer)
argv = commandArgs(trailingOnly=T)
dist.file = paste0(argv[1], ".dist")

base = basename(argv[1])
dist = as.matrix(read.delim(dist.file)[, -1])
n = ncol(dist)

pdf(paste0(base, ".pdf"))
image(1:n, 1:n, dist,
      axes=FALSE,
      xlab="",
      ylab="",
      col=brewer.pal(9, "Blues"),
      main=paste("Distance of", base))
axis(1, 1:n, colnames(dist), srt=45, xpd=T, tick=F)
axis(2, 1:n, colnames(dist), srt=90, tick=F)

plot(hclust(as.dist(dist)), main=base)
dev.off()
