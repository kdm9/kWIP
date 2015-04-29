#!/usr/bin/env Rscript
t = read.delim('stdin')
m = as.matrix(t[,-1])
#m = as.matrix(t)
image(m)
plot(hclust(as.dist(m)))
