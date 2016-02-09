library(R.matlab)
library(RColorBrewer)
source('utilities.R')
load('dictFitDataNLS.RData')
mat.data <- readMat('embTemplate.mat')
template <- mat.data$template[,,1]
cols <- rev(brewer.pal(11, 'Spectral'))


genes <- c('Kr', 'gt', 'hb', 'CG13894')

plotGene <- function(g) {

  g.idx <- which(geneNames == g)
  g.vals <- X[, g.idx]
  g.img <- generateImage(template, g.vals)
  im2plot <- t(g.img)
  im2plot <- im2plot[,seq(ncol(im2plot), 1, by=-1)]
  par(mar=c(0, 0, 1, 0))
  image(im2plot, col=cols, main=g, xaxt='n', yaxt='n')
}

par(mfrow=c(2, 2))
sapply(genes, plotGene)
