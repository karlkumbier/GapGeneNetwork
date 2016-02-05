library(R.matlab)
library(RColorBrewer)
source('utilities.R')
load('dictFitDataNLS.RData')
mat.data <- readMat('embTemplate.mat')
template <- mat.data$template[,,1]
cols <- rev(brewer.pal(11, 'Spectral'))

n.dict <- ncol(Dstd) - 1
dict.mat <- array(0, c(16, 32, n.dict))

for (i in 1:n.dict) {
  dict.mat[,,i] <- generateImage(template, Dstd[,i])
  dev.new()
  im2plot <- t(dict.mat[,,i])
  im2plot <- im2plot[,seq(ncol(im2plot), 1, by=-1)]
  image(im2plot)
}

# Local network analysis for transcription factors expressed in pp1 and pp2
tf.data <- X[, tfInd]
tf.alphas <- alpha[, tfInd]
tf.names <- geneNames[tfInd]

pp.idx <- 7
pp.region <- 6:8 #first three segmentation stripes
expressed <- function(x) any(x > 0.1)
local.tf.idcs <- which(apply(tf.alphas[pp.region,], 2, expressed))
local.tf.data <- tf.data[, local.tf.idcs]
local.tf.names <- tf.names[local.tf.idcs]

pp.weights <- Dstd[, pp.idx]
rescale <- function(x) return(x / sum(x))
pp.weights <- rescale(pp.weights)

local.correlation <- weightedCor(local.tf.data, pp.weights)
local.correlation <- mergeDuplicates(local.correlation, local.tf.names)

cex.lab <- 0.9
pdf('heatmap.pdf')
heatmap(local.correlation, cexRow=cex.lab, cexCol=cex.lab)
dev.off()
