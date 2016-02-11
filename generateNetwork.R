library(R.matlab)
library(RColorBrewer)
library(pcalg)
library(igraph)
source('utilities.R')
load('dictFitDataNLS.RData')
mat.data <- readMat('embTemplate.mat')
template <- mat.data$template[,,1]
cols <- brewer.pal(11, 'RdYlBu')

n.dict <- ncol(Dstd) - 1
dict.mat <- array(0, c(16, 32, n.dict))

#for (i in 1:n.dict) {
#  dict.mat[,,i] <- generateImage(template, Dstd[,i])
#  dev.new()
#  im2plot <- t(dict.mat[,,i])
#  im2plot <- im2plot[,seq(ncol(im2plot), 1, by=-1)]
#  image(im2plot)
#}

# Local network analysis for transcription factors expressed in pp1 and pp2
tf.data <- X[, tfInd]
tf.alphas <- alpha[, tfInd]
tf.names <- geneNames[tfInd]

#pp.centers <- c(5:9, 17, 20)
#pp.neighbors <- list(c(4, 6), c(4, 7), c(6, 8), c(7, 9), c(8, 17), c(9, 20), c(17, 20))
# hindgut pp as determined by annotations
pp.centers <- c(17:19)
pp.neighbors <- list(c(16, 18), c(10, 19, 20), c(18, 21))

local.correlation.mats <- list()
local.weighted.expression <- list()
for (i in 1:length(pp.centers)) {
  pp.idx <- pp.centers[i]
  pp.region <- pp.neighbors[[i]]
  expressed <- function(x) any(x > 0.1)
  local.tf.idcs <- which(apply(tf.alphas[pp.region,], 2, expressed))
  local.tf.data <- tf.data[, local.tf.idcs]
  local.tf.names <- tf.names[local.tf.idcs]

  pp.weights <- Dstd[, pp.idx]
  rescale <- function(x) return(x / sum(x))
  pp.weights <- rescale(pp.weights)
  local.weighted.expression[[i]] <- apply(local.tf.data, 2, '*', pp.weights)

  local.correlation <- weightedCor(local.tf.data, pp.weights)
  local.correlation <- mergeDuplicates(local.correlation, local.tf.names)
  local.correlation.mats[[i]] <- local.correlation
}

# Look at spectral clustering for pp7 (based on CG13894 experiment)
# Look at hindgut pp networks
cor.list <- local.correlation.mats[1:3]
threshold <- sapply(cor.list, getQtThreshold, qt.threshold=0.05)
threshold <- c(min(thresholds[1,]), max(thresholds[2,]))
modules <- getLocalModules(cor.list)

network.subsets <- lapply(cor.list, subsetNetwork, genes=modules$genes)
# only consider observations with the same sign across all pp
#network.subsets <- forceSignAgreement(networkSubsets) 


setwd('./plots/hindgut')
for (i in 1:length(network.subsets)) {

  pdf(paste0('localNetwork', i, '.pdf'))
  plotCorGraph(network.subsets[[i]], qt.threshold=0.025)
  dev.off()
}

for (i in 1:length(modules$mod)) {
  mod <- modules$mod[[i]]
  if (length(mod) <= 1) next
  module.subsets <- lapply(network.subsets, subsetNetwork, genes=mod)
  for (j in 1:length(module.subsets)) {

    pdf(paste0('moduleNetworkConstantThresh_m', i, '_pp', j, '.pdf'))
    plotCorGraph(module.subsets[[j]], qt.threshold=0.5, scale=1)
    dev.off()
  }
} 

# consider the module with gap genes
subset.network <- lapply(cor.list[1:2], subsetNetwork, genes=modules$mod[[4]]) 
subset.modules <- getLocalModules(subset.network)

