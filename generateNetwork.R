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

pp.centers <- c(5:9, 17, 20)
pp.neighbors <- list(c(4, 6), c(4, 7), c(6, 8), c(7, 9), c(8, 17), c(9, 20), c(17, 20))

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
cor.list <- local.correlation.mats[1:3]
threshold <- sapply(cor.list, getQtThreshold, qt.threshold=0.05)
threshold <- c(min(thresholds[1,]), max(thresholds[2,]))
modules <- getLocalModules(cor.list)

network.subsets <- lapply(cor.list, subsetNetwork, genes=modules$genes)
setwd('./plots')
for (i in 1:length(network.subsets)) {

  pdf(paste0('localNetwork', i, '.pdf'))
  plotCorGraph(network.subsets[[i]], threshold=threshold)
  dev.off()
}

for (i in 1:length(modules$mod)) {
  mod <- modules$mod[[i]]
  print(i)
  if (length(mod) <= 1) next
  module.subsets <- lapply(cor.list, subsetNetwork, genes=mod)
  for (j in 1:length(module.subsets)) {

    print(j)
    pdf(paste0('moduleNetworkConstantThresh_m', i, '_pp', j, '.pdf'))
    plotCorGraph(module.subsets[[j]], threshold=threshold, scale=2)
    dev.off()
  }
} 















# PC algorithm analysis
# don't want n=405 since pixels are highly dependent and many are 0
local.pc.fit <- list()
for (i in 1:length(pp.centers)) {
  labs <- rownames(local.correlation.mats[[i]])
  local.pc.fit[[i]] <- pc(suffStat=list(C=local.correlation.mats[[i]], n=405), indepTest=gaussCItest, alpha=0.1, labels=labs)
}

# analysis of graph data
pdf('testGraph.pdf')
plot(fit)
dev.off()

# TODO: for all nodes in any local correlation newtorks, compare edges for all
# local networks where the node occurs
local.network.genes <- sapply(local.correlation.mats, rownames)
local.network.genes <- unique(unlist(local.genes))

all.children <- list()
for (g in local.network.genes) {

  # check each network to determine if gene is expressed in the spatial region
  g.children <- list()
  for (i in 1:length(local.pc.fit)) {
    graph <- local.pc.fit[[i]]@graph
    graph.edges <- graph@edgeL
    graph.nodes <- graph@nodes

    g.idx <- which(graph.nodes == g)
    #TODO: differentiate between case where node is in graph but empty and where
    #node is not in graph
    if (length(g.idx) > 0) {
      if (length(unlist(graph.edges[[g.idx]])) > 0) { 
        g.children[[i]] <- graph.nodes[unlist(graph.edges[[g.idx]])]
      } else {
        g.children[[i]] <- 'NONE'
      }
    } else {
      g.children[[i]] <- numeric(0)
    }
  }
  all.children[[g]] <- g.children
}





# Determine genes that are in all local correlation networks for segmentation
# stripes
local.genes <- sapply(local.correlation.mats, rownames)
global.genes <- unique(unlist(local.genes))
for (i in 1:length(local.genes)) {
  global.genes <- intersect(global.genes, local.genes[[i]])
}

local.correlation.subset <- array(0, c(length(global.genes), length(global.genes), length(local.genes)))
for (i in 1:length(local.genes)) {
  local.correlation.subset[,,i] <- local.correlation.mats[[i]][global.genes, global.genes]
}

temp <- apply(local.correlation.subset, MARGIN=c(1, 2), mean)
rownames(temp) <- global.genes
colnames(temp) <- global.genes
heatmap(temp)

cex.lab <- 0.9
pdf('heatmap.pdf')
heatmap(local.correlation, cexRow=cex.lab, cexCol=cex.lab)
dev.off()
