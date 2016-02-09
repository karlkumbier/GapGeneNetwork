library(R.matlab)
library(RColorBrewer)
library(pcalg)
source('utilities.R')
load('dictFitDataNLS.RData')
mat.data <- readMat('embTemplate.mat')
template <- mat.data$template[,,1]
cols <- brewer.pal(11, 'RdYlGn')

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
cors <- local.correlation.mats[[1]]
spectral.output5 <- spectralVectors(cors, threshold=0.25, n.eigen=1)$sv
spectral.cluster5 <- kmeans(spectral.output5, centers=2)
clusters5 <- spectral.cluster5$cluster

cors <- local.correlation.mats[[2]]
spectral.output6 <- spectralVectors(cors, threshold=0.25, n.eigen=1)$sv
spectral.cluster6 <- kmeans(spectral.output6, centers=2)
clusters6 <- spectral.cluster6$cluster

cors <- local.correlation.mats[[3]]
spectral.output7 <- spectralVectors(cors, threshold=0.25, n.eigen=1)$sv
spectral.cluster7 <- kmeans(spectral.output7, centers=2)
clusters7 <- spectral.cluster7$cluster

jointly.expressed <- intersect(rownames(spectral.output5), rownames(spectral.output6))
jointly.expressed <- intersect(rownames(spectral.output7), jointly.expressed)
s5.idcs <- which(rownames(spectral.output5) %in% jointly.expressed)
s6.idcs <- which(rownames(spectral.output6) %in% jointly.expressed)
s7.idcs <- which(rownames(spectral.output7) %in% jointly.expressed)
v5 <- spectral.output5[s5.idcs,]
v6 <- spectral.output6[s6.idcs,]
v7 <- spectral.output7[s7.idcs,]
v <- cbind(v5, v6, v7)


cols <- brewer.pal(11, 'RdYlBu')
heatmap(v, labRow=rownames(pp.cors), col=cols)
spectral.clustering <- kmeans(spectral.vectors, centers=3)
clusters <- spectral.clustering$cluster
genes <- rownames(pp.cors)[-idcs2remove]










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
