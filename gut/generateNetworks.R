library(R.matlab)
source('utilities.R')
load('dictFitDataNLS.RData')

tf.idcs <- which(geneNames %in% tfNames)
tf.data <- X[, tf.idcs]
tf.alphas <- alpha[, tf.idcs]
tf.names <- geneNames[tf.idcs]

pp.centers <- c(12, 2)
pp.neighbors <- list(c(4, 7, 9), c(6, 8))
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

cors <- local.correlation.mats[[1]]
genes <- rownames(cors)
eps <- 0.1 
n.tree <- 50
noisy.cors <- replicate(n.tree, list(generateNoisyCor(cors, eps)), simplify=FALSE)
trees <- lapply(noisy.cors, spectralSplit, qt.thresh=0.25)
gene.sims <- geneSimilarity(trees, genes)
gene.sims.mat <- generateSimilarityMatrix(gene.sims, genes)
sorted.spectral.pairs <- sortSets(gene.sims)

gene.pairs <- combn(genes, 2, simplify=FALSE)
pair.cors <- sapply(gene.pairs, function(g) cors[g[1], g[2]])
sorted.cor.pairs <- sortSets(list(gene.set=gene.pairs, set.prop=pair.cors))

cor.adjacency <- generateAdjacency(cors, qt.thresh=0.9)
cor.adjacency <- cor.adjacency * sign(cors)
spectral.adjacency <- gene.sims.mat
spectral.adjacency <- spectral.adjacency * sign(cors)
plotCorGraph(cor.adjacency, qt.thres=0.9)
dev.new()
plotCorGraph(spectral.adjacency, qt.thresh = 0.9)
