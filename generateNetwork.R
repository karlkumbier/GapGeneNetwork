library(R.matlab)
library(RColorBrewer)
library(igraph)
source('utilities.R')
load('dictFitDataNLS.RData')
mat.data <- readMat('embTemplate.mat')
template <- mat.data$template[,,1]
cols <- brewer.pal(11, 'RdYlBu')

n.dict <- ncol(Dstd) - 1
dict.mat <- array(0, c(16, 32, n.dict))

# Local network analysis for transcription factors expressed in pp1 and pp2
tf.data <- X[, tfInd]
tf.alphas <- alpha[, tfInd]
tf.names <- geneNames[tfInd]

pp.centers <- c(5:9, 17, 20)
pp.neighbors <- list(c(4, 6), c(4, 7), c(6, 8), c(7, 9), c(8, 17), c(9, 20), c(17, 20))
# hindgut pp as determined by annotations
#pp.centers <- c(17:19)
#pp.neighbors <- list(c(16, 18), c(10, 19, 20), c(18, 21))

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

cors.pp7 <- local.correlation.mats[[3]]
n.nodes <- nrow(cors.pp7)
set.seed(47)

# stability based analysis of spectral clustering: take subsets of nodes and
# split resulting graphs using spectral clustering. How often do genes appear in
# the same leaf nodes
n.trees <- 100

# stability through subsampling
node.samples <- replicate(n.trees, unique(sample(1:n.nodes, replace=TRUE)), simplify=FALSE)
cor.samples <- lapply(node.samples, function(s) list(cors.pp7[s, s])) # subsample

#stability through noise
#set.seed(47)
#eps <- 0.5
#cor.samples <- lapply(1:n.trees, function(s) list(generateNoisyCor(cors.pp7, eps)))

# TODO: get rid of n.splits
cluster.subsamples <- lapply(cor.samples, spectralSplit, qt.thresh=0.25)
gene.similarity.output <- geneSimilarity(cluster.subsamples, rownames(cors.pp7))
gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, rownames(cors.pp7))

gene.names <- rownames(cors.pp7)
gene.pairs <- combn(gene.names, 2, simplify=FALSE)
cor.vals <- sapply(gene.pairs, function(p) cors.pp7[p[1], p[2]])
cors.ordered <- sortSets(list(gene.sets=gene.pairs, set.prop=cor.vals))
sim.ordered <- sortSets(gene.similarity.output) 

gene.names <- rownames(cors.pp7)
gene.pairs <- combn(gene.names, 2, simplify=FALSE)
cor.vals <- sapply(gene.pairs, function(p) cors.pp7[p[1], p[2]])
sorted.spectral.pairs <- sortSets(gene.similarity.output)
sorted.cor.pairs <- sortSets(list(gene.set=gene.pairs, set.prop=cor.vals))


h.cluster <- hclust(dist(1-gene.sim.mat))
clusters <- cutree(h.cluster, k=4)
gene.clusters <- rownames(cors.pp7)[clusters == 2]
cluster.subsamples <- lapply(cor.samples, spectralSplit, n.splits=2, qt.threshold=0.1)
gene.sim.4 <- geneSimilarity(cluster.subsamples, gene.clusters,  set.size=4) 

# Comparison of spectral clustering and correlation based results
library(gplots)
diag(gene.sim.mat) <- 1
plot.dir <- './corrlationComparison/'
dir.create(plot.dir)

pdf(paste0(plot.dir, 'localSpectral.pdf')) 
plotHeatmap(gene.sim.mat)
dev.off()

pdf(paste0(plot.dir, 'localCor.pdf')) 
plotHeatmap(cors.pp7)
dev.off()

pdf(paste0(plot.dir, 'localCor_90.pdf')) 
thrsh <- getQtThreshold(cors.pp7, 0.9)
plotHeatmap(cors.pp7, threshold=thrsh)
dev.off()

pdf(paste0(plot.dir, 'localCor_70.pdf'))  
thrsh <- getQtThreshold(cors.pp7, 0.7)
plotHeatmap(cors.pp7, threshold=thrsh)
dev.off()

# Comparison of global and local results
global.cors <- cor(tf.data)
global.cors <- mergeDuplicates(global.cors, tf.names)

pdf(paste0(plot.dir, 'globalCor.pdf')) 
plotHeatmap(global.cors, cex=0.4)
dev.off()

pdf(paste0(plot.dir, 'globalCor_90.pdf')) 
thrsh <- getQtThreshold(global.cors, 0.9)
plotHeatmap(global.cors, threshold=thrsh, cex=0.4)
dev.off()

global.cors.thrsh <- global.cors
global.cors.thrsh[global.cors > thrsh[1] & global.cors < thrsh[2]] <- 0


set.seed(47)
eps <- 0.05
cor.samples <- lapply(1:n.trees, function(s) list(generateNoisyCor(global.cors, eps)))

cluster.subsamples <- lapply(cor.samples, spectralSplit, n.splits=2, qt.threshold=0.25)
gene.similarity.output <- geneSimilarity(cluster.subsamples, rownames(global.cors))
gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, rownames(global.cors))
