library(R.matlab)
library(RColorBrewer)
library(igraph)
source('../utilities.R')
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
  pp.region <- c(pp.idx, pp.neighbors[[i]])
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

# stabiliti y based analysis of spectral clustering: take subsets of nodes and
# split resulting graphs using spectral clustering. How often do genes appear in
# the same leaf nodes

root.dir <- '/Users/Karl/Desktop/LBL/NMF/gapGeneNetwork/networkPlotsOrgan/'
for (2 in 1:length(local.correlation.mats)) {
  print(i)
  cur.dir <- paste0(root.dir, 'pp', pp.centers[i], '/')
  dir.create(cur.dir, recursive=TRUE)
  setwd(cur.dir)
  
  gap.genes <- c('Kr', 'kni', 'hb', 'tll', 'hkb', 'gt')

  set.seed(47)
  n.trees <- 100

  # stability through subsampling
  cors <- local.correlation.mats[[i]]
  n.nodes <- nrow(cors)
  gene.names <- rownames(cors)
  node.samples <- replicate(n.trees, unique(sample(1:n.nodes, replace=TRUE)), simplify=FALSE)
  cor.samples <- lapply(node.samples, function(s) list(cors[s, s])) # subsample
  #cor.samples <- c(cor.samples, lapply(node.samples, function(s) list(cors[s, s])))

  #stability through noise
  #set.seed(47)
  #eps <- 0.25
  #cor.samples <- lapply(1:n.trees, function(s) list(generateNoisyCor(cors, eps)))
  cluster.subsamples <- lapply(cor.samples, spectralSplit, qt.thresh=0.25)
  gene.similarity.output <- geneSimilarity(cluster.subsamples, gene.names)
  gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, gene.names)
  gene.sim.mat <- gene.sim.mat / max(gene.sim.mat)
  diag(gene.sim.mat) <- 1

  pdf(paste0('moduleImagePP', pp.centers[i], '.pdf'))
  h.cluster <- hclust(dist(1-gene.sim.mat))
  o <- h.cluster$order
  plotModuleImage(gene.sim.mat, o, col=colorRampPalette(c('white', 'orange', 'red'))(10))
  dev.off()

  pdf(paste0('moduleNetworkPP', pp.centers[i], '.pdf'))
  par(mar=rep(0, 4))
  plotCorGraph(gene.sim.mat * sign(cors), qt.thresh=0.9, seed=10)
  dev.off()

  # local networks for inverse and marginal
  pdf(paste0('corNetworkPP', pp.centers[i], '.pdf'))
  cor.adj <- generateAdjacency(cors, qt.thresh=0.95)
  plotNetwork(cor.adj * cors, gap.genes)
  dev.off()

  pdf(paste0('gapNetwork_090_PP', pp.centers[i], '.pdf'))
  cor.adj <- generateAdjacency(cors, qt.thresh=0.9)
  plotNetworkSubset(cor.adj * cors, gap.genes)
  dev.off()

  pdf(paste0('corNetworkInvPP', pp.centers[i], '.pdf'))
  cor.adj <- generateAdjacency(solve(cors), qt.thresh=0.95)
  plotNetwork(cor.adj * solve(cors), gap.genes)
  dev.off()

  pdf(paste0('corNetworkInvi_gap_PP', pp.centers[i], '.pdf'))
  cor.adj <- generateAdjacency(solve(cors), qt.thresh=0.9)
  plotNetworkSubset(cor.adj * solve(cors), gap.genes)
  dev.off()
}

plotCorGraph(gene.sim.mat * cors, qt.thresh=0.95, emph.nodes=gap.genes, seed=50)
