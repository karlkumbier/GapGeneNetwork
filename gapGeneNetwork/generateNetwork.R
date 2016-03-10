library(abind)
library(R.matlab)
library(RColorBrewer)
library(igraph)
library(fields)
library(Rgraphviz)
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

# stability based analysis of spectral clustering: take subsets of nodes and
# split resulting graphs using spectral clustering. How often do genes appear in
# the same leaf nodes
root.dir <- '/Users/Karl/Desktop/LBL/NMF/gapGeneNetwork/networkPlotsPermPW/'
for (i in 1:length(pp.centers)) {
  print(i)
  cur.dir <- paste0(root.dir, 'pp', pp.centers[i], '/')
  dir.create(cur.dir, recursive=TRUE)
  setwd(cur.dir)
  
  gap.genes <- c('Kr', 'kni', 'hb', 'tll', 'hkb', 'gt')
  set.seed(47)
  n.trees <- 100

  # permutation based testing for correlation significance and fdr thresholding
  n.rep <- 250
  alpha <- 0.01
  cors <- localCor(pp.centers[i], pp.neighbors[[i]], tf.data, tf.names) 
  
  #cor.permutation <- replicate(n.rep, localCor(pp.centers[i],
  # stability through subsampling
  n.nodes <- nrow(cors)
  gene.names <- rownames(cors)
  #eps <- 25
  #tf.data.noisy <- replicate(n.trees, tf.data + rexp(length(tf.data), eps), simplify=FALSE)
  #cor.samples <- lapply(tf.data.noisy, function(m) list(localCor(pp.centers[i], pp.neighbors[[i]], m, tf.names)))
  
  node.samples <- replicate(n.trees, unique(sample(1:n.nodes, replace=TRUE)), simplify=FALSE)
  cor.samples <- lapply(node.samples, function(s) list(cors[s, s])) # subsample

  #stability through noise
  #set.seed(47)
  #eps <- 0.25
  #cor.samples <- lapply(1:n.trees, function(s) list(generateNoisyCor(cors, eps)))
  cor.samples <- lapply(cor.samples, function(m) return(list(m)))
  cluster.subsamples <- lapply(cor.samples, spectralSplit, qt.thresh=0.5)
  gene.similarity.output <- geneSimilarity(cluster.subsamples, gene.names)
  gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, gene.names)
 
  pdf(paste0('moduleImagePP', pp.centers[i], '.pdf'))
  gene.sim.mat <- gene.sim.mat / max(gene.sim.mat)
  diag(gene.sim.mat) <- 1
  h.cluster <- hclust(dist(1-gene.sim.mat))
  o <- h.cluster$order
  plotModuleImage(gene.sim.mat, o, col=colorRampPalette(c('white', 'orange', 'red'))(10))
  dev.off()


  pdf(paste0('moduleNetworkPP', pp.centers[i], '.pdf'))
  par(mar=rep(0, 4))
  diag(gene.sim.mat) <- 0
  plotCorGraph(gene.sim.mat * sign(cors), qt.thresh=0.95, emph.nodes=gap.genes)
  dev.off()

  pdf(paste0('corImagePP', pp.centers[i], '.pdf'))
  diag(cors.t) <- 1
  h.cluster <- hclust(dist(1-abs(cors.t)))
  o <- h.cluster$order
  plotModuleImage(abs(cors.t), o, col=colorRampPalette(c('white', 'orange', 'red'))(10))
  dev.off()

  pdf(paste0('corNetworkPP', pp.centers[i], '.pdf'))
  par(mar=rep(0, 4))
  diag(cors.t) <- 0
  plotCorGraph(cors.t, qt.thresh=0.95, emph.nodes=gap.genes)
  dev.off()

}

