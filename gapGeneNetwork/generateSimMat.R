library(abind)
library(RColorBrewer)
library(fields)
library(Rgraphviz)
library(glasso)
setwd('~/Desktop/LBL/NMF/gapGeneNetwork/')
source('./utilities.R')
load('dictFitDataNLS.RData')

n.dict <- ncol(Dstd) - 1
dict.mat <- array(0, c(16, 32, n.dict))

# Local network analysis for transcription factors expressed in pp1 and pp2
tf.data <- X[, tfInd]
tf.alphas <- alpha[, tfInd]
tf.names <- geneNames[tfInd]

pp.centers <- c(1:21)
pp.neighbors <- list(c(2:4), c(1, 4, 5), c(4, 14), c(1, 2, 3, 13), c(4, 6), 
                     c(4, 5, 7), c(6, 8), c(7, 9), c(8, 17), c(11), c(10, 12), 
                     c(11, 13), c(12, 14),c(13, 15), c(14, 16), c(15), 
                     c(9, 18, 20), c(17), c(20), c(17, 19, 21), c(20))


for (i in 1:length(pp.centers)) { 
  print(paste('PP', i))
  gap.genes <- c('Kr', 'kni', 'hb', 'tll', 'hkb', 'gt')
  set.seed(47)
  n.trees <- 50
  
  cors <- localCor(pp.centers[i], pp.neighbors[[i]], tf.data, tf.alphas, tf.names, rank=TRUE)
  #cors <- localCor(pp.centers[i], pp.neighbors[[i]], X, alpha, geneNames)
  n.nodes <- nrow(cors)
  gene.names <- rownames(cors)
  
  node.samples <- replicate(n.trees, unique(sample(1:n.nodes, replace=TRUE)), simplify=FALSE)
  cor.samples <- lapply(node.samples, function(s) list(cors[s, s])) # subsample
  
  
  #cluster.subsamples <- lapply(cor.samples, spectralSplit, qt.thresh=0.75)
  #gene.similarity.output <- geneSimilarity(cluster.subsamples, gene.names)
  #gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, gene.names)
  
  #par(mar=rep(0, 4))
  #diag(gene.sim.mat) <- 0
  #plotCorGraph(gene.sim.mat * sign(cors), qt.thresh=0.95, emph.nodes=gap.genes)
  
  root.dir <- paste0('simMatsInvTF4-6/pp', pp.centers[i], '/')
  dir.create(root.dir, recursive=TRUE)
  #thresholdMats <- function(mat, thresh) mat * generateAdjacency(mat, thresh)
  #thresh <- seq(.6, 0.99, by=0.01)
  rhos <- seq(0.5, 0.95, by=0.05)
  
  a <- Sys.time()
  inv.mats <- lapply(rhos, function(r) {
    inv.mat <- glasso(cors, rho=r, approx=TRUE)$wi
    rownames(inv.mat) <- gene.names
    colnames(inv.mat) <- gene.names
    diag(inv.mat) <- 0
    return(inv.mat)
  })
  b <- Sys.time()
  print(b-a)
  
  #thresh.cor.mats <- lapply(thresh, thresholdMats, mat=cors)
  #thresh.cor.inv <- lapply(thresh, thresholdMats, mat=cors.inv)
  
  saveMat <- function(mat, thrsh) save(file=paste0(root.dir, 'simMat', 10*thrsh, '.Rdata'), mat)
  for (i in 1:length(rhos)) {
    saveMat(inv.mats[[i]], rhos[i])
  }
}
