library(abind)
library(RColorBrewer)
library(fields)
library(Rgraphviz)

root.dir <- '/Users/Karl/Desktop/LBL/NMF/gapGeneNetwork/'
setwd(root.dir)
source('./utilities.R')
load('dictFitDataNLS.RData')

n.dict <- ncol(Dstd) - 1
dict.mat <- array(0, c(16, 32, n.dict))

# Local network analysis for transcription factors expressed in segmentation stripes
tf.data <- X[, tfInd]
tf.alphas <- alpha[, tfInd]
tf.names <- geneNames[tfInd]

pp.centers <- c(5:9, 17, 20)
pp.neighbors <- list(c(4, 6), c(4, 7), c(6, 8), c(7, 9), c(8, 17), c(9, 20), c(17, 20))


for (i in 1:length(pp.centers)) {
  print(i)
  #cur.dir <- paste0(root.dir, 'pp', pp.centers[i], '/')
  #dir.create(cur.dir, recursive=TRUE)
  #setwd(cur.dir)
  
  gap.genes <- c('Kr', 'kni', 'hb', 'tll', 'hkb', 'gt')
  set.seed(47)
  n.trees <- 100
  
  cors <- localCor(pp.centers[i], pp.neighbors[[i]], tf.data, tf.names) 
  adj.mat <- generateAdjacency(cors, qt.thresh=0.9)
  zero.idcs <- which(rowSums(adj.mat) == 0)
  
  if (length(zero.idcs)) {
    adj.mat <- adj.mat[-zero.idcs, -zero.idcs]
    cors <- cors[-zero.idcs, -zero.idcs]
  }
  
  transition.mat <- generateTransitionMatrix(adj.mat)
  diffusion.mat <- diffusionDist(transition.mat, k=25)
  
  diffusion.mat <- diffusion.mat / max(diffusion.mat )
  threshold <- getQtThreshold(diffusion.mat, qt.thresh=0.5, remove.zero=TRUE)[1]
  diffusion.mat[diffusion.mat > threshold] <- 1
  diffusion.mat <- 1 - diffusion.mat
  heatmap(diffusion.mat)
  
  
  
  # stability through subsampling
  n.nodes <- nrow(cors)
  gene.names <- rownames(cors)
  node.samples <- replicate(n.trees, unique(sample(1:n.nodes, replace=TRUE)), simplify=FALSE)
  cor.samples <- lapply(node.samples, function(s) list(cors[s,s])) # subsample
  
  cluster.subsamples <- lapply(cor.samples, spectralSplit, qt.thresh=0.7)
  gene.similarity.output <- geneSimilarity(cluster.subsamples, gene.names)
  
  if(false) {
    gene.sim.mat <- matrix(0, nrow=n.nodes, ncol=n.nodes)
    rownames(gene.sim.mat) <- gene.names
    colnames(gene.sim.mat) <- gene.names
    sim.idcs <- which(gene.similarity.output$set.prop != 0)
    for (i in 1:length(sim.idcs)) {
      ii <- sim.idcs[i]
      temp.genes <- gene.similarity.output$gene.set[[ii]]
      temp.sim <- gene.similarity.output$set.prop[[ii]]
      gene.comb <- combn(temp.genes, 2)
      for (j in 1:ncol(gene.comb)) {
        gene.sim.mat[gene.comb[1, j], gene.comb[2, j]] <- temp.sim
        gene.sim.mat[gene.comb[2, j], gene.comb[1, j]] <- temp.sim
      }
    }
  }
  gene.sim.mat <- generateSimilarityMatrix(gene.similarity.output, gene.names)
  
  pdf(paste0('moduleImagePP', pp.centers[i], '.pdf'))
  gene.sim.mat <- 1 - diffusion.mat
  diag(gene.sim.mat) <- 1
  h.cluster <- hclust(dist(1-gene.sim.mat))
  o <- h.cluster$order
  plotModuleImage(gene.sim.mat, o, col=colorRampPalette(c('white', 'orange', 'red'))(10))
  dev.off()
  
  
  pdf(paste0('moduleNetworkPP', pp.centers[i], '.pdf'))
  #gene.sim.mat <- 1 - diffusion.mat
  par(mar=rep(0, 4))
  diag(gene.sim.mat) <- 0
  plotCorGraph(gene.sim.mat * sign(cors), qt.thresh=0.9, emph.nodes=gap.genes)
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

