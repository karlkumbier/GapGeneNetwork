# Plotting functions for embryo images and correlation networks

plotEmb <- function(D, width, height) {
  # Convert a column vector of expression values for plotting
  # args:
  #  D: a column vector of expression levels
  #  width: image width
  #  height: image height 
  img <- matrix(D, nrow=height, ncol=width)
  return(img)

} 
 
generateImage <- function(template, data) {

  # plot expression levels onto a template that defines the embryo ellipse
  # args:
  #  template: a numeric matrix with the same size as the desired image. An
  #   entry of 1 indicates embryo pixels and 0 indicates non-embryo
  #  data: a column vector of values to plot, with length equal to the number of
  #   1s in template
  embryo.idcs <- which(template == 1, arr.ind=TRUE)
  outer.idcs <- which(template == 0, arr.ind=TRUE)
  n.row <- nrow(template)
  n.col <- ncol(template) 
  
  full.img <- matrix(0, nrow=n.row, ncol=n.col)
  for (j in 1:nrow(embryo.idcs)) {    
      full.img[embryo.idcs[j,1], embryo.idcs[j,2]] <- data[j]
  }   
  for (j in 1:nrow(outer.idcs)) {
      full.img[outer.idcs[j,1], outer.idcs[j,2]] <- (-1)
  }   
  
  return(full.img)
} 

plotModuleImage <- function(sim.mat, order, ...) {

  sim.mat <- sim.mat[order, order]
  image(sim.mat, axes=F, ...)
  mtext(text=rownames(sim.mat), side=2, line=0.3, at=seq(0,1,1/nrow(sim.mat)), las=1, cex=25/nrow(sim.mat))
  mtext(text=rownames(sim.mat), side=1, line=0.3, at=seq(0,1,1/nrow(sim.mat)), las=2, cex=25/nrow(sim.mat))
  image.plot(sim.mat, legend.only=TRUE, horizontal=TRUE, smallplot=c(0.05, 0.45, 0.06, 0.09), ...) 

}

getWeight <- function(cor.mat, edge.name) {
  # return graph edge weight associated with a pair of nodes
  nodes <- unlist(strsplit(edge.name, '~'))
  wt <- cor.mat[nodes[1], nodes[2]]
  return(wt)
} 

getWeights <- function(cor.mat, edge.names) {
  # return graph edge weights associated with a list of node pairs
  wts <- sapply(edge.names, getWeight, cor.mat=cor.mat)
  return(wts)
} 

plotCorGraph <- function(cor.mat, qt.thresh=0.5, emph.nodes=NULL, seed=47) {
  
  # plot graph corresponding to a correlation matrix after threshing at
  # specified value or quantile
  # args:
  #  cor.mat: correlation matrix
  #  qt.thresh: indicates a quantile to thresh at
  adj.mat <- generateAdjacency(cor.mat, qt.thresh=qt.thresh)
  graph <- new("graphAM", adjMat=adj.mat, edgemode="undirected")
  node.names <- nodes(graph)
  edge.names <- edgeNames(graph)

  edge.weight <- getWeights(adj.mat * cor.mat, edge.names)
  edge.color <- ifelse(edge.weight > 0, 'blue', 'red') 
  node.color <- rep('green3', length(emph.nodes))
  names(node.color) <- emph.nodes

  eAttrs <- list(color=edge.color)
  nAttrs <- list(fontcolor=node.color)
  #attrs <- list(node=list(fontsize=20, shape="ellipse", fixedsize=FALSE))
  attrs <- list(node=list(fontsize=20, shape="plaintext", fixedsize=FALSE))
  plot(graph, "fdp", attrs=attrs, nodeAttrs=nAttrs, edgeAttrs=eAttrs)
}


plotHeatmap <- function(sim.mat, thresh=NULL, cex=1) {
  
  # plot nicer heatmaps using heatmaps.2 with adjusted margins
  # args:
  #  sim.mat: a similarity matrix obtained through e.g. spectral split
  #  thresh: if specified, all values withing thresh will be set to 0
  #  cex: graphical parameter for label text size
  if (!is.null(thresh)) {
    sim.mat[sim.mat > thresh[1] & sim.mat < thresh[2]] <- 0
  }
  heatmap.2(abs(sim.mat), dendrogram="none", trace="none", key=FALSE, 
                lhei=c(1,7), lwid=c(1, 9), cexRow=cex, cexCol=cex)
}

# Functions used for generating locally weighted correlation matrices
weightedCorVector <- function(x, y, w) {

  # Calculate the weighted correlation between two vectors
  # args:
  #  x: a numeric vector
  #  y: a numeric vector
  #  w: a weight vector that sums to 1, giving the weight of each component in
  #   the correlation calculation
  weightedVar <- function(x, w) {
    weighted.mean <- sum(w * x)
    weighted.var <- sum( w * (x - weighted.mean) ^ 2)
    return(weighted.var)
  }

  weightedCov <- function(x, y, w) {

    weighted.mean.x <- sum(w * x)
    weighted.mean.y <- sum(w * y)
    weighted.cov <- sum(w * (x - weighted.mean.x) * (y - weighted.mean.y))
    return(weighted.cov)
  }

  weighted.cor <- weightedCov(x, y, w) / (sqrt(weightedVar(x, w)) * sqrt(weightedVar(y, w)))
  return(weighted.cor)
}

weightedCor <- function(data.mat, w) {

  # A wraper to compute weighted correlation between all vectors in a matrix
  # args:
  #  data.mat: the data matrix
  #   w: a weight vector that sums to 1, giving the weight of each component in
  #   the correlation calculation
  n.obs <- ncol(data.mat)
  cor.mat <- matrix(0, nrow=n.obs, ncol=n.obs)
  for (i in 1:n.obs) {
    for (j in i:n.obs) {
      weighted.cor <- weightedCorVector(data.mat[,i], data.mat[,j], w)
      cor.mat[i, j] <- weighted.cor
      cor.mat[j, i] <- weighted.cor
    }
  }
  return(cor.mat)
}

localCor <- function(pp.center, pp.neighbors, tf.data, tf.names) {

  pp.region <- c(pp.center, pp.neighbors)
  expressed <- function(x) any(x > 0.1)
  local.tf.idcs <- which(apply(tf.alphas[pp.region,], 2, expressed))
  local.tf.data <- tf.data[, local.tf.idcs]
  local.tf.names <- tf.names[local.tf.idcs]

  pp.weights <- Dstd[, pp.center]
  rescale <- function(x) return(x / sum(x))
  pp.weights <- rescale(pp.weights)
  
  local.cor <- weightedCor(local.tf.data, pp.weights)
  local.cor <- mergeDuplicates(local.cor, local.tf.names)
  return(local.cor)
}
  
#weightedRankCorVector <- function(x, y, w) {
#
#  n <- length(x)
#  ranks.x <- order(x)
#  ranks.y <- order(y)
#  weighted.rank.diff <- sum(w * (ranks.x - ranks.y)^2)
#  weighted.cor <- 1 - 6 * weighted.rank.diff / (n ^ 2 - 1)
#  return(weighted.cor)
#}
#
#weightedRankCor <- function(data.mat, w) {
#  
#  n.obs <- ncol(data.mat)
#  cor.mat <- matrix(0, nrow=n.obs, ncol=n.obs)
#  for (i in 1:n.obs) {
#    for (j in i:n.obs) {
#      weighted.cor <- weightedRankCorVector(data.mat[,i], data.mat[,j], w)
#      cor.mat[i, j] <- weighted.cor
#      cor.mat[j, i] <- weighted.cor
#    }
#  }
#  return(cor.mat)
#}

mergeMax <- function(x) {

  # Return the maximum absolute value of each row in a matrix
  # args:
  #  x: a numeric matrix
  maxCor <- function(x) {
    max.idx <- which.max(abs(x))
    return(x[max.idx])
  }
  row.max <- t(apply(x, 1, maxCor))
  return(c(row.max))
}

mergeDuplicates <- function(cor.mat, node.names) {

  # Merge duplicated observations in a correlation matrix by taking the maximum
  # absolute correlation between duplicated variables and all other variables
  # args:
  #  cor.mat: a numeric matrix of correlations
  #  node.names: a character vector specifying the name of each variable in
  #   cor.mat
  
  if (nrow(cor.mat) != length(node.names)) stop('dimension mismatch')
  
  # determine the indices of all duplicated observations 
  duplicate.idcs <- which(duplicated(node.names))
  duplicate.names <- node.names[duplicate.idcs]
  matchName <- function(target, x) return(which(x == target))
  duplicates <- lapply(duplicate.names, matchName, x=node.names)


  # take maximum absolute correlation across all replicates
  mergeMat <- function(x, idcs) {
    merged <- mergeMax(x[, idcs])
    x[idcs, ] <- matrix(rep(merged, length(idcs)), ncol=ncol(x), byrow=TRUE)
    x[, idcs] <- matrix(rep(merged, length(idcs)), nrow=nrow(x))
    return(x)
  }
  for (i in 1:length(duplicates)) {
    cor.mat <- mergeMat(cor.mat, duplicates[[i]])
  }
  
  # remove duplicated observations
  idcs2remove <- unlist(sapply(duplicates, '[', -1))
  cor.mat <- cor.mat[-idcs2remove, -idcs2remove] 
  rownames(cor.mat) <- node.names[-idcs2remove]
  colnames(cor.mat) <- node.names[-idcs2remove]
  return(cor.mat)
}

# Functions for spectral clustering of gene networks
laplacianSVD <- function(cor.mat, qt.thresh, n.eigen=1) {

  # Calculate the spectral decomposition of the normalized Laplacian for a
  # correlation matrix, where adjacency is determined by threshing
  # correlations at a specific quantile
  # args:
  #  cor.mat: a matrix of correlations
  #  thresh: a numeric value between 0 and 1 specifying the quantile to
  #   set adjacency at
  #  n.eigen: the number of eigenvectors to return, starting with n-1 and
  #  working backwards

  # determine threhsold value based on quantiles of correlation matrix
  if (nrow(cor.mat) == 1) {
    return(1)
  } else{

    # calculate adjacency and normalizing matrices
    adjacency <- generateAdjacency(cor.mat, qt.thresh)
    
    n <- nrow(adjacency)
    if (all(adjacency == 0)) return(rep(0, n))

    col.sums <- colSums(adjacency)
    col.sums.sqrt <- sqrt(col.sums)
    d.sqrt <- col.sums.sqrt * diag(n)

    # remove indices with no neighbors in the adjacency matrix
    idcs2remove <- which(colSums(d.sqrt == 0) == n)
    if (length(idcs2remove) != 0) {
      d.sqrt <- d.sqrt[-idcs2remove, -idcs2remove]
      adjacency <- adjacency[-idcs2remove, -idcs2remove]
      n <- n - length(idcs2remove)
    }                          
                               
    L <- diag(n) - solve(d.sqrt) %*% adjacency %*% solve(d.sqrt)
    e <- eigen(L)              
    e.vectors <- matrix(e$vectors[,(n-n.eigen):(n-1)], nrow=n)
    
    # add back removed observations with value of 0
    if (length(idcs2remove) != 0) {
      temp <- numeric(nrow(cor.mat))
      temp[-idcs2remove] <- e.vectors
      e.vectors <- matrix(temp, nrow=nrow(cor.mat))
    }
    rownames(e.vectors) <- rownames(cor.mat)
    return(e.vectors)
    }
}

spectralSplit <- function(cor.list, qt.thresh=0.25) { 

  # Recursively split graph until cluster sizes are below a specified value.
  # Stop splitting clusters if they fall below a certain size
  # args:
  #  cor.list: a list of correlation matrices, should be length one if trying to
  #   split a single network
  #  min.size: integer that determines when to stop splitting a given cluster
  #  max.size: integer that determines when to exit the function. Cuts will
  #   continue until all clusters are below max.size 
  # qt.thresh: a numeric value between 0 and 1 specifying the quantile to
  #   set adjacency at
  set.seed(47)

  # split laplacian eigenvector into two groups and determine modularity of the
  # clustering
  vec <- lapply(cor.list, laplacianSVD, qt.thresh=qt.thresh)  
  clusterVectors <- function(v) {
    if(length(v) > 2){ 
      return(kmeans(v, center=2, nstart=10)$cluster)
    } else if (length(v) == 2) { 
      return(c(1, 2))
    } else {
      return(1)
    }
  }
  clusters <- lapply(vec, clusterVectors)  
  adj.list <- lapply(cor.list, function(c) generateAdjacency(c, qt.thresh=qt.thresh))
  mod <- mapply(function(a, c) modularity(a, c), adj.list, clusters)
  # TODO: gene sometimes not appearing in cluster tree after selection, fix this
  # sif modularity is positive, split genes according to spectral clustering
  genes <- lapply(cor.list, rownames)
  splitGenes <- function(clusters, gene.names, mod) {

    n.clusters <- length(unique(clusters))
    if (mod > 0) {
      lapply(1:n.clusters, function(c) gene.names[clusters==c])
    } else {
      list(gene.names)
    }
  }
  gene.clusters <- mapply(splitGenes, clusters, genes, mod, SIMPLIFY=FALSE)
  gene.clusters <- unlist(gene.clusters, recursive=FALSE)

  # if modularity is positive, split correlation matrices, else return the
  # original matrix
  splitCors <- function(clusters, cors, mod) {
    
    n.clusters <- length(unique(clusters))
    if (mod > 0) {
      cors <- lapply(1:n.clusters, function(c) as.matrix(cors[clusters==c, clusters==c]))   
    } else {
      cors <- list(cors)
    }
    return(cors)    
  }
  cor.list <- mapply(splitCors, clusters, cor.list, mod, SIMPLIFY=FALSE)
  cor.list <- unlist(cor.list, recursive=FALSE)

  if (any(mod > 0)) {
    return(spectralSplit(cor.list, qt.thresh))
  } else {
    return(gene.clusters)
  }
}

geneSimilarity <- function(cluster.tree, genes, set.size=2) {

  gene.set <- combn(genes, set.size, simplify=FALSE)

  # determine which trees contain both genes in a given pair
  setInTree <- function(set, tree) return(all(set %in% unlist(tree)))
  containsSet <- function(set) which(sapply(cluster.tree, setInTree, set=set))
  set.trees <- lapply(gene.set, containsSet) 
  
  # for a given pair, determine the proportion of times the pairs occur in the
  # same cluster
  sameLeaf <- function(set, tree) {
    idcs <- sapply(set, function(s) {
              which(sapply(tree, function(t) s %in% t))
    })
    return(length(unique(idcs)) == 1)
  }
  sameLeafProportion <- function(set, shared.idcs) {
    joint.trees <- cluster.tree[shared.idcs]
    if (length(joint.trees) == 0) return(0) # genes in set never co-occur
    proportion <- mean(sapply(joint.trees, sameLeaf, set=set))
    return(proportion)
  }
  set.proportions <- mapply(sameLeafProportion, set=gene.set, shared.idcs=set.trees) 

  return(list(gene.set=gene.set, set.prop=set.proportions))
}

 
generateSimilarityMatrix <- function(gene.sim, genes) {
  # function for generating a similarity matrix for gene pairs based on the
  # output of gene similarity
  # args:
  #  gene.sim: output from geneSimilarity
  #  genes: a vector of gene names for which to determine pairwise similarity 
  pair.list <- gene.sim$gene.set
  pair.similarity <- gene.sim$set.prop 
  similarity.matrix <- matrix(0, nrow=length(genes), ncol=length(genes)) 
  rownames(similarity.matrix) <- genes         
  colnames(similarity.matrix) <- genes
  for (i in 1:length(pair.list)) {
    sim <- pair.similarity[i]
    g1 <- pair.list[[i]][1]
    g2 <- pair.list[[i]][2]
    similarity.matrix[g1, g2] <- sim
    similarity.matrix[g2, g1] <- sim
  }
  return(similarity.matrix)
}

# General network utility functions for dealing with correlation matrices
subsetNetwork <- function(cor.mat, genes) {

  # subset a correlation matrix to indicated genes
  cor.names <- rownames(cor.mat)
  gene.idcs <- cor.names %in% genes
  return(cor.mat[gene.idcs, gene.idcs])
}

getQtThreshold <- function(cor.mat, qt.thresh, remove.zero=FALSE) {

  # Calculate a specified quantile for a correlation matrix
  # args:
  #  cor.mat: a matrix of correlations
  #  qt.thresh: a numeric value between 0 and 1 specifying the quantile

  qt.thresh <- max(qt.thresh, 1-qt.thresh)
  cor.vector <- cor.mat[upper.tri(cor.mat, diag=FALSE)]
  if (remove.zero) cor.vector <- cor.vector[cor.vector != 0]
  thresh <- quantile(cor.vector, c(1-qt.thresh, qt.thresh))
  return(thresh)
}


generateNoisyCor <- function(cor.mat, eps) {

  # generate a correlation matrix based on iid gaussian variables that can be
  # used to add noise to an estimated correlation matrix
  n <- nrow(cor.mat)
  noise <- matrix(rnorm(n ^ 2), nrow=n)
  noise.cor <- cor(noise) * eps
  diag(noise.cor) <- 0
  return(cor.mat + noise.cor)
}

generateAdjacency <- function(cor.mat, qt.thresh) {

  # threshold a correlation matrix at specified quantile and generate
  # corresponding adjacency matrix
  thresh <- getQtThreshold(cor.mat, qt.thresh)
  cor.mat[cor.mat >= thresh[1] & cor.mat <= thresh[2]] <- 0
  cor.mat[cor.mat != 0] <- 1
  diag(cor.mat) <- 0
  return(cor.mat)
}

modularity <- function(adjacency, label) {
 
  # calculate the modularity associated with a specific labeling of the
  # adjacency matrix 
  if (all(adjacency == 0)) return(0)

  lab.vals <- unique(label)
  if (length(lab.vals) > 2) stop('maximum of two classes')
  lab1.idcs <- which(label == lab.vals[1])
  lab2.idcs <- which(label == lab.vals[2])
  label[lab1.idcs] <- 1
  label[lab2.idcs] <- -1
  
  m <- sum(adjacency) 
  degree <- rowSums(adjacency)
  joint.degree <- degree %*% t(degree) / (2 * m)
  modularity <- (1 / (4*m)) * (t(label) %*% (adjacency - joint.degree) %*% label)
  return(modularity)
}

sortSets <- function(gene.sims) {

  set <- gene.sims$gene.set
  sort.vals <- gene.sims$set.prop
  order.idcs <- order(abs(sort.vals), decreasing=TRUE)
  return(list(set=set[order.idcs], vals=sort.vals[order.idcs]))
}

getGeneNetwork <- function(cor.mat, thresh, genes) {
    
  adjacency <- generateAdjacency(cor.mat, thresh=thresh)
  adjacency <- adjacency * sign(cor.mat)
  gene.idcs <- which(rownames(adjacency) %in% genes)
  gene.network <- adjacency[gene.idcs, gene.idcs]
  return(gene.network)
}

setModules <- function(sim.mat, clusters) {

  if (length(clusters) != nrow(sim.mat)) stop('dimension mismatch')

  n.clusters <- length(unique(clusters))
  cutMat <- function(sim.mat, c, clusters) {
   cluster.idcs <- which(clusters==c)
   sim.mat[cluster.idcs, -cluster.idcs] <- 0
   sim.mat[-cluster.idcs, cluster.idcs] <- 0
   return(sim.mat)
  }

  for (i in 1:n.clusters) {
    sim.mat <- cutMat(sim.mat, i, clusters)
  }
  return(sim.mat)
}
