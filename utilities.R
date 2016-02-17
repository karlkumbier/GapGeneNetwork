
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

spectralPreprocess <- function(cor.mat, threshold, n.eigen) {

  # Calculate the spectral decomposition of the normalized Laplacian for a
  # correlation matrix, where adjacency is determined by thresholding
  # correlations at a specific quantile
  # args:
  #  cor.mat: a matrix of correlations
  #  threshold: a numeric value between 0 and 1 specifying the quantile to
  #   set adjacency at
  #  n.eigen: the number of eigenvectors to return, starting with n-1 and
  #  working backwards

  # determine threhsold value based on wuantiles of correlation matrix
  if (nrow(cor.mat) == 1) {
    return(list(sv=1, e=1))
  } else{
    qt.threshold <- getQtThreshold(cor.mat, threshold) 

    # calculate adjacency and normalizing matrices
    adjacency <- cor.mat
    diag(adjacency) <- 0
    adjacency[adjacency <= qt.threshold[1] | adjacency >= qt.threshold[2]] <- 1
    adjacency[adjacency != 1] <- 0

    n <- nrow(adjacency)
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
    spectral.vectors <- matrix(e$vectors[,(n-n.eigen):(n-1)], nrow=n)
    
    # add back removed observations with value of 0
    if (length(idcs2remove) != 0) {
      temp <- numeric(nrow(cor.mat))
      temp[-idcs2remove] <- spectral.vectors
      spectral.vectors <- matrix(temp, nrow=nrow(cor.mat))
    }
    rownames(spectral.vectors) <- rownames(cor.mat)
    return(list(sv=spectral.vectors, e=e$values))
    }
}

spectralSplit <- function(cor.list, n.splits=1, qt.threshold=0.25) { 

  # Recursively split graph until cluster sizes are below a specified value.
  # Stop splitting clusters if they fall below a certain size
  # args:
  #  cor.list: a list of correlation matrices, should be length one if trying to
  #   split a single network
  #  min.size: integer that determines when to stop splitting a given cluster
  #  max.size: integer that determines when to exit the function. Cuts will
  #   continue until all clusters are below max.size 
  # qt.threshold: a numeric value between 0 and 1 specifying the quantile to
  #   set adjacency at
  set.seed(47)
  sv <- lapply(cor.list, function(c) spectralPreprocess(c, threshold=qt.threshold, 1)$sv)
   
  clusterVectors <- function(v) {
    ifelse(length(v) > 2, 
           return(kmeans(v, center=2, nstart=10)$cluster), 
           return(rep(1, length(v))))
  }
  clusters <- lapply(sv, clusterVectors)
  genes <- lapply(cor.list, rownames)
  
  splitGenes <- function(clusters, gene.names) {
    lapply(1:length(unique(clusters)), function(c) gene.names[clusters==c])
  }
  gene.clusters <- mapply(splitGenes, clusters, genes, SIMPLIFY=FALSE)
  gene.clusters <- unlist(gene.clusters, recursive=FALSE)

  splitCors <- function(clusters, cors) {
    lapply(1:length(unique(clusters)), function(c) as.matrix(cors[clusters==c, clusters==c]))
  }
  cor.list <- mapply(splitCors, clusters, cor.list, SIMPLIFY=FALSE)
  cor.list <- unlist(cor.list, recursive=FALSE)

  if (length(cor.list) < 2 ^ n.splits) {
    return(spectralSplit(cor.list, n.splits, qt.threshold))
  } else {
    return(gene.clusters)
  }
}

geneSimilarity <- function(cluster.tree, genes, set.size=2) {

  # TODO: gene sometimes not appearing in cluster tree after selection, fix this
  gene.set <- combn(genes, set.size, simplify=FALSE)

  # determine which trees contain both genes in a given pair
  setInTree <- function(set, tree) return(all(set %in% unlist(tree))
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
    proportion <- mean(sapply(joint.trees, sameLeaf, set=set))
    return(proportion)
  }
  set.proportions <- mapply(sameLeafProportion, set=gene.set, shared.idcs=pair.trees) 

  return(list(gene.set=gene.set, set.prop=set.proportions))
}

 
generateSimilarityMatrix <- function(pair.list, pair.similarity, genes) {
   
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
  
getLocalModules <- function(cor.list, threshold=0.25) {

  # Determine modules that are consistent across multiple principle patterns
  # args:
  #  cor.list: a list of correlation matrices, typically determined from
  #   multiple pp
  # threhold: a numeric value between 0 and 1 specifying the quantile to
  #   set adjacency at
  set.seed(47)
  spectral.outputs <- lapply(cor.list, spectralPreprocess, threshold=threshold, n.eigen=1)
  spectral.vectors <- sapply(spectral.outputs, '[', 'sv')
  spectral.clustering <- lapply(spectral.vectors, kmeans, center=2, nstart=10)
  spectral.clusters <- sapply(spectral.clustering, '[', 'cluster')

  # find all genes that are expressed in all local networks
  gene.names <- lapply(cor.list, rownames)
  jointly.expressed <- Reduce(intersect, gene.names)

  # for each jointly expressed gene, determine its cluster in the different
  # local networks
  cluster.mat <- matrix(0, nrow=length(jointly.expressed), ncol=length(cor.list))
  for (i in 1:length(jointly.expressed)) {

    gene <- jointly.expressed[i]
    gene.idcs <- sapply(gene.names, function(g) which(g==gene))
    cluster.mat[i,] <- mapply(function(clusters, idx) clusters[idx], spectral.clusters, gene.idcs)
  }

  module.mat <- matrix(FALSE, nrow=length(jointly.expressed), ncol=length(jointly.expressed))
  for (i in 1:nrow(module.mat)) {
    for (j in i:nrow(module.mat)) {
      mod <- all(cluster.mat[i,] == cluster.mat[j,])
      module.mat[i,j] <- mod
      module.mat[j,i] <- mod
    }
  }
  rownames(module.mat) <- jointly.expressed
  colnames(module.mat) <- jointly.expressed

  modules <- unique(apply(module.mat, MARGIN=1, which))
  modules <- sapply(modules, function(m) jointly.expressed[m])
  return(list(mod=modules, mat=module.mat, genes=jointly.expressed))
}

subsetNetwork <- function(cor.mat, genes) {

  cor.names <- rownames(cor.mat)
  gene.idcs <- cor.names %in% genes
  return(cor.mat[gene.idcs, gene.idcs])
}

plotCorGraph <- function(cor.mat, threshold=NULL, qt.threshold=0.5, scale=0.5) {

  
  if (is.null(threshold)) {
    threshold <- getQtThreshold(cor.mat, qt.threshold) 
  }
  # remove edges for correlations below quantile threshold
  diag(cor.mat) <- 0
  cor.mat[cors > threshold[1] & cors < threshold[2]] <- 0

  if (! all(cor.mat == 0)) {
    graph <- graph.adjacency(cor.mat, mode='lower', weighted=TRUE)
    edge.sign <- sign(E(graph)$weight)
    edge.weight <- abs(E(graph)$weight)
    E(graph)$color <- ifelse(edge.sign > 0, 'green', 'red')

    graph.layout <- layout.circle(graph)
    plot(graph, edge.width=exp(scale * edge.weight) / scale, vertex.label.cex=0.5,
         layout=graph.layout)
  } else{
    warning('no interactions at specified threshold')
  }
}

getQtThreshold <- function(cor.mat, qt.threshold) {

  # Calculate a specified quantile for a correlation matrix
  # args:
  #  cor.mat: a matrix of correlations
  #  qt.threshold: a numeric value between 0 and 1 specifying the quantile

  qt.threshold <- max(qt.threshold, 1-qt.threshold)
  cor.vector <- cor.mat[upper.tri(cor.mat, diag=FALSE)]
  threshold <- quantile(cor.vector, c(1-qt.threshold, qt.threshold))
  return(threshold)
}

forceSignAgreement <- function(networks) {

  # ensure that sign agrees for all pairwise interactions in all correlation
  # matrices of networks by setting those that disagree to 0
  cor.signs <- sapply(networks, sign)
  same.sign <- apply(cor.signs, 1, function(r) length(unique(r)) == 1)
  matchedSigns <- function(mat, same.sign) {
      mat[!same.sign] <- 0
    return(mat)
  } 
  networks <- lapply(network.subsets, matchedSigns, same.sign=same.sign)
  return(networks)
}
