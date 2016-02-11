
plotEmb <- function(D, width, height) {
  
  img <- matrix(D, nrow=height, ncol=width)
  return(img)

} 
 

generateImage <- function(template, data) {

  # function for plotting expression data on embryo template
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

  # calculate the weighted correlation between two vectors x and y
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

  # wraper to compute weighted correlation matrix
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

  # given a matrix with multiple column observations, take the maximum across
  # rows
  maxCor <- function(x) {
    max.idx <- which.max(abs(x))
    return(x[max.idx])
  }
  row.max <- t(apply(x, 1, maxCor))
  return(c(row.max))
}

mergeDuplicates <- function(cor.mat, node.names) {

  if (nrow(cor.mat) != length(node.names)) stop('dimension mismatch')
  # take maximum correlations for duplicated genes and delete 
  duplicate.idcs <- which(duplicated(node.names))
  duplicate.names <- node.names[duplicate.idcs]
  matchName <- function(target, x) return(which(x == target))
  duplicates <- lapply(duplicate.names, matchName, x=node.names)

  mergeMat <- function(x, idcs) {
    merged <- mergeMax(x[, idcs])
    x[idcs, ] <- matrix(rep(merged, length(idcs)), ncol=ncol(x), byrow=TRUE)
    x[, idcs] <- matrix(rep(merged, length(idcs)), nrow=nrow(x))
    return(x)
  }
  for (i in 1:length(duplicates)) {
    cor.mat <- mergeMat(cor.mat, duplicates[[i]])
  }
  idcs2remove <- unlist(sapply(duplicates, '[', -1))
  cor.mat <- cor.mat[-idcs2remove, -idcs2remove] 
  rownames(cor.mat) <- node.names[-idcs2remove]
  colnames(cor.mat) <- node.names[-idcs2remove]
  return(cor.mat)
}

spectralVectors <- function(cors, threshold, n.eigen) {

  # determine threhsold value based on wuantiles of correlation matrix
  qt.threshold <- getQtThreshold(cors, threshold) 

  # calculate adjacency and normalizing matrices
  adjacency <- cors
  diag(adjacency) <- 0
  adjacency[adjacency <= qt.threshold[1] | adjacency >= qt.threshold[2]] <- 1
  adjacency[adjacency != 1] <- 0

  n <- nrow(adjacency)
  col.sums <- colSums(adjacency)
  col.sums.sqrt <- sqrt(col.sums)
  d.sqrt <- col.sums.sqrt * diag(n)

  # remove indices with no neighbors
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
    temp <- numeric(nrow(cors))
    temp[-idcs2remove] <- spectral.vectors
    spectral.vectors <- matrix(temp, nrow=nrow(cors))
  }
  rownames(spectral.vectors) <- rownames(cors)
  return(list(sv=spectral.vectors, e=e$values))
}

getLocalModules <- function(cor.list, thrsh=0.25, n.eigen=1) {

  set.seed(47)
  spectral.outputs <- lapply(cor.list, spectralVectors, threshold=thrsh, n.eigen=n.eigen)
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

subsetNetwork <- function(cors, genes) {

  # subset network to include only observations specified by genes
  cor.names <- rownames(cors)
  gene.idcs <- cor.names %in% genes
  return(cors[gene.idcs, gene.idcs])
}

plotCorGraph <- function(cors, threshold=NULL, qt.threshold=0.5, scale=0.5) {

  
  if (is.null(threshold)) {
    threshold <- getQtThreshold(cors, qt.threshold) 
  }
  # remove edges for correlations below quantile threshold
  diag(cors) <- 0
  cors[cors > threshold[1] & cors < threshold[2]] <- 0

  if (! all(cors == 0)) {
    graph <- graph.adjacency(cors, mode='lower', weighted=TRUE)
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

getQtThreshold <- function(cors, qt.threshold) {

  qt.threshold <- max(qt.threshold, 1-qt.threshold)
  cor.vector <- cors[upper.tri(cors, diag=FALSE)]
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
