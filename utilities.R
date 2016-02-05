
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
  return(row.max)
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
    x[idcs, ] <- merged
    x[, idcs] <- merged
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




