library(R.matlab)
load('geneExpSym.RData')
load('dictFitDataNLS.Rdata')
mat.data <- readMat('osMat.mat')


annot.genes <- unlist(mat.data$geneNames)
gut.gene.idcs <- match(geneNames, annot.genes)
annot.genes <- which(!is.na(gut.gene.idcs))
annot.alpha <- alpha[, annot.genes]
gut.gene.idcs <- gut.gene.idcs[annot.genes]
gut.gene.mat <- mat.data$osMat[gut.gene.idcs,]
colnames(gut.gene.mat) <- colnames(geneExpSym)

permTest <- function(x, labels, n.perms=100) {

  
  permDiff <- function(x.p, l.p) abs(mean(x.p[l.p == 1]) - mean(x.p[l.p == 0]))
  observed.diff <- permDiff(x, labels)
  permutations <- replicate(n.perms, sample(labels, length(labels)))
  perm.diffs <- apply(permutations, 2, permDiff, x.p=x)
  perm.p <- mean(perm.diffs > observed.diff)
  return(perm.p)
}

n.dict <- nrow(alpha)
n.annot <- ncol(gut.gene.mat)
p.val.dict <- matrix(0, nrow=n.dict, ncol=ncol(gut.gene.mat))
rownames(p.val.dict) <- paste0('PP', 1:12)
colnames(p.val.dict) <- colnames(gut.gene.mat)

for (i in 1:n.annot) {

  print(paste0('annot atom: ', i, '...'))
  annot <- gut.gene.mat[,i]
  
  for (j in 1:n.dict) {
    alpha.temp <- alpha[j,]
    exp.idcs <- which(annot==1)
    unexp.idcs <- which(annot==0)
    if (length(exp.idcs) == 0 | length(unexp.idcs) == 0) {
      p.val.dict[j,i] <- 1
    } else {
      p.val.dict[j,i] <- wilcox.test(alpha.temp[exp.idcs], alpha.temp[unexp.idcs], 'greater')$p.val
    }#p.val.dict[j,i] <- permTest(alpha.temp, annot)
  }
}

annot.names <- colnames(p.val.dict)
idcs <- sapply(c('gut', 'Gut', 'Pole'), grep, x=annot.names)
idcs <- unlist(idcs)
#nan.idcs <- apply(p.val.dict, 2, function(col) any(is.nan(col)))
sim.mat <- p.val.dict[,idcs]

library(fields)
pdf('fateMap.pdf')
cols <- colorRampPalette(c('white', 'orange', 'red'))(25)
image(sim.mat^2, axes=F, col=cols)
mtext(text=rownames(sim.mat), side=1, line=0.3, at=seq(0,1,length.out=nrow(sim.mat)), las=1, cex=10/nrow(sim.mat))
mtext(text=colnames(sim.mat), side=2, line=0.3, at=seq(0,1,length.out=ncol(sim.mat)), las=1, cex=10/ncol(sim.mat))
image.plot(sim.mat, legend.only=TRUE, horizontal=TRUE, smallplot=c(0.05, 0.45, 0.06, 0.09), col=cols) 
dev.off()
