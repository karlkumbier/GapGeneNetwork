organ.dir <- Sys.getenv('ORGAN_DIR')
library(GOstats)
library(GeneAnswers)
library(org.Dm.eg.db)
library(GO.db)
setwd('~/Desktop/LBL/NMF/gapGeneNetwork/')
source('./utilities.R')
load('./dictFitDataNLS.RData')
load(paste0(organ.dir, '/files/annotMap/goData/flybaseGO.Rdata'))

pp.centers <- 5:8
pp.neighbors <- list(c(4, 6), c(4, 5, 7), c(6, 8), c(7, 9))
test.over <- list()
expressed.genes <- c()
universe.genes <- c()
adj.thresh <- seq(0.8, 0.95, by=0.05)
diff.k <- 2:8
qt.thresh <- seq(0.8, 0.95, by=0.05)

for (a.t in adj.thresh) {
  for (d.k in diff.k) {
    #for (q.t in qt.thresh) {
      iter <- paste('processing:', a.t, d.k, '...')
      print(iter)
      cg.sim <- list()
      genes <- c()
      thresh <- 0
      for (i in 1:length(pp.centers)) {
        
        cors <- localCor(pp.centers[i], pp.neighbors[[i]], X, alpha, geneNames)
        adj.mat <- generateAdjacency(cors, qt.thresh=a.t)
        diffusion.mat <- diffusionDist(adj.mat, k=d.k)
        #thresh <- max(0, getQtThreshold(diffusion.mat, q.t, remove.zero=T)[1])
        genes <- c(genes, rownames(cors))
        cg.sim[[i]] <- diffusion.mat['CG13894',]
      }
      
      genes <- unique(genes)
      combinePP <- function(l, genes) {
        gene.mat <- matrix(nrow=4, ncol=length(genes))
        colnames(gene.mat) <- genes
        for (i in 1:length(l)) {
          temp <- l[[i]]
          gene.mat[i,names(temp)] <- temp
        }
        return(gene.mat)
      }
      
      gene.mat <- combinePP(cg.sim, genes)
      all.mean <- colMeans(gene.mat)
      universe.genes <- names(all.mean)
      temp <- all.mean[!is.na(all.mean)]
      n.expressed <- 50
      expressed.genes <- names(temp)[order(temp, decreasing=F)[1:n.expressed]]
     
       #TODO UNIVERSE GENES ARE WRONG FOR THIS!!!
      #universe.genes <- c(universe.genes, rownames(cors))
      #expressed.genes <- c(expressed.genes, universe.genes[which(cg.sim < thresh[1])])
      
      universe.genes <- unique(universe.genes)
      expressed.genes <- unique(expressed.genes)
      
      # Map gene names to entrez indices
      flybase.idcs <- match(expressed.genes, flybase.data$SYM)
      expressed.fb <- as.character(flybase.data$ID[na.omit(flybase.idcs)])
      
      flybase.idcs <- match(universe.genes, flybase.data$SYM)
      universe.fb <- as.character(flybase.data$ID[na.omit(flybase.idcs)])
      
      fb.entrez <- unlist(as.list(org.Dm.egFLYBASE2EG))
      expressed.entrez <- as.character(na.omit(fb.entrez[expressed.fb]))
      universe.entrez <- as.character(na.omit(fb.entrez[universe.fb]))
      
      # test for enrichment using conditional hypergeometric test
      p.cutoff <- 1
      params.over <- new('GOHyperGParams', geneIds=expressed.entrez, 
                         universeGeneIds=universe.entrez, annotation='org.Dm.eg.db', 
                         ontology='BP', pvalueCutoff=p.cutoff, conditional=TRUE, 
                         testDirection='over')
      test.over.full <- hyperGTest(params.over)
      file.name <- paste0('./cgGOEnrichment/diffusion_', a.t, '_', d.k, 'Full.Rdata')
      save(file=file.name, test.over.full)
    }
  }
#}

