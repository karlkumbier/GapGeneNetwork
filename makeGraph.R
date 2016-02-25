source('genCircle.R')

# Entire network
gapGenes <- c('Kr', 'kni', 'hb', 'gt', 'tll', 'hkb')
nameStrUniq <- rownames(edgeMatTemp)#unique(tf.names) 
edgeMatTemp <- generateAdjacency(cors, qt.thresh=0.9) * cors
nameStrUniq <- rownames(edgeMatTemp)#unique(tf.names) 
numStrTF = length(nameStrUniq);
theta = seq(from = 360 , to = 0, length.out = numStrTF);
lowerIdx = floor(numStrTF/4);
upperIdx = ceiling(numStrTF/4*3);
theta[lowerIdx:upperIdx] = 180 + theta[lowerIdx:upperIdx];
coord = genCircle(numStrTF+1); 
coord$x = coord$x[1:numStrTF];
coord$y = coord$y[1:numStrTF];
idxGap = NULL;
for (i in 1:length(gapGenes)){
      idxGap = c(idxGap, which(nameStrUniq==gapGenes[i]));
}


distX = as.dist(1 - edgeMatTemp);
hr <- hclust(distX, method="average"); # h-clustering
idx0 = hr$order;

nameStrUniq1 = nameStrUniq[idx0];
edgeMatTemp = edgeMatTemp[idx0,idx0];

nb = colSums(abs(edgeMatTemp)>0.001);
nb[nb == 1] = 0;
cex = 2*sqrt(nb)/sqrt(max(nb));

plot(coord$x[-idxGap],coord$y[-idxGap],pch = 16,cex = cex[-idxGap],
     xaxt = 'n', yaxt = 'n',xlim = c(-1.1,1.1),ylim = c(-1.1,1.1),
     xlab = '',ylab = '', col = '#FF000080',bty = 'n');
points(coord$x[idxGap],coord$y[idxGap],pch = 16,cex = cex[idxGap], col = '#FF000080');
for (k in 1:(numStrTF-1)){
  for (j in (k+1):numStrTF){
    if (edgeMatTemp[k,j]<=-0.001){
      segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#FF000099',lwd = 1);
    }
    if (edgeMatTemp[k,j]>0.001){
      segments(coord$x[k],coord$y[k],coord$x[j],coord$y[j],col = '#0000FF99',lwd = 1);
    }
  }
}
for (k in 1:length(nameStrUniq)){
  if (sum(gapGenes == nameStrUniq1[k])){
    text(coord$x[k]*1.09,coord$y[k]*1.09,nameStrUniq1[k],srt = theta[k],cex = 0.5,col ='red');
  }else{
    text(coord$x[k]*1.09,coord$y[k]*1.09,nameStrUniq1[k],srt = theta[k],cex = 0.55);
  }
}
