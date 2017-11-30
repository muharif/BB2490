count=read.table('../Data/Done_deal.txt',header = 1, row.names = 1)
gene_length=read.table('../Data/gene_length.txt',header = 1)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
tpm=countToTpm(count,gene_length$Length)

temp=colnames(count)
conds=c()
for(i in temp){
  conds=append(conds,substr(i,1, 1))
}

stages<-unique(conds)
coldef.stage<-c("black","red","green","blue","cyan","magenta")
names(coldef.stage) <- stages
col.stage <- coldef.stage[conds]

library('Rtsne')
library(rgl)
PC<-prcomp(t(tpm))
plot3d(PC$x[,1:3],col=col.stage,main="first PCA")


tsne.out <- Rtsne(t(log2(tpm+1)),perplexity=10)

plot(tsne.out$Y,col=col.stage,main="tSNE")
legend("topleft",stages,col=coldef.stage,pch=16,cex=0.5,bty='n')
