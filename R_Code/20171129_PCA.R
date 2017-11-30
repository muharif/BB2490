count=read.table('../Data/raw_counts.txt',header = 1, row.names = 1)
count=read.table('../Data/Done_deal.txt',header = 1, row.names = 1)

metadata=read.table('../Data/metadata_patient_code.txt',header = 1,stringsAsFactors = F)
metadata$status=metadata$NASH_IDENTIFIER
for(i in 1:dim(metadata)[1]){
  if(grepl("^([A-Z])", metadata[i,1], perl = TRUE)){
    metadata[i,1] = metadata[i,1]
  } else {
    metadata[i,1] = paste0('X',metadata[i,1])
  }
  metadata[i,3]=substr(metadata[i,2], 1, 1)
}

library(geneLenDataBase)
data(hg19.ensGene.LENGTH)
gene_length=read.table('../Data/gene_length.txt',header = 1)
transcript=intersect(rownames(count),gene_length[,1])


subjects=intersect(colnames(count),metadata[,2])
metadata2=metadata[metadata[,1] %in% subjects,]
count2=count[transcript,subjects]
count2=count2[order(rownames(count2)),,drop=FALSE]

trans_length=hg19.ensGene.LENGTH[hg19.ensGene.LENGTH$Gene %in% transcript,]
trans_length=trans_length[order(trans_length$Gene),]

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
tpm=countToTpm(count2,trans_length$Length)

tpm=tpm[,metadata2$Sample_ID]

write.table(metadata2,'../Data/Metadata_Processed.txt',sep='\t')
write.table(count2,'../Data/Count_Processed.txt',sep='\t')
write.table(trans_length,'../Data/gene2transcript.txt',sep='\t')

stages<-unique(metadata2$status)
coldef.stage<-c("black","red","green","blue","cyan","magenta")
names(coldef.stage) <- stages
col.stage <- coldef.stage[metadata2$status]


PC<-prcomp(log2(t(tpm)+1))

library(rgl)
plot3d(PC$x[,1:3],col=col.stage,main="first PCA")
# text3d(pc$scores[,1:3],texts=rownames(iris))
text3d(PC$loadings[,1:3], texts=rownames(PC$loadings), col="red")
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}
lines3d(coords, col="red", lwd=4)
# legend("topleft",stages,pch=coldef.stage,cex=0.5,bty='n')
