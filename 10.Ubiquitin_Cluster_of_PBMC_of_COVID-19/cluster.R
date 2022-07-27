
rm(list = ls())
library(ConsensusClusterPlus)
getwd()
workDir=getwd()
setwd(workDir)
rt=read.table("Ubi_genes.txt",sep="\t",header=T,check.names=F,row.names = 1)           
rt2=as.matrix(rt)

rt1=rt2
rt1=log2(rt1+1)
rt1=rt1[rowMeans(rt1)>0,]
data=rt1

data[is.na(data)] <- 0

maxK=9
results = ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
             plot="pdf" )#


icl<-calcICL(results, plot="pdf")


clusterNum=2                 
cluster=results[[clusterNum]][["consensusClass"]]
outTab=cbind(t(rt2),cluster)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="cluster.txt",sep="\t",quote=F,row.names=F)
