rm(list=ls())
library(pheatmap)
library(Rtsne)
library(ggfortify)
library(mvtnorm)
 
rt2=read.table("ubi_virus.txt",sep="\t",header=T,check.names=F) 
rownames(rt2)<-rt2[,1]
rt2<-rt2[,2:ncol(rt2)]

rt2=rt2[(order(rt2[,1])),]
summary(factor(rt2[,1]))



rt1=rt2[,4:ncol(rt2)]
rt=as.matrix(log2(rt1+1))  

group=rt2[,2]

df=cbind(as.data.frame(rt),group=rt2[,2])
head(df)

pca_dat3 <- prcomp( df[,1:(ncol(df)-1)], rank. = 3)
library(factoextra)
fviz_pca_ind(pca_dat3,geom.ind="point",pointsize=4,
             pointshape=21,fill.ind=group, ellipse.type = "confidence",  
             addEllipses = TRUE,
             legend.titl="Groups",title="PCA of ubiquitin proteins")+theme_grey()


