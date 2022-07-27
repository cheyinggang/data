
rm(list = ls())
library(limma)


rt=read.table("cluser_Ubi_gene.txt",sep = "\t")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
colnames(rt)=rt[1,]
data1=rt[2:nrow(rt),2:ncol(rt)]
data=data1[-1,]
summary(factor(data1[1,]))
p <- data1[1,]
group <- data1[1,]
group <- factor(group,levels = c("cluster1","cluster2"))
dimnames=list(rownames(data),colnames(data))
exprSet=matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
exprSet=log2(exprSet+1)
exprSet_symbol1 <- aggregate(x = exprSet,
                             by = list(rownames(exprSet)),
                             FUN = mean)

rownames(exprSet_symbol1)=exprSet_symbol1[,1]
exprSet_symbol1=exprSet_symbol1[,2:ncol(exprSet_symbol1)]
dimnames=list(rownames(exprSet_symbol1),colnames(exprSet_symbol1))
exprSet=matrix(as.numeric(as.matrix(exprSet_symbol1)),nrow=nrow(exprSet_symbol1),dimnames=dimnames)
exprSet=exprSet[rowMeans(exprSet)>0,]
res.pca <- prcomp(t(exprSet), scale = TRUE)
library(factoextra)
fviz_pca_ind(res.pca, geom.ind = "point",
             col.ind = group,
             palette = c( "tomato3","deepskyblue4"),
           
             addEllipses = TRUE, ellipse.type = "convex",
             legend.title = "Cluster",title="PCA of COVID-19 patients"
)
