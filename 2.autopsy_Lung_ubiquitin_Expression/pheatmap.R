

rm(list = ls())

rt=read.table("lung_Ubi.txt",sep="\t",header=T,row.names=1,check.names=F)    
rt1=rt[,-1]
rt1=t(rt1)
rt1=log2(rt1+1)
rt2=rt1[rowMeans(rt1)>0,]
library(ggpubr)
library(pheatmap)
Type=as.data.frame(rt[,1])
colnames(Type)=c("group")
rownames(Type)=colnames(rt1)
summary(factor(Type[,1]))


pheatmap(rt2, annotation_col=Type,gaps_col = c(2),
         color =colorRampPalette(c("#00468B99", "white", "#ED0000FF"))(50),
         cluster_rows = T,
         cluster_cols =F,
         fontsize=12,
         fontsize_row=8,
          scale="row",
         show_colnames=F,
        show_rownames=F,
         fontsize_col=3)


 