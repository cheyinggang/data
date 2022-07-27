

rm(list = ls())

rt=read.table("A549_SARS-CoV-2_ubi.txt",sep="\t",header=T,row.names=1,check.names=F)   
rt1=t(rt[,-(1:2)])
rt1=log2(rt1+1)
rt1=rt1[rowMeans(rt1)>0,]
library(pheatmap)
Type=as.data.frame(rt[,1:2])
rownames(Type)=colnames(rt1)
pheatmap(rt1, annotation_col=Type,gaps_col = c(18),
         color =colorRampPalette(c("#0099B4FF", "white", "#ED0000FF"))(50),
         cluster_rows = T,
         cluster_cols =F,
         fontsize=12,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
        show_rownames=F,
         fontsize_col=3)

