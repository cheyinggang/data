

rm(list = ls())
rt=read.table("xCell.txt",sep="\t",header=T,row.names=1,check.names=F)    



library(pheatmap)
Type=read.table("clinical_SIG.sig.txt",sep="\t",header=T,row.names=1,check.names=F)

row.names(Type)=colnames(rt)

pheatmap(rt, annotation_col=Type,
         color = colorRampPalette(c("blue","white", "red"))(50),
         
         cluster_cols =F,
         fontsize=12,
         fontsize_row=15,
         scale="row",
         show_colnames=F,gaps_col = c(47,100),
         fontsize_col=3)


