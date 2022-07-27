


rm(list = ls())

rt<- read.table("GO_ubiquitin.txt",sep = "\t")
colnames(rt)=rt[1,]
row.names(rt)=rt[,1]
rt=t(rt[-1,-1])
expre=rt[2:nrow(rt),2:ncol(rt)]
dimnames=list(rownames(expre),colnames(expre))
heatdata=matrix(as.numeric(as.matrix(expre)),nrow=nrow(expre),dimnames=dimnames)

Type=as.data.frame(rt[1,])
Type=as.data.frame(Type[-1,])
colnames(Type)=c("group")
rownames(Type)=colnames(heatdata)


GO= as.data.frame(rt[,1])
GO<-as.data.frame(GO[-1,])
colnames(GO)=c("GO")
rownames(GO)=rownames(heatdata)
library(pheatmap)
library(viridisLite)
library(colorspace)
pdf("Figure_D.pdf",height=22,width=30)
pheatmap(heatdata, 
         cluster_rows = T,
         fontsize=20,
         fontsize_row=14,
         cluster_cols = F,
         annotation_col =Type,
         annotation_row =GO,
         annotation_legend=TRUE, 
         show_rownames = T,
         show_colnames = F,gaps_col = c(100),
         scale = "row",
         color = colorRampPalette(c("blue4", "white", "red4"))(30)
)
dev.off()







