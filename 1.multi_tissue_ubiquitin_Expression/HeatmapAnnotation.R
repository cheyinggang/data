

rm(list = ls())
rt=read.table("UbiGene_Anno.txt",sep="\t",header=T,row.names=1,check.names=F)   
rt1=rt[,-(1:12)]
rt1=t(rt1)
rt1=log2(rt1+1)
rt2=rt1[rowMeans(rt1)>0,]
library(ggpubr)
library(pheatmap)
Type=as.data.frame(rt[,1:12])

summary(factor(Type[,10]))
library(ComplexHeatmap)


grid.newpage()
pushViewport(viewport(y = 0.7, height = 0.6, width = 0.8))

ha = HeatmapAnnotation(df=Type)

draw(ha, 1:10)

rt2=as.data.frame(rt2)
Heatmap(rt2, name = "expression",show_row_names = F,
        show_column_names = F,
        column_title = NULL,bottom_annotation =  ha
)


