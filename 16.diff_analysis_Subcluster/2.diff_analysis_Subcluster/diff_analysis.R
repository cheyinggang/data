
rm(list = ls())
library(limma)


rt=read.table("cluser_gene.txt",sep = "\t")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
colnames(rt)=rt[1,]
data1=rt[2:nrow(rt),2:ncol(rt)]
data=data1[-1,]
summary(factor(data1[1,]))
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
save(exprSet,file = "exprSet.Rdata")
exprSet=exprSet[rowMeans(exprSet)>0,]
res.pca <- prcomp(t(exprSet), scale = TRUE)
library(factoextra)
fviz_pca_ind(res.pca,col.ind = group)

design <- model.matrix(~group)
colnames(design) <- levels(group)
design
fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "allDiff.Rdata")

library(dplyr)
diffgene <- allDiff %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)
diffgene <- subset(allDiff,abs(logFC) >1 & adj.P.Val < 0.05)
save(diffgene,group,file = "diffgene.Rdata")




rm(list = ls())
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "allDiff.Rdata")

df <-allDiff
head(df)
df$threshold = factor(ifelse(df$P.Value  < 0.05 & abs(df$logFC) >= 1, ifelse(df$logFC >= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) 

for_label <- df %>% 
  filter(abs(logFC) >2& adj.P.Val< 0.05)

p1=ggplot(df,aes(x=logFC,y= -log10(P.Value),size = abs(logFC),fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_size(limits  = c(0, 4.5))+
  scale_fill_manual(values=c("#fe0000","#13fc00","#bdbdbd"))+
  geom_text_repel(
    data = df[df$P.Value<0.05&abs(df$logFC)>0.5,],
    aes(label = gene),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  ylab('-log10 (P value)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
  theme_classic( 
    base_line_size = 1
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+ 
  labs(title = "The difference of non-COVID-19 and COVID-19 patients")+
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", 
                                   size = 7, 
                                   face = "bold")
  ) 

p1
ggsave("Figure_B.pdf",width = 8,height = 6)
