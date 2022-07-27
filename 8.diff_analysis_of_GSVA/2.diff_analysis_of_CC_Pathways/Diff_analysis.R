
rm(list = ls())
library(limma)

rt=read.table("UBI_CC_pathways.txt",sep = "\t",header = T,row.names = 1)
data=rt

group <- c(rep("COVID_19",100),rep("non_COVID_19",26)) 

group <- factor(group,levels = c("non_COVID_19","COVID_19"))

dimnames=list(rownames(data),colnames(data))
exprSet=matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)

exprSet_symbol1 <- aggregate(x = exprSet,
                             by = list(rownames(exprSet)),
                             FUN = mean)

rownames(exprSet_symbol1)=exprSet_symbol1[,1]
exprSet_symbol1=exprSet_symbol1[,2:ncol(exprSet_symbol1)]
dimnames=list(rownames(exprSet_symbol1),colnames(exprSet_symbol1))
exprSet=matrix(as.numeric(as.matrix(exprSet_symbol1)),nrow=nrow(exprSet_symbol1),dimnames=dimnames)



save(exprSet,file = "exprSet.Rdata")

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
  filter(abs(logFC) >0.5)

diffgene <- subset(allDiff,abs(logFC) >0.5 & adj.P.Val < 0.05)
test <- allDiff[allDiff$adj.P.Val < 0.05& abs(allDiff$logFC)>1,]
save(diffgene,group,file = "diffgene.Rdata")


rm(list = ls())
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
load(file = "allDiff.Rdata")
data=allDiff
data$gene <- rownames(data)
df <- data

head(df)
df$threshold = factor(ifelse(df$P.Value  < 0.05 & abs(df$logFC) >= 0.5, ifelse(df$logFC >= 0.5 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df)
p1=ggplot(df,aes(x=logFC,y= -log10(P.Value),size = abs(logFC),fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_size(limits  = c(0, 1))+ 
  scale_fill_manual(values=c("#fe0000","#13fc00","#bdbdbd"))+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
  theme_classic(  
    base_line_size = 1
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+ 
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
for_label <- df %>% 
  filter(abs(logFC) >0.5& adj.P.Val< 0.05)

p1 +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,alpha = 0.6,fontface="bold", 
    color="black"
  )+ theme_classic(base_size = 8)+
  theme(axis.title.y =  element_text(size = 20,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 20,angle = 45,vjust = 0.5,hjust = 0.5))+
  theme(axis.title.x =  element_text(size = 20,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20,angle = 45,vjust = 0.5,hjust = 0.5))+
  
  theme(legend.position = "NA")


ggsave("Figure_B", height = 12, width = 15)


