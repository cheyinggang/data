

rm(list = ls())

rt=read.table("ubi_anno.txt",sep="\t",header=T,row.names=1,check.names=F)   
rt1=t(rt[,-1])

rt2=rt1[rowMeans(rt1)>0,]

library(ggpubr)
library(pheatmap)
Type=as.data.frame(rt[,1])
colnames(Type)=c("COVID-19")
rownames(Type)=colnames(rt1)
summary(factor(Type[,1]))


pheatmap(rt2, annotation_col=Type,gaps_col = c(67),
         color =colorRampPalette(c( "white","#ED0000FF"))(50),
         cluster_rows = T,
         cluster_cols =T,
         fontsize=12,
         fontsize_row=10,
         show_colnames=F,
        show_rownames=T,
         fontsize_col=3)

library(tidyverse)
rt1=as.data.frame(t(rt1))
dat <- rt1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Complement_cell,value = Expression,-Sample)
dat$COVID_19 = rt[,1]

ggplot(dat,aes(Complement_cell,Expression,fill = COVID_19 )) + 
  geom_boxplot(outlier.shape = 21,color = "white") + 
  labs(x = "", y = "Enrichment scores") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5))+
  theme_classic()  + 
  theme(axis.title.y =  element_text(size = 10,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 12,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "right") + 
  scale_fill_manual(values = c("deepskyblue1", "tomato2"))+ stat_compare_means(aes(group = COVID_19,label = ..p.signif..),method = "t.test")
ggsave("Figure_F.pdf",width =20,height =5)



