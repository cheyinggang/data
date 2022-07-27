
rm(list = ls())

rt=read.table("TNF_genes.txt",sep="\t",header=T,row.names=1,check.names=F)  
rt1=rt[,-1]
rt1=t(rt1)
rt1=log2(rt1+1)
library(tidyverse)
library(ggpubr)
rt1=as.data.frame(t(rt1))
dat <- rt1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Complement_Gene,value = Expression,-Sample)
dat$Group = rt[,1]

ggplot(dat,aes(Complement_Gene,Expression,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "white") + 

  labs(title = "TNF family genes", 
       x="",y = "Expression") +
  theme(legend.position = "top") + theme_classic()  + 
  theme(axis.title.y =  element_text(size = 10,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 12,angle = 45,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  
  theme(axis.text.x = element_text(angle=60,vjust = 0.5))+
  scale_fill_manual(values = c("tomato2","#0099B4FF"))+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
ggsave("Figure_H.pdf",width = 10,height = 6.1)

