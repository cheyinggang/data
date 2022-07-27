

rm(list = ls())
library(tidyverse)
rt=read.table("immune_group.txt",sep="\t",header=T,row.names=1,check.names=F)
rt2=t(rt[,-1])
group=rt[,1]
phylum=as.data.frame(rt2)
phylum_per <- as.data.frame(lapply(phylum, function(x) x / sum(x)))
row.names(phylum_per) <- row.names(phylum)

phylum.ave <- apply(phylum_per, 1, FUN=mean)
phylum.2 <- cbind(phylum_per, phylum.ave)[order(-phylum.ave),]


phylum.2 <- subset(phylum.2[1:10,], select=-phylum.ave)
phylum.2 <- rbind(phylum.2, others=apply(phylum.2, 2, function(x){1-sum(x)}))
library(reshape2)

PhylumID=row.names(phylum.2)
phylum.2=cbind(PhylumID,phylum.2)
phylum.2$PhylumID <- factor(phylum.2$PhylumID, levels = rev(phylum.2$PhylumID))
phylum.gg <- melt(phylum.2, id.vars="PhylumID", variable.name="SampleID", value.name="Abundance")
head(phylum.gg)
phylum.gg$group <- group



library(wesanderson)
library(colortools)
library(ggpubr)

ggbarplot(phylum.gg, x = "SampleID", y="Abundance", color="black", fill="PhylumID",
          legend="right", 
          legend.title="Top immune cells", main="Relative immune cells per patient",
          font.main = c(12,"bold", "black"), font.x = c(12, "bold"), 
          font.y=c(12,"bold")) + 
  theme_classic()  +
  rotate_x_text() + 
  scale_fill_manual(values=c("gray",rev(wheel("skyblue3")[1:10]))) +
  facet_grid(~ group, scales = "free_x", space='free',) + 
  labs(x = "Sample", y = "Relative proportion") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"), 
        plot.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"),
        axis.title.y =  element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"),
        legend.text =  element_text(size = 10,face = "bold")
            ) 
ggsave(filename = "Figure_B.pdf", device="pdf", width=30, height=6)


