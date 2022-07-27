

library(ggpubr)
rm(list = ls())
rt2=read.table("violin.txt",sep="\t",header=T,check.names=F) 
COLS <- c("tomato3", "deepskyblue3")
head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$HFDP45, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of HFD45",x="cluster", y = "HFD45")+ 
  stat_compare_means() # Add global p-value
p1=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$ventilator_free_days, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of ventilator_free_days",x="cluster", y = "ventilator_free_days")+ 
  stat_compare_means() 
p2=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")





head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$`ferritin(ng/ml)`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Ferritin",x="cluster", y = "Ferritin(ng/ml)")+ 
  stat_compare_means() 
p3=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")



head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$`crp (mg/l)`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of CRP",x="cluster", y = "CRP (mg/l)")+ 
  stat_compare_means() 
p4=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$`ddimer (mg/l_feu)`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of ddimer",x="cluster", y = "ddimer (mg/l_feu)")+ 
  stat_compare_means()
p5=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=rt2$`procalcitonin (ng/ml)`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Procalcitonin",x="cluster", y = "Procalcitonin (ng/ml)")+ 
  stat_compare_means()
p6=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")
library(cowplot)

plot_grid(p1,p2,p3,p4,p5,p6,align = "h",nrow =3)



library(ggstatsplot)
p7=ggpiestats(rt2,  cluster, Sex, palette = 'Set2',title = 'Gender')+ 
  scale_fill_manual(values = COLS)

head(rt2)

p8=ggpiestats(rt2,  cluster, ICU, palette = 'Set2',title = 'ICU') + 
  scale_fill_manual(values = COLS) 


head(rt2)

p9=ggpiestats(rt2,  cluster, mechanical_ventilation, palette = 'Set2',title = 'mechanical_ventilation')+ 
  scale_fill_manual(values = COLS)  


library(cowplot)
plot_grid(p7,p8,p9,align = "h",nrow = 3)



