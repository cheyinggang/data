rm(list = ls())
housing.df=read.table(file = "data/UBI_mechanism.txt",header = T,sep = "\t",row.names = 1)

clinical=housing.df[,1:2]
exp=housing.df[,3:357]
expr=exp

housing.df=cbind(clinical,expr)

set.seed(1) ## to get the same sequence of numbers
# randomly sample 60% of the row IDs for training; the remaining 40% serve as validation
train.rows <- sample(rownames(housing.df), dim(housing.df)[1]*0.6)

# collect all the columns with training row ID into training set:
housing.train <- housing.df[train.rows, ]

# assign row IDs that are not already in the training set, into validation
valid.rows <- setdiff(rownames(housing.df), train.rows)
housing.valid <- housing.df[valid.rows, ]

write.table(housing.valid,"housing.valid.txt",sep = "\t",quote = F)
glmnet_input=housing.train
glmnet_input<-glmnet_input[order(glmnet_input[,1]),]

x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

library(survival)
y <- data.matrix(Surv(glmnet_input[,1],glmnet_input[,2]))
## lasso
library(glmnet)
cv.fit <- cv.glmnet(x, y, family="cox") 
plot(cv.fit)

if(T){
  library(ggplot2)
  library(broom)
  tidied_cv = tidy(cv.fit)
  ggplot(tidied_cv, aes(lambda, estimate)) +
    geom_point(color = "red",size=2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    scale_x_log10()+
    geom_vline(xintercept = cv.fit$lambda.min, lty=2,color = "red",size=1) +
    geom_vline(xintercept = cv.fit$lambda.1se, lty=2) +
    ylab("Partial Likelihood Deviance")+
    annotate("text",x=cv.fit$lambda.min/5,size=5,
             y =max(tidied_cv$estimate),
             label =paste0("lambda.min = ",
                           round(cv.fit$lambda.min,3),
                           "\n",
                           "number = ",tidied_cv$nzero[tidied_cv$lambda==cv.fit$lambda.min]))+
    ggtitle("")+
    theme_bw()
}


fit <- glmnet(x, y, family = "cox")

fit_predict <- fit

plot(fit, label = TRUE)
if(T){
  (Coefficients <- coef(fit, s = cv.fit$lambda.min))
  (Active.Index <- which(as.matrix(Coefficients) != 0))
  (Active.Coefficients  <- Coefficients[Active.Index])
  
  (gene_found <- row.names(Coefficients)[Active.Index])
  save(gene_found,file = "gene_found.Rdata")
}

if(T){
  tidy_covariates = tidy(cv.fit$glmnet.fit)
  index <- tidy_covariates$term %in% gene_found
  library(dplyr)
  data_point <- tidy_covariates %>% 
    filter(index) %>% 
    arrange(lambda) %>% 
    distinct(term,.keep_all = T)
  library(ggrepel)
  library(ggplot2)
  ggplot(tidy_covariates)+
    geom_line(data = subset(tidy_covariates,index), aes(x=lambda, y=estimate, color=as.factor(term)),size=1) +
    geom_line(data = subset(tidy_covariates,!index), aes(x=lambda, y=estimate, color=as.factor(term)),size=0.8,alpha=0.1) +
    guides(color=FALSE) +
    geom_vline(xintercept = cv.fit$lambda.min, lty=2,color='red') +
    scale_x_log10()+
    ylab("Coefficients")+
    geom_point(data =data_point,aes(lambda,estimate,color=as.factor(term)),size=3)+
    geom_text_repel(data =data_point,aes(lambda,estimate,label=term),size=3)+
    theme_classic()
}
ggsave("Figure_D_Coefficients.pdf",width =2,height = 2)
if(T){
  coxdata <- cbind(glmnet_input[,c(1,2)],glmnet_input[,gene_found])
  res <- data.frame()
  genes = gene_found
  for (i in 1:length(genes)) {
    print(i)
    surv =as.formula(paste('Surv(ventilator_free_days, mechanical_ventilation)~', genes[i]))
    cox = coxph(surv, data = coxdata)
    cox = summary(cox)
    p.value=signif(cox$wald["pvalue"], digits=2)
    HR =signif(cox$coef[2], digits=2);#exp(beta)
    HR.confint.lower = signif(cox$conf.int[,"lower .95"], 2)
    HR.confint.upper = signif(cox$conf.int[,"upper .95"],2)
    CI <- paste0("(", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res[i,1] = genes[i]
    res[i,2] = HR
    res[i,3] = CI
    res[i,4] = p.value
  }
  names(res) <- c("ID","HR","95% CI","p.value")
  save(res,file = "res_HR.Rdata")
}

  
  coxresult <- cbind(res[,],Active.Coefficients)
  
  colnames(coxresult) <- c("Symbol", 
                           "HR", "95% CI", "p.value", "Lasso Coefficient")
 write.table(coxresult,file = "coefficience.txt",sep = "\t",row.names = F,quote = F)
  coxresult$HR <- as.numeric(as.character(coxresult$HR))
 coxresult <-  dplyr::arrange(coxresult,desc(HR))


  x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])
  
  rt <- glmnet_input  
  rt <- rt[,1:2]
  order =rownames(rt)
  library(glmnet)
  rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")
  # best AUC cutoff
  rt <- as.data.frame(apply(rt,2,as.numeric))
  rownames(rt) <- rownames(glmnet_input)

  

  source("survivalROC_NEW.R")
  cutoffROC <- function(t,rt){
    rt=rt
    ROC_rt <- survivalROC_NEW(Stime=rt[,1], status=rt[,2], marker = rt[,3], 
                              predict.time =t, method="KM")

    youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)

    index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1

    result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
  }

  t <- 30
  (res.cut <- cutoffROC(t,rt))
  cutoff <- res.cut[1]

  

  library(dplyr)
  if(T){
    rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
    rt$order <- order
    rt$riskScore = as.numeric(rt$riskScore)
    rt = rt %>% arrange(riskScore)
    rt$num =seq(1,length(rownames(rt)))
  }
  

  if(T){
    cutnum <- sum(rt$risk=="low")
    library(ggplot2)
    (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
        geom_point(aes(col=risk),size=3)+
        geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
        geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
        geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size = 8,alpha=0.8)+
        theme_classic()+
        theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
              axis.text.x=element_text(size = 20),
              #axis.ticks.x=element_blank(),
              axis.title.y =  element_text(size = 20),
              axis.text.y = element_text(size = 15))+
        scale_color_manual(values=c("#FC4E07","#00AFBB"))+
        scale_x_continuous(limits = c(0,NA),expand = c(0,0))
      )
  }
  

  if(T){
    event = factor(rt$mechanical_ventilation)
    (plot.sur=ggplot(rt,aes(x=num,y=ventilator_free_days))+
        geom_point(aes(col=event),size=3)+
        geom_vline(aes(xintercept = cutnum),linetype="dashed")+
        theme_classic()+
        theme(axis.title.x=element_text(size = 15),
              axis.text.x=element_text(size = 20),
              axis.ticks.x=element_blank(),
              axis.title.y =  element_text(size = 20),
              axis.text.y = element_text(size = 15))+      
        scale_color_manual(values=c("#00AFBB","#FC4E07"))+
        scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
  }
  

  if(T){
    tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
    library(data.table)
    library(plyr)
    library(scales)
    tmp.m <- reshape2::melt(tmp)
    colnames(tmp.m)=c("Name", "variable", "value")
    tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
    (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
        geom_tile(aes(fill = rescale), colour = "white") + 
        scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y = element_text(size = 15))
      )
  }

  library(cowplot)
  plot_grid(plot.point, plot.sur, plot.h,
            align = 'v',ncol = 1,axis="b")
  ggsave("Figure_EFG_riskscore.pdf",width =8,height = 10)
  ##########################################################################################
  
  

  if(T){
    library("survival")
    library("survminer")
    my.surv <- Surv(rt$ventilator_free_days, rt$mechanical_ventilation)
    group <- rt$risk
    survival_dat <- data.frame(group = group)
    fit <- survfit(my.surv ~ group)
    

    group <- factor(group, levels = c("low", "high"))
    data.survdiff <- survdiff(my.surv ~ group)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    x = summary(coxph(Surv(ventilator_free_days, mechanical_ventilation)~riskScore, data = rt))
    HR = signif(x$coef[2], digits=2)
    up95 = signif(x$conf.int[,"upper .95"],2)
    low95 = signif(x$conf.int[,"lower .95"], 2)
    HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
    CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
    rt <- rt[order(rt[,"riskScore"],decreasing = T),]
    ggsurvplot(fit, data = survival_dat ,
             
               conf.int = F, 
              
               censor = F,
               palette = c("#D95F02","#1B9E77"), 
               
               font.legend = 11,
            
               legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                             paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
               ylab=c('mechanical ventilation probability '),
      
               pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                          paste("p = ",round(p.val,3), sep = "")),
                            HR, CI, sep = "\n"))
  }
  ggsave("Figure_H_me_probability.pdf",width = 5.5,height = 5.5)
  
  ##ROC
  if(T){
    library(tidyr)
    library(purrr)
    library(survivalROC)
    ## Define a helper function to evaluate at various t
    survivalROC_helper <- function(t) {
      survivalROC(Stime=rt$ventilator_free_days, status=rt$mechanical_ventilation, marker = rt$riskScore, 
                  predict.time =t, method="KM")
    }
    ## Evaluate 3,5,10
    survivalROC_data <- data_frame(t = c(10,20,30)) %>%
      mutate(survivalROC = map(t, survivalROC_helper),
             ## Extract scalar AUC
             auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
             ## Put cut off dependent values in a data_frame
             df_survivalROC = map(survivalROC, function(obj) {
               dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
             })) %>%
      dplyr::select(-survivalROC) %>%
      unnest() %>%
      arrange(t, FP, TP)
    ## Plot
    survivalROC_data1 <- survivalROC_data %>% 
      mutate(auc =sprintf("%.3f",auc))%>% 
      unite(year, t,auc,sep = " days AUC: ")
    
    AUC =factor(survivalROC_data1$year)
    survivalROC_data1 %>%
      ggplot(mapping = aes(x = FP, y = TP)) +
      geom_path(aes(color= AUC),size=2)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      theme_classic() +
      theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y =  element_text(size = 15),
                        axis.text.y = element_text(size = 15),
                        legend.text = element_text(size = 12))
  }
  
  ggsave("Figure_I_ROC.pdf",width = 5.5,height = 5.0)
  
  
  
  
  
  
  
  
  
  
  
  
  #***validation
  glmnet_input=read.table(file = "housing.valid.txt",header = T,sep = "\t",row.names = 1)
  x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])
  
  rt <- glmnet_input  
  rt <- rt[,1:2]
  order =rownames(rt)

  library(glmnet)
  rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")

  rt <- as.data.frame(apply(rt,2,as.numeric))
  rownames(rt) <- rownames(glmnet_input)

  rt=read.table(file = "data/rt.txt",header = T,sep = "\t")

  source("survivalROC_NEW.R")
  cutoffROC <- function(t,rt){
    rt=rt
    ROC_rt <- survivalROC_NEW(Stime=rt[,1], status=rt[,2], marker = rt[,3], 
                              predict.time =t, method="KM")
    youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)
    index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1
    result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
  }

  cutoff <- 0.97

  
  # 准备数据
  library(dplyr)
  if(T){
    rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
    rt$order <- order
    rt$riskScore = as.numeric(rt$riskScore)
    rt = rt %>% arrange(riskScore)
    rt$num =seq(1,length(rownames(rt)))
  }
  

  if(T){
    cutnum <- sum(rt$risk=="low")
    library(ggplot2)
    (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
        geom_point(aes(col=risk),size=3)+
        geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
        geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
        geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size =8,alpha=0.8)+
        theme_classic()+
        theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
              axis.text.x=element_text(size = 20),
              #axis.ticks.x=element_blank(),
              axis.title.y =  element_text(size = 20),
              axis.text.y = element_text(size = 15))+
        scale_color_manual(values=c("#FC4E07","#00AFBB"))+
        scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
  }
  

  if(T){
    event = factor(rt$mechanical_ventilation)
    (plot.sur=ggplot(rt,aes(x=num,y=ventilator_free_days))+
        geom_point(aes(col=event),size=3)+
        geom_vline(aes(xintercept = cutnum),linetype="dashed")+
        theme_classic()+
        theme(axis.title.x=element_text(size = 15),
              axis.text.x=element_text(size = 20),
              axis.ticks.x=element_blank(),
              axis.title.y =  element_text(size = 20),
              axis.text.y = element_text(size = 15))+
        scale_color_manual(values=c("#00AFBB","#FC4E07"))+
        scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
  }
  

  if(T){
    tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
    library(data.table)
    library(plyr)
    library(scales)
    tmp.m <- reshape2::melt(tmp)
    colnames(tmp.m)=c("Name", "variable", "value")
    tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
    (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
        geom_tile(aes(fill = rescale), colour = "white") + 
        scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y = element_text(size = 15)))
  }

  library(cowplot)
  plot_grid(plot.point, plot.sur, plot.h,

            align = 'v',ncol = 1,axis="b")
  ggsave("Figure_JKL_riskscore_validation.pdf",width =8,height = 10)
  ##########################################################################################
  
  

  if(T){
    library("survival")
    library("survminer")
    my.surv <- Surv(rt$ventilator_free_days, rt$mechanical_ventilation)
    group <- rt$risk
    survival_dat <- data.frame(group = group)
    fit <- survfit(my.surv ~ group)
    

    group <- factor(group, levels = c("low", "high"))
    data.survdiff <- survdiff(my.surv ~ group)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    x = summary(coxph(Surv(ventilator_free_days, mechanical_ventilation)~riskScore, data = rt))
    HR = signif(x$coef[2], digits=2)
    up95 = signif(x$conf.int[,"upper .95"],2)
    low95 = signif(x$conf.int[,"lower .95"], 2)
    HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
    CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
    
    rt <- rt[order(rt[,"riskScore"],decreasing = T),]
    ggsurvplot(fit, data = survival_dat ,
               
               conf.int = F, 
              
               censor = F, 
               palette = c("#D95F02","#1B9E77"), 
               
               font.legend = 11,
              
               legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                             paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
               ylab=c('mechanical ventilation probability '),
               
               pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                          paste("p = ",round(p.val,3), sep = "")),
                            HR, CI, sep = "\n"))
  }
  ggsave("Figure_M_me_probability_validation.pdf",width = 5.5,height = 5.5)
  
  ##ROC
  if(T){
    library(tidyr)
    library(purrr)
    library(survivalROC)
    ## Define a helper function to evaluate at various t
    survivalROC_helper <- function(t) {
      survivalROC(Stime=rt$ventilator_free_days, status=rt$mechanical_ventilation, marker = rt$riskScore, 
                  predict.time =t, method="KM")
    }
    ## Evaluate 3,5,10
    survivalROC_data <- data_frame(t = c(10,20,30)) %>%
      mutate(survivalROC = map(t, survivalROC_helper),
             ## Extract scalar AUC
             auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
             ## Put cut off dependent values in a data_frame
             df_survivalROC = map(survivalROC, function(obj) {
               dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
             })) %>%
      dplyr::select(-survivalROC) %>%
      unnest() %>%
      arrange(t, FP, TP)
    ## Plot
    survivalROC_data1 <- survivalROC_data %>% 
      mutate(auc =sprintf("%.3f",auc))%>% 
      unite(year, t,auc,sep = " days AUC: ")
    
    AUC =factor(survivalROC_data1$year)
    survivalROC_data1 %>%
      ggplot(mapping = aes(x = FP, y = TP)) +
      geom_path(aes(color= AUC),size=2)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      theme_classic() +
      theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),
                                                axis.title.y =  element_text(size = 15),
                                                axis.text.y = element_text(size = 15),
                                                legend.text = element_text(size = 12))
  }
  
  ggsave("Figure_N_ROC_validation.pdf",width = 5.5,height = 5.0)
  
  