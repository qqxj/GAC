#训练集，内部验证集，外部验证集的km

rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
library(survival)
library(survminer)
load(file = "./12.2RSF11.11/RSF_genes.Rdata")

#训练集的km-----------
{
  load(file = "./12.2RSF11.11/训练集.Rdata")
  head(datset1)
  datset1 <- datset1[, c("ID","OS.time","OS","risk.score", RSF_genes)]
  datset1 <- subset(datset1,OS.time > 0)
  datset1 <- subset(datset1,OS.time < 3000)
  
  datset1$risk.group <- ifelse(datset1$risk.score > median(datset1$risk.score), 
                               "High", "Low")
  datset1$risk.group <- factor(datset1$risk.group, levels = c("Low","High"))
  surv_obj <- Surv(time = datset1$OS.time, event = datset1$OS)
  fit <- survfit(surv_obj ~ risk.group, data = datset1)
  ggsurvplot(fit,
             data = datset1,
             risk.table = F,
             pval = TRUE,
             conf.int = TRUE,
             palette = c("#27FF42", "#FF890B"),   # 蓝=Low, 红=High
             legend.title = "Risk Group",
             legend.labs = c("Low Risk","High Risk"),  # 与 factor 顺序对应
             xlab = "Time (days)",
             ylab = "Overall Survival Probability",
             title = "Training Set Risk Score Survival Curve")
  #pdf6*6 
  #美化
  {
    datset1$risk <- factor(datset1$risk.group, levels = c("Low", "High"))
    fit <- survfit(Surv(OS.time, OS) ~ risk, data = datset1)
    
    ggsurvplot(
      fit,
      data = datset1,
      palette = c("#33FF57", "#FF5733"),   # Low=绿色, High=红色
      conf.int = TRUE, conf.int.style = "step",
      pval = TRUE, pval.method = TRUE,
      risk.table = FALSE, risk.table.pos = "in",
      legend = c(0.85, 0.85),
      legend.title = "Risk Group",
      legend.labs = c("Low", "High"),
      title = "Survival Curve for Training Set",
      xlab = "Time (Years)",
      surv.median.line = "hv",
      ggtheme = theme_bw(base_size = 20)
    )
  }
  
}
#师兄的代码（训练集km散点图热图）--------
##### 生存分析
#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
#table(TG_phe$Id %in% dat1$Id)
#phe <- TG_phe[TG_phe$Id %in% dat1$Id,]
rm(list=ls())
load(file = "./12.2RSF11.11/RSF_genes.Rdata")
load(file = "./12.2RSF11.11/训练集.Rdata")
dat1 <- datset1
phe <- datset1
fac_mix <- datset1
{
  # factGene <- c("CYGB","CLVS1","CBY3","CPNE6")
  library(survival)
  library(survminer)
  library(ggplotify)
  library(cowplot)
  library(Hmisc)
  library(pheatmap)
  library(gridExtra)
  #s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(factGene,collapse = ' + '))))
  #model <- coxph(s, data = fac_mix )
  #summary(model,data=fac_mix)
  # RiskScore <- predict(model,type = "risk")
  RiskScore <- fac_mix$risk.score
  #names(RiskScore) = rownames(fac_mix)
  fp <- RiskScore
  phe <- fac_mix
  phe$id <- phe$ID
  
  # 生存图 
  {
    #最佳节点
    # res.cut <- surv_cutpoint(dat1, #数据集
    #                          time = "OS.time", #生存时间
    #                          event = "OS", #生存状态
    #                          variables = "exp") #需要计算的数据列名
    # 
    # summary(res.cut) #查看数据最佳截断点及统计量
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    dat1$OS.time=dat1$OS.time/365
    sfit <- survfit(Surv(OS.time, OS)~RiskGroup, data=dat1)
    # ggsurvplot(sfit, pval=TRUE,xlab ="Time(Years)",surv.median.line = "hv",pval.method = T)
    # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
    #            risk.table =TRUE,pval =TRUE,
    #            conf.int =TRUE,xlab ="Time in years",
    #            ggtheme =theme_light(),
    #            ncensor.plot = TRUE)
    ggsurvplot(sfit,
               palette = c("#FF8300", "#27FF42"),#"#FF5733", "#33FF57"
               conf.int = T,conf.int.style='step', 
               pval = T,pval.method = T,
               risk.table = F,risk.table.pos='in',
               legend=c(0.85,0.85),
               legend.title="Risk Group",
               legend.labs=c("High","Low"),
               title="Survival Curve for Training Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 20))
    #PDF10*10
  }
  
  # 散点+热图 中位数 #pdf14*10
  {
    fp_dat=data.frame(patientid=phe$id,fp=phe$risk.score)
    fp_dat$RiskGroup= ifelse(fp_dat$fp>= median(fp_dat$fp),'High','Low')
    library(dplyr)
    fp_dat=arrange(fp_dat,fp)
    
    sur_dat=data.frame(patientid=phe$id,time=phe$OS.time/365,Status=phe$OS)
    sur_dat$Status=ifelse(sur_dat$Status==0,'Alive','Dead')
    sur_dat$Status=factor(sur_dat$Status,levels = c("Dead","Alive"))
    sur_dat$time <- sur_dat$time
    rownames(sur_dat)=sur_dat$patientid
    sur_dat=sur_dat[fp_dat$patientid,]
    
    exp_dat=dat1[,RSF_genes]
    rownames(exp_dat)=phe$id
    exp_dat=exp_dat[fp_dat$patientid,]
    fp_dat$patientid=1:length(fp)
    sur_dat$patientid=1:length(fp)
    rownames(exp_dat)=1:length(fp)
    rownames(sur_dat)=1:length(fp)
    
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp_dat$fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    
    
    # p1 <- ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=RiskGroup))+
    #   scale_colour_manual(values = c("#FF6666","#00CCCC"))+
    #   theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
    #   geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    #   # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
    #   geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
    #              colour="black", linetype="dotted",size=0.8) +
    #   theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
    #         axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
    #         legend.text=element_text(size=14),legend.title =element_text(size=14))
    # p1
    
    library(ggsci)
    library(scales)
    
    fp_dat$catpo <- rep(median(fp_dat$fp),length(fp_dat$patientid))
    p1 <- ggplot(fp_dat,aes(patientid))+
      geom_ribbon(aes(ymin = catpo, ymax = fp, fill = RiskGroup), alpha = 1)+
      scale_fill_manual(values = c("#FF8300", "#27FF42")) +
      #scale_fill_jco()+
      #scale_colour_manual(values = c("#FF6666","#00CCCC"))+#"#FF5733", "#33FF57"
      theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
      geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
      # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
                 colour="black", linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p1
    
    pal <- pal_jco('default')(10)
    p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=Status),size = 3)+theme_bw()+
      scale_colour_manual(values = c("#FF8300", "#27FF42"))+
      labs(x="Patient ID(Increasing Risk Score)",y="Survival Time(Years)")+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),colour="black", 
                 linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p2
    
    # mycolors <- colorRampPalette(c("#64b5f6", "#fffde7", "#ff5252"), bias = 1.2)(100)
    # mycolors <- colorRampPalette(c("#FF6666", "White", "#00CCCC"), bias = 1.2)(100)
    mycolors <- colorRampPalette(c("#27FF42", "#FFFFFF","#FF8300"), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
                show_rownames = T,cluster_rows = F,
                fontsize_row = 14,fontsize = 14)
    p3
    plots = list(p1,p2,as.ggplot(as.grob(p3)))
    lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
    grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
    # plots = list(p1,as.ggplot(as.grob(p3)))
    # lay1 = rbind(c(rep(1,7)),c(rep(2,7)))
    # grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10))
  }
}





#内部训练集km------
{
  rm(list=ls())
  load(file = "./12.2RSF11.11/RSF_genes.Rdata")
  load(file = "./12.2RSF11.11/内部验证集.Rdata")
  head(datset2)
  datset2 <- datset2[, c("ID","OS.time","OS","risk.score", RSF_genes)]
  datset2 <- subset(datset2,OS.time > 0)
  datset2 <- subset(datset2,OS.time < 3000)
  
  datset2$risk.group <- ifelse(datset2$risk.score > median(datset2$risk.score), 
                               "High", "Low")
  datset2$risk.group <- factor(datset2$risk.group, levels = c("Low","High"))
  surv_obj <- Surv(time = datset2$OS.time, event = datset2$OS)
  fit <- survfit(surv_obj ~ risk.group, data = datset2)
  ggsurvplot(fit,
             data = datset2,
             risk.table = F,
             pval = TRUE,
             conf.int = TRUE,
             palette = c("#27FF42", "#FF890B"),   # 蓝=Low, 红=High
             legend.title = "Risk Group",
             legend.labs = c("Low Risk","High Risk"),  # 与 factor 顺序对应
             xlab = "Time (days)",
             ylab = "Overall Survival Probability",
             title = "Internal validation set risk score survival curve")
  #pdf6*6
}
#师兄的代码（内部验证集km散点图热图）--------
##### 生存分析
#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
#table(TG_phe$Id %in% dat1$Id)
#phe <- TG_phe[TG_phe$Id %in% dat1$Id,]
#rm(list=ls())
#load(file = "./内部验证集.Rdata")
dat1 <- datset2
phe <- datset2
fac_mix <- datset2
{
  # factGene <- c("CYGB","CLVS1","CBY3","CPNE6")
  library(survival)
  library(survminer)
  library(ggplotify)
  library(cowplot)
  library(Hmisc)
  library(pheatmap)
  library(gridExtra)
  #s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(factGene,collapse = ' + '))))
  #model <- coxph(s, data = fac_mix )
  #summary(model,data=fac_mix)
  # RiskScore <- predict(model,type = "risk")
  RiskScore <- fac_mix$risk.score
  #names(RiskScore) = rownames(fac_mix)
  fp <- RiskScore
  phe <- fac_mix
  phe$id <- phe$ID
  
  # 生存图 
  {
    #最佳节点
    # res.cut <- surv_cutpoint(dat1, #数据集
    #                          time = "OS.time", #生存时间
    #                          event = "OS", #生存状态
    #                          variables = "exp") #需要计算的数据列名
    # 
    # summary(res.cut) #查看数据最佳截断点及统计量
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    dat1$OS.time=dat1$OS.time/365
    sfit <- survfit(Surv(OS.time, OS)~RiskGroup, data=dat1)
    # ggsurvplot(sfit, pval=TRUE,xlab ="Time(Years)",surv.median.line = "hv",pval.method = T)
    # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
    #            risk.table =TRUE,pval =TRUE,
    #            conf.int =TRUE,xlab ="Time in years",
    #            ggtheme =theme_light(),
    #            ncensor.plot = TRUE)
    ggsurvplot(sfit,
               palette = c( "#FF890B","#27FF42"),
               conf.int = T,conf.int.style='step', 
               pval = T,pval.method = T,
               risk.table = F,risk.table.pos='in',
               legend=c(0.85,0.85),
               legend.title="Risk Group",
               legend.labs=c("High","Low"),
               title="Survival Curve for Training Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 20))
    #pdf10*10
  }
  
  # 散点+热图 中位数#pdf14*10
  {
    fp_dat=data.frame(patientid=phe$id,fp=phe$risk.score)
    fp_dat$RiskGroup= ifelse(fp_dat$fp >= median(fp_dat$fp),'High','Low')
    library(dplyr)
    fp_dat=arrange(fp_dat,fp)
    
    sur_dat=data.frame(patientid=phe$id,time=phe$OS.time/365,Status=phe$OS)
    sur_dat$Status=ifelse(sur_dat$Status == 0,'Alive','Dead')
    sur_dat$Status=factor(sur_dat$Status,levels = c("Dead","Alive"))
    sur_dat$time <- sur_dat$time
    rownames(sur_dat)=sur_dat$patientid
    sur_dat=sur_dat[fp_dat$patientid,]
    
    exp_dat=dat1[,RSF_genes]
    rownames(exp_dat)=phe$id
    exp_dat=exp_dat[fp_dat$patientid,]
    fp_dat$patientid=1:length(fp)
    sur_dat$patientid=1:length(fp)
    rownames(exp_dat)=1:length(fp)
    rownames(sur_dat)=1:length(fp)
    
    RiskGroup = ifelse(fp<median(fp_dat$fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    
    library(ggsci)
    library(scales)
    
    fp_dat$catpo <- rep(median(fp_dat$fp),length(fp_dat$patientid))
    p1 <- ggplot(fp_dat,aes(patientid))+
      geom_ribbon(aes(ymin = catpo, ymax = fp, fill = RiskGroup), alpha = 1)+
      scale_fill_manual(values = c("#FF890B","#27FF42")) +
      #scale_fill_jco()+
      #scale_colour_manual(values = c("#FF6666","#00CCCC"))+#"#FF5733", "#33FF57"
      theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
      geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
      # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
                 colour="black", linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p1
    
    pal <- pal_jco('default')(10)
    p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=Status),size = 3)+theme_bw()+
      scale_colour_manual(values = c("#FF890B","#27FF42"))+
      labs(x="Patient ID(Increasing Risk Score)",y="Survival Time(Years)")+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),colour="black", 
                 linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p2
    
    mycolors <- colorRampPalette(c("#27FF42","#FFFFFF","#FF890B"), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
                show_rownames = T,cluster_rows = F,
                fontsize_row = 14,fontsize = 14)
    p3
    plots = list(p1,p2,as.ggplot(as.grob(p3)))
    lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
    grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
    #pdf14*10
  }
}





#外部验证集km-----
{
  rm(list=ls())
  load(file = "./12.2RSF11.11/RSF_genes.Rdata")
  load(file = "./12.2RSF11.11/外部验证集.Rdata")
  head(datset3)
  datset3 <- datset3[, c("ID","OS.time","OS","risk.score", RSF_genes)]
  datset3 <- subset(datset3,OS.time > 0)
  #datset3 <- subset(datset3,OS.time < 3000)
  
  datset3$risk.group <- ifelse(datset3$risk.score > median(datset3$risk.score), 
                               "High", "Low")
  datset3$risk.group <- factor(datset3$risk.group, levels = c("Low","High"))
  surv_obj <- Surv(time = datset3$OS.time, event = datset3$OS)
  fit <- survfit(surv_obj ~ risk.group, data = datset3)
  ggsurvplot(fit,
             data = datset3,
             risk.table = F,
             pval = TRUE,
             conf.int = TRUE,
             palette = c("#FF8300", "#27FF42"),   # 蓝=Low, 红=High
             legend.title = "Risk Group",
             legend.labs = c("Low Risk","High Risk"),  # 与 factor 顺序对应
             xlab = "Time (days)",
             ylab = "Overall Survival Probability",
             title = "External validation set risk score survival curve")
  #pdf6*6
}
#师兄的代码（外部验证集km散点图热图）--------
##### 生存分析
#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
#table(TG_phe$Id %in% dat1$Id)
#phe <- TG_phe[TG_phe$Id %in% dat1$Id,]
# rm(list=ls())
# load(file = "./外部验证集.Rdata")
dat1 <- datset3
phe <- datset3
fac_mix <- datset3
{
  # factGene <- c("CYGB","CLVS1","CBY3","CPNE6")
  library(survival)
  library(survminer)
  library(ggplotify)
  library(cowplot)
  library(Hmisc)
  library(pheatmap)
  library(gridExtra)
  #s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(factGene,collapse = ' + '))))
  #model <- coxph(s, data = fac_mix )
  #summary(model,data=fac_mix)
  # RiskScore <- predict(model,type = "risk")
  RiskScore <- fac_mix$risk.score
  #names(RiskScore) = rownames(fac_mix)
  fp <- RiskScore
  phe <- fac_mix
  phe$id <- phe$ID
  
  # 生存图 
  {
    #最佳节点
    # res.cut <- surv_cutpoint(dat1, #数据集
    #                          time = "OS.time", #生存时间
    #                          event = "OS", #生存状态
    #                          variables = "exp") #需要计算的数据列名
    # 
    # summary(res.cut) #查看数据最佳截断点及统计量
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    dat1$OS.time=dat1$OS.time/365
    sfit <- survfit(Surv(OS.time, OS)~RiskGroup, data=dat1)
    # ggsurvplot(sfit, pval=TRUE,xlab ="Time(Years)",surv.median.line = "hv",pval.method = T)
    # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
    #            risk.table =TRUE,pval =TRUE,
    #            conf.int =TRUE,xlab ="Time in years",
    #            ggtheme =theme_light(),
    #            ncensor.plot = TRUE)
    ggsurvplot(sfit,
               palette = c("#FF8300", "#27FF42"),
               conf.int = T,conf.int.style='step', 
               pval = T,pval.method = T,
               risk.table = F,risk.table.pos='in',
               legend=c(0.85,0.85),
               legend.title="Risk Group",
               legend.labs=c("High","Low"),
               title="Survival Curve for Training Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 20))
    #pdf10*10
  }
  
  # 散点+热图 中位数
  {
    fp_dat=data.frame(patientid=phe$id,fp=phe$risk.score)
    fp_dat$RiskGroup= ifelse(fp_dat$fp>= median(fp_dat$fp),'High','Low')
    library(dplyr)
    fp_dat=arrange(fp_dat,fp)
    
    sur_dat=data.frame(patientid=phe$id,time=phe$OS.time/365,Status=phe$OS)
    sur_dat$Status=ifelse(sur_dat$Status==0,'Alive','Dead')
    sur_dat$Status=factor(sur_dat$Status,levels = c("Dead","Alive"))
    sur_dat$time <- sur_dat$time
    rownames(sur_dat)=sur_dat$patientid
    sur_dat=sur_dat[fp_dat$patientid,]
    
    exp_dat=dat1[,RSF_genes]
    rownames(exp_dat)=phe$id
    exp_dat=exp_dat[fp_dat$patientid,]
    fp_dat$patientid=1:length(fp)
    sur_dat$patientid=1:length(fp)
    rownames(exp_dat)=1:length(fp)
    rownames(sur_dat)=1:length(fp)
    
    RiskGroup = ifelse(fp<median(fp_dat$fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    
    library(ggsci)
    library(scales)
    
    fp_dat$catpo <- rep(median(fp_dat$fp),length(fp_dat$patientid))
    p1 <- ggplot(fp_dat,aes(patientid))+
      geom_ribbon(aes(ymin = catpo, ymax = fp, fill = RiskGroup), alpha = 1)+
      scale_fill_manual(values = c("#FF8300", "#27FF42")) +
      #scale_fill_jco()+
      #scale_colour_manual(values = c("#FF6666","#00CCCC"))+#"#FF5733", "#33FF57"
      theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
      geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
      # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
                 colour="black", linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p1
    
    pal <- pal_jco('default')(10)
    p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=Status),size = 3)+theme_bw()+
      scale_colour_manual(values = c("#FF8300", "#27FF42"))+
      labs(x="Patient ID(Increasing Risk Score)",y="Survival Time(Years)")+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),colour="black", 
                 linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p2
  
    
    mycolors <- colorRampPalette(c("#27FF42", "#FFFFFF", "#FF8300"), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
                show_rownames = T,cluster_rows = F,
                fontsize_row = 14,fontsize = 14)
    p3
    plots = list(p1,p2,as.ggplot(as.grob(p3)))
    lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
    grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
    #pdf14*10
  }
}

