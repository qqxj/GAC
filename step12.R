
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
#模型基因
#RSF_genes <- c("IFNG","RAC3","AGT","SLC11A1","INHA","CARD11","LGR5","TYMP","FGFR4")
load(file = "./12.2RSF11.11/RSF_genes.Rdata")

#列线图----------
# nomaogram(诺莫图)
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2) 
library(survival)
library(rms)

# load(file = "/home/datahup/syj/ESCC/机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
# load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
load(file = "./1.download_bulk/dat.output.Rdata")
load(file = "./1.download_bulk/phe.output.Rdata")
load(file = "./12.2RSF11.11/训练集.Rdata")

phe <- subset(phe,Id %in% datset1$ID)
dat <- dat[,phe$Id]

colnames(phe)
train_phe=phe[,c(1,3,4,5,6,7,8)]
train_phe$id = phe$Id
risk = data.frame(row.names = datset1$ID, 
                  RiskScore = datset1$risk.score, 
                  Id = datset1$ID)
rownames(train_phe) <- train_phe$Id
train_phe <- train_phe[rownames(risk),]
train_phe = merge(train_phe, risk , by = 'Id')
rownames(train_phe) = train_phe$Id
train_phe = train_phe[,-c(1,8)]
colnames(train_phe)
# 自定义排序顺序
#custom_order <- c("Age", "Gender", "tnm_stage", "OS", "OS.time", "RiskScore")
custom_order <- c("Age", "Gender", "T","N", "OS", "OS.time", "RiskScore")
train_phe <- train_phe[, custom_order]
train_phe$OS <- ifelse(train_phe$OS == "Alive",0,1)
colnames(train_phe)[colnames(train_phe) == "T"] <- "T_stage"
colnames(train_phe)[colnames(train_phe) == "N"] <- "N_stage"
#save(train_phe,file = "./列线图数据.RData")

{
  Nmdat <- train_phe
  
  # Nmdat <- Nmdat[,c('RiskScore','OS','OS.time')]
  # Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste('RiskScore', sep = '', collapse = '+')))
  
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,
                   fun = list(surv1,surv3,surv5),
                   lp = T,
                   funlabel = c('1-year survival Probability',
                                '3-year survival Probability',
                                '5-year survival Probability'),
                   maxscale = 100,
                   fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  
  #pdf(file = '../../Fig/Step4-1/nomaogram.pdf',width = 8,height = 6)
  plot(nomo,
       lplabel = 'Linear Preadictor',
       xfrac = .35,
       varname.label = T,
       varname.label.sep = '=',
       ia.space = .2,
       tck = NA,
       tcl = 0.2,
       lmgp = 0.3,
       points.label = 'Points',
       total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  #dev.off()
  
  #列线图优化添加进入其他的临床数据#pdf8*6
  {
    library(rms)
    library(survival)
    
    #1 准备数据
    Nmdat <- train_phe
    
    Nmdat <- subset(Nmdat,Nmdat$T_stage  != "X")
    Nmdat <- subset(Nmdat,Nmdat$N_stage  != "X")
    Nmdat$Gender <- ifelse(Nmdat$Gender == "male", "Male", "Female")
    # 若 OS.time 单位是月，可转为天（比如×30），根据你数据实际情况决定
    # Nmdat$OS.time <- Nmdat$OS.time * 30
    # 转换类型
    Nmdat$Age <- as.numeric(Nmdat$Age)        # 连续变量
    Nmdat$Gender <- factor(Nmdat$Gender)      # 因子变量
    Nmdat$T_stage <- factor(Nmdat$T_stage)
    Nmdat$N_stage <- factor(Nmdat$N_stage)
    
    # 2 设置数据分布信息（rms模型必需）
    dd <- datadist(Nmdat)
    options(datadist = "dd")
    
    # 3 构建多因素 Cox 模型
    multivarl <- as.formula(Surv(OS.time, OS) ~ Age + Gender + T_stage + N_stage + RiskScore)
    cox_model <- cph(multivarl,data = Nmdat,
      x = TRUE,y = TRUE,surv = TRUE,time.inc = 365)
    
    # 4 构建生存函数
    surv <- Survival(cox_model)
    surv1 <- function(x) surv(1*365, x)
    surv3 <- function(x) surv(3*365, x)
    surv5 <- function(x) surv(5*365, x)
    
    # 5 绘制列线图
    nomo <- nomogram(
      cox_model,
      fun = list(surv1, surv3, surv5),
      lp = TRUE,
      funlabel = c('1-year survival Probability', 
                   '3-year survival Probability',
                   '5-year survival Probability'),
      maxscale = 100,
      fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
    )
    
    plot(
      nomo,
      lplabel = 'Linear Predictor',
      xfrac = 0.35,
      varname.label = TRUE,
      varname.label.sep = '=',
      ia.space = 0.2,
      tck = NA,
      tcl = 0.2,
      lmgp = 0.3,
      points.label = 'Points',
      total.points.label = 'Total Points',
      total.sep.page = FALSE,
      cap.labels = FALSE,
      cex.var = 0.6,
      cex.axis = 0.6,
      lwd = 0.8,
      label.every = 1,
      col.grid = gray(c(0.8, 0.95))
    )
  }
  
  ## 校准曲线
  ## 参数说明：
  ## 1、绘制校正曲线前需要在模型函数中添加参数x=T, y=T，详细内容参考帮助
  ## 2、u需要与之前模型中定义好的time.inc一致，即365或730；
  ## 3、m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）
  ## 而m代表每组的样本量数，因此m*3应该等于或近似等于样本量；
  ## 4、b代表最大再抽样的样本量
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=139, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=139, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=139, B=1000)
  
  #pdf(file = '../Figer/nomaogram_calibrate.pdf',width = 8,height = 6)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = 'blue',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = 'blue',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = 'green',
       col = 'green',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  abline(0,1,lty=1,lwd=1)
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("blue","green","#FF0033"),
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  #dev.off()
}


#美化列线图
{
  library(survival)
  colnames(train_phe)
  library(dplyr)
  train_phe <- train_phe %>% filter(T_stage != "X")
  train_phe <- train_phe %>% filter(N_stage != "X")
  train_phe$Age <- ifelse(train_phe$Age < 30, "Young",ifelse(train_phe$Age <= 60, "Middle-aged", "Old"))
  
  Coxfit<-coxph(Surv(OS.time,OS==1)~Age+Gender+T_stage+N_stage+RiskScore,data=train_phe)  
  library(regplot)
  regplot(Coxfit, plots=c("density","boxes") ,
          observation=FALSE,points=TRUE,
          title="Survival Nomogram",
          failtime=c(366,731,1827),
          dencol="#ADD8E6",boxcol="#FFCCCB",droplines=TRUE)
  regplot(Coxfit, plots=c("bean","boxes"), 
          observation=FALSE,points=TRUE,#observation=survdata[20,],
          title="Survival Nomogram", failtime=c(365,730,1826),
          prfail=F, clickable=TRUE,  
          dencol="green", boxcol="yellow",droplines=TRUE)
}



























































































