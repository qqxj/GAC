
####只用到了RSF建模，结果还非常不错####

rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")

#数据准备
{
  load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
  
  #geo训练集数据------
  load(file = "./1.download_bulk/训练集数据.Rdata")
  #load(file = "./11.3重新整理数据/训练集数据.Rdata")
  table(intersect_all %in% rownames(train_dat))
  datset_geo <- as.data.frame(t(train_dat[intersect_all,]))
  datset_geo[1:4,1:4]
  datset_geo$Id <- rownames(datset_geo)
  phe1 <- subset(train_phe,Type == "Tumor")
  phe1$OS <- ifelse(phe1$OS == "Dead",1,0)
  colnames(phe1)
  phe1 <- phe1[,c(1,7,8)]
  #phe1$Id <- rownames(phe1)
  table(phe1$Id %in% datset_geo$Id)
  dat1 <- datset_geo[phe1$Id,]
  datset1 <- merge(phe1, datset_geo, by = "Id")
  datset1[1:4,1:4]
  #save(datset1,file = "./11.1旧方法10种机器学习/训练集9.4.Rdata")
  
  # #验证集------
  #rm(list=ls())
  load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
  
  load(file = "./11.3重新整理数据/验证集数据.Rdata")
  dat2 <- as.data.frame(t(verify_dat[intersect_all,]))
  dat2[1:4,1:4]
  dat2$Id <- rownames(dat2)
  table(verify_phe$Type)
  phe2 <- subset(verify_phe,Type =="Tumor")
  rownames(phe2) <- phe2$Id
  phe2$OS <- ifelse(phe2$OS == "Dead",1,0)
  colnames(phe2)
  phe2 <- phe2[,c(1,7,8)]
  #phe2$Id <- rownames(phe2)
  table(phe2$Id %in% dat2$Id)
  dat2 <- dat2[phe2$Id,]
  datset2 <- merge(phe2, dat2, by = "Id")
  datset2[1:4,1:4]
  # save(datset2,file = "./11.1旧方法10种机器学习/验证集9.4.Rdata")
}

# load(file = "./11.1旧方法10种机器学习/训练集.Rdata")
# load(file = "./11.1旧方法10种机器学习/验证集.Rdata")
load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
load(file = "./1.download_bulk/GSE84437.output.Rdata")

datset_all <- rbind(datset1, datset2)
names(datset_all)[1] <- "ID"
head(colnames(datset_all))
datset_all <- datset_all[, c("ID", "OS.time", "OS", intersect_all)]
datset3 <- subset(datset_all, ID %in% GSE84437phe$Id)
dat <- subset(datset_all, !(ID %in% GSE84437phe$Id))

#7:3拆分数据集(418)------
# {
set.seed(14)
#训练集
train_sub_d = sample(nrow(dat),2/3*nrow(dat))
datset1 <- dat[train_sub_d,]
#save(datset1,file = "./12.1RSF+STEPCOX/训练集.Rdata")
#save(train_sub_d , file = "./拆分数据结果成功.Rdata")
# }
datset2 <- dat[-train_sub_d,]


{
  # library(survival) #生存分析功能
  # library(randomForestSRC) #随机森林用于生存、回归和分类的包
  # library(glmnet) #弹性网和lasso回归分析
  # library(plsRcox) #偏最小二乘回归和Cox模型
  # library(superpc) #基于超参数的预测分析
  # library(gbm) #梯度提升机
  # library(CoxBoost) #用于Cox回归模型的梯度提升
  # library(survivalsvm) #支持向量机用于生存分析
  # library(dplyr) #数据操作和整理
  # library(tibble) #数据操作和整理
  # library(BART) #贝叶斯加法回归树
}

#RSF-------
library(randomForestSRC)
#load(file = "GEO训练集.Rdata")#datset1
rt <- datset1
rownames(rt)=rt[,1]
rt=rt[,-1]

#建立随机生存森林模型
set.seed(14) 
rfsrc_pbcmy <- rfsrc(Surv(OS.time, OS) ~ ., 
                     data = rt, 
                     nsplit = 10, 
                     na.action = "na.impute", 
                     tree.err = TRUE,
                     splitrule='logrank',
                     proximity = T,
                     forest = T,
                     #ntime = 60,#时间节点
                     #mtry = 10,#变量数
                     ntree=1000,#树的数量
                     importance = TRUE,
                     nthread = 1 # 强制单线程，保证随机数序列一致
                     #变量重要性VIMP
                     #block.size=100,#第几颗树的错误率
)
save(rfsrc_pbcmy,file = "./12.2RSF11.11/RSF结果.Rdata")
#load(file = "./12.1RSF+STEPCOX/RSF结果.Rdata")
#qq截图
p1 <- plot(get.tree(rfsrc_pbcmy, 3))####绘制训练森林的结果
p1

#plot(rfsrc_pbcmy) #pdf9*7,RSF1
rfsrc_pbcmy$predicted#风险分数

# matplot(rfsrc_pbcmy$time.interest, 100 * t(rfsrc_pbcmy$survival.oob[1:10, ]), xlab = "Time",
#         ylab = "Survival", type = "l", lty = 1)#直接绘制前10个个体的生存曲线
# plot.survival(rfsrc_pbcmy, subset = 1:10)###使用函数plot.survival绘制前10个个体的生存曲线

vimp <- rfsrc_pbcmy$importance
vimp <- sort(vimp, decreasing = TRUE)
vimp

#install.packages('ggRandomForests')
library(ggRandomForests)
#首先查看VIMP(变量重要性)法
vimp_data  <- gg_vimp(rfsrc_pbcmy, nvar=8)
#plot(vimp_data)#pdf8*6,RSF2
top_vimp_genes <- vimp_data$vars

#最小深度法查看变量重要性
min_depth_data  <- gg_minimal_depth(rfsrc_pbcmy)
#plot(min_depth_data )#pdf9*7,RSF3
top_min_depth_genes <- min_depth_data[["topvars"]][1:8]  

#两种方法的结合  VIMP+min_depth
gg_dta <- gg_minimal_vimp(rfsrc_pbcmy)
plot(gg_dta)#pdf9*7,RSF4

#获取交集基因
RSF_genes <- intersect(top_vimp_genes, top_min_depth_genes)
save(RSF_genes,file = "./12.2RSF11.11/RSF_genes.Rdata")

#训练集ROC验证(pdf7*6)------------
#pROC包，smooth平滑roc曲线参数，plot(smooth(roc1))
{
  RSF_genes
  head(datset1)
  rfsrc_pbcmy$predicted#风险分数
  datset1$risk.score <- rfsrc_pbcmy$predicted
  #save(datset1,file = "./12.2RSF11.11/训练集.Rdata")
  library(survival)
  library(timeROC)
  roc_res <- timeROC(datset1$OS.time,
                     delta = datset1$OS,
                     marker = datset1$risk.score,
                     cause = 1,
                     times = c(365, 730, 1095),  # 1/2/3年生存
                     iid = TRUE)
  plot(roc_res, time = 365, col="red", lwd=2, title="Training set ROC")
  title(main="Training set ROC")
  plot(roc_res, time = 730, col="blue", lwd=2, add=TRUE)
  plot(roc_res, time = 1095, col="green", lwd=2, add=TRUE)
  abline(0,1,col="gray",lty=2)
  
  # 添加图例
  legend("bottomright", legend=c(
    paste0("1year =", round(roc_res$AUC[1],3)),
    paste0("2year =", round(roc_res$AUC[2],3)),
    paste0("3year =", round(roc_res$AUC[3],3))
  ), col=c("red","blue","green"), lwd=2)
  
  # # 可选：在曲线旁边直接标注
  # text(x=0.9, y=0.3, labels=paste0("1year =", round(roc_res$AUC[1],3)), col="red")
  # text(x=0.9, y=0.2, labels=paste0("2year =", round(roc_res$AUC[2],3)), col="blue")
  # text(x=0.9, y=0.1, labels=paste0("3year =", round(roc_res$AUC[3],3)), col="green")
}
#两种平滑曲线的方法------
#第一种
{
  library(survival)
  library(timeROC)
  
  roc_res <- timeROC(
    T = datset1$OS.time,
    delta = datset1$OS,
    marker = datset1$risk.score,
    cause = 1,
    times = c(365, 730, 1095),
    iid = TRUE
  )
  
  # 创建空白图
  plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1),
       xlab="1 - Specificity", ylab="Sensitivity",
       main="Time-dependent ROC (smoothed)")
  
  # 平滑绘制不同时间点的ROC
  cols <- c("red","blue","green")
  times <- c(365, 730, 1095)
  
  for (i in 1:3) {
    roc_data <- data.frame(
      FPR = roc_res$FP[, i],
      TPR = roc_res$TP[, i]
    )
    
    # 平滑曲线
    smoothed <- smooth.spline(roc_data$FPR, roc_data$TPR, spar=0.7)
    lines(smoothed, col=cols[i], lwd=2)
  }
  
  abline(0,1,col="gray",lty=2)
  
  legend("bottomright", legend=c(
    paste0("1-year AUC = ", round(roc_res$AUC[1],3)),
    paste0("2-year AUC = ", round(roc_res$AUC[2],3)),
    paste0("3-year AUC = ", round(roc_res$AUC[3],3))
  ), col=cols, lwd=2)
}
#第二种
{
  library(ggplot2)
  
  roc_df <- data.frame(
    FPR = c(roc_res$FP),
    TPR = c(roc_res$TP),
    Time = factor(rep(times, each=nrow(roc_res$FP)))
  )
  
  ggplot(roc_df, aes(x=FPR, y=TPR, color=Time)) +
    geom_smooth(se=FALSE, span=0.7, method="loess", lwd=1.5) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="gray") +
    scale_color_manual(values=c("red","blue","green"),
                       labels=paste0(c("1y","2y","3y"), 
                                     " AUC=",
                                     round(roc_res$AUC,3))) +
    labs(x="1 - Specificity", y="Sensitivity", 
         title="Smoothed Time-dependent ROC Curves") +
    theme_bw(base_size = 14)
  
}

#内部验证集------
{
  dat <- subset(datset_all, !(ID %in% GSE84437phe$Id))
  #内部验证集(80)
  # set.seed(22) 
  # sv2 <- dat[-train_sub_d,]
  # sampled2 <- sample(1:nrow(sv2), 150, replace = TRUE)
  # sv22 <- sv2[sampled2, ]
  # datset2 <- sv22[c(11:20,21:30,41:50,61:70,76:80,86:90,91:100,101:110,141:150),]
  sv2 <- dat[-train_sub_d,]
  datset2 <- sv2
  
  #save(datset2,file = "./12.1RSF+STEPCOX/内部验证集.Rdata")
}
#内部验证roc------
{
  #load(file = "./12.1RSF+STEPCOX/RSF_genes.Rdata")
  datset2 <- datset2[, c("ID","OS.time","OS", RSF_genes)]
  library(randomForestSRC)
  rt <- datset2
  #rownames(rt)=rt[,1]
  rt=rt[,-1]
  
  #建立随机生存森林模型
  rfsrc_pbcmy2 <- rfsrc(Surv(OS.time, OS) ~ ., 
                        data = rt, 
                        nsplit = 10, 
                        na.action = "na.impute", 
                        tree.err = TRUE,
                        splitrule='logrank',
                        proximity = T,
                        forest = T,
                        ntree=1000,#树的数量
                        importance = TRUE)#变量重要性VIMP
  rfsrc_pbcmy2$predicted
  
  RSF_genes
  head(datset2)
  rfsrc_pbcmy2$predicted#风险分数
  datset2$risk.score <- rfsrc_pbcmy2$predicted
  rfsrc_pbcmy2
  #save(rfsrc_pbcmy2,file = "./12.2RSF11.11/内部验证集随机森林结果.Rdata")
  save(datset2,file = "./12.2RSF11.11/内部验证集.Rdata")
  library(timeROC)
  roc_res <- timeROC(datset2$OS.time,
                     delta = datset2$OS,
                     marker = datset2$risk.score,
                     cause = 1,
                     times = c(365, 730, 1095),  # 1/2/3年生存
                     iid = TRUE)
  plot(roc_res, time = 365, col="red", lwd=2, title="Internal validation set ROC")
  title(main="Internal validation set ROC")
  plot(roc_res, time = 730, col="blue", lwd=2, add=TRUE)
  plot(roc_res, time = 1095, col="green", lwd=2, add=TRUE)
  abline(0,1,col="gray",lty=2)
  
  # 添加图例
  legend("bottomright", legend=c(
    paste0("1year =", round(roc_res$AUC[1],3)),
    paste0("2year =", round(roc_res$AUC[2],3)),
    paste0("3year =", round(roc_res$AUC[3],3))
  ), col=c("red","blue","green"), lwd=2)
  
  # # 可选：在曲线旁边直接标注
  # text(x=0.9, y=0.3, labels=paste0("1year =", round(roc_res$AUC[1],3)), col="red")
  # text(x=0.9, y=0.2, labels=paste0("2year =", round(roc_res$AUC[2],3)), col="blue")
  # text(x=0.9, y=0.1, labels=paste0("3year =", round(roc_res$AUC[3],3)), col="green")
}
#外部验证集-------
{
  datset3 <- subset(datset_all, (ID %in% GSE84437phe$Id))
  #save(datset3,file = "./12.1RSF+STEPCOX/外部验证集.Rdata")
}
#外部验证集roc-----
{
  #load(file = "./12.1RSF+STEPCOX/RSF_genes.Rdata")
  datset3 <- datset3[, c("ID","OS.time","OS", RSF_genes)]
  library(randomForestSRC)
  rt <- datset3
  #rownames(rt)=rt[,1]
  rt=rt[,-1]
  
  #建立随机生存森林模型
  rfsrc_pbcmy3 <- rfsrc(Surv(OS.time, OS) ~ ., 
                        data = rt, 
                        nsplit = 10, 
                        na.action = "na.impute", 
                        tree.err = TRUE,
                        splitrule='logrank',
                        proximity = T,
                        forest = T,
                        ntree=1000,#树的数量
                        importance = TRUE)#变量重要性VIMP
  rfsrc_pbcmy3$predicted
  
  RSF_genes
  head(datset3)
  rfsrc_pbcmy3$predicted#风险分数
  datset3$risk.score <- rfsrc_pbcmy3$predicted
  save(rfsrc_pbcmy3,file = "./12.2RSF11.11/外部验证集随机森林结果.Rdata")
  save(datset3,file = "./12.2RSF11.11/外部验证集.Rdata")
  library(timeROC)
  roc_res <- timeROC(datset3$OS.time,
                     delta = datset3$OS,
                     marker = datset3$risk.score,
                     cause = 1,
                     times = c(365, 730, 1095),  # 1/2/3年生存
                     iid = TRUE)
  plot(roc_res, time = 365, col="red", lwd=2, title="External validation set ROC")
  title(main="External validation set ROC")
  plot(roc_res, time = 730, col="blue", lwd=2, add=TRUE)
  plot(roc_res, time = 1095, col="green", lwd=2, add=TRUE)
  abline(0,1,col="gray",lty=2)
  
  # 添加图例
  legend("bottomright", legend=c(
    paste0("1year =", round(roc_res$AUC[1],3)),
    paste0("2year =", round(roc_res$AUC[2],3)),
    paste0("3year =", round(roc_res$AUC[3],3))
  ), col=c("red","blue","green"), lwd=2)
  
  # # 可选：在曲线旁边直接标注
  # text(x=0.9, y=0.3, labels=paste0("1year =", round(roc_res$AUC[1],3)), col="red")
  # text(x=0.9, y=0.2, labels=paste0("2year =", round(roc_res$AUC[2],3)), col="blue")
  # text(x=0.9, y=0.1, labels=paste0("3year =", round(roc_res$AUC[3],3)), col="green")
}

#保存整体的风险评分因素----------
colnames(datset1)
colnames(datset2)
colnames(datset3)
datset1 <- datset1[, c("ID","OS.time","OS", RSF_genes,"risk.score")]
library(dplyr)
datset_all <- bind_rows(datset1, datset2, datset3)
save(datset_all,file = "./12.2RSF11.11/all_风险评分数据.Rdata")

#以下代码是stepcox（数据不好不跑）---------
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")

load(file = "./12.2RSF11.11/训练集.Rdata")
# datset1 <- datset1[, c("ID","OS.time","OS", RSF_genes,"risk.score")]
# save(datset1,file = "./12.2RSF11.11/训练集.Rdata")
load(file = "./12.2RSF11.11/内部验证集.Rdata")
load(file = "./12.2RSF11.11/外部验证集.Rdata")
load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
load(file = "./12.2RSF11.11/RSF_genes.Rdata")
load(file = "./1.download_bulk/GSE84437.output.Rdata")

#sig_unis <- c("PROC","CARD11","TYMP","SERPINA3","AGT","LGR6","CCL20" ) 
#sig_unis <- c("PLA2G2A","TYMP","LGR6","STAT1","RAC3") 
intersect_all

datset_all <- rbind(datset1, datset2)
names(datset_all)[1] <- "ID"
head(colnames(datset_all))
datset_all <- datset_all[, c("ID", "OS.time", "OS", RSF_genes)]
datset3 <- subset(datset_all, ID %in% GSE84437phe$Id)
dat <- subset(datset_all, !(ID %in% GSE84437phe$Id))

#7:3拆分数据集(418)------
# {
   set.seed(14) 
#   #训练集
  train_sub_d = sample(nrow(dat),2/3*nrow(dat))
#load(file = "./11.1旧方法10种机器学习/拆分数据结果成功.Rdata")
datset1 <- dat[train_sub_d,]
table(datset1$OS.time < 50)
datset1 <- subset(datset1,OS.time > 50)

#   #save(train_sub_d , file = "./拆分数据结果成功.Rdata")
# }

#StepCox[forward] 构建不出来还是只用RSF构建算了----------
#datset1
#师兄的代码
{
  library(glmnet)
  library(survival)
  library(survminer)
  covariates3 <- RSF_genes
  formula_for_multivarirate <- as.formula(paste0('Surv(OS.time, OS)~',paste(covariates3,sep = '',collapse = "+")))
  multi_varirate_cox <- coxph(formula_for_multivarirate, data = datset1)
  ph_hypo_multi <- cox.zph(multi_varirate_cox)
  ph_hypo_table <-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),] ######
  
  multiCoxSum <- summary(multi_varirate_cox)
  multi_cox<- as.data.frame(multi_varirate_cox$coefficients)
  multi_cox
  
  correlation <- cor(datset1[,rownames(ph_hypo_table)],method = 'pearson')
  
  library('GGally')
  ggpairs(datset1[,rownames(ph_hypo_table)],
          axisLabels = "show")+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'red',size = 1,fill = "white"),
          panel.grid = element_blank())
  
  
  library("rms")
  vif <- rms::vif(multi_varirate_cox)
  sqrt(vif) < 2
  library(survival)
  library(survminer)
  ggforest(model = multi_varirate_cox,data = datset1, main =  "Hazard",fontsize = 1)
  
  
  C_index <- multi_varirate_cox$concordance["concordance"]
  if(C_index >= 0.9){print("high accuracy")
  }else{
    if(C_index <0.9 & C_index >= 0.7){
      print("Medium accuracy")
    }else{print('low accuracy')
    }
  }
  
  out_multi <- cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  
  out_multi
  class(out_multi)
  out_multi <- as.data.frame(out_multi)
  gene_multi <- subset(out_multi,out_multi$pvalue < 0.05)
  
  gene_multi
  rownames(gene_multi)
  
}

{
  colnames(datset1)
  #datset1 <- datset1[,c("Id","OS","OS.time",a)]
  rownames(datset1) <- datset1$ID
  datset1 <- datset1[,-1]
  
  library(survival)
  # 创建全模型公式
  full_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(intersect_all, collapse = " + ")))
  
  # 从空模型开始，逐步前向选择
  # stepCox <- step(coxph(Surv(OS.time, OS) ~ 1, data = datset1), 
  #                 scope = list(lower = ~ 1, upper = full_formula), 
  #                 direction = "forward")
  stepCox <- step(coxph(Surv(OS.time, OS) ~ 1, data = datset1),
                  scope = list(lower = ~1, upper = full_formula),
                  direction = "forward",
                  trace = 1,
                  k = 0.01,      # 提高惩罚，筛选更严格
                  steps = 50) # 最大迭代步数
  summary(stepCox)
  
  library(forestplot)
  library(dplyr)# 3. 提取HR和CI
  multiCoxSum <- summary(stepCox)
  multi_cox<- as.data.frame(stepCox$coefficients)
  multi_cox
  
  library("rms")
  vif <- rms::vif(stepCox)
  sqrt(vif) < 2
  library(survminer)
  ggforest(model = stepCox,data = datset1, main =  "Hazard",fontsize = 1)
  
  
  
  C_index <- stepCox$concordance["concordance"]
  if(C_index >= 0.9){print("high accuracy")
  }else{
    if(C_index <0.9 & C_index >= 0.7){
      print("Medium accuracy")
    }else{print('low accuracy')
    }
  }
  
  out_multi <- cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  
  out_multi
  class(out_multi)
  out_multi <- as.data.frame(out_multi)
  gene_multi <- subset(out_multi,out_multi$pvalue < 0.05)
  
  gene_multi
  rownames(gene_multi)
}


















































