#100种机器学习算法
#https://blog.csdn.net/zfyyzhys/article/details/141048663#:~:text=Golang%20%E4%BD%9C%E4%B8%BA%E4%B8%80

rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")

#数据准备
{
  load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
  
  #geo训练集数据------
  load(file = "./1.download_bulk/训练集数据.Rdata")
  
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
  #save(datset1,file = "./11.5.101次机器学习1028/训练集1028.Rdata")
  
  # #验证集------
  rm(list=ls())
  load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
  
  load(file = "./1.download_bulk/验证集数据.Rdata")
  #load(file = "./11.3重新整理数据/验证集数据.Rdata")
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
  #save(datset2,file = "./11.5.101次机器学习1028/验证集1028.Rdata")
}

library(CoxBoost)
library(fastAdaboost)
library(Mime1)

rm(list=ls())
setwd("/home/datahup/syj/GAC/")

load(file = "./11.5.101次机器学习1028/训练集1028.Rdata")
load(file = "./11.5.101次机器学习1028/验证集1028.Rdata")
load(file = "/home/datahup/syj/GAC/9.韦恩图/交集基因1028.Rdata")
#外部验证集
load(file = "./1.download_bulk/GSE84437.output.Rdata")

genelist <- data.frame(gene = intersect_all)
genelist <- as.character(genelist$gene)

library(dplyr)

head(colnames(datset1))#顺序id、time、os、基因
head(colnames(datset2))#顺序id、time、os、基因
datset_all <- rbind(datset1, datset2)
#datset1 <- datset_all

#datset1 <- datset1[, c(1, 3, 2, 4:ncol(datset1))]
names(datset_all)[1] <- "ID"
#datset2 <- datset2[, c(1, 3, 2, 4:ncol(datset2))]
#names(datset2)[1] <- "ID"

# datset2$ID
# is_tcga <- grepl("^TCGA", datset2$ID)#排序数据
# datset2_sorted <- datset2[order(!is_tcga), ]
# range(datset2_sorted$OS.time)
# range(datset1$OS.time)
# datset2_sorted <- subset(datset2_sorted,datset2_sorted$OS.time > 0)
datset_all <- subset(datset_all,datset_all$OS.time > 0)
datset_all$OS.time <- datset_all$OS.time/365
#datset2_sorted$OS.time <- datset2_sorted$OS.time/365
# datset2 <- subset(datset2,datset2$OS.time > 0)
# datset2 <- subset(datset2,datset2$OS.time < 4000)


datset3 <- subset(datset_all, ID %in% GSE84437phe$Id)
dat <- subset(datset_all, !(ID %in% GSE84437phe$Id))
#dat <- datset1



list_train_vali_Data <- list(Training_set = datset1, 
                             Internal_validation_set = datset2,
                             External_validation_set = datset3)

# library(rms)
# dd <- datadist(list_train_vali_Data$Training_set)  # 使用 Training_set 定义数据分布对象
# options(datadist = "dd")  # 设置全局选项
# library(RSpectra)
# setwd("/home/datahup/syj/GAC/11.5.101次机器学习1028/")

#开始建模------- 
{
  res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Training_set,
                         list_train_vali_Data = list_train_vali_Data,
                         unicox.filter.for.candi = T,
                         unicox_p_cutoff = 0.5,
                         candidate_genes = genelist,
                         mode = 'all',nodesize =5,seed = 5201314 )
  
  # 只跑一个
  #?ML.Dev.Prog.Sig
  # res1 <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Training_set,
  #                         list_train_vali_Data = list_train_vali_Data,
  #                         candidate_genes = genelist,#筛选的基因集
  #                         unicox.filter.for.candi = TRUE,
  #                         unicox_p_cutoff = 0.5,
  #                         mode = 'single',
  #                         single_ml = "RSF",#Lasso RSF
  #                         nodesize = 5,
  #                         seed = 5201314 )
  #res1
  #save(res,file = "./100种机器学习10.14.Rdata")
}

#?cindex_dis_all
p1 <- cindex_dis_all(res,
                     #color = c("#1f77b4", "#ff7f0e", "#2ca02c","black","yellow"),#五个颜色
                     validate_set = names(list_train_vali_Data)[-1],
                     order = names(list_train_vali_Data),
                     width = 0.5)#,height = 0.5
#保存pdf
ggsave("./11.5.101次机器学习1028/cindex_plot11.11.pdf", plot = p1, width = 6, height = 12)

##提取风险指数
# train <- res1$riskscore$RSF$Training_set
# a[1:4,]
# validation <- res1$riskscore$RSF$External_validation_set


#km验证
# {
#   # 加载包
#   library(survival)
#   library(survminer)
#   
#   # 提取风险评分数据
#   risk_data <- res1[["riskscore"]][["RSF"]][["Training_set"]]
#   
#   # 按中位数分组
#   median_rs <- median(risk_data$RS, na.rm = TRUE)
#   risk_data$group <- ifelse(risk_data$RS > median_rs, "High", "Low")
#   
#   # 构建生存对象
#   surv_obj <- Surv(time = risk_data$OS.time, event = risk_data$OS)
#   
#   # 拟合 Kaplan-Meier 曲线
#   fit <- survfit(surv_obj ~ group, data = risk_data)
#   
#   # 绘图
#   ggsurvplot(fit, data = risk_data,
#              pval = TRUE,              # 显示 p 值
#              risk.table = TRUE,        # 显示风险表
#              surv.median.line = "hv",  # 添加中位生存线
#              palette = c("#E7B800", "#2E9FDF"))
#   
# }











