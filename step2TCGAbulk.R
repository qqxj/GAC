#下载bulk数据：胃癌
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
gc()

#TCGA数据----
 #407例
{
  setwd("/opt/disease/TCGA/STAD/")
  count <- read.table(file = "../STAD/TCGA-STAD.htseq_counts.tsv.gz",
                          header = T,row.names = 1)
  phe <- read.delim(file = "../STAD/TCGA-STAD.GDC_phenotype.tsv.gz")
  dim(count)
  dim(phe)
  setwd("/home/datahup/syj/GAC/1.download_bulk/")
  # count[1:4,1:4]
  # range(count)
  colnames(phe)
  # table(phe$primary_diagnosis.diagnoses)
  # #不需要做亚型划分，胃癌中95%的癌症为胃腺癌，胃腺癌包括了腺型、肠型、弥漫型、印戒细胞癌、黏液腺癌等
  # #5%非上皮源性肿瘤（罕见）：胃淋巴瘤、胃肠间质瘤（GIST）、神经内分泌肿瘤（如胃类癌）
  # table(phe$sample_type.samples)#443例癌症、101例正常组
  # phe$days_to_death.demographic
  # table(phe$vital_status.demographic)
  
  # 处理表达信息
  # ID转换
  ids <- read.delim(file = "./gencode.v22.annotation.gene.probeMap")
  ids1 <- data.frame(ids$id,ids$gene) #取出ID跟Gene对应关系
  names(ids1) <- c("id","gene")
  dat = count
  dat$id <- rownames(count) #取出Ensenmble Id这列
  
  head(ids1)
  table(ids1$id %in% dat$id)
  library(dplyr)
  dat <- inner_join(ids1,dat, by="id")#内连接函数
  length(dat$gene)
  length(unique(dat$gene)) # 有重复Gene
  dat <- dat[!duplicated(dat$gene),] # 去除重复gene
  rownames(dat) <- dat$gene
  dat <- dat[,-c(1,2)]
  
  #除去表达矩阵中没有的临床样本
  phe$submitter_id.samples <- gsub("-",".",phe$submitter_id.samples)#将-修改为.，与表达矩阵中的点与临床数据相对应，好进行取交集
  dat <- dat[,colnames(dat) %in% phe$submitter_id.samples]#取交集,改变表达矩阵的列于临床数据对应
  phe <- phe[phe$submitter_id.samples %in% colnames(dat),]#取交集改变临床的列与表达矩阵对应
  
  # 处理临床信息
  # 添加分组信息
  rownames(phe) <- phe[,1]
  group_text <- phe
  group_text <- group_text[c(-1)]   # #改变行名
  pd <- data.frame(group_text$submitter_id,
                   group_text$sample_type_id.samples)
  table(pd$group_text.sample_type_id.samples)
  pd[pd==6] <- 1
  table(pd$group_text.sample_type_id.samples)
  
  Group<-factor(pd$group_text.sample_type_id.samples,levels=c('1','11'))  ## 11 -normal,  1-toumal
  design<-model.matrix(~0+Group)
  colnames(design)<-c('Tumor','Normal')
  rownames(design)<-rownames(group_text)
  
  group_list = design[colnames(dat),]# 按dat列名排序
  
  library(stringr)
  phe <- group_text
  group = str_split(as.character(phe$sample_type.samples),' ',simplify = T)[,1]  
  table(group)
  
  phe2 <- data.frame(rownames(phe),
                     phe$age_at_initial_pathologic_diagnosis,
                     phe$pathologic_M,
                     phe$pathologic_N,
                     phe$pathologic_T,
                     phe$gender.demographic,
                     phe$days_to_death.demographic,
                     phe$days_to_last_follow_up.diagnoses,
                     phe$vital_status.demographic,
                     phe$tumor_stage.diagnoses,
                     phe$sample_type.samples)
  names(phe2) <- c("Id","Age","M","N",
                   "T","Gender","Days to Death","Days to Last Follow",
                   "OS","Stage","Type")
  
  # # 处理临床信息   暂时不需要进行亚型分类
  # # 根据临床分子亚型对样本进行区分
  # # if (!require("BiocManager", quietly = TRUE))
  # #   install.packages("BiocManager")
  # # 
  # # BiocManager::install("TCGAbiolinks")
  # library(TCGAbiolinks)
  # subtypes <- PanCancerAtlas_subtypes()
  # BRCA <- subset(subtypes,subtypes$cancer.type == "BRCA")
  # table(BRCA$Subtype_mRNA)
  # BRCA$pan.samplesID <- gsub("-",".",BRCA$pan.samplesID) # gsub()替换掉字符串中所有查找到的指定字符
  # BRCA$ID = substr(BRCA$pan.samplesID, 1,16) #substr() 能对给定的字符串对象取出子集，其参数是子集所处的起始和终止位置
  # SubType <- data.frame(SubType = BRCA$Subtype_mRNA, Id = BRCA$ID,  row.names = BRCA$ID)
  # phe2 <- merge(phe2,SubType,by = "Id")
  TCGAdat <- dat
  TCGAphe <- phe2
  save(group_list,TCGAdat,TCGAphe,
       file = "./Step1 TCGA.output.Rdata")
}
load(file = "./1.download_bulk/Step1 TCGA.output.Rdata")
#GEO------
 #GSE84437  433例
{ rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.download_bulk/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
   #eSet <- getGEO(filename = "./GSE84437_series_matrix.txt.gz")
  #不清楚什么情况，手动下载好矩阵你在重新使用getGEO下载就好了
  eSet <- getGEO('GSE84437', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)

  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  colnames(phe)
  phe <- data.frame(Id = rownames(phe),
                    Age = phe$`age:ch1`,
                    Gender = phe$`Sex:ch1`,
                    OS = phe$`death:ch1` ,
                    #Grade = phe$`grade:ch1`,
                    OS.time = phe$`duration overall survival:ch1`,
                    N = phe$`pnstage:ch1`,
                    T = phe$`ptstage:ch1`,
                    Type = phe$`tissue:ch1`)
  str(phe)
  phe <- na.omit(phe)
  phe$OS.time <- as.numeric(as.character(phe$OS.time)) 
  phe$OS.time <- phe$OS.time * 30
  phe$OS <- ifelse(phe$OS == "1","Dead","Alive")
  phe$Type <- ifelse(phe$Type == "gastric cancer","Tumor","NA")
  table(phe$Type)
  
  
  # GSE42568_phe <- GSE42568_phe[GSE42568_phe$Tissue == "breast cancer",]
  expr <- as.data.frame(raw_exprSet[,phe$Id])
  dim(expr)
  expr <- na.omit(expr)
  range(expr)
  #[1]    30.39947 38105.52035
  expr <- log2(expr + 1)
  
  head(rownames(expr))
  # https://blog.csdn.net/weixin_40739969/article/details/103186027
  #BiocManager::install("illuminaHumanv3.db")
  library(illuminaHumanv3.db)
  ls("package:illuminaHumanv3.db")
  ids=toTable(illuminaHumanv3SYMBOL) #toTable这个函数：提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
  head(ids) #head为查看前六行
  
  dat = expr
  table(ids$probe_id %in%  rownames(dat))
  # FALSE  TRUE 
  # 64 29367 
  ids=ids[ids$probe_id %in%  rownames(dat),]
  dat=dat[rownames(dat) %in% ids$probe_id,] 
  
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
  
  GSE84437dat <- dat
  #19192基因 * 433样本
  GSE84437phe <- phe
  save(GSE84437dat,GSE84437phe,file = "./GSE84437.output.Rdata")
}
load(file = "./1.download_bulk/GSE84437.output.Rdata")
 #GSE66229:ACRG 370例
{ 
  rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.download_bulk/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
  eSet <- getGEO('GSE66229', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  
  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  colnames(phe)
  #没有患者临床数据，PMID: 25894828有临床数据，但还是只有肿瘤的，
  table(phe$characteristics_ch1)
  
  phe_normal <- subset(phe,characteristics_ch1 == "tissue: normal gastric tissue")
  phe_normal$characteristics_ch1 <- ifelse(phe_normal$characteristics_ch1 == "tissue: normal gastric tissue","Normal","Tumor")
  phe_normal <- data.frame(Id = rownames(phe_normal),
                           Type = phe_normal$characteristics_ch1,
                           Tumor_ID = phe_normal$`patient:ch1`)
  
  phe_tumor <- subset(phe,characteristics_ch1 == "tissue: Gastric tumor")
  phe_tumor$characteristics_ch1 <- ifelse(phe_tumor$characteristics_ch1 == "tissue: Gastric tumor","Tumor","Normal")
  phe_tumor <- data.frame(Id = rownames(phe_tumor),
                          Type = phe_tumor$characteristics_ch1,
                          Tumor_ID = phe_tumor$`patient:ch1`)
  phe_tumor$Tumor_ID <- as.numeric(phe_tumor$Tumor_ID)
  
  #文章pmid25894828有临床数据
  library(readxl)  
  library(openxlsx)  
  pmid25894828 <- read_excel("./pmid25894828.xls",col_names = TRUE)
  colnames(pmid25894828 )
  phe25894828 <- data.frame(Tumor_ID = pmid25894828$`Tumor ID`,
                            Age = pmid25894828$age,
                            Gender = pmid25894828$sex,
                            Stage = pmid25894828$pStage,
                            T = pmid25894828$T,
                            N = pmid25894828$N,
                            M = pmid25894828$M,
                            OS.time = pmid25894828$"OS\n(months)" ,
                            OS = pmid25894828$`FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss`)
  table(phe25894828$Gender)
  phe25894828$Gender <- ifelse(phe25894828$Gender == "M","male","female")
  class(phe25894828$OS.time)
  phe25894828$OS.time <- phe25894828$OS.time * 30
  table(phe25894828$OS)
  phe25894828 <- subset(phe25894828, OS != 4 & OS != 5)
  phe25894828$OS <- ifelse(phe25894828$OS %in% c(0, 1), 0,
                                  ifelse(phe25894828$OS %in% c(2, 3), 1, NA))
  save(phe25894828,file = "./phe25894828.Rdata")
  
  table(phe25894828$Tumor_ID %in% phe_tumor$Tumor_ID)
  phe_tumor <- phe_tumor[phe_tumor$Tumor_ID %in% phe25894828$Tumor_ID,]
  str(phe_tumor$Tumor_ID)
  str(phe25894828$Tumor_ID)
  phe2 <- merge(phe_tumor, phe25894828, by = "Tumor_ID")
  rownames(phe2) <- phe2$Id
  
  table(phe_normal$Tumor_ID %in% phe25894828$Tumor_ID)
  phe_normal <- phe_normal[phe_normal$Tumor_ID %in% phe25894828$Tumor_ID,]
  phe25894828_1 <- phe25894828[phe25894828$Tumor_ID %in% phe_normal$Tumor_ID,]
  phe3 <- merge(phe_normal, phe25894828_1, by = "Tumor_ID")
  rownames(phe3) <- phe3$Id
  
  #合并
  str(phe3)
  str(phe2)
  phe4 <- rbind(phe2, phe3)
  phe4 <- subset(phe4, !(Id %in% c("GSM1523926", "GSM1524080", "GSM1523942")))
  
  #raw_exprSet <- as.data.frame(raw_exprSet[,phe$Id])
  raw_exprSet[1:4,1:4]
  # ID转换 GPL570 对应R包：hgu133plus2
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db")
  ids=toTable(hgu133plus2SYMBOL) #toTable这个函数：提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
  head(ids) #head为查看前六行
  
  dat = raw_exprSet
  table(ids$probe_id %in%  rownames(dat))
  # TRUE 
  # 43101
  ids=ids[ids$probe_id %in%  rownames(dat),]
  dat[1:4,1:4]   
  dat=dat[rownames(dat) %in% ids$probe_id,] 
  
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
  
  dat <- as.data.frame(dat)
  dat <- dat[,phe4$Id]
  
  GSE66229dat <- dat
  #20824基因 * 100正常样本
  GSE66229phe <- phe4
  save(GSE66229dat,GSE66229phe,file = "./GSE66229.output.Rdata")
}
load(file = "./1.download_bulk/GSE66229.output.Rdata")
#GSE54129
#GSE26942
#GSE62254  300
#GSE15459 
#GSE84433  373
#GSE26901

#去除批次效应-----
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
gc()
load(file = "./1.download_bulk/Step1 TCGA.output.Rdata")
load(file = "./1.download_bulk/GSE84437.output.Rdata")
load(file = "./1.download_bulk/GSE66229.output.Rdata")
{
  #首先合并表达矩阵
  common_genes <- Reduce(intersect, list(
    rownames(GSE66229dat),
    rownames(GSE84437dat),
    rownames(TCGAdat)
  ))
  length(common_genes)
  
  GSE66229dat_1 <- GSE66229dat[common_genes, ]
  GSE84437dat_1 <- GSE84437dat[common_genes, ]
  TCGAdat_1 <- TCGAdat[common_genes, ]
  
  #合并
  merged_expr <- cbind(GSE66229dat_1, GSE84437dat_1, TCGAdat_1)
  table(duplicated(colnames(merged_expr)))#没有重复的样本
  
  #查看拼接点表达量高低#保存pdf格式8*6
  pdf("./1.download_bulk/GSE66229和GSE84437拼接点数据表达情况.pdf", width = 10, height = 8)
  par(mar = c(8, 4, 3, 2))  
  boxplot(merged_expr[, 360:380], las = 2,
          main = "Expression of splice junction data for GSE66229 and GSE84437")
  dev.off()
  
  boxplot(merged_expr[,795:815],las="2")
  pdf("./1.download_bulk/TCGA和GSE84437拼接点数据表达情况.pdf", width = 10, height = 8)
  par(mar = c(10, 4, 3, 2))  
  boxplot(merged_expr[,795:815],las="2",
          main = "Expression of splice junction data for GSE84437 and TCGA")
  dev.off()
  
  #创建分组信息
  group <- data.frame(row.names = colnames(merged_expr),
                      Sample=colnames(merged_expr))
  # 按样本数量指定分组标签
  group$Dataset <- c(
    rep("GSE66229", 370),
    rep("GSE84437", 803 - 370),
    rep("TCGA", nrow(group) - 803)
  )
  
  library(factoextra)
  library(FactoMineR)
  #dat_group$DataSet[dat_group$DataSet == "GTEx"] <- "GSE53624"
  before_exp_pca <- PCA(t(merged_expr),
                        scale.unit=T,ncp=5,graph=F)
  before_exp_pca.plot <- fviz_pca_ind(before_exp_pca,
                                     axes=c(1,2),
                                     label="none",
                                     addEllipses = T,
                                     ellipse.level=0.9,
                                     habillage = factor(group$Dataset),
                                     palette = "aaas",
                                     mean.point=F,
                                     title="")
  #保存pdf为8*6
  before_exp_pca.plot
  
  library(limma)
  dat_exp <- removeBatchEffect(merged_expr[,rownames(group)],
                               batch = group$Dataset)
  
  #查看拼接点表达量高低#保存pdf格式8*6
  pdf("./1.download_bulk/GSE66229和GSE84437拼接点数据表达情况(去批次后).pdf", width = 10, height = 8)
  par(mar = c(8, 4, 3, 2))  
  boxplot(dat_exp[, 360:380], las = 2,
          main = "ESplicing point expression of GSE66229 and GSE84437 after de-batching")
  dev.off()
  
  pdf("./1.download_bulk/TCGA和GSE84437拼接点数据表达情况(去批次后).pdf", width = 10, height = 8)
  par(mar = c(10, 4, 3, 2))  
  boxplot(dat_exp[,795:815],las="2",
          main = "ESplicing point expression of GSE84437 and TCGA after de-batching")
  dev.off()
  
   #dat_exp <- as.matrix(dat)
  #去批次后PCA
  library(factoextra)
  library(FactoMineR)
  #dat_group$DataSet[dat_group$DataSet == "GTEx"] <- "GSE53624"
  before_exp_pca <- PCA(t(dat_exp),
                        scale.unit=T,ncp=5,graph=F)
  before_exp_pca.plot <- fviz_pca_ind(before_exp_pca,
                                      axes=c(1,2),
                                      label="none",
                                      addEllipses = T,
                                      ellipse.level=0.9,
                                      habillage = factor(group$Dataset),
                                      palette = "aaas",
                                      mean.point=F,
                                      title="")
  #保存pdf为8*6
  before_exp_pca.plot
  
  dat <- as.data.frame(dat_exp)
  save(dat,file = "./1.download_bulk/dat.output.Rdata")
}
load(file = "./1.download_bulk/dat.output.Rdata")

#处理临床信息--------
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
gc()
load(file = "./1.download_bulk/Step1 TCGA.output.Rdata")
load(file = "./1.download_bulk/GSE84437.output.Rdata")
load(file = "./1.download_bulk/GSE66229.output.Rdata")

{
  #顺序不要变
  #  merged_expr <- cbind(GSE66229dat_1, GSE84437dat_1, TCGAdat_1)
  colnames(GSE66229phe)
  colnames(GSE84437phe)
  colnames(TCGAphe)
  #处理GSE66229phe
  {
    GSE66229phe1 <- GSE66229phe[,-c(1,6,9)]
    table(GSE66229phe1$OS)
    GSE66229phe1$OS <- ifelse(GSE66229phe1$OS == "1","Dead","Alive")
  }
  #处理GSE84437phe
  {
    table(GSE84437phe$T)
    GSE84437phe$T <- ifelse(GSE84437phe$T == "T1", 1,
                            ifelse(GSE84437phe$T == "T2", 2,
                                   ifelse(GSE84437phe$T == "T3", 3, 4)))
    
    table(GSE84437phe$N)
    GSE84437phe$N <- ifelse(GSE84437phe$N == "N0", 0,
                            ifelse(GSE84437phe$N == "N1", 1,
                                   ifelse(GSE84437phe$N == "N2", 2, 3)))
  }
  #处理TCGAphe
  {
    TCGAphe$OS.time <- ifelse(is.na(TCGAphe$`Days to Death`),TCGAphe$`Days to Last Follow`,TCGAphe$`Days to Death`)
    TCGAphe1 <- TCGAphe[,-c(3,7,8,10)]
    table(TCGAphe1$N)
    TCGAphe1$N <- ifelse(TCGAphe1$N %in% c("N0"), 0,
                         ifelse(TCGAphe1$N %in% c("N1"), 1,
                                ifelse(TCGAphe1$N %in% c("N2"), 2,
                                       ifelse(TCGAphe1$N %in% c("N3", "N3a", "N3b"), 3, "X"))))
    table(TCGAphe1$T)
    TCGAphe1$T <- ifelse(TCGAphe1$T %in% c("T1", "T1a", "T1b"), 1,
                         ifelse(TCGAphe1$T %in% c("T2", "T2a", "T2b"), 2,
                                ifelse(TCGAphe1$T %in% c("T3"), 3,
                                       ifelse(TCGAphe1$T %in% c("T4", "T4a", "T4b"), 4, "X"))))
    table(TCGAphe1$Type)
    TCGAphe1$Type <- ifelse(TCGAphe1$Type == "Primary Tumor","Tumor","Normal")
    TCGAphe1 <- na.omit(TCGAphe1)
    table(TCGAphe1$OS)
    TCGAphe1 <- TCGAphe1[TCGAphe1$OS != "Not Reported", ]
    
  }
  
  colnames(GSE66229phe1)
  colnames(GSE84437phe)
  colnames(TCGAphe1)
  
  # 统一列顺序（按 GSE66229phe1 的列顺序来排列其他两个数据框）
  GSE84437phe <- GSE84437phe[, colnames(GSE66229phe1)]
  TCGAphe1 <- TCGAphe1[, colnames(GSE66229phe1)]
  # 合并三个数据框
  phe <- rbind(GSE66229phe1, GSE84437phe, TCGAphe1)
  save(phe,file = "./1.download_bulk/phe.output.Rdata")
}
load(file = "./1.download_bulk/phe.output.Rdata")

#整理出训练集和验证集-----
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
gc()

load(file = "./1.download_bulk/dat.output.Rdata")
load(file = "./1.download_bulk/phe.output.Rdata")

#验证集geo
{
  head(phe$Id)
  tail(phe$Id)
  verify_phe <- phe[grepl("^GSM", phe$Id), ]
  verify_dat <- dat[,verify_phe$Id]
  verify_dat[1:4,1:4]
  save(verify_dat,verify_phe,file = "./1.download_bulk/验证集数据.Rdata")
}
load(file = "./1.download_bulk/验证集数据.Rdata")
#训练集TCGA
{
  train_phe <- phe[grepl("^TCGA", phe$Id), ]
  train_dat <- dat[,train_phe$Id]
  save(train_dat,train_phe,file = "./1.download_bulk/训练集数据.Rdata")
}
load(file = "./1.download_bulk/训练集数据.Rdata")




