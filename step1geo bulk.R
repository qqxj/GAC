#GSE15459 没有临床数据
{ 
  rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.1下载其他的geo数据/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
  #eSet <- getGEO(filename = "./GSE84437_series_matrix.txt.gz")
  #不清楚什么情况，手动下载好矩阵你在重新使用getGEO下载就好了
  eSet <- getGEO('GSE15459', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  
  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  boxplot(raw_exprSet[,1:50],las="2")
  
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

#GSE15460 没有临床数据
{ 
  rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.1下载其他的geo数据/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
  #eSet <- getGEO(filename = "./GSE84437_series_matrix.txt.gz")
  #不清楚什么情况，手动下载好矩阵你在重新使用getGEO下载就好了
  eSet <- getGEO('GSE15460', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  
  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  boxplot(raw_exprSet[,1:50],las="2")
  
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

#GSE62254 没有临床数据
{ 
  rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.1下载其他的geo数据/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
  #eSet <- getGEO(filename = "./GSE84437_series_matrix.txt.gz")
  #不清楚什么情况，手动下载好矩阵你在重新使用getGEO下载就好了
  eSet <- getGEO('GSE62254', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  
  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  boxplot(raw_exprSet[,1:50],las="2")
  
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

#GSE102556 没有临床数据
{ 
  rm(list=ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/1.1下载其他的geo数据/")
  
  library(GEOquery)
  #超时选择官网手动下载（点击 "Series Matrix File" 下载），并手动导入
  #eSet <- getGEO(filename = "./GSE84437_series_matrix.txt.gz")
  #不清楚什么情况，手动下载好矩阵你在重新使用getGEO下载就好了
  eSet <- getGEO('GSE102556', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  
  b = eSet[[1]]
  raw_exprSet= exprs(b) 
  phe=pData(b)
  
  boxplot(raw_exprSet[,1:50],las="2")
  
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