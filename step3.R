## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/GAC/2.download_singlecell/")
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

#先试着读一个看看能不能读进来------
# # 读入 10x Genomics 格式数据
# seurat_obj <- Read10X(data.dir = "./GSE183904_RAW.tar")
# 
# # 创建 Seurat 对象
# seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = "GSE183904")
# 
# a = read.table("./GSM5573466_sample1.csv.gz" ,row.names = T,col.names = T)#, quote = ''

a <- read.table(gzfile("./2.download_singlecell/GSE183904/GSM5573467_sample2.csv.gz"), 
                header = TRUE, row.names = 1, sep = ",")
# #这一句也可以用
# a <- read.table("./GSM5573466_sample1.csv.gz",
#                 header = TRUE, row.names = 1, sep = ",")
NSCLC7 <- CreateSeuratObject(counts = a, project = "STAD",
                             min.cells = 10, min.features = 500)

#GSE183904:26个癌症样本119862个细胞-----
{
  setwd("/home/datahup/syj/GAC/2.download_singlecell/GSE183904")
  
  file_paths <- list.files(pattern = "\\.csv\\.gz$", full.names = TRUE)
  file_paths  
  
  # 创建一个空列表来存储Seurat对象
  seurat_list <- list()
  # 循环读取每个文件并创建Seurat对象
  for (file in file_paths) {
    # 生成样本名（去除路径和扩展名）
    sample_name <- tools::file_path_sans_ext(basename(file)) # 去掉.csv.gz
    sample_name <- sub("\\.csv$", "", sample_name) # 再去掉.csv
    
    # 读取数据
    counts <- read.table(gzfile(file), header = TRUE, row.names = 1, sep = ",")
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(counts = counts,
                                     project = "GSE183904",
                                     min.cells = 10,
                                     min.features = 500)
    
    # 将对象存入列表
    seurat_list[[sample_name]] <- seurat_obj
  }
  
  # 查看结果
  names(seurat_list)
  
  # 取出第一个样本作为起始对象
  GSE183904 <- seurat_list[[1]]
  
  # 依次合并后续样本
  for (i in 2:length(seurat_list)) {
    GSE183904 <- merge(GSE183904, y = seurat_list[[i]])
  }
  
  # # 修改 orig.ident 为统一值 "GSE183904"
  # GSE183904$orig.ident <- "GSE183904"
  
  # 查看meta信息
  table(GSE183904$orig.ident)
  
  save(GSE183904,file = "./GSE183904_RAW.Rdata")
  
}
load(file = "./GSE183904_RAW.Rdata")

#先试试能不能读进去(太少了，先尝试在国家生信中心寻找样本ing)---------
#GSE163558 3个
#GSE184198  1个
#GSE134520  1个
#GSE112302 2（每个测了三组一共6个）
#GSE167297 5例

#GSE163558(太少了不用了) 3例:13705个------
{
  setwd("/home/datahup/syj/GAC/2.download_singlecell/")
  # 获取所有子文件夹（10x格式数据）
  folders <- list.dirs(path = "./GSE163558/", recursive = FALSE, full.names = TRUE)
  # 创建空列表保存每个Seurat对象
  seurat_list <- list()
  # 批量读取并创建 Seurat 对象
  for (folder in folders) {
    sample_name <- basename(folder)  # 提取文件夹名作为样本名
    
    # 读取10x格式数据
    seurat_counts <- Read10X(data.dir = folder)
    
    # 创建Seurat对象
    seu <- CreateSeuratObject(counts = seurat_counts,
                              project = "GSE163558",
                              min.cells = 10,
                              min.features = 500)
    
    # 添加样本标识
    seu$orig.ident <- sample_name
    
    # 保存到列表
    seurat_list[[sample_name]] <- seu
  }
  # 查看
  names(seurat_list)
  
  # # 合并所有 Seurat 对象
  # GSE163558 <- merge(x = seurat_list[[1]],
  #                     y = seurat_list[-1],
  #                     add.cell.ids = names(seurat_list),
  #                     project = "GSE163558")
  
  # 取出第一个样本作为起始对象
  GSE163558 <- seurat_list[[1]]
  
  # 依次合并后续样本
  for (i in 2:length(seurat_list)) {
    GSE163558 <- merge(GSE163558, y = seurat_list[[i]])
  }
  # 查看
  table(GSE163558$orig.ident)
  GSE163558$orig.ident <- "GSE163558"
  
  #保存
  save(GSE163558, file = "./GSE163558_RAW.Rdata")
  
}
load(file = "./GSE163558_RAW.Rdata")

#GSE268238 38例 248520个细胞-----
{
  setwd("/home/datahup/syj/GAC/2.download_singlecell/")
  scRNA <- readRDS("/home/datahup/syj/GAC/2.download_singlecell/GSE268238/GSE268238_GC.PNI.rds")
  table(scRNA@meta.data[["orig.ident"]])
  table(scRNA@meta.data[["sample"]]) #38例患者，16例癌旁样本
  #scRNA@meta.data[["orig.ident"]] <- "GSE268238"
  
  # 1. 提取原始计数矩阵
  counts<- scRNA@assays[["RNA"]]@counts
  # 2. 使用原始矩阵创建新的 Seurat 对象（可以设定新的过滤条件）
  GSE268238 <- CreateSeuratObject(counts = counts,
                                  project = "GSE268238",
                                  min.cells = 10,
                                  min.features = 500)
  
  save(GSE268238,file = "./GSE268238_RAW.Rdata")
}
load(file = "./GSE268238_RAW.Rdata")

#NGDC(OMIX)国家生物信息中心------
#数据库介绍，以及使用ftp下载的方法： https://www.jianshu.com/p/555e42062233
#                                   https://ngdc.cncb.ac.cn/omix/download-guide/Submission-Guidance-FAQ-OMIX-CN.pdf
# https://ngdc.cncb.ac.cn/omix/release/OMIX001073
#OMIX001073 包括 23 个原发性胃癌样本。
#PMID：40124374 
#备注：所有细胞得矩阵无法下载，里面格式错误，我们依次下载好注释好的细胞矩阵

#读入单个数据成功----
{
  seurat_obj <- Read10X(data.dir = "./OMIX/OMIX001073-20-01/")
  
  OMIX001073 <- CreateSeuratObject(counts = seurat_obj,
                                   project = "OMIX001073",
                                   min.cells = 10,
                                   min.features = 500)
  table(OMIX001073$orig.ident)
  # OMIX001073 
  # 19864 
  
  #save(OMIX001073,file = "./OMIX001073_RAW.Rdata")
  load(file = "./OMIX001073_RAW.Rdata")
}
#编写循环依次读入 OMIX:129177个细胞----
{
  setwd("/home/datahup/syj/GAC/2.download_singlecell/")
  
  # 获取所有子文件夹（10x格式数据）
  folders <- list.dirs(path = "./OMIX/", recursive = FALSE, full.names = TRUE)
  
  # 创建空列表保存每个Seurat对象
  seurat_list <- list()
  
  # 批量读取并创建 Seurat 对象
  for (folder in folders) {
    sample_name <- basename(folder)  # 提取文件夹名作为样本名
    
    # 读取10x格式数据
    seurat_counts <- Read10X(data.dir = folder)
    
    # 创建Seurat对象
    seu <- CreateSeuratObject(counts = seurat_counts,
                              project = "OMIX001073",
                              min.cells = 10,
                              min.features = 500)
    
    # 添加样本标识
    seu$orig.ident <- sample_name
    
    # 保存到列表
    seurat_list[[sample_name]] <- seu
  }
  # 查看
  names(seurat_list)
  
  # 合并所有 Seurat 对象
  OMIX001073 <- merge(x = seurat_list[[1]],
                          y = seurat_list[-1],
                          add.cell.ids = names(seurat_list),
                          project = "OMIX001073")
  
  # 取出第一个样本作为起始对象
  OMIX001073 <- seurat_list[[1]]
  
  # 依次合并后续样本
  for (i in 2:length(seurat_list)) {
    OMIX001073 <- merge(OMIX001073, y = seurat_list[[i]])
  }
  # 查看
  table(OMIX001073$orig.ident)
  OMIX001073$orig.ident <- "OMIX001073"
  
  #保存
  save(OMIX001073, file = "./OMIX001073_RAW.Rdata")
}
load(file = "./OMIX001073_RAW.Rdata")


#harmony质控，合并数据----
 #OMIX001073、GSE183904、GSE268238
 #OMIX001073 26例 129177个细胞
 #GSE268238 38例 248520个细胞
 #GSE183904:26例癌症样本119862个细胞
{
  rm(list = ls())
  options(stringsAsFactors = F)
  gc()
  setwd("/home/datahup/syj/GAC/2.download_singlecell/")
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(stringr)
  
  load(file = "./GSE183904_RAW.Rdata")
  load(file = "./GSE268238_RAW.Rdata")
  load(file = "./OMIX001073_RAW.Rdata")
  
  # Step 1: 确保每个对象的 orig.ident 统一设置
  table(GSE183904$orig.ident) 
  table(GSE268238$orig.ident)
  table(OMIX001073$orig.ident)
  
  # Step 2: 合并三个对象
  combined <- merge(GSE183904,
                    y = list(GSE268238, OMIX001073),
                    add.cell.ids = c("GSE183904", "GSE268238", "OMIX001073"))#,project = "GAC_Combined"
  table(combined@meta.data[["orig.ident"]])
  
  # Step 3: 标准单细胞预处理流程
  combined = NormalizeData(combined) %>%
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose=FALSE)
  # combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  # combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  # combined <- ScaleData(combined, verbose = FALSE)
  # combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  
  # Step 4: PCA 可视化（去除批次效应前）
  options(repr.plot.height = 4,repr.plot.width = 7)
  p1 <- DimPlot(combined, reduction = "pca",group.by = "orig.ident")
  p2 <- VlnPlot(combined, features = "PC_1", pt.size = 0, group.by = "orig.ident")
  library(cowplot)
  library(ggplot2)
  combined_plot <- plot_grid(p1, p2, ncol = 2)
  ggsave(filename = "./PCA_去除批次效应前.pdf", plot = combined_plot,
         width = 12, height = 6, units = "in")
  
  # Harmony整合
  library(harmony)
  scRNA_harmony <- RunHarmony(combined, group.by.vars = "orig.ident")
  table(scRNA_harmony@meta.data[["orig.ident"]])
  
  #去除批次效应后
  p3 <- DimPlot(scRNA_harmony, reduction = "harmony",group.by = "orig.ident")
  p4 <- VlnPlot(scRNA_harmony, features = "harmony_1", pt.size = 0, group.by = "orig.ident")
  library(cowplot)
  library(ggplot2)
  combined_plot <- plot_grid(p3, p4, ncol = 2)
  ggsave(filename = "./PCA_去除批次效应后.pdf", plot = combined_plot,
         width = 12, height = 6, units = "in")
  
  save(scRNA_harmony,file = "./scRNA_harmony.Rdata")
}
load(file = "./scRNA_harmony.Rdata")






