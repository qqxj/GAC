#后台进程号[1] 520136

## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
setwd("/home/datahup/syj/GAC/")

load(file = '/home/datahup/syj/GAC/3.Single_cell_annotation/KNN_scRNA_harmony.Rdata')
table(scRNA_harmony@active.ident)
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.6"#选择稍微多一点的

Epi <- scRNA_harmony[, scRNA_harmony@active.ident %in% c("6", "7", "9", "22")]
table(Epi@active.ident)
#33947 个上皮细胞
rm(scRNA_harmony)

###细胞重命名
{
  Idents(object = Epi)
  table(Idents(object = Epi))
  #旧名字
  Epi[["old.ident"]] <- Idents(object = Epi)
  #新名字
  Idents(object = Epi) <- "Epi_cell"
  #再次查看细胞名字
  Idents(object = Epi)
  levels(Epi)
  table(Idents(Epi))
  #保存重命名细胞
  Epi[["orig.ident"]] <- Idents(object = Epi)
}
table(Epi@active.ident)

#单细胞流程
{
  Epi <- NormalizeData(Epi)
  Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
  # 查看最高变的10个基因
  {
    top10 <- head(VariableFeatures(Epi), 10)
    # 画出不带标签或带标签基因点图
    plot1 <- VariableFeaturePlot(Epi)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2
    #ggsave(filename = "Epi-Top10Gene.png",width = 20,height = 10,path = "../../Fig/Step1-4/")
  }
  # 数据归一化 + 线形降维
  all.genes <- rownames(Epi)
  Epi <- ScaleData(Epi, features = all.genes)
  Epi <- RunPCA(Epi, features = VariableFeatures(object = Epi))
  # 查看PCA结果
  {
    print(Epi[["pca"]], dims = 1:5, nfeatures = 5)
    
    VizDimLoadings(Epi, dims = 1:2, reduction = "pca")
    #ggsave(filename = "Epi-PCA1.png",width = 16,height = 10,path = "../../Fig/Step1-4/")
    
    DimPlot(Epi, reduction = "pca", raster=FALSE)
    #ggsave(filename = "Epi-PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-4/") 
    
    DimHeatmap(Epi, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
    #ggsave(filename = "Epi-PCA3_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")
    
    DimHeatmap(Epi, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
    #ggsave(filename = "Epi-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")
    
    #save(Epi,file = "res0.4_Epi.Rdata")
  }
  ElbowPlot(Epi, ndims = 50)#肘部图 
  # 综合以上方法，选择---25---个主成成分作为参数用于后续分析。
  #ggsave(filename = "res0.4-Epi-Sub-ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-4/")
  #细胞聚类 KNN算法
  library(clustree)
  Epi <- FindNeighbors(Epi, dims = 1:25)#即选取前25个主成分来分类细胞。
  Epi <- FindClusters(object = Epi,
                      resolution = c(seq(0,1,by = 0.1)))
  clustree(Epi@meta.data, prefix = "RNA_snn_res.") 
  #ggsave(filename = "Epi-Sub-resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step1-4/")
  
  #选取resolution = 0.2 作为后续分析参数
  Idents(object = Epi) <- "RNA_snn_res.0.2"
  Epi@meta.data$seurat_clusters = Epi@meta.data$RNA_snn_res.0.2
  head(Idents(Epi), 5)#查看前5个细胞的分类ID
  table(Epi@active.ident)
  
  # UMAP 8*6
  Epi <- RunUMAP(Epi, dims = 1:25)
  DimPlot(Epi, reduction = "umap", label = T,label.box = T,raster=T)
  #ggsave(filename = "Epi-Sub-UMAP-label.png",width = 18,height = 12,path = "../../Fig/Step1-4/")
  
  # TSNE 8*6
  Epi <- RunTSNE(Epi, dims = 1:25,check_duplicates = FALSE)
  DimPlot(Epi, reduction = "tsne", label = T,label.box = T,raster=T)
  #Epi-Sub-TSNE
  
  save(Epi,file = "./7.copykat/Epi.Rdata")
}
load(file = "./7.copykat/Epi.Rdata")






rm(list=ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/GAC/7.copykat/")
library(Seurat)
library(infercnv)

load(file = "/home/datahup/syj/GAC/7.copykat/Epi.Rdata")
table(Epi@active.ident)
#开始
counts <- as.matrix(Epi@assays$RNA@counts)
##运行时间较长 
#时间可能要两天   Time difference of 2.45014 days
library(copykat)
cnv <- copykat(rawmat=counts,
               ngene.chr=5,
               sam.name="ESCC",#忘记改参数了
               n.cores=8)

# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV(拷贝数变异）预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
save(list = ls(), "/home/datahup/syj/GAC/7.copykat/cnv.rdata")

rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/GAC/7.copykat/")
mallignant <- read.delim("./ESCC_copykat_prediction.txt")

mall <- mallignant
rownames(mall) <- mall[,1]
length(unique(mall$cell.names)) # 有重复
mall <- mall[!duplicated(mall$cell.names),] # 去除重复
rownames(mall) <- mall[,1]
names(mall)

table(mall$copykat.pred)
#aneuploid     非整倍体(突变染色体)
#diploid       正常染色体
#not.defined   未识别(随便：想归类到哪归类到哪，也可以归纳为非整倍体)

#归纳数据、修改名字
mall$copykat.pred <- ifelse(mall$copykat.pred == "diploid","Epithelial cell","Cancer cell")
table(mall$copykat.pred)
 # Cancer cell Epithelial cell 
 # 13549           20398 

b <- data.frame('copykat.pred' = mall[,-1])
rownames(b) <- rownames(mall)
table(b$copykat.pred)

#加载上皮细胞-----
load(file = "/home/datahup/syj/GAC/7.copykat/Epi.Rdata")
table(Epi@active.ident)
scRNA <- Epi

# 把细胞的良恶性信息加入metadata
library(Seurat)
scRNA <- AddMetaData(scRNA, metadata = b)
table(scRNA@meta.data[["copykat.pred"]])
save(scRNA,file = "/home/datahup/syj/GAC/7.copykat/copykat_Epi.Rdata")

# UMAP
scRNA <- RunUMAP(scRNA, dims = 1:25)
p1 <- DimPlot(scRNA, group.by = "RNA_snn_res.0.2", 
              label = T,label.box = T,
              raster = T) + NoLegend()
p2 <- DimPlot(scRNA, group.by = "copykat.pred",
              label=T,label.box = T,
              #cols = pal,
              raster = T) + NoLegend()
pc <- p1 + p2
pc

# TSNE
scRNA <- RunTSNE(scRNA, dims = 1:25,check_duplicates = FALSE)
p3 <- DimPlot(scRNA,group.by = "RNA_snn_res.0.2",
              reduction = "tsne",
              label = T,label.box = T,raster=T)+NoLegend()
p4 <- DimPlot(scRNA,group.by = "copykat.pred",
              reduction = "tsne",
              label = T,label.box = T,raster=T)+NoLegend()
pc1 <- p3 + p4
pc1






























