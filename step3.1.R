rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/GAC/")
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

load(file = "./2.download_singlecell/scRNA_harmony.Rdata")
table(scRNA_harmony@meta.data[["orig.ident"]])
# GSE183904  GSE268238 OMIX001073 
# 119862     248520     129177 

###细胞重命名
{
  Idents(object = scRNA_harmony)
  table(Idents(object = scRNA_harmony))
  #旧名字
  scRNA_harmony[["old.ident"]] <- Idents(object = scRNA_harmony)
  #新名字
  Idents(object = scRNA_harmony) <- "GAC_cell"
  #再次查看细胞名字
  Idents(object = scRNA_harmony)
  levels(scRNA_harmony)
  table(Idents(scRNA_harmony))
  #保存重命名细胞
  scRNA_harmony[["orig.ident"]] <- Idents(object = scRNA_harmony)
  
}

#细胞质控QC------
{
  # 计算线粒体基因比例 线粒体是独立遗传的，不是染色体上基因控制的 默认<10%
  scRNA_harmony[["percent.mt"]] <- PercentageFeatureSet(scRNA_harmony, assay = "RNA", pattern = "^MT-")
  head(scRNA_harmony@meta.data)
  # nFeature_RNA代表每个细胞测到的基因数目。
  # nCount_RNA代表每个细胞测到所有基因的表达量之和。
  # percent.mt代表测到的线粒体基因的比例。
  VlnPlot(scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3, pt.size = 0)  & 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # 设置x轴标签不倾斜
                                       axis.title.x = element_blank())  # 隐藏x轴标题
  #保存为10*7
  plot1 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)+
    theme( # 图例
          legend.position = c(0.95, 0.95),  # 右上方
          legend.justification = c(1, 1),
          legend.text = element_text(size = 12, face = "bold"),   
          legend.title = element_text(size = 15, face = "bold"), 
          plot.title = element_text(size = 15, face = "bold"), #标题
          axis.title.x = element_text(size = 13, face = "bold"),  # X轴标题
          axis.title.y = element_text(size = 13, face = "bold"))+  # Y轴标题
            labs(title = "Pre-QC mitochondrial gene percentage",
            x = "nCount_RNA",
            y = "percent.mt")     
  plot2 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)+
    theme( # 图例
      legend.position = c(0.95, 0.95),  # 右上方
      legend.justification = c(1, 1),
      legend.text = element_text(size = 12, face = "bold"),   
      legend.title = element_text(size = 15, face = "bold"), #标题
      plot.title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold"),  # X轴标题
      axis.title.y = element_text(size = 13, face = "bold"))+  # Y轴标题
    labs(title = "Pre-QC: Distribution of Gene Counts vs UMI Totals",
         x = "nCount_RNA",
         y = "nFeature_RNA")    

  plot1 + plot2
  #保存为10*6散点图文件太大，建议直接绘制好用png格式保存

  #根据上图进行质控
  table(scRNA_harmony[["percent.mt"]] < 10)
  table(scRNA_harmony@meta.data[["nCount_RNA"]] < 10000)
  scRNA_harmony <- subset(scRNA_harmony, subset = nFeature_RNA > 500 & nFeature_RNA < 6000
                            & percent.mt < 10 &
                            nCount_RNA > 500 & nCount_RNA < 10000)
  
  
  table(scRNA_harmony@active.ident)
  # GAC_cell 
  # 273924 
  #计算红血细胞基因比例 红细胞没有细胞核，没有转录组 默认<3%
  rownames(scRNA_harmony)[grep("^HB[^(p)]", rownames(scRNA_harmony))]
  scRNA_harmony <- PercentageFeatureSet(scRNA_harmony, "^HB[^(p)]", col.name = "percent_hb")
  table(scRNA_harmony[["percent_hb"]] < 1)
  scRNA_harmony = subset(scRNA_harmony, subset = percent_hb < 1)
  
  VlnPlot(scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3 ,pt.size = 0) &
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # 设置x轴标签不倾斜
        axis.title.x = element_blank())  # 隐藏x轴标题
  #保存为10*7
  
  plot1 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)+
    theme( # 图例
      legend.position = c(0.95, 0.95),  # 右上方
      legend.justification = c(1, 1),
      legend.text = element_text(size = 12, face = "bold"),   
      legend.title = element_text(size = 15, face = "bold"), 
      plot.title = element_text(size = 15, face = "bold"), #标题
      axis.title.x = element_text(size = 13, face = "bold"),  # X轴标题
      axis.title.y = element_text(size = 13, face = "bold"))+  # Y轴标题
    labs(title = "After-QC mitochondrial gene percentage",
         x = "nCount_RNA",
         y = "percent.mt")     
  plot2 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)+
    theme( # 图例
      legend.position = c(0.95, 0.95),  # 右上方
      legend.justification = c(1, 1),
      legend.text = element_text(size = 12, face = "bold"),   
      legend.title = element_text(size = 15, face = "bold"), #标题
      plot.title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold"),  # X轴标题
      axis.title.y = element_text(size = 13, face = "bold"))+  # Y轴标题
    labs(title = "After-QC: Distribution of Gene Counts vs UMI Totals",
         x = "nCount_RNA",
         y = "nFeature_RNA")    
  
  plot1 + plot2
  #保存为10*6散点图文件太大，建议直接绘制好用png格式保存
  
  save(scRNA_harmony,file = "./3.Single_cell_annotation/QC_scRNA_harmony.Rdata")
  
}
load(file = "./3.Single_cell_annotation/QC_scRNA_harmony.Rdata")

#标准工作流程------
{
  scRNA_harmony <- NormalizeData(scRNA_harmony)
  scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scRNA_harmony)
  scRNA_harmony <- ScaleData(scRNA_harmony, features = all.genes)
  scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(object = scRNA_harmony))
  
  {
    # 查看最高变的10个基因
    top10 <- head(VariableFeatures(scRNA_harmony), 10)
    # 画出不带标签或带标签基因点图
    plot1 <- VariableFeaturePlot(scRNA_harmony)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2
    #保存pdf12*6
  }
  
  ElbowPlot(scRNA_harmony, ndims = 50)
  #pdf8*5
  # ggsave(filename = "ElbowPlot.pdf",width = 12,height = 10,
  #        path = "./Fig/step 2/")
  
  #降维聚类 knn
  ###根据肘部图参数选择25###
  scRNA_harmony = FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:12) ###选择25###
  scRNA_harmony = FindClusters(object = scRNA_harmony, 
                               resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
  library(clustree)
  clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.") 
  #pdf20*14
  # ggsave(filename = "resolution(0.3).pdf",width = 20,height = 14,
  #        path = "/home/syj/NSCC/Fig/step 2/")
  
  save(scRNA_harmony,file = '/home/datahup/syj/GAC/3.Single_cell_annotation/KNN_scRNA_harmony.Rdata')
}

load(file = '/home/datahup/syj/GAC/3.Single_cell_annotation/KNN_scRNA_harmony.Rdata')

#可视化----
{
  Idents(object = scRNA_harmony) <- "RNA_snn_res.0.6"#选择稍微多一点的
  head(Idents(scRNA_harmony), 5)#查看前5个细胞的分类ID
  table(Idents(scRNA_harmony))
  
  # UMAP
  scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:12)
  plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T, raster=T, label.box = T) 
  plot1
  # ggsave(filename = "scRNA_harmony_umap.pdf", plot = plot1, width = 8,height = 6,
  #        path = "/home/syj/NSCC/Fig/step 2/")
  
  # TSNE
  scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:12,check_duplicates = FALSE)
  plot2 = DimPlot(scRNA_harmony, reduction = "tsne", label = T,label.box = T,raster=T)
  plot2
  # ggsave(filename = "scRNA_harmony_tsne.png",
  #        width = 8,height = 6,path = "../../Fig/Step1-1/")
  
}
##找每个簇的差异基因------
{
  library(dplyr)
  Idents(object = scRNA_harmony) <- "RNA_snn_res.0.6"
  table(Idents(scRNA_harmony))
  markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
  top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  #DoHeatmap(scRNA_harmony, features = top10$gene) + NoLegend() #参考，保存12*6
  #ggsave(filename = "Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-1/")
  table(top10$cluster)
  write.csv(markers,file = './3.Single_cell_annotation/First_markers.csv')
  write.csv(top10,file = './3.Single_cell_annotation/First_markers_top10.csv')
  
}

## SingleR(跳过自己注释)--------
{
  sce = scRNA_harmony
  library(Seurat)
  library(celldex)
  library(SingleR)
  sce_for_SingleR <- GetAssayData(sce, slot="data")
  sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.6
  clusters=sce@meta.data$seurat_clusters
  
  #加载注释包
  #load("/home/syj/NSCC/Data/Annotation packages/singleRref.Rdata")
  
  Blue.ref <- celldex::BlueprintEncodeData()
  pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                           method = "cluster", clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
  pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                           method = "cluster", clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
  pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                           method = "cluster", clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  Mona.ref <- celldex::MonacoImmuneData()
  pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                           method = "cluster", clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  Nover.ref <- celldex::NovershternHematopoieticData()
  pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                            method = "cluster", clusters = clusters, 
                            assay.type.test = "logcounts", assay.type.ref = "logcounts")
  cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                      Blue=pred.Blue.ref$labels,
                      DICE=pred.DICE.ref$labels,
                      HPCA=pred.HPCA.ref$labels,
                      Mona=pred.Mona.ref$labels,
                      Nover=pred.Nover.ref$labels )
  head(cellType)
  sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
  sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
  sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
  sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
  sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']
  
  scRNA_harmony = sce
  
  #保存scRNA_harmony和机器注释结果
  #save(scRNA_harmony, file = "./save/step 2/Step2_output_1.Rdata")
  #write.csv(cellType,file = './save/step 2/Step2_celltype_1.csv')
}


#手工注释(cellmarker2.0、chatgpt、文献)以及可视化-------------
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
library(ggplot2)
setwd("/home/datahup/syj/GAC/3.Single_cell_annotation/")

load(file = '/home/datahup/syj/GAC/3.Single_cell_annotation/KNN_scRNA_harmony.Rdata')
table(Idents(scRNA_harmony)) #25dim
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.6"


#先更换注释信息------
{
  new.cluster.ids <- c(
    "0" = "T cell",
    "1" = "T cell",
    "2" = "B cell",
    "3" = "T cell",
    "4" = "Plasma cell",
    "5" = "Monocyte/Macrophage",
    "6" = "Gastric epithelial cell",
    "7" = "Gastric epithelial cell",
    "8" = "Endothelial cell",
    "9" = "Gastric epithelial cell",  
    "10" = "T cell",
    "11" = "T cell",
    "12" = "Dendritic cell",
    "13" = "Fibroblast",
    "14" = "T cell",
    "15" = "Mast cell",
    "16" = "Dendritic cell",
    "17" = "T cell",
    "18" = "Fibroblast",
    "19" = "T cell",
    "20" = "Smooth muscle cell",
    "21" = "T cell",
    "22" = "Gastric epithelial cell",
    "23" = "B cell",
    "24" = "Endothelial cell",
    "25" = "Smooth muscle cell")
  
  Idents(scRNA_harmony) <- "RNA_snn_res.0.6"
  scRNA_harmony$celltype <- plyr::mapvalues(
    x = Idents(scRNA_harmony),
    from = names(new.cluster.ids),
    to = new.cluster.ids)
  table(Idents(scRNA_harmony))
  Idents(scRNA_harmony) <- scRNA_harmony@meta.data[["celltype"]]
    
  # UMAP
   #scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:12)
  plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T, raster=T, label.box = T) 
  plot1
  #"UMAP_Cellmarker.pdf" 8*6
  
  # TSNE
   #scRNA_harmony <- RunTSNE(scRNA_harmony, dims = 1:12,check_duplicates = FALSE)
  plot2 = DimPlot(scRNA_harmony, reduction = "tsne", label = T,label.box = T,raster=T)
  plot2
  #"TSNE_Cellmarker.pdf" 8*6
  save(scRNA_harmony,file= "./KNN_scRNA_harmony.Rdata")
}

#检查每个簇中的基因表达情况(可跳过)-----------
{
  top10 <- read.csv('./First_markers_top10.csv',row.names = 1)
  Gene <- subset(top10, top10$cluster==26) # 依次从0开始检查每个簇的top10marker基因
  p_all_markers=DotPlot(scRNA_harmony,
                        features = Gene$gene,
                        cols = c("lightgrey", "purple"),
                        scale = T,assay='RNA' )+
    theme(axis.text.x=element_text(angle=45,hjust = 1))
  p_all_markers
  Gene$gene
  
  #热图(先看一下结果）
  library(pheatmap)
  library(ggplot2)
  load(file = '/home/syj/NSCC/save/step 2/Step2_output.Rdata')
  table(Idents(scRNA_harmony))
  top10 <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_celltype_0.4.csv')#,row.names = 1
  # #1
  # DoHeatmap(Lymph, features = top10$gene) + NoLegend()
  # ggsave(filename = "Top10-MarkerGene.pdf", plot = plot, width = 12,height = 8,
  #        path = "/home/syj/NSCC/Fig/step 2/")
  #2
  library(ComplexHeatmap)
  library(ggunchull)
  library(scRNAtoolVis)
  library(ggplot2)
  #show_col(pal_d3('category20')(20))
  library(RColorBrewer)
  pal <- brewer.pal(20, "Set3")
  #pal = pal_d3('category20')(20)
  library(viridis)
  pal <- viridis(30, option = "D", alpha = 0.9)
  pdf(file = './Fig/step 2/scRNA_harmony_pheatmap_0.4.pdf',width = 25,height = 25)
  AverageHeatmap(object = scRNA_harmony,
                 markerGene = top10$gene,
                 column_names_rot = 0, 
                 htCol = c("blue", "yellow", "red")
  )
  dev.off()
  
  #调整top_10行的排列顺序，重新画上方热图
  table(top10$cluster)
  custom_order  <- factor(top10$cluster, levels = c("0", "2", "1", "3","4","9","5","6","7","8"))  # 按照需要定义排序顺序
  top10_1 <- top10[order(custom_order), ]
  top10_1 <- top10_1[!(top10_1$gene %in% c("LTB", "IGFBP7","CALD1","HLA−DRA","HLA−DQB1","HLA−DQA1","SPARCL1")), ]
  
  top10 <- top10_1
}

#手工注释------
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
library(ggplot2)
setwd("/home/datahup/syj/GAC/3.Single_cell_annotation/")

load(file= "./KNN_scRNA_harmony.Rdata")

#热图(先看一下结果）
library(pheatmap) 
library(ggplot2)
#需要数据单细胞数据、gene名
table(Idents(scRNA_harmony))
gene <- c("CD3E","TRBC2","IL7R","CD8A",#T
          "CD74","HLA-DRB1","TYROBP","HLA-DMA",#DC
          "COL3A1","COL1A2","DCN","LUM",#成纤维
          "TPSAB1","CPA3","TPSD1",#肥大
          "MS4A1","CD79A","BANK1",#B
          "TAGLN","ACTA2","TPM2",#"MYL9","CALD1","RGS5", #平滑肌
          "KRT8","AGR2","MUC5AC","ANXA10","PGA3","LIPF",#"GKN1",  #腺上皮
          "PLVAP","RAMP2",#内皮
          "JCHAIN","MZB1","IGHA1",#"DERL3","XBP1" #浆细胞
          "S100A8","IL1RN","FCN1")#单核巨噬

library(ComplexHeatmap)
library(ggunchull)
library(scRNAtoolVis)
library(ggplot2)
#show_col(pal_d3('category20')(20))
library(RColorBrewer)
pal <- brewer.pal(20, "Set3")
#pal = pal_d3('category20')(20)
library(viridis)
pal <- viridis(30, option = "D", alpha = 0.9)
pdf(file = './scRNA_harmony_pheatmap_0.4.pdf',width = 25,height = 25)
AverageHeatmap(object = scRNA_harmony,
               markerGene = gene,
               column_names_rot = 0, 
               htCol = c("blue", "yellow", "red")
)
dev.off()

# marker的气泡图
library(Seurat)
library(ggplot2)
#需要数据：单细胞数据集、gene名
 #Idents(object = scRNA_harmony) <- scRNA_harmony@meta.data[["FirstAnnotation"]]
table(Idents(scRNA_harmony))

gene
pdf(file = './scRNA_harmony_DotPlot_0.4.pdf',width = 15,height = 15)
DotPlot(scRNA_harmony,
        features = gene,
        cols = c("lightgrey", "purple"),
        scale = T,
        assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
dev.off()

p <- DotPlot(scRNA_harmony, features = gene,
             assay='RNA' ) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p

#美化图片----
{
  #load(file = "/home/datahup/syj/NSCC/step1_1NoteCell/Step2_output.Rdata")
  #cell_cluster <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_celltype_0.4.csv')#,row.names = 1
  #top10 <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_markers_0.4.csv',row.names = 1)

  gene_cell_pairs <- data.frame(
    gene = c("CD3E","TRBC2","IL7R","CD8A",#T
              "CD74","HLA-DRB1","TYROBP","HLA-DMA",#DC
              "COL3A1","COL1A2","DCN","LUM",#成纤维
              "TPSAB1","CPA3","TPSD1",#肥大
              "MS4A1","CD79A","BANK1",#B
              "TAGLN","ACTA2","TPM2",#"MYL9","CALD1","RGS5", #平滑肌
              "KRT8","AGR2","MUC5AC","ANXA10","PGA3","LIPF",#"GKN1",  #腺上皮
              "PLVAP","RAMP2",#内皮
              "JCHAIN","MZB1","IGHA1",#"DERL3","XBP1" #浆细胞
              "S100A8","IL1RN","FCN1"),#单核巨噬
    cell_type = c(rep("T cell", 4),
                  rep("Dendritic cell", 4),
                  rep("Fibroblast", 4),
                  rep("Mast cell", 3),
                  rep("B cell", 3),
                  rep("Smooth muscle cell",3),
                  rep("Gastric epithelial cell", 6),
                  rep("Endothelial cell", 2),
                  rep("Plasma cell", 3),
                  rep("Monocyte/Macrophage", 3)))
  #定义顺序
  gene_cell_pairs$cell_type <- factor(gene_cell_pairs$cell_type,
                                      levels = c("Monocyte/Macrophage","Plasma cell","Endothelial cell",
                                                 "Gastric epithelial cell","Smooth muscle cell","B cell",
                                                 "Mast cell","Fibroblast","Dendritic cell","T cell"))
  #可视化
  pdf(file = './scRNA_DotPlot_美化.pdf', width = 16, height = 16, paper = "special", onefile = TRUE)
  DotPlot(scRNA_harmony,
                features = split(gene_cell_pairs$gene, gene_cell_pairs$cell_type),) +
    RotatedAxis() + 
    theme(
      panel.border = element_rect(color = "black"),
      panel.spacing = unit(1, "mm"),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
    )+
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
  dev.off()
}



























