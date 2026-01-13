## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
library(ggplot2)
library(monocle)
library(patchwork)
library(dplyr)
setwd("/home/datahup/syj/GAC/")

load(file = "./7.copykat/copykat_Epi.Rdata")
table(scRNA@active.ident)
Idents(object = scRNA) <- scRNA@meta.data[["copykat.pred"]]
table(scRNA@active.ident)
# Epithelial cell     Cancer cell 
# 20398               13549 

# Cancer_cell  = scRNA[,scRNA@active.ident %in% "Cancer cell"]
# table(Cancer_cell@active.ident)
# save(Cancer_cell,file = "./8.拟时序/上皮癌细胞.Rdata")
#单细胞流程(跳过)--------
 #scRNA <- Cancer_cell
rm(Cancer_cell)
{
  library(Seurat)
  scRNA <- NormalizeData(scRNA)
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scRNA)
  scRNA <- ScaleData(scRNA, features = all.genes)
  scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"))
  scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
  ElbowPlot(scRNA)
  #细胞聚类 KNN算法
  library(clustree)
  scRNA <- FindNeighbors(scRNA, dims = 1:11)#查看肘部图,即选取前**个来分类细胞。
  scRNA <- FindClusters(object = scRNA,
                        resolution = c(seq(0,1,by = 0.1)))
  clustree(scRNA@meta.data, prefix = "RNA_snn_res.") 
  #选择Cluster = 0.1
  Idents(object = scRNA) <- "RNA_snn_res.0.1"
  scRNA@meta.data$seurat_clusters = scRNA@meta.data$RNA_snn_res.0.1
  
  Cancer_cell <- scRNA
  save(Cancer_cell,file = "./8.拟时序/上皮癌细胞.Rdata")
}
#load(file = "./8.拟时序/上皮癌细胞.Rdata")

# monocle  -----------
Neutrophil <- scRNA
rm(scRNA)
{ # Step1 Seurat -> monocle
  {
    #转化为稀疏矩阵
    data <- as(as.matrix(Neutrophil@assays$RNA@counts), 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = Neutrophil@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
  }

  # Step 2 估计size factor和离散度
  {
    monocle_cds <- estimateSizeFactors(monocle_cds)#估计每个细胞中mRNA的大小 
    monocle_cds <- estimateDispersions(monocle_cds,cores = 8)#估计数据中基因表达的离散度
    # 过滤低质量的细胞
    monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
    print(head(fData(monocle_cds)))#显示有多少个基因
  }
  save(monocle_cds,file = "./8.拟时序/monocle_cds.Rdata")
  
  HSMM=monocle_cds

  # Step 4 构建轨迹
  # Step 4.1: 选择定义过程的基因 三种方法均为无监督
  {
    # #使用clusters差异表达基因
    deg.cluster <- FindAllMarkers(Neutrophil)
    deg.cluster <- FindAllMarkers(
      Neutrophil,
      min.pct = 0.1,              # 不要求表达比例
      logfc.threshold = 0.1,      # 不设 fold change 阈值
      only.pos = FALSE          # 包括上调和下调基因
    )
    save(deg.cluster,file = "./8.拟时序/deg.cluster.Rdata")
    diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
    HSMM <- setOrderingFilter(HSMM, diff.genes)
    plot_ordering_genes(HSMM)
    ##使用seurat选择的高变基因
    # var.seurat <- VariableFeatures(Neutrophil)
    # var.seurat <- Neutrophil@assays[["RNA"]]@var.features
    # HSMM <- setOrderingFilter(HSMM, var.seurat)
    # plot_ordering_genes(HSMM)
    # #使用monocle选择的高变基因
    # disp_table <- dispersionTable(HSMM)
    # disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    # unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1*dispersion_fit)
    # %>% pull(gene_id) %>% as.character()
    
    # HSMM <- setOrderingFilter(HSMM, disp.genes)
    # plot_ordering_genes(HSMM)
  }
  
  # 可视化 tSNE
  {
    #降维
    HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                            reduction_method = 'tSNE', verbose = T,
                            check_duplicates = FALSE)
    #聚类
    HSMM <- clusterCells(HSMM, num_clusters = 2)
    plot_cell_clusters(HSMM)
    # plot_cell_clusters(HSMM, 1, 2, color = "CellType",
    #                    markers = c("MYF5", "ANPEP"))
  }
 

  
  # Step 4.2 降维
  {
    HSMM <- reduceDimension(HSMM, max_components = 3, ncenter = 50,
                            method = 'DDRTree')
  }

  # Step 4.3 按照轨迹排序细胞
    HSMM <- orderCells(HSMM)
    {
    # 每次重新加载monocle包都会报以下错误：
    # Error in if (class(projection) != "matrix") projection <- as.matrix(projection) : 
    #   the condition has length > 1
    ##### 解决办法：
    #运行
    #trace('project2MST', edit = T, where = asNamespace("monocle")) 
    # 在弹出窗口中把if (elass(projection) != "matrix") 注释，并保存，重新运行即可
    
    # (逆向调节）使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
    # cds <- orderCells(cds, root_state = 3) #把State3设成拟时间轴的起始点
    ## !!!若时序得分和预期得分颠倒，可将Pseudotime的值进行反向转换，将Pseudotime列的值更新为最大值减去原始Pseudotime值。即原始值越大，转换后的值越小，原始值越小，转换后的值越大。
    # pData(cds)$Pseudotime <- max(pData(cds)$Pseudotime) - pData(cds)$Pseudotime
    }
    ls()
  save(data,deg.cluster,diff.genes,fd,fData,HSMM,monocle_cds,Neutrophil,pd,file = "./8.拟时序/HSMM.Rdata")
}

#可视化结果-----
{
  #先展示差异基因的排序图
  plot_cell_trajectory(HSMM, color_by = "copykat.pred") +
    geom_point(aes(color = copykat.pred), size = 0.5) +
    scale_color_manual(values = c(
      "Epithelial cell" = "#1f78b4",     # 蓝色
      "Cancer cell" = "#e31a1c"    # 红色
    )) +
    theme_classic()
  
  plot_cell_trajectory(HSMM, color_by = "copykat.pred")+
    facet_wrap(~seurat_clusters, nrow = 1)
  plot_cell_trajectory(HSMM, color_by = "State")
 
  plot_cell_trajectory(HSMM, color_by = "Pseudotime")
  plot_cell_trajectory(HSMM, color_by = "Pseudotime",
                       show_branch_points = F,  # show_branch_points 是否展示分枝节点
                       show_tree = T)  +        # show_tree 是否展示连线
    scale_color_steps(high = 'green', low = 'purple')
  #拟时序轨迹错了，纠正
  test=orderCells(HSMM,reverse = T) 
  save(test,file = "./8.拟时序/test.Rdata")
  plot_cell_trajectory(test, color_by = "Pseudotime")
  plot_cell_trajectory(test, color_by = "State")
  #ggsave(filename = 'State.pdf',width = 10,height = 5,path = './photo/step 3-1/')
  
  #伪时间热图
  #deg.cluster差异基因热图
  library(dplyr)
  Marker <- unique(deg.cluster$gene)
  #差异基因分析
  Time_diff <- differentialGeneTest(HSMM[Marker,], cores = 8, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  #pull 函数提取列
  Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
  # num_clusters为人为设置的聚类
  
  #自定义筛掉一些基因
  
  library(ComplexHeatmap)
  library(circlize)  # 用于 colorRamp2
  
  Time_gene <- Time_genes
  a <- c("PLA2G2A", "TFF3", "GAPDH", "IL32", "COL9A3", "PRAP1", "PSMA7", "CDKN2A", "PKDCC", "AC103702.2", "SPINK1", "MT1G", "MT2A", "PGC", "TFF2", "TFF1", "LIPF", "PGA3",#1:100
         "KLF6", "RPS21", "S100A11", "S100A10", "S100A4", "CLDN4", "TGM2", "KRT7", "CAMK2N1", "LAMA5", "C3", "KRT17", "MSLN", "HOXB6", "PMEPA1", "COL16A1", "UBE2C",#101:200
         "ACTB", "FTL", "FABP1", "GSTP1", "HAGHL", "TXN", "PFDN4", "RAMP1", "GRINA",#201:300
         "GSTA1", "RPS19", "S100A6", "LYZ", "SELENOM", "TRAF2", "TOP2A",#301：400
         "MRPL15", "NGRN", "CCT3", "SPTSSA", "LTBP3", "NTMT1", "PRPF3", "CNPY3", "GID8", "MKKS", "COL1A1", "MSMB",#401:500
         "PI3", "TMSB10", "RPL8", "NDUFA4", "UBD", "SF3B4", "SLC12A7", "EDN1", "MEA1", "DGUOK", "SNRPC", "CDC26", "MRPS26", "MRPL3", "VWA1",#501:600
         "SMARCE1", "NUSAP1", "PTTG1", "HSPE1", "RPL30", "RPL28", "RPS27", "BCAP31", "GTF2A2", "NAXE", "BUD31", "COPS6", "STMP1", "POLR2K", "POLD2", "BANF1", "MDH2", "C1orf35", "AGRN", "GRB7", "TNRC18", "C7orf50", "PSMA2", "HMGB1", "GPX4", "S100A13", "COX6A1",
         "NCOR1", "LIPH", "GALNT6", "GCNT1", "PSMD4", "LSM2", "SIVA1", "SSB", "SF3B6", "TRMT112", "PSMB5", "RNF7", "ATP6V1F", "PDCD5", "PSMA6", "HSPD1", "HSP90AB1", "KRT8", "REG4", "MIF", "RPLP1", "RPS18", "TNFAIP2", "SERPINA1", "RPL39",
         "H2AFZ", "DBI", "HMGN2", "TUBB", "SUB1", "UQCRFS1", "ERBB2", "LLGL2", "TUBA1C", "KLF5", "CBR1", "PROM1", "CIB1", "SLPI", "SYNGR2", "ERGIC3", "DNMBP", "VPS13A", "RETSAT", "HECTD1", "MYL12B", "RPS7", "RPS5", "RPL36A", "RPLP2", "RPLP0", "RPL7A", "HIST1H1C", "IFI30", "RPS2", "RPS12", "RPS15", "FTH1", "RPL19", "RPL32", "RPS3"
         )
  Time_gene <- Time_gene[!(Time_gene %in% a)]
  pdf("./8.拟时序/pseudotime_heatmap1.pdf", width = 8, height = 12) 
  p=plot_pseudotime_heatmap(HSMM[Time_gene,], 
                            num_clusters=2,
                            show_rownames=F, 
                            return_heatmap=T)
  dev.off()
  #ggsave("Time_Marker.pdf", p, width = 10, height = 12,path = './photo/step 3-1/')
  
  #Pseudotime的热图
  # 提取每个簇的基因
  p$tree_row
  clusters <- cutree(p$tree_row, k = 3)
  clustering <- data.frame(clusters)
  clustering[,1] <- as.character(clustering[,1])
  colnames(clustering) <- "Gene_Clusters"
  table(clustering)
  
  #修改拟时序结果
  library(ggpubr)
  library(Seurat)
  df <- pData(test) 
  table(df$copykat.pred)
  # df$NewState <- ifelse(df$State == '1', 'State1',
  #                       ifelse(df$State == '2', 'State2', 'State3'))
  df$NewState <- df$copykat.pred
  table(df$NewState)
  
  State <- data.frame(State = df$NewState,
                      row.names = rownames(df),
                      Cell = rownames(df))
  State <- State[order(State$State),] 
  table(State$State)
  Neutrophil <- AddMetaData(Neutrophil, metadata = State)
  Idents(object = Neutrophil) <- "State"
  table(Idents(Neutrophil))
  
  Statemarker <- FindAllMarkers(Neutrophil)
  #这个结果上面有deg.cluster.Rdata中，但是我们要去重最后剩余991个基因
  #save(Statemarker,file = './8.拟时序/癌细胞_Statemarker.Rdata')
}
save(list = ls(), file = "./8.拟时序/HSMM.RData")

#修改后的可视化----------
{
  #热图
  table(Statemarker$cluster)
  top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  pdf("./8.拟时序/DoHeatmap.pdf", width = 8, height = 10) 
  DoHeatmap(Neutrophil,
            features = top50$gene,
            #group.by = "State",
            group.colors = c("red","blue"),#,"green"
            size = 3,
            assay = "RNA",
            angle = 0,
            label = F) #+
  dev.off()
  #scale_fill_gradientn(colors = c("white","grey","firebrick3"))#颜色
  #theme(axis.text.y = element_blank())#隐藏基因名
  # ggsave(filename = '基因热图.png',width = 8,height = 10,
  #        path = './')
  # ggsave(filename = '基因热图.pdf', width = 8, height = 10, path = './')#直接糊掉了
  #拟时序热图
  library(dplyr)
  table(Statemarker$cluster)
  top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
  gene0 <- unique(top50$gene)
  Time_diff <- differentialGeneTest(HSMM[gene0,], cores = 8, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
  # num_clusters为人为设置的聚类
  library(viridis)
  p <- plot_pseudotime_heatmap(HSMM[Time_genes,],
                               num_clusters = 3,
                               show_rownames=F,
                               return_heatmap=T,
                               hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
  # 默认的颜色就很好看
  # hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
  # hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))
  # hmcols = viridis(256))
  ggsave(filename = '拟时序热图.png',
         plot = p,width = 8,height = 6,path = './')
  ggsave(filename = '拟时序热图.pdf',
         plot = p, width = 8, height = 12, path = './')
}

#go分析每条通路的结果-----------
rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/GAC/8.拟时序/")
load(file = './HSMM.RData')
save(Statemarker,file = "./Statemarker.Rdata")
rm(list = ls())
load(file = './Statemarker.Rdata')
table(Statemarker$cluster)
duplicated_genes <- Statemarker$gene[duplicated(Statemarker$gene)]
table(duplicated_genes)

Statemarker$cluster <- ifelse(Statemarker$cluster == "Epithelial cell", "State1", "State2")

# 通路富集
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
GO_GSEA_list <- list()
KEGG_GSEA_list <- list()

for (i in 1:2) {
  name <- paste0('State',i)
  dat <- subset(Statemarker,Statemarker$cluster == name)
  dat$SYMBOL <- dat$gene
  
  ENTREZID1 <- bitr(unique(dat$SYMBOL), fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
  
  dat <- subset(dat,dat$SYMBOL %in% ENTREZID1$SYMBOL)
  dat <- merge(dat,ENTREZID1)
  geneList=dat$avg_log2FC
  names(geneList)=dat$ENTREZID
  geneList=sort(geneList,decreasing = T)
  
  #GO GSEA
  GO <- gseGO(
    geneList, #gene_fc
    ont = "ALL",# "BP"、"MF"和"CC"或"ALL"
    OrgDb = org.Hs.eg.db,#人类注释基因
    keyType = "ENTREZID",
    pvalueCutoff = 0.9,
    pAdjustMethod = "BH")#p值校正方法
  
  sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
  head(sortGO)
  dim(sortGO)
  # write.table(sortGO,paste0(paste0('State',i),'_','gsea_sortGO.txt')) #保存结果
  GO_GSEA_list[[i]] <- sortGO
  GO_GSEA_list[[i+3]] <- GO
  
  #KEGG GSEA
  KEGG <- gseKEGG(
    geneList,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.9,
    pAdjustMethod = "BH")
  
  
  sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
  head(sortKEGG)
  dim(sortKEGG)
  # write.table(sortKEGG,paste0(paste0('State',i),'_','gsea_sortKEGG.txt')) #保存结果
  KEGG_GSEA_list[[i]] <- sortKEGG
  KEGG_GSEA_list[[i+3]] <- KEGG
  
}
save(KEGG_GSEA_list,GO_GSEA_list,file = './State_GSEA.Rdata')
load(file = './State_GSEA.Rdata')

# 查看富集情况-KEGG
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggsci)
pal = pal_ucscgb()(20)

#State1 正常
State1_KEGG <- as.data.frame(KEGG_GSEA_list[1]) %>% filter(pvalue < 0.05)
write.csv(State1_KEGG,file = './State1_KEGG.csv')
State1_KEGG <- read.csv('./State1_KEGG.csv', row.names = 1)
#paths <- State1_KEGG$ID#选取你需要展示的通路ID
 #State1_KEGG <- State1_KEGG[c(4,),] 
paths <- c("hsa04668","hsa03050","hsa04670","hsa03010","hsa04514","hsa04530","hsa03040")
gseaplot2(KEGG_GSEA_list[[4]],paths, 
          pvalue_table = TRUE,
          title = 'Gastric epithelial cell KEGG Enrichment',
          color = ggsci::pal_aaas()(7),#color = pal[1:9]
          rel_heights = c(2, 0.3,0.3),
          subplots = c(1:3))
# ggsave(filename = 'State1 KEGG Enrichment.png',
#        width = 12,height = 8,
#        path = './')

#State2
State2_KEGG <- as.data.frame(KEGG_GSEA_list[2]) %>% filter(pvalue < 0.05)
write.csv(State2_KEGG,file = './State2_KEGG.csv')
State2_KEGG <- read.csv('./State2_KEGG.csv', row.names = 1)
#paths <- State2_KEGG$ID#选取你需要展示的通路ID
paths <- c("hsa04514","hsa03050","hsa04670","hsa03040","hsa03010","hsa04530","hsa04064","hsa04971")
gseaplot2(KEGG_GSEA_list[[5]],paths,
          pvalue_table = TRUE,
          title = 'Cancer cell KEGG Enrichment',
          color = ggsci::pal_aaas()(8),#color = pal[1:2]
          rel_heights = c(2, 0.3,0.3),
          subplots = c(1:3))
# paths <- c("hsa05171", "hsa03010")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment')
# ggsave(filename = 'State2 KEGG Enrichment.png',
#        width = 12,height = 8,
#        path = './')

















































