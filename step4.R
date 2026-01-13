#下载bulk数据：胃癌
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")
gc()
library(dplyr)
library(DESeq2)

load(file = "./1.download_bulk/训练集数据.Rdata")
table(train_phe$Type)

{
  dat <- train_dat
  dat <- na.omit(dat)
  dat <- 2^dat #-1 # FPKM转位Count
  All_dat <- round(dat,digits = 0) # 取整
  
  #构建dds矩阵
  countData <- All_dat
  condition <- factor(train_phe$Type)
  condition#注意:正常组在前、癌症在后
  #Levels:normal cancer
  head(condition)
  countData[1:4,1:4]
  
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
  head(dds)
  dim(dds)
  dds <- DESeq(dds)#时间较长
  resultsNames(dds)
  res <- results(dds)
  summary(res)
  
  DEG_INFO <- as.data.frame(res)
  
  # 提取差异分析结果
  # 获取padj（p值经过多重校验校正后的值）小于0.01，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
  table(res$padj<0.05) #取P值小于0.01的结果
  res <- res[order(res$padj),]
  diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因
  
  Up_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > 1)) 
  Up_gene_deseq2 <- row.names(Up_gene_deseq2) # 上调差异基因
  
  Down_gene_deseq2 <-  subset(res,padj < 0.05 & (log2FoldChange <  -1))
  Down_gene_deseq2 <- row.names(Down_gene_deseq2) # 下调差异基因
  
  resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  
  save(resdata,diff_gene_deseq2,Up_gene_deseq2,Down_gene_deseq2,
       file = "./4.DEG/DEseq2.Rdata")
}
load(file = "./4.DEG/DEseq2.Rdata")

#火山图----
{
  res = resdata
  rownames(res) = res$Row.names
  res = res[,-1]
  library(ggplot2)
  for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                            'padj' = res$padj,
                            'State' = rep('No', length(res$log2FoldChange)),
                            row.names = rownames(res))
  up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05))
  down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05))
  for_volcano[up_sig_indices,'State'] <- 'Up'
  for_volcano[down_sig_indices,'State'] <- 'Down'
  for_volcano$State <- as.factor(for_volcano$State)
  for_volcano$padj <- -log10(for_volcano$padj)
  
  this_tile <- paste0('Cutoff for logFC is 1',
                      '\nThe number of Up gene is ',nrow(for_volcano[for_volcano$State =='Up',]) ,
                      '\nThe number of Down gene is ',nrow(for_volcano[for_volcano$State =='Down',]))
  
  p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
    geom_point(size = I(1))+
    scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
    geom_vline(xintercept = c(1, -1), lty=2, size=I(0.4), colour = 'grey11')+
    geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
          panel.grid = element_blank())+
    labs(x='log2FoldChange', y = '-log10Pvalue')+
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=15,hjust = 0.5)) 
  
  p#保存pdf8*6
}

#随便选出一个上调基因进行箱线图验证----
{
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  library(stringr)
  library(cowplot)
  
  rm(list = ls())
  options(stringsAsFactors = F)
  setwd('/home/datahup/syj/GAC/')
  
  load(file = "./1.download_bulk/训练集数据.Rdata")
  load(file = "./4.DEG/DEseq2.Rdata")
  
  #上调基因验证
  #Gene <- Up_gene_deseq2[1:5]
  expr1 <- as.data.frame(t(train_dat[intersect_all,]))
  expr1$Group <- train_phe$Type
  #expr1$Group <- ifelse(expr1$Group == "cancer","Tumor","Normal")
  rownames(expr1)
  # 箱图
  { # [1] "COL10A1" "CST1"    "ESM1"    "MTHFD1L" "COL11A1"
    #[1] "CXCL8"   "FGFR4"   "CXCL16"  "CCL20"   "MUC4"    "PLA2G2A"
    p1 <- ggplot(data=expr1)+ 
      geom_boxplot(mapping=aes(x=Group,y=PLA2G2A,colour = Group ), 
                   alpha = 0.5,
                   size=1.2,
                   width = 0.6)+ 
      geom_jitter(mapping=aes(x=Group,y=PLA2G2A,colour = Group), #散点
                  alpha = 0.3,size=2)+
      scale_color_manual(limits=c("Tumor","Normal"), 
                         values=c("red","blue"))+ #颜色
      geom_signif(mapping=aes(x=Group,y=PLA2G2A), # 不同组别的显著性
                  comparisons = list(c("Tumor", "Normal")), 
                  map_signif_level=T, # T显示显著性，F显示p value
                  # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
                  # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
                  size=1, # 修改线的粗细
                  textsize = 6, # 修改显著性标记的大小
                  test = "t.test")+ # 检验的类型
      theme_classic(  # 主题设置，这个是无线条主题
        base_line_size = 1 # 坐标轴的粗细
      )+
      guides(colour="none")+ # 删除图注
      labs(title='PLA2G2A',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
      theme(plot.title = element_text(size = 18,
                                      colour = "black",
                                      face = "bold",
                                      vjust = 1.9,
                                      hjust = 0.5),
            axis.title.y = element_text(size = 12, 
                                        # family = "myFont", 
                                        color = "black",
                                        face = "bold.italic", 
                                        vjust = 1.9, 
                                        hjust = 0.5, 
                                        angle = 90),
            axis.title.x = element_blank(),#删除X轴标签
            legend.title = element_text(color="black", # 修改图例的标题
                                        size=15, 
                                        face="bold"),
            legend.text = element_text(color="black", # 设置图例标签文字
                                       size = 10, 
                                       face = "bold"),
            axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                       color = "black", # 颜色
                                       face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                       vjust = 0.5, # 位置
                                       hjust = 0.5, 
                                       angle = 0), #角度
            axis.text.y = element_text(size = 13,  
                                       color = "black",
                                       face = "bold.italic", 
                                       vjust = 0.5, 
                                       hjust = 0.5, 
                                       angle = 0) 
      )
    
    p1#pdf8*6
  }
  
  #下调基因验证
  Gene <- Down_gene_deseq2[1:5]
  expr1 <- as.data.frame(t(train_dat[Gene,]))
  expr1$Group <- train_phe$Type
  #expr1$Group <- ifelse(expr1$Group == "cancer","Tumor","Normal")
  #rownames(expr1)
  Gene
  # 箱图
  { #[1] "ENPP7"   （"SYNPO2L" "S100G"   "TFAP2B" 这三个没有显著型） "FLG"  
    p1 <- ggplot(data=expr1)+ 
      geom_boxplot(mapping=aes(x=Group,y= FLG,colour = Group ), 
                   alpha = 0.5,
                   size=1.2,
                   width = 0.6)+ 
      geom_jitter(mapping=aes(x=Group,y= FLG,colour = Group), #散点
                  alpha = 0.3,size=2)+
      scale_color_manual(limits=c("Tumor","Normal"), 
                         values=c("red","blue"))+ #颜色
      geom_signif(mapping=aes(x=Group,y= FLG), # 不同组别的显著性
                  comparisons = list(c("Tumor", "Normal")), 
                  map_signif_level=T, # T显示显著性，F显示p value
                  # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
                  # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
                  size=1, # 修改线的粗细
                  textsize = 6, # 修改显著性标记的大小
                  test = "t.test")+ # 检验的类型
      theme_classic(  # 主题设置，这个是无线条主题
        base_line_size = 1 # 坐标轴的粗细
      )+
      guides(colour="none")+ # 删除图注
      labs(title='FLG',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
      theme(plot.title = element_text(size = 18,
                                      colour = "black",
                                      face = "bold",
                                      vjust = 1.9,
                                      hjust = 0.5),
            axis.title.y = element_text(size = 12, 
                                        # family = "myFont", 
                                        color = "black",
                                        face = "bold.italic", 
                                        vjust = 1.9, 
                                        hjust = 0.5, 
                                        angle = 90),
            axis.title.x = element_blank(),#删除X轴标签
            legend.title = element_text(color="black", # 修改图例的标题
                                        size=15, 
                                        face="bold"),
            legend.text = element_text(color="black", # 设置图例标签文字
                                       size = 10, 
                                       face = "bold"),
            axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                       color = "black", # 颜色
                                       face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                       vjust = 0.5, # 位置
                                       hjust = 0.5, 
                                       angle = 0), #角度
            axis.text.y = element_text(size = 13,  
                                       color = "black",
                                       face = "bold.italic", 
                                       vjust = 0.5, 
                                       hjust = 0.5, 
                                       angle = 0) 
      )
    
    p1#pdf8*6
  }
  
}
#美化----
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/')

load(file = "./1.download_bulk/训练集数据.Rdata")
load(file = "./4.DEG/DEseq2.Rdata")

#第一种
{
  library(ggplot2)
  data <- resdata
  data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))
  
  res = resdata
  rownames(res) = res$Row.names
  res = res[,-1]
  for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                            'padj' = res$padj,
                            'State' = rep('No', length(res$log2FoldChange)),
                            row.names = rownames(res))
  up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05))
  down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05))
  for_volcano[up_sig_indices,'State'] <- 'Up'
  for_volcano[down_sig_indices,'State'] <- 'Down'
  for_volcano$State <- as.factor(for_volcano$State)
  for_volcano$padj <- -log10(for_volcano$padj)
  
  this_tile <- paste0('Cutoff for logFC is 1',
                      '\nThe number of Up gene is ',nrow(for_volcano[for_volcano$State =='Up',]) ,
                      '\nThe number of Down gene is ',nrow(for_volcano[for_volcano$State =='Down',]))
  
  p <- ggplot(data,aes(log2FoldChange, -log10(padj)))+
    # 横向水平参考线：
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
    # 纵向垂直参考线：
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
    # 散点图:
    geom_point(aes(size=-log10(padj), color= -log10(padj)))+
    # 指定颜色渐变模式：
    scale_color_gradientn(values = seq(0,1,0.2),
                          colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
    # 指定散点大小渐变模式：
    scale_size_continuous(range = c(1,1.5))+
    # 主题调整：
    theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=15,hjust = 0.5)) 
  
  p#pdf为8*6
}

#第二种
{
  # data <- resdata
  # data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))
  
  ggplot(data,aes(log2FoldChange, -log10(padj)))+
    # 横向水平参考线：
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
    # 纵向垂直参考线：
    geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
    # 散点图:
    geom_point(aes(size=-log10(padj), color= -log10(padj)))+
    # 指定颜色渐变模式：
    scale_color_gradientn(values = seq(0,1,0.2),
                          colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
    # 指定散点大小渐变模式：
    scale_size_continuous(range = c(1,3))+
    # 主题调整：
    theme_bw()+
    # 调整主题和图例位置：
    theme(panel.grid = element_blank(),
          legend.position = c(0.01,0.7),
          legend.justification = c(0,1)
    )+
    # 设置部分图例不显示：
    guides(col = guide_colourbar(title = "-Log10_q-value"),
           size = "none")+
    # 添加标签：
    geom_text(aes(label=label, color = -log10(padj)), size = 3, vjust = 1.5, hjust=1)+
    # 修改坐标轴：
    xlab("Log2FC")+
    ylab("-Log10(FDR q-value)")+
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=15,hjust = 0.5)) 
}


#差异分析热图-----
{
  rm(list = ls())
  options(stringsAsFactors = F)
  setwd('/home/datahup/syj/GAC/')
  
  load(file = "./1.download_bulk/训练集数据.Rdata")
  load(file = "./4.DEG/DEseq2.Rdata")
  
  #Up_gene_deseq2 <- subset(resdata,padj < 0.05 & (log2FoldChange > 1)) 
  Up_gene_deseq2 <- Up_gene_deseq2[1:100]#<- row.names(Up_gene_deseq2) # 上调差异基因
  #Down_gene_deseq2 <-  subset(resdata,padj < 0.05 & (log2FoldChange <  -1))
  Down_gene_deseq2 <- Down_gene_deseq2[1:100]# 下调差异基因
  
  # 合并上下调基因
  deg_genes <- c(Up_gene_deseq2, Down_gene_deseq2)
  
  # 2. 提取表达矩阵
  mat <- train_dat[rownames(train_dat) %in% deg_genes, ]
  
  # 3. 行标准化（Z-score）
  mat <- t(scale(t(as.matrix(mat))))
  
  # 4. 基因分组信息（上调/下调）
  gene_group <- data.frame(Regulation = ifelse(rownames(mat) %in% Up_gene_deseq2, "Up", "Down"))
  rownames(gene_group) <- rownames(mat)   # 把行名对齐到 mat
  table(gene_group$Regulation)
  
  # 5. 画热图
  library(pheatmap)
  pheatmap(mat,
           annotation_row = gene_group,  # 标记上下调
           show_rownames = FALSE,
           show_colnames = F,
           cluster_rows = T,
           cluster_cols = F,
           main = "DEG Heatmap (Up & Down)")
  
}





























































