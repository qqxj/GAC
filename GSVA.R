
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/')

load(file = "./1.download_bulk/dat.output.Rdata")
load(file = "./1.download_bulk/phe.output.Rdata")
load(file = "./12.2RSF11.11/RSF_genes.Rdata")

gene <- RSF_genes
#基因矩阵
NSCLCcount <- dat
dat = as.data.frame(t(NSCLCcount))
phe <- subset(phe,phe$Type == "Tumor")
dat0 <- dat[phe$Id,]
df = data.frame(gene = dat0[,gene],
                row.names = rownames(dat0)) 
colnames(df)
colnames(df) <- sub("^gene\\.", "", colnames(df))

#Step1 高低分组做差异分析----------
library(DESeq2)
setwd('/home/datahup/syj/GAC/15.探索基因功能/')
for (i in 1:length(gene)) {
  df = df[order(df[,i]),]
  a = quantile(df[,i],c(0.4,0.6))
  df$group = ifelse(df[,i] <= a[1],'Low',
                    ifelse(df[,i] >= a[2],'High', 'no'))
  table(df$group)  
  df1 = subset(df,df$group == 'High' | df$group == 'Low')
  
  #构建dds矩阵
  dat1 = as.data.frame(t(dat0[rownames(df1),]))
  range(dat1)
  dat1 <- 2^dat1 -1 #转Count
  dat1 <- round(dat1,digits = 0) # 取整
  dat1[dat1 < 0 ] = 0
  countData <- dat1
  countData[1:6,1:6]
  Group <- data.frame(group = df1$group, row.names = rownames(df1))
  condition <- factor(Group$group)
  table(condition)
  head(condition)
  # 差异分析
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
  head(dds)
  dim(dds)
  dds <- DESeq(dds) 
  resultsNames(dds)
  res <- results(dds)
  summary(res)
  DEG_INFO <- as.data.frame(res)
  # 提取差异分析结果
  table(res$padj<0.05) #取P值小于0.05的结果
  res <- res[order(res$padj),]
  resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  rownames(resdata) <- resdata$Row.names
  resdata <- resdata[,-1]
  resdata <- na.omit(resdata)
  diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因
  diff_gene_deseq2_df <- resdata[diff_gene_deseq2,1:6]
  diff_gene_deseq2_df$state <- ifelse(diff_gene_deseq2_df$log2FoldChange > 1, 'Up','Down')
  table(diff_gene_deseq2_df$state)
  
  write.csv(diff_gene_deseq2_df,file = paste0(gene[i],'_DEG.csv'))
  
  # 火山图
  library(ggplot2)
  for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                            'padj' = res$padj,
                            'State' = rep('No', length(res$log2FoldChange)))
  up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05))
  down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05))
  for_volcano[up_sig_indices,'State'] <- 'Up'
  for_volcano[down_sig_indices,'State'] <- 'Down'
  for_volcano$State <- as.factor(for_volcano$State)
  for_volcano$padj <- -log10(for_volcano$padj)
  
  this_tile <- paste0('Cutoff for logFC is 1',
                      '\nThe number of Up gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Up',]) ,
                      '\nThe number of Down gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Down',]))
  
  
  p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
    geom_point(size = I(0.7))+
    scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
    geom_vline(xintercept = c(1, -1), lty=2, size=I(0.4), colour = 'grey11')+
    geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
          panel.grid = element_blank())+
    labs(x='log2FoldChange', y = '-log10Pvalue')+
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=15,hjust = 0.5)) 
  
  p
  
  ggsave(filename = paste0(gene[i],'_vol.png'),plot = p,width = 8,height = 6,
         path = "./")
  
}

#  go|KEGG-GSVA--------
#rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/15.探索基因功能/')
library(stringr)
#kegg
kk_list = list()
kegg_gsva_list = list()
#go
go_list = list()
go_gsva_list = list()

# 加载GO通路数据库
library(clusterProfiler)
{
  source("getGoTerm.R")
  GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
  save(GO_DATA, file = "GO_DATA.RData")
  findGO <- function(pattern, method = "key"){
    
    if(!exists("GO_DATA"))
      load("GO_DATA.RData")
    if(method == "key"){
      pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
    } else if(method == "gene"){
      pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
    }
    
    colnames(pathways) = "pathway"
    
    if(length(pathways) == 0){
      cat("No results!\n")
    } else{
      return(pathways)
    }
  } # 用于寻找 GO ID
  getGO <- function(ID){
    
    if(!exists("GO_DATA"))
      load("GO_DATA.RData")
    allNAME = names(GO_DATA$PATHID2EXTID)
    if(ID %in% allNAME){
      geneSet = GO_DATA$PATHID2EXTID[ID]
      names(geneSet) = GO_DATA$PATHID2NAME[ID]
      return(geneSet)     
    } else{
      cat("No results!\n")
    }
  } # 获取 GO geneSet
  load("GO_DATA.RData") # 载入数据 GO_DATA
}

# TYMP ----
colnames(df)
{
  df = df[order(df[,1]),]
  a = quantile(df[,1],c(0.4,0.6))
  df[,paste0(gene[1],'_group')] = ifelse(df[,1] <= a[1],'Low',
                                         ifelse(df[,1] >= a[2],'High', 'no'))
  table(df$TYMP_group)
  df1 = subset(df,df[,paste0(gene[1],'_group')] == 'High' | df[,paste0(gene[1],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[1],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  #deg_entr$ENTREZID <- as.numeric(deg_entr$ENTREZID)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.5)
  kk_list[[gene[1]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  #kegg-GSVA
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[1],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[1]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[1]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  
  #go-gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) #gsva
  go_gsva_list[[gene[1]]] <- go_gsva
  
  #print(paste0(gene[1],' is over'))
  save(kk_list,kegg_gsva_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

#IFNG----
colnames(df)
{
  df = df[order(df[,2]),]
  a = quantile(df[,2],c(0.4,0.6))
  df[,paste0(gene[2],'_group')] = ifelse(df[,2] <= a[1],'Low',
                                         ifelse(df[,2] >= a[2],'High', 'no'))

  table(df[,paste0(gene[2],'_group')])  
  df1 = subset(df,df[,paste0(gene[2],'_group')] == 'High' | df[,paste0(gene[2],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[2],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.5)
  kk_list[[gene[2]]] <- kk_deg
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[2],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[2]]] <- kegg_gsva
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.5, 
                     qvalueCutoff   = 1, 
                     readable       = TRUE)
  go_list[[gene[2]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  go_gsva <- gsva(df3, go_genelist, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) #gsva
  go_gsva_list[[gene[2]]] <- go_gsva
  
  #print(paste0(gene[6],' is over'))
  save(kk_list,kegg_gsva_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

#ITGAX----
colnames(df)
{
  df = df[order(df[,3]),]
  a = quantile(df[,3],c(0.4,0.6))
  df[,paste0(gene[3],'_group')] = ifelse(df[,3] <= a[1],'Low',
                                         ifelse(df[,3] >= a[2],'High', 'no'))
  
  table(df[,paste0(gene[3],'_group')])  
  df1 = subset(df,df[,paste0(gene[3],'_group')] == 'High' | df[,paste0(gene[3],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[3],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.5,
                       qvalueCutoff =0.5)
  kk_list[[gene[3]]] <- kk_deg
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[3],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist, 
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[3]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[3]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[3]]] <- go_gsva
  save(kk_list,kegg_gsva_list,go_gsva_list,file = "./基因富集结果.Rdata")

  #print(paste0(gene[11],' is over'))
}

#GBP5----
colnames(df)
{
  df = df[order(df[,4]),]
  a = quantile(df[,4],c(0.4,0.6))
  df[,paste0(gene[4],'_group')] = ifelse(df[,4] <= a[1],'Low',
                                         ifelse(df[,4] >= a[2],'High', 'no'))
  table(df[,paste0(gene[4],'_group')])  
  df1 = subset(df,df[,paste0(gene[4],'_group')] == 'High' | df[,paste0(gene[4],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[4],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
  kk_list[[gene[4]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[4],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[4]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[4]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #go-gsva
  go_gsva <- gsva(df3, go_genelist, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[4]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
  save(kk_list,kegg_gsva_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

#GBP4----
colnames(df)
{
  df = df[order(df[,5]),]
  a = quantile(df[,5],c(0.4,0.6))
  df[,paste0(gene[5],'_group')] = ifelse(df[,5] <= a[1],'Low',
                                         ifelse(df[,5] >= a[2],'High', 'no'))
  table(df[,paste0(gene[5],'_group')])  
  df1 = subset(df,df[,paste0(gene[5],'_group')] == 'High' | df[,paste0(gene[5],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[5],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
  kk_list[[gene[5]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[5],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[5]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[5]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #go-gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[5]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
  save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

#STAT1----
colnames(df)
{
  df = df[order(df[,6]),]
  a = quantile(df[,6],c(0.4,0.6))
  df[,paste0(gene[6],'_group')] = ifelse(df[,6] <= a[1],'Low',
                                         ifelse(df[,6] >= a[2],'High', 'no'))
  table(df[,paste0(gene[6],'_group')])  
  df1 = subset(df,df[,paste0(gene[6],'_group')] == 'High' | df[,paste0(gene[6],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[6],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.5)
  kk_list[[gene[6]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[6],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[6]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[6]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #go-gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[6]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
  save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

#CD84----
colnames(df)
{
  df = df[order(df[,7]),]
  a = quantile(df[,7],c(0.4,0.6))
  df[,paste0(gene[7],'_group')] = ifelse(df[,7] <= a[1],'Low',
                                         ifelse(df[,7] >= a[2],'High', 'no'))
  table(df[,paste0(gene[7],'_group')])  
  df1 = subset(df,df[,paste0(gene[7],'_group')] == 'High' | df[,paste0(gene[7],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[7],'_DEG.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
  kk_list[[gene[7]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[7],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[7]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[7]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #go-gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[7]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
  save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,file = "./基因富集结果.Rdata")
}

save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
     file = './Step6-1 Result_1.Rdata')

# 6.28 GSVA得分差异分析 找到关键通路--------
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/15.探索基因功能/')

load(file = 'Step6-1 Result_1.Rdata')
load(file = "../12.2RSF11.11/RSF_genes.Rdata")
gene = RSF_genes

KEGG_logFC_list = list()
GO_logFC_list = list()

#go和kegg的gsva差异分析(循环,跳过挨个跑)-------
for (i in gene) {
  
  kk_re = kk_list[[i]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[i]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  
  go_re = go_list[[i]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[i]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  
  group = data.frame(row.names = rownames(df), group = df[,paste0(i,'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[i]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[i]] = GO_logFC
  
  
}
# save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
#      KEGG_logFC_list,GO_logFC_list,
#      file = 'Step6-1 Result.Rdata')

library(stringr)
#TYMP------
gene
{
  kk_re = kk_list[[1]]@result
  # kk_re_1 = as.data.frame(kk_list[["PMEPA1"]])
  # table(kk_re$Description %in% kk_re_1$Description)
  
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[1]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[1]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[1]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("TYMP",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[1]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[1]] = GO_logFC
  
}

#IFNG-------
gene
{
  kk_re = kk_list[[2]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[2]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[2]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[2]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  go_gsva_df <- na.omit(go_gsva_df)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("IFNG",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[2]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[2]] = GO_logFC
  
}

#ITGAX------
gene
{
  kk_re = kk_list[[3]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[3]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[3]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[3]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("ITGAX",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[3]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[3]] = GO_logFC
  
}

#GBP5-------
gene
{
  kk_re = kk_list[[4]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[4]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[4]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[4]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("GBP5",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[4]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[4]] = GO_logFC
  
}

#GBP4------
gene
{
  kk_re = kk_list[[5]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[5]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[5]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[5]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("GBP4",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[5]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[5]] = GO_logFC
  
}

#STAT1------
gene
{
  kk_re = kk_list[[6]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[6]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[6]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[6]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("STAT1",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[6]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[6]] = GO_logFC
  
}
#CD84------
gene
{
  kk_re = kk_list[[7]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[7]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[7]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[7]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("CD84",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[7]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[7]] = GO_logFC
  
}
save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
     KEGG_logFC_list,GO_logFC_list,
     file = 'Step6-1 Result.Rdata')

# 作图--------
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/15.探索基因功能/')

load(file = 'Step6-1 Result.Rdata')
load(file = "../12.2RSF11.11/RSF_genes.Rdata")
library(dplyr)
library(tidyr)
library(tibble)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(stringr)

gene <- RSF_genes
# TYMP-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[1]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["TYMP"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("NOD-like receptor signaling pathway",
                              "Th1 and Th2 cell differentiation",
                              "Th17 cell differentiation",
                              "Toll-like receptor signaling pathway",
                              "Chemokine signaling pathway",
                              "PPAR signaling pathway",
                              "Glycerolipid metabolism",
                              "Hormone signaling",
                              "Wnt signaling pathway",
                              "Fat digestion and absorption",
                              "Linoleic acid metabolism",
                              "Vascular smooth muscle contraction",
                              "cGMP-PKG signaling pathway")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("TYMP",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in TYMP Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[1]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.25)
  
  go_gsva_df = as.data.frame(go_gsva_list[["TYMP"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c( "lymphocyte chemotaxis",
                               "chemokine activity",
                               "cellular response to type II interferon",
                               "pyroptosis",
                               "chemokine-mediated signaling pathway",
                               "MHC protein complex binding",
                               "CXCR chemokine receptor binding",
                               "peptide antigen assembly with MHC class II protein complex",
                               "peptide antigen assembly with MHC protein complex")]
   
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in TYMP Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  #mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in PMEPA1 Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])
  #   + stat_compare_means(aes(group = group),
  #   + method = "wilcox.test",
  #   + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','TYMP',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./TYMP总.pdf", 
         plot = p, 
         width = 8, 
         height = 12)
}

# IFNG-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[2]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["IFNG"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("Natural killer cell mediated cytotoxicity",
                              "Cytokine-cytokine receptor interaction",
                              "Chemokine signaling pathway",
                              "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                              "Toll-like receptor signaling pathway",
                              "Th1 and Th2 cell differentiation",
                              "Th17 cell differentiation",
                              "Inflammatory bowel disease",
                              "cAMP signaling pathway",
                              "cGMP-PKG signaling pathway")]
                              #"Vascular smooth muscle contraction",
  #"Gastric acid secretion"
                              
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("IFNG",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in IFNG Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[2]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.25)
  
  go_gsva_df = as.data.frame(go_gsva_list[["IFNG"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c("regulatory T cell differentiation",
                              "cytokine receptor binding",
                              "alpha-beta T cell activation",
                              "response to type II interferon",
                              "natural killer cell mediated immunity",
                              "regulation of immune effector process",
                              "pyroptosis",
                              "lymphocyte mediated immunity",
                              "interleukin-1 beta production"
                              )]#"cell killing"#"positive regulation of tyrosine phosphorylation of STAT protein",
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in IFNG Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in IFNG Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','IFNG',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./IFNG总.pdf", 
         plot = p, 
         width = 8, 
         height = 12)
}

# ITGAX-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[3]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["ITGAX"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("JAK-STAT signaling pathway",
                              "Chemokine signaling pathway",
                              "Cytokine-cytokine receptor interaction",
                              "IgSF CAM signaling",
                              "Toll-like receptor signaling pathway",
                              "Inflammatory bowel disease",
                              "Th17 cell differentiation",
                              "Rap1 signaling pathway",
                              "Th1 and Th2 cell differentiation")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("ITGAX",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in ITGAX Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[3]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.4)
  
  go_gsva_df = as.data.frame(go_gsva_list[["ITGAX"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c( "integrin-mediated signaling pathway",
                               "Fc receptor signaling pathway",
                               "immune response-regulating cell surface receptor signaling pathway involved in phagocytosis",
                               "regulation of phagocytosis",
                               #"myeloid leukocyte activation",
                               "positive regulation of interleukin-4 production",
                               "mast cell activation",
                               #"regulation of podosome assembly",
                               "interleukin-2 production",
                               #"inhibitory MHC class I receptor activity",
                               "regulation of macrophage chemotaxis",
                               "T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell",
                               "regulation of type 2 immune response")]
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in ITGAX Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in ITGAX Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','ITGAX',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./ITGAX总.pdf", 
         plot = p, 
         width = 10.5, 
         height = 12)
}

# GBP5-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[4]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["GBP5"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("NOD-like receptor signaling pathway",
                              "Toll-like receptor signaling pathway",
                              "Cytokine-cytokine receptor interaction",
                              "Chemokine signaling pathway",
                              "Natural killer cell mediated cytotoxicity",
                              "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                              "Th17 cell differentiation",
                              "Th1 and Th2 cell differentiation",
                              "T cell receptor signaling pathway",
                              "IgSF CAM signaling",
                              "Phagosome")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("GBP5",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in GBP5 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[4]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.35)
  
  go_gsva_df = as.data.frame(go_gsva_list[["GBP5"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c("regulation of CD4-positive, alpha-beta T cell activation",
                              "positive regulation of alpha-beta T cell differentiation",
                              "response to type II interferon",
                              "positive regulation of leukocyte differentiation",
                              "pyroptosis",
                              "natural killer cell activation",
                              "positive regulation of lymphocyte differentiation",
                              "cytokine receptor activity")]
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP5 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP5 Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','GBP5',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./GBP5总.pdf", 
         plot = p, 
         width = 8, 
         height = 12)
}

# GBP4-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[5]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["GBP4"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("NOD-like receptor signaling pathway",
                              "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                              "Chemokine signaling pathway",
                              "Toll-like receptor signaling pathway",
                              "Cytokine-cytokine receptor interaction",
                              "T cell receptor signaling pathway",
                              "Th1 and Th2 cell differentiation",
                              "Th17 cell differentiation",
                              "Inflammatory bowel disease",
                              "IgSF CAM signaling",
                              "Phagosome")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("GBP4",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in GBP4 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[5]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.4)
  
  go_gsva_df = as.data.frame(go_gsva_list[["GBP4"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c("response to type II interferon",
                              "pyroptosis",
                              "regulation of alpha-beta T cell activation",
                              "positive regulation of mononuclear cell proliferation",
                              "lymphocyte mediated immunity",
                              "CD8-positive, alpha-beta T cell activation",
                              "antigen receptor-mediated signaling pathway",
                              "immune response to tumor cell")]
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP4 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP5 Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','GBP4',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./GBP4总.pdf", 
         plot = p, 
         width = 8, 
         height = 12)
}

# STAT1-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[6]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["STAT1"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("Toll-like receptor signaling pathway",
                              "NOD-like receptor signaling pathway",
                              "Chemokine signaling pathway",
                              "Vascular smooth muscle contraction",
                              "cGMP-PKG signaling pathway",
                              "Gastric acid secretion")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("STAT1",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in STAT1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[6]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.15)
  
  go_gsva_df = as.data.frame(go_gsva_list[["STAT1"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c("cellular response to type II interferon",
                              "response to type II interferon",
                              "pyroptosis",
                              "chemokine receptor binding",
                              "cell killing",
                              "leukocyte mediated cytotoxicity",
                              "leukocyte mediated immunity",
                              "CXCR chemokine receptor binding",
                              "cytokine activity",
                              "T cell chemotaxis",
                              "neutrophil chemotaxis",
                              "lymphocyte chemotaxis")]
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in STAT1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP5 Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','STAT1',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./STAT1总.pdf", 
         plot = p, 
         width = 6.9, 
         height = 12)
}

# CD84-----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[7]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["CD84"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  colnames(kk_gsva_df)
  kk_gsva_df <- kk_gsva_df[,c("IgSF CAM signaling",
                              "Leukocyte transendothelial migration",
                              "Phagosome",
                              "Chemokine signaling pathway",
                              "Natural killer cell mediated cytotoxicity",
                              "Cytokine-cytokine receptor interaction",
                              "Th17 cell differentiation",
                              "Th1 and Th2 cell differentiation",
                              "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                              "T cell receptor signaling pathway",
                              "Toll-like receptor signaling pathway",
                              "NF-kappa B signaling pathway")]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("CD84",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in CD84 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  # plot1 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  #   # split violin
  #   geom_violin(alpha = 1, trim = FALSE, color = NA, width = 1.2) +
  #   #geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   #geom_boxplot(outlier.shape = 21,color = "black") + 
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  #   coord_flip() 
  # plot1
  # 保存为PDF文件
  # ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
  #        plot = plot1, width = 8, height = 12)
  # GO
  # go_re = go_list[["PMEPA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[7]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.65)
  
  go_gsva_df = as.data.frame(go_gsva_list[["CD84"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c("positive regulation of tumor necrosis factor superfamily cytokine production",
                              "regulation of macrophage activation",
                              "interleukin-10 production",
                              "endolysosomal toll-like receptor signaling pathway",
                              "positive regulation of macrophage chemotaxis",
                              "CD40 signaling pathway",
                              "neutrophil activation involved in immune response")]
  
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in CD84 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  #格式800*1200
  # plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  #   # split violin
  #   geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  #   # mean point
  #   stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  #   # errorbar
  #   stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
  #                size = 0.3,
  #                position = position_dodge(0.2)) +
  #   theme_bw() + 
  #   labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  #   stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  #   coord_flip() 
  # ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)
  
  
  #将两者合并
  colnames(kk_gsva_df)
  colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
  colnames(go_gsva_df)
  colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
  kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
  merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
  rownames(merged_df) <- merged_df$Row.names
  merged_df <- merged_df[,-1]
  
  merged_df2 <- merged_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, 
                            group$group, 'unkown')
  table(merged_df2$group)
  
  mypalette = pal_ucscgb()(10)
  # ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  #   geom_boxplot(outlier.shape = 21,color = "black") + 
  #   theme_bw() + 
  #   labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GBP5 Grouping') +
  #   theme(legend.position = "right") + 
  #   theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
  #         plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  #   scale_fill_manual(values = mypalette[c(1,6,3)])+
  #   stat_compare_means(aes(group = group),
  #                      + method = "wilcox.test",
  #                      + label.y = 1.1,label = 'p.signif')
  
  #ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')
  
  library(ggplot2)
  library(ggunchained) 
  # library(devtools)  
  # install_github("JanCoUnchained/ggunchained")
  
  mypalette = pal_ucscgb()(26)
  #格式800*1200
  p <- ggplot(merged_df2,aes(x = Pathway, y = `GSVA_Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','CD84',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#953984","#FFC82D"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + 
    coord_flip() 
  p
  ggsave("./CD84总.pdf", 
         plot = p, 
         width = 8.8, 
         height = 12)
}
