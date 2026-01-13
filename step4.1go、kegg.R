
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")

load(file = "./4.DEG/DEseq2.Rdata")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	
library(GOplot)
#转换基因id------
input_gene=unique(as.vector(diff_gene_deseq2))
entrezIDs=mget(input_gene, 
               org.Hs.egSYMBOL2EG, 
               ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO----------
kk=enrichGO(gene=gene,
            keyType = "ENTREZID",
            OrgDb=org.Hs.eg.db,
            pvalueCutoff=0.5,
            qvalueCutoff=0.5,
            ont="ALL", 
            readable=T)
save(kk,file = "./4.2.癌症组差异基因的go、kegg分析/go.Rdata")
GO=as.data.frame(kk)
GO=GO[(GO$pvalue < 0.05 & GO$qvalue < 0.05),]
GO <- subset(GO,Count > 10)

#挑选与课题相关的活性通路
#GO <- GO[c(1:20),]
GO$Description[1201:1247]
GO <- subset(GO,GO$ID %in% c("GO:0070374","GO:0030198","GO:0097529","GO:0002065","GO:0030858","GO:0043114","GO:0010463",#胃癌相关的
                             "GO:0006959","GO:0002526","GO:0002688","GO:0006910","GO:0050727","GO:0090025","GO:0019933",
                             "GO:0045580","GO:0010463","GO:005086","GO:0006898","GO:0002698","GO:0045807"))#免疫相关的
GO$Description                             
#在此筛选自己想要的通路
selected_GO_ids <- c("GO:0070374","GO:0002065","GO:0030858",#"GO:0030198","GO:0097529",
                     "GO:0010463", #胃癌相关的#"GO:0043114",
                     "GO:0006959","GO:0006910","GO:0050727",#"GO:0002526","GO:0002688",
                     "GO:0090025","GO:0019933","GO:0045580","GO:0050863","GO:0006898",
                     "GO:0045807") #免疫相关的  #"GO:0002698",
# 筛选GO数据框
GO_filtered <- subset(GO, ID %in% selected_GO_ids)
GO_filtered$Description   
# 创建新的enrichResult对象只包含筛选的通路
kk_filtered <- kk  
kk_filtered@result <- kk_filtered@result[kk_filtered@result$ID %in% selected_GO_ids, ]


#柱状图
#pdf(file="GObarplot.pdf", width=10, height=7)
barplot(kk_filtered, 
        drop=TRUE, 
        showCategory=20,
        label_format=50, 
        split="ONTOLOGY") +  # 移除 color 参数
  facet_grid(ONTOLOGY~., scale='free') +
  ggtitle("Selected GO Terms Related to Gastric Cancer and Immunity")
print(bar)
#dev.off()
barplot(kk_filtered, 
        drop=TRUE, 
        showCategory=20,
        label_format=50, 
        split="ONTOLOGY",
        color = "pvalue") +  # 使用p值作为颜色深浅
  facet_grid(ONTOLOGY~., scale='free') +
  ggtitle("Selected GO Terms Related to Gastric Cancer and Immunity") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

#气泡图
#pdf(file="GObubble.pdf", width=10, height=7)
dotplot(kk_filtered, showCategory=20,
            orderBy="GeneRatio", 
            label_format=130,
            split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

#dev.off()

cnetplot(kk_filtered) #网络图展示富集功能和基因的包含关系
emapplot(kk_filtered) #网络图展示各富集功能之间共有基因关系
heatplot(kk_filtered) #热图展示富集功能和基因的包含关系

#GO美化-------
library(GOplot)
library(ggplot2)
library(tidydr)

# load(file ="GO.Rdata" )
# load(file = "DEseq2.Rdata")

#重新选择通路-----
rownames(resdata) <- resdata$Row.names
resdata1 <- resdata[diff_gene_deseq2,]

{
  table(GO_filtered$ONTOLOGY)
  # CC <- subset(GO,ONTOLOGY == "CC")#没有
  # MF <- subset(GO,ONTOLOGY == "MF")#没有
  # BP <- subset(GO,ONTOLOGY == "BP")
  #table(GO$Count)
  BP <- subset(GO_filtered, Count > 5)
  #BP <- BP[c(3, 5, 12, 19, 26, 27, 30, 32, 51, 52, 59, 60, 62, 63, 64, 67, 69, 70, 71, 79, 83, 91, 92, 93, 101, 102, 107, 109, 121, 122, 123, 125, 137),]
  #BP1 <- BP[c(5,6,13,14,21,29,37,40,48,52,83,90,102),]
  #BP2 <- BP[c(52,83,90,102),]
}

rt <- data.frame(row.names = resdata$Row.names,
                 gene = resdata$Row.names,
                 logFC = resdata$log2FoldChange)
rt <- rt[diff_gene_deseq2,]

kegg1 <- strsplit(BP$geneID,"/") 
# a1 <- data.frame(Term = BP$Description[1],Genes = kegg1[1])
# a2 <- data.frame(Term = BP$Description[2],Genes = kegg1[2])
# a3 <- data.frame(Term = BP$Description[3],Genes = kegg1[3])
# a4 <- data.frame(Term = BP$Description[4],Genes = kegg1[4])
# a5 <- data.frame(Term = BP$Description[5],Genes = kegg1[5])
# a6 <- data.frame(Term = BP$Description[6],Genes = kegg1[6])
# a7 <- data.frame(Term = BP$Description[7],Genes = kegg1[7])
# a8 <- data.frame(Term = BP$Description[8],Genes = kegg1[8])
# a9 <- data.frame(Term = BP$Description[9],Genes = kegg1[9])
# a10 <- data.frame(Term = BP$Description[10],Genes = kegg1[10])
# a11 <- data.frame(Term = BP$Description[11],Genes = kegg1[11])
# a12 <- data.frame(Term = BP$Description[12],Genes = kegg1[12])
# 
# 
# colnames(a1)[2] <- "gene"
# colnames(a2)[2] <- "gene"
# colnames(a3)[2] <- "gene"
# colnames(a4)[2] <- "gene"
# colnames(a5)[2] <- "gene"
# colnames(a6)[2] <- "gene"
# colnames(a7)[2] <- "gene"
# colnames(a8)[2] <- "gene"
# colnames(a9)[2] <- "gene"
# colnames(a10)[2] <- "gene"
# colnames(a11)[2] <- "gene"
# colnames(a12)[2] <- "gene"
# 
# 
# all1 <- rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
# all1$value <- 1
all1 <- do.call(rbind, lapply(seq_along(kegg1), function(i) {
  data.frame(
    Term = BP$Description[i],
    gene = kegg1[[i]],
    value = 1
  )
}))

library(tidyr)
all2 <- spread(all1,key = "Term",value = "value",fill = 0)
rownames(all2) <- all2$gene
all2 <- all2[-1]

id <- intersect(rownames(all2),rt$gene)
rt <- rt[id,]
rt <- rt[order(rt$logFC,decreasing = T),]
all2 <- all2[rownames(rt),]
all2$logFC <- rt$logFC
#colnames(all2)[10] <- "logFC"

GOChord(all2,
        space = 0,
        gene.order = "logFC",
        gene.space = 0.25,
        gene.size = 0,
        border.size = 0.1,
        process.label = 7)
#美化
# GOChord(all2, space = 0.02, 
#         gene.order = 'logFC', 
#         gene.space = 0.25, gene.size = 0)
col <- c("#FF0000","#FF7F00","#FFFF00","#00FF00","#00FFFF","#8B00FF",
         "#FF0000","#FF7F00","#FFFF00","#00FF00","#00FFFF","#8B00FF")
GOChord(all2, 
        space = 0, #GO term 紧贴排列
        gene.order = "logFC", 
        gene.space = 0.25,#基因标签贴着 logFC 环
        gene.size = 0,
        lfc.col = c("#FF0000","#00FF00","blue"),
        process.label = 10,
        border.size = 0.1,
        ribbon.col = col)
#pdf12*12,图例用18*12提取用ai处理

#KEGG--------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	

load(file = "./4.DEG/DEseq2.Rdata")
input_gene <- diff_gene_deseq2
input_gene=unique(as.vector(input_gene))

# entrezIDs=mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
# entrezIDs=as.character(entrezIDs)
# gene=entrezIDs[entrezIDs!="NA"]
# gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
entrezIDs <- bitr(input_gene, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
#kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk_deg <- enrichKEGG(gene = entrezIDs$ENTREZID,
                     organism = 'hsa',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1)
#kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk_deg)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(input_gene[match(strsplit(x,"/")[[1]],as.character(entrezIDs))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue < 0.05 & KEGG$qvalue < 0.1),]

#保存
save(kk_deg,KEGG,file = "./4.2.癌症组差异基因的go、kegg分析/kk.Rdata")
#write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)
#KEGG_data <- read.table("KEGG.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

KEGG$Description[50:61]
a <- c("hsa04024","hsa04110","hsa04060","hsa03320","hsa04610",
       "hsa04657","hsa04512","hsa04514","hsa04115","hsa04022")
KEGG_1 <- KEGG[a,]

#柱状图
#pdf(file="KEGGbarplot.pdf", width=9, height=7)
# barplot(kk_deg, drop=TRUE, 
#         showCategory= 10,
#         label_format=130, 
#         color = col)
#dev.off()

#气泡图
#pdf(file="KEGGbubble.pdf", width = 9, height = 7)
# dotplot(KEGG_1, showCategory = 10,
#         orderBy="GeneRatio", 
#         label_format=130, color=colorSel)
#dev.off()

#美化------
ggplot(KEGG_1,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust)
           ,stat='identity')+
  scale_fill_gradient(low="#FFCC33",high="#CC6666")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x="Count",
       y="",
       title="KEGG enrichment")

# library(ggplot2)
# KEGG_1$Description <- factor(KEGG_1$Description, levels = rev(KEGG_1$Description))
# ggplot(KEGG_1[1:10,], aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
#   geom_point() +
#   scale_color_gradient(low="blue", high="red") +
#   theme_bw() +
#   labs(x="GeneRatio", y="", color="Adjusted P", size="Gene Count")
mytheme1 <- theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 14, hjust = 0.5, face ="bold"),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 11))
#x轴为基因数量
p3 <- ggplot(data = KEGG_1,
             aes(x = Count, y = reorder(Description, Count))) +
  geom_point(aes(size = Count, color = -log10(pvalue)), shape = 18) +
  scale_color_distiller(palette ="Spectral", direction = 1) +
  scale_size_continuous(range = c(4, 8)) +  # 反转大小
  labs(x ="Gene Number",
       y="",
       title="Dotplot of Enriched KEGG Pathways",
       size="gene number") +
  theme_bw() +
  mytheme1
p3

#x轴为基因占比
KEGG_1$GeneRatio <- sapply(KEGG_1$GeneRatio, function(x) {
  nums <- as.numeric(unlist(strsplit(x, "/")))
  nums[1] / nums[2]
})
KEGG_1$GeneRatio <- KEGG_1$GeneRatio *10
p3 <- ggplot(data = KEGG_1,
             aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +  
  geom_point(aes(size = Count, color = -log10(pvalue)), shape = 18) +
  scale_color_distiller(palette ="Spectral", direction = 1) +
  scale_size_continuous(range = c(4, 9)) +
  labs(x ="Gene Ratio",
       y="",
       title="Bubble chart of KEGG pathways for differential genes in cancer groups",
       size="gene number") +
  theme_bw() +
  mytheme1
p3
#pdf
ggsave("./4.2.癌症组差异基因的go、kegg分析/癌症组差异基因KEGG.pdf",
       plot = p3,
       width = 9,   # 宽更大
       height = 8,  # 高更大
       units = "in")
































