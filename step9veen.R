#7.9

rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/GAC/")

#过表达基因(上调基因)
load(file = "./4.DEG/DEseq2.Rdata")
table(duplicated(Up_gene_deseq2))
up <- unique(Up_gene_deseq2)

#snv突变基因
#lusc.snv = read.delim(gzfile('/home/datahup/syj/NSCC/step8_1TCGA_snv/TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "./5.snv/snv_data.Rdata")
table(duplicated(snv_data$Hugo_Symbol))
table(duplicated(mut_gene)) #这个
snv <- unique(mut_gene)

#单细胞拟时序相关基因集
#load(file = './8.拟时序/Statemarker.Rdata')
#load(file = './8.拟时序/deg.cluster10.Rdata')
deg_cluster <- read_excel("./8.1拟时序基因/deg_cluster.xlsx")
deg.cluster <- deg_cluster
table(duplicated(deg.cluster$gene)) 
cnv <- unique(deg.cluster$gene)
# Statemarker1 <- subset(Statemarker,cluster == "Epithelial cell")
# Statemarker2 <- subset(Statemarker,cluster != "Epithelial cell")
# table(Statemarker1$gene == Statemarker2$gene)
#intersect_genes <- intersect(Statemarker$gene, deg.cluster$gene)


#癌症免疫相关基因集
#imm_genes <- read.table("./9.韦恩图/GeneList副本.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library(readxl)
imm_genes <- read_excel("./9.韦恩图/GeneList副本.xlsx")
table(duplicated(imm_genes$Symbol)) 

#Mut_genes <- unique(Statemarker$gene)[unique(Statemarker$gene) %in% unique(lusc.snv$gene)]

#交集
{
  # 1. 上调基因
  up_genes <- Up_gene_deseq2#unique()
  # 2. SNV突变基因
  snv_genes <- unique(mut_gene)#
  # 3. 拟时序相关基因
  pseudo_genes <- unique(deg.cluster$gene)#
  # 4. 免疫相关基因
  immune_genes <- imm_genes$Symbol#unique()
  
  intersect_all <- Reduce(intersect, list(up_genes, snv_genes, pseudo_genes, immune_genes))#,
  length(intersect_all)
  a <- unique(intersect_all)  # 可查看具体有哪些公共基因
  
  save(intersect_all,file = "./9.韦恩图/交集基因1028.Rdata")
  # 存储在列表中
  gene_sets <- list(
    Up = up_genes,
    SNV = snv_genes,
    CNV = pseudo_genes,
    Immune = immune_genes
  )
  


library(ggvenn)
ggvenn(list(Up_gene = up_genes,
            SNV = snv_genes,
            CNV = pseudo_genes,
            Cor_Immune = immune_genes),
  fill_color = c("red", "blue", "green", "purple"),
  text_size = 5)
#pdf8*10
}


























































