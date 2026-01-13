#高低风险与ICD和ICP调节之间的关联
rm(list = ls())
options(stringsAsFactors = F)

setwd("/home/datahup/syj/GAC/")
load(file = "./1.download_bulk/dat.output.Rdata")
load(file = "./1.download_bulk/phe.output.Rdata")

load(file = "./12.2RSF11.11/训练集.Rdata")
load(file = "./12.2RSF11.11/RSF_genes.Rdata")

# ICD <- c("CALR","CXCL10","FPR1","HGF","IFNAR1","IFNAR2","IFNB1","IFNE","IFNK","IFNW1","MET",
#          "P2RX7","TLR4","LRP1","EIF2A","EIF2AL3","EIF2AK2","EIF2AK4","EIF2AK1","HMGB1","ANXA1","PANX1",
#          "P2RY2","IFNA1","IFNA2","TLR3")
# table(rownames(dat) %in% ICD)

ICP <-  c("ADORA2A","BTLA","BTNL2","CD160","CD200","CD200R1","CD244","CD27","CD274","CD276",
          "CD28","CD40","CD40LG","CD44","CD48","CD70","CD80","CD86","CTLA4","HAVCR2","HHLA2",
          "ICOS","ICOSLG","IDO1","IDO2","KIR3DL1","LAG3","LAIR1","LGALS9","NRP1","PDCD1","PDCD1LG2",
          "TIGIT","TMIGD2","TNFRSF14","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14","TNFSF15",
          "TNFSF18","TNFSF4","TNFSF9","VSIR","VTCN1")
table(rownames(dat) %in% ICP)

phe53625 <- phe
dat53625 <- dat

# rownames(phe53625) <- phe53625$Id
# phe53625 <- phe53625[datset1$ID,]
# dat53625 <- as.data.frame(t(dat53625[c(RSF_genes,ICD),]))
# colnames(dat53625)
# dat53625 <- dat53625[, !grepl("^NA", colnames(dat53625))]
# dat53625 <- dat53625[datset1$ID,]
 
rownames(phe53625) <- phe53625$Id
phe53625 <- phe53625[datset1$ID,]
dat53625 <- as.data.frame(t(dat53625[c(RSF_genes,ICP),]))
dat53625 <- dat53625[phe53625$Id,]
colnames(dat53625)
dat53625 <- dat53625[, !grepl("^NA", colnames(dat53625))]

exprSet <- dat53625
# exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
#   exprSet[,rownames(gene)[2]]*gene[2,2]+
#   exprSet[,rownames(gene)[3]]*gene[3,2]+
#   exprSet[,rownames(gene)[4]]*gene[4,2]+
#   exprSet[,rownames(gene)[5]]*gene[5,2]
exprSet$RiskScore <- datset1$risk.score
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
#exprSet <- exprSet[, !(colnames(exprSet) %in% gene$sig_gene_multi_cox)]

ImmSet <- exprSet
colnames(ImmSet)
ImmSet2 <- ImmSet[,8:52]

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == rownames(ImmSet) ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)

library(ggpubr)
#ICP:18*12
ggplot(ImmSet3, aes(x = reorder(`Immune Cell Type`, -`Estimating Score`), 
                    y = `Estimating Score`, fill = RiskGroup)) +
  # 小提琴图层
  #geom_violin(position = position_dodge(0.9), alpha = 0.7, width = 1.25, trim = TRUE) +
  # 箱线图图层
  geom_boxplot(width = 0.5, show.legend = T, 
               position = position_dodge(0.9), 
               color = 'black', alpha = 0.8, 
               outlier.shape = 21, outlier.colour = "black") +
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = " ", y = 'Gene Expression') +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4,
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")) +  # 增加 x 轴标题的上边距
  # 颜色设置
  scale_fill_manual(values = c("#FF0000",'#00CCFF'), name = 'Risk Group') +
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup, label = ..p.signif..), method = "wilcox.test") +
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5),
        panel.grid.minor = element_blank()) 































































