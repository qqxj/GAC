
# 28种免疫细胞的浸润情况
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/14.免疫微环境/')

# 获取28种免疫细胞基因集 ---------
library(dplyr)
library(tidyverse)
{
  geneSet <- read.csv("./CellReports.txt",header = F,sep = "\t",) # 用EXCEL打开删除NA列
  class(geneSet)
  geneSet <- geneSet %>%
    column_to_rownames("V1")%>%t()
  a <- geneSet
  a <- a[1:nrow(a),]
  set <- colnames(a)
  l <- list()
  for (i in set) {
    x <-  as.character(a[,i])
    x <- x[nchar(x)!=0]
    x <-  as.character(x)
    l[[i]] <-x
  }
  #save(l,file = "./CellReports.Rdata")
}
load(file = "./CellReports.Rdata")
#准备表达矩阵、临床信息-----
load(file = "../1.download_bulk/dat.output.Rdata")
load(file = "../1.download_bulk/phe.output.Rdata")

df <- scale(dat)
df <- as.data.frame(df)

# 运行ssGSEA--------
library(GSVA)
library(limma)
{
  exprSet <- as.matrix(dat)
  ssgsea<- gsva(exprSet, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  #ssgsea<- gsva(exprSet, l,method='ssgsea',mx.diff=F,verbose = F)
  
  Immune2 = as.data.frame(t(ssgsea))
  Immune2$id <- rownames(Immune2)
  # save(Immune2,file = './FPKM_ssGSEA.Rdata')
}
load(file = './FPKM_ssGSEA.Rdata')

# 加载风险评分指数 -------------
#load(file = "/home/datahup/syj/GAC/12.1RSF+STEPCOX/风险评分数据.Rdata")
load(file = "../12.1RSF+STEPCOX/训练集.Rdata")
load(file = "../12.1RSF+STEPCOX/RSF_genes.Rdata")
#df <- read.csv(file = 'df_scale.csv')
#rownames(df) <- df$X
#df <- df[,-1]
df[1:6,1:6]
exprSet <- as.data.frame(t(df[RSF_genes,]))
# dim(exprSet)
# exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
#   exprSet[,rownames(gene)[2]]*gene[2,2]+
#   exprSet[,rownames(gene)[3]]*gene[3,2]+
#   exprSet[,rownames(gene)[4]]*gene[4,2]+
#   exprSet[,rownames(gene)[5]]*gene[5,2]
table(rownames(exprSet) %in% datset1$ID)
exprSet <- exprSet[datset1$ID,]
all(rownames(exprSet) == datset1$ID)#检查顺序和信息完全一致
exprSet$RiskScore <- datset1$risk.score

exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
# table(phe$Type)
# phe <- subset(phe,phe$Type == 'Tumor')
# exprSet <- exprSet[phe$Id,]

ImmSet <- merge(Immune2,exprSet,by = "id")
rownames(ImmSet) <- ImmSet$id
colnames(ImmSet)
ImmSet2 <- ImmSet[,2:29]
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
dim(ImmSet2)
dim(ImmSet3)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == ImmSet$id ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)
#save(ImmSet,ImmSet3,Immune2,file = './ssGSEA_output.Rdata')

# 小提琴图----------
#备注：主要负责抗原呈递细胞（APCs）：
#Activated dendritic cell、Plasmacytoid dendritic cell、Macrophage、
#B cell(包括 Activated、Immature、Memory)
library(ggplot2)
library(tidyverse)
library(ggpubr)
colnames(ImmSet3)
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  # # 小提琴图层
  # geom_violin(position = position_dodge(0.9),alpha = 1.2,
  #             width = 1.8,trim = T,
  #             color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  # outlier.shape = 21,
  # outlier.shape = NA
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  # 颜色设置
  # scale_fill_manual(values = c('Low'='#398AB9','High'='red'),
  #                   name = '') +
  scale_fill_manual(values = c("#FF0000",'#00CCFF'))+ 
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),#..p.format..   ..p.signif..
                     method = "wilcox.test")#kruskal.test  p.signif wilcox.test
# #添加显著性标记
# stat_compare_means(aes(group=RiskGroup),
#                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
# 
#                                     symbols = c("***", "**", "*", "NS")),label = "p.signif",
#                    label.y = 1.7,size = 4) + ylim(0.3,1.7)
ggsave(filename = './11.14结果/ssGSEA.pdf',width = 12,height = 7)















































































