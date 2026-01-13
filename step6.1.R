#设置工作路径
rm(list = ls())
options(stringsAsFactors = F)
gc()

#加载这个数据
setwd("/home/datahup/syj/GAC/")
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

load(file = "/home/datahup/syj/GAC/7.copykat/copykat_Epi.Rdata")
NSCC_list <- scRNA
table(Idents(NSCC_list))
Idents(NSCC_list) <- NSCC_list@meta.data[["copykat.pred"]]
table(Idents(NSCC_list))
# Epithelial cell     Cancer cell 
# 20398           13549 

# 获取当前的Idents
current_idents <- Idents(NSCC_list)
# 创建一个新的标识向量
new_idents <- ifelse(current_idents == "Epithelial cell", "Epithelial cell", "CNV cancer cell")
# 将新的标识赋值给Seurat对象
Idents(NSCC_list) <- new_idents
# 检查结果
table(Idents(NSCC_list))

#NSCLC.Integrate@active.ident = NSCLC.Integrate@active.ident
genes_to_check = c("S100A4","ERBB2","CEACAM5",#"EPCAM",
                   "GKN1")#"MUC5AC","TFF1","PGA3",
#普通的热图或者气泡图(可视化)
DotPlot(NSCC_list,features = unique(genes_to_check)) + RotatedAxis()
#保存pdf8*6

#两类基因的对比箱线图，显示p值
# 提取基因表达数据
expression_data <- FetchData(NSCC_list, vars = c("S100A4", "ERBB2","CEACAM5", "GKN1","ident"))
head(expression_data)
table(expression_data$ident)
expression_data$group <- expression_data$ident

# 加载必要的包
library(ggplot2)
library(ggpubr)

# 绘制 S100A4 的箱线图
p1 <- ggplot(expression_data, aes(x = group, y = GKN1, fill = group)) +
  geom_boxplot() +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",  # 显示 p 值
    comparisons = list(c("Epithelial cell", "CNV cancer cell")),  # 指定比较的组
    label.y = max(expression_data$GKN1) * 1.1  # 将 p 值标签放在箱线图上方
  ) +
  theme_classic() +
  labs(title = "GKN1 Expression", x = "Cell Type", y = "Expression Level")
p1













































