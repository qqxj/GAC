library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(cowplot)

rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/GAC/')

load(file = "./1.download_bulk/dat.output.Rdata")
load(file = "./1.download_bulk/phe.output.Rdata")
#load(file = "./12.2RSF11.11/训练集.Rdata")
load(file = "./12.2RSF11.11/外部验证集.Rdata")

colnames(datset3)
Gene <- colnames(datset3)[4:10]
datset3$Group <- ifelse(datset3$risk.score > median(datset3$risk.score), "High", "Low")
rownames(datset3) <- datset3$ID 

expr1 <- datset3
# 检查分组情况
#rownames(dat) <- dat$Id
# a <- dat[dat$type == "Normal",]
# table(substr(rownames(a), 14,15))
# expr1 <- as.data.frame(t(NSCLCcount[Gene,rownames(sv1)]))
# # expr1 <- as.data.frame(t(NSCLCcount[Gene,]))
# table(substr(rownames(expr1), 14,15))
# expr1$Group <- ifelse(substr(rownames(expr1), 14,15) == 11 ,"Normal","Tumor")
# table(expr1$Group)
# Gene

# expr1 <- as.data.frame(t(dat[Gene,datset1$ID]))
# expr1 <- expr1[datset1,]
# expr1$Group <- datset1$risk.score
# expr1$Group <- ifelse(expr1$Group == "cancer","Tumor","Normal")
# rownames(expr1)
#expr1 <- expr1[grep("^GSM129", rownames(expr1)), ]



# 箱图------
{#TYMP
  p1 <- ggplot(data=expr1)+ 
    geom_boxplot(mapping=aes(x=Group,y=TYMP,colour = Group ), 
                 alpha = 0.5,
                 size=1.2,
                 width = 0.6)+ 
    geom_jitter(mapping=aes(x=Group,y=TYMP,colour = Group), #散点
                alpha = 0.3,size=2)+
    scale_color_manual(limits=c("High","Low"), 
                       values=c("red","blue"))+ #颜色
    geom_signif(mapping=aes(x=Group,y=TYMP), # 不同组别的显著性
                comparisons = list(c("High","Low")), # 哪些组进行比较
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
    labs(title='TYMP',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
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
  p1
}

# 分半小提琴图------
library(devtools)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)

dat53625 <- dat[Gene,datset3$ID]
{
  a <- Gene
  gene1 <- as.data.frame(t(dat53625[a[1],]))
  names(gene1) <- "Expression"
  gene1$Gene <- rep(a[1],length(rownames(gene1)))
  #gene1$Type <- phe53625$group
  gene1$Type <- datset3$Group
  
  gene2 <- as.data.frame(t(dat53625[a[2],]))
  names(gene2) <- "Expression"
  gene2$Gene <- rep(a[2],length(rownames(gene2)))
  #gene2$Type <- ifelse(substr(rownames(gene2), 14,15) == 11 ,"Normal","Tumor")
  #gene2$Type <- phe53625$group
  gene2$Type <- datset3$Group
  
}

a <- Gene
for (i in 1:7) {
  # 创建临时数据框
  gene_temp <- as.data.frame(t(dat53625[a[i], ]))
  names(gene_temp) <- "Expression"
  gene_temp$Gene <- rep(a[i], nrow(gene_temp))  # 修正为 a[i]
  gene_temp$Type <- datset3$Group
  
  # 使用 assign 创建动态变量名
  assign(paste0("gene", i), gene_temp)
}

df <- rbind(gene1,gene2,gene3,gene4,gene5,gene6,gene7)#这一步要添加基因数量
df$ID <- rownames(df)
rownames(df) <- NULL
#可视化
{
  #mypalette = pal_ucscgb()(20)
  ggplot(df,aes(x = Gene, y = Expression, fill = Type)) +
    # split violin
    geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 0,color = 'black',hjust = 0.5),
          legend.position = 'top') + # 图例的位置
    # scale_fill_brewer(palette = 'Set1') +
    # scale_fill_jco(name = '') +
    labs(y="Expression(log2(Count+1))",x=NULL)+ # 添加标题，x轴，y轴内容
    # scale_color_lancet(name = '')+
    #scale_fill_manual(values = mypalette[c(6,1)]) + 
    scale_fill_manual(values = c("#32067A","#B8E54E")) + 
    # scale_fill_manual(values = c(pal[1],pal[2]),name = '') +
    # ylim(1,8) + #限定y轴范围
    # 添加显著性标记
    stat_compare_means(aes(group=Type),
                       # symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "NS")),
                       label = "p.signif",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "NS")),
                       label.y = 13.5,size = 5)
  #保存为pdf12*8
}
#美化癌症和癌旁可视化--------
ggplot(df, aes(x = Gene, y = Expression, fill = Type)) +
  # 箱线图
  geom_boxplot(alpha = 0.5, color = "black", width = 0.6, outlier.shape = NA, 
               position = position_dodge(0.8)) +
  # 散点图（jittered散点）
  geom_jitter(aes(color = Type), width = 0.2, height = 0, size = 2, alpha = 0.7) +
  # mean point
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.8), size = 3, color = "red") +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = 0.2, 
               size = 0.3, position = position_dodge(0.8)) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 0, color = 'black', hjust = 0.5),
        legend.position = 'top') + # 图例的位置
  labs(y = "Expression(log2(Count+1))", x = NULL) + # 添加标题，x轴，y轴内容
  scale_fill_manual(values = c("#FF0000FF", "#999999FF")) + 
  # 添加显著性标记
  stat_compare_means(aes(group = Type),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "NS")),
                     label.y = 22, size = 5)



#km-----
load(file = "./12.2RSF11.11/外部验证集.Rdata")
colnames(datset3)
#range(datset1$OS.time)
#sv2 <- subset(sv2,OS.time > "")
dat3 <- datset3[,c("ID","OS","OS.time",Gene)]
dat3$OS.time=dat3$OS.time/365
dat4 = data.frame(OS = dat3$OS,
                  OS.time = dat3$OS.time,
                  row.names = dat3$ID)

# 最佳生存节点
library(survival)
library(survminer)
library(ggsci)
#load(file = "../100种机器学习算法-8.1/geo数据格式.Rdata")
#dat3 <- datset1
# TYMP PDF8*8 --------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "TYMP") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  TYMP = ifelse(dat3$TYMP<res.cut$cutpoint[,1],"Low","High")
  TYMP = factor(TYMP)
  sfit <- survfit(Surv(OS.time, OS)~TYMP, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="TYMP",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for TYMP", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# IFNG --------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "IFNG") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  IFNG = ifelse(dat3$IFNG<res.cut$cutpoint[,1],"Low","High")
  IFNG = factor(IFNG)
  sfit <- survfit(Surv(OS.time, OS)~IFNG, data=dat4)
  p2 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="IFNG",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for IFNG", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p2
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# ITGAX ----------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "ITGAX") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  ITGAX = ifelse(dat3$ITGAX<res.cut$cutpoint[,1],"Low","High")
  ITGAX = factor(ITGAX)
  sfit <- survfit(Surv(OS.time, OS)~ITGAX, data=dat4)
  p3 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="ITGAX",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for ITGAX", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p3
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# GBP5 -------------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "GBP5") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  GBP5 = ifelse(dat3$GBP5<res.cut$cutpoint[,1],"Low","High")
  GBP5 = factor(GBP5)
  sfit <- survfit(Surv(OS.time, OS)~GBP5, data=dat4)
  p4 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="GBP5",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for GBP5", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p4
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# GBP4 ---------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "GBP4") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  GBP4 = ifelse(dat3$GBP4<res.cut$cutpoint[,1],"Low","High")
  GBP4 = factor(GBP4)
  sfit <- survfit(Surv(OS.time, OS)~GBP4, data=dat4)
  p5 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="GBP4",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for GBP4", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p5
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# STAT1 ---------
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "STAT1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  STAT1 = ifelse(dat3$STAT1<res.cut$cutpoint[,1],"Low","High")
  STAT1 = factor(STAT1)
  sfit <- survfit(Surv(OS.time, OS)~STAT1, data=dat4)
  #fit_test <- survdiff(Surv(OS.time, OS) ~ STAT1, data = dat4)
  #p_raw <- 1 - pchisq(fit_test$chisq, df = 1)
  # p_fmt <- floor(p_raw * 100) / 100
  # p_fmt <- format(p_fmt, nsmall = 2)
  # p6 = ggsurvplot(
  #   sfit,
  #   palette = c("#32067A", "#B8E54E"),
  #   conf.int = TRUE, conf.int.style = 'step',
  #   pval = paste0("p = ", p_fmt),     # ← 手动写入 p 值
  #   pval.method = FALSE,              # 不用自动
  #   legend = c(0.85, 0.85),
  #   legend.title = "STAT1",
  #   legend.labs = c("High", "Low"),
  #   title = "Survival Curve for STAT1",
  #   xlab = "Time(Years)",
  #   surv.median.line = "hv",
  #   ggtheme = theme_bw(base_size = 12))
  p6 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="STAT1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for STAT1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p6
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# CD84 ---------
{
  # dat3 <- datset3[,c("ID","OS","OS.time",Gene)]
  # dat3$OS.time=dat3$OS.time/365
  # dat4 = data.frame(OS = dat3$OS,
  #                   OS.time = dat3$OS.time,
  #                   row.names = dat3$ID)
  res.cut <- surv_cutpoint(dat3, #数据集
                          time = "OS.time", #生存时间
                          event = "OS", #生存状态
                          variables = "CD84") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  CD84 = ifelse(dat3$CD84<res.cut$cutpoint[,1],"Low","High")
  CD84 = factor(CD84)
  sfit <- survfit(Surv(OS.time, OS)~CD84, data=dat4)
  p7 = ggsurvplot(sfit,
                  palette = c("#32067A","#B8E54E"), 
                  #palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="CD84",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for CD84", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p7
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}

#合并
library(patchwork)
final_plot1 <- (p1$plot | p2$plot | p3$plot | p4$plot) 
final_plot2 <-  (p5$plot | p6$plot | p7$plot)
final_plot1#pdf32*8
final_plot2#pdf24*8

