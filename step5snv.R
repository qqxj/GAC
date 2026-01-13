
rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd('/home/datahup/syj/GAC/')
library(maftools)

#snv突变数据
gac.snv = read.delim(gzfile('/opt/disease/TCGA/STAD/TCGA-STAD.muse_snv.tsv.gz'))
#cnv突变数据
gac.cnv = read.delim(gzfile('/opt/disease/TCGA/STAD/TCGA-STAD.cnv.tsv.gz'))#没有参考意义



#临床数据
load(file = "./1.download_bulk/训练集数据.Rdata")
gac.snv$Sample_ID <- gsub("-", ".", gac.snv$Sample_ID)
table(train_phe$Id %in% gac.snv$Sample_ID)
table(rownames(train_dat) %in% gac.snv$gene)

#cnv数据——————这里不加载cnv数据是想后面韦恩图用上
 # setwd("/home/datahup/syj/ESCC/SNP-TCGA-6.1/")
 # load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
 # gene = unique(Statemarker$gene)
 # table(gene %in% lusc.snv$gene)
 # snv <- lusc.snv[lusc.snv$gene %in% gene, ]
 # table(duplicated(snv$gene))
snv <- gac.snv
table(duplicated(snv$gene))


#数据清洗(用不了、用NSCC文件里的就行)-----
{
  dat <- train_phe
  tmp <- snv
  
  #dat$Id2 = gsub('[.]','-',dat$Id)
  dat$Id2 <- dat$Id
  rownames(dat) = dat$Id2
  length(unique(tmp$gene))
  length(unique(tmp$Sample_ID))
  phe = subset(dat,dat$Id2 %in% unique(tmp$Sample_ID))
  tmp <- subset(tmp,unique(tmp$Sample_ID) %in% dat$Id2)
  head(tmp)   
  colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                    "Chromosome", "Start_Position", 
                    "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                    "HGVSp_Short" , 'effect' ,"Consequence",
                    "vaf" ) 
  tmp$Entrez_Gene_Id =1
  tmp$Center ='ucsc'
  tmp$NCBI_Build ='GRCh38'
  tmp$Strand ='+'
  tmp$Variant_Classification = tmp$effect
  tail(sort(table(tmp$Variant_Classification )))
  tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
  tmp$Variant_Type = ifelse(
    tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
    'SNP','INDEL'
  )
  table(tmp$Variant_Type )
  tmp2 = subset(tmp, tmp$Tumor_Sample_Barcode  %in% dat$Id2 ) # 添加临床信息
  phe$age = ifelse(phe$Age <= 60,'<=60','>60')
  tmp2$age = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$age, 'Unknow')
  table(tmp2$age)
  tmp2$Gender = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$Gender, 'Unknow')
  table(tmp2$Gender)
  tmp2$`Survive State` = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$OS, 'Unknow')
  table(tmp2$`Survive State`)
  neu.state = subset(tmp2,tmp2$Hugo_Symbol %in% tmp$Hugo_Symbol)
  
  tcga.neu.state = read.maf(maf = tmp,
                            vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))
  
  save(neu.state,tcga.neu.state,tmp2,file = "./5.snv/snv.Rdata")
  
}
load(file = "./5.snv/snv.Rdata")
# 添加临床信息-----
{
  a = tcga.neu.state@clinical.data
  a$age = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$age,'Unkown')
  a$age[is.na(a$age)] = 'Unkown'
  table(a$age)
  a$Gender = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$Gender,'Unkown')
  table(a$Gender)
  a$`Survive State` = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$`Survive State`,'Unkown')
  table(a$`Survive State`)
  tcga.neu.state@clinical.data = a
}

####瀑布图-----
{
  oncoplot(maf = tcga.neu.state,top = 10) # 高频突变的前10个基因
  plot <- oncoplot(maf = tcga.neu.state, top = 10)
  # ggsave(file = "./5.snv/前十个基因突变——瀑布图.png", 
  #        plot, width = 10, height = 6, dpi = 300)
  
  getFields(tcga.neu.state)
  getClinicalData(tcga.neu.state) #查看临床信息
  getSampleSummary(tcga.neu.state)#查看每个样品发生突变的情况，此处就可以计算tumor mutation load,TML=Missense_Mutation/外显子数。
  
  plotmafSummary(maf = tcga.neu.state, rmOutlier = TRUE, 
                 addStat = 'median', dashboard = TRUE,
                 titvRaw=FALSE)#绘制整体的突变情况
  
  library(ggsci)
  col = pal_lancet()(9)#6种
  [1] "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF"
  [7] "#AD002AFF" "#ADB6B6FF" "#1B1919FF"  
  col <- c("#00468BFF", "#ED0000FF",  "#0099B4FF" ,
           "#925E9FFF", "#FDAF91FF","#ADB6B6FF", "#42B540FF")
  table(tcga.neu.state@data[["effect"]])
  names(col) = c('frameshift_variant','missense_variant','synonymous_variant',
                 'intron_variant', '3_prime_UTR_variant' , 'stop_gained')
  
  tcga.neu.state
  getClinicalData(x = tcga.neu.state)
  
  # 修改列名
  colnames(tcga.neu.state@clinical.data)[colnames(tcga.neu.state@clinical.data) == "Survive State"] <- "Survive_State"
  anocol <- list(
    age = c("<=60" =  "#42B540FF", ">60" = "#ED0000FF", "Unkown" = "#BDBDBD"),
    Gender = c("male" = "#0099B4FF", "female" = "#ED0000FF", "Unkown" = "#ADB6B6FF"),
    `Survive_State` = c("Alive" = "#42B540FF", "Dead" = "#ED0000FF", "Unkown" = "#ADB6B6FF")
  )
  
  pdf(file  = './5.snv/neu_oncoplot.pdf',width = 12,height = 8)
  oncoplot(maf = tcga.neu.state, top = 30, colors = col,
           draw_titv = T,
           clinicalFeatures = c('age','Gender','Survive_State'),
           annotationColor = anocol,
           sortByAnnotation = TRUE )
  dev.off()
  
  # 转换 颠倒 
  pdf(file  = './5.snv/neu_titv.pdf',width = 12,height = 8)
  titv(tcga.neu.state, useSyn = FALSE, plot = TRUE, file = NULL)
  dev.off()
}

####提取snv突变数据####
{
  snv_data <- tcga.neu.state@data[tcga.neu.state@data$Variant_Type == "SNP", ]
  selected_columns <- c("Tumor_Sample_Barcode","Hugo_Symbol", "Chromosome", "Start_Position", 
                        "End_Position", "Reference_Allele", 
                        "Tumor_Seq_Allele2", "Variant_Classification")
  
  #只提取出错义突变基因（这个是最后筛出来突变基因，先用这个）
  missense_data <- snv_data[grepl("missense_variant", snv_data$Variant_Classification), ]
  length(unique(missense_data$Hugo_Symbol))
  #[1] 15768
  mut_gene1 <- unique(missense_data$Hugo_Symbol)
  save(missense_data,mut_gene1,file = "./5.snv/snv_data_RAW.Rdata")
  
  #很奇怪得到的结果是不一样的（这个是提前筛出来突变基因）
  snv_data <- snv_data[, selected_columns, with = FALSE]
  length(unique(snv_data$Hugo_Symbol))
  #[1] 15920
  mut_gene <- unique(snv_data$Hugo_Symbol)
  save(mut_gene,snv_data,file = "./5.snv/snv_data.Rdata")
}





