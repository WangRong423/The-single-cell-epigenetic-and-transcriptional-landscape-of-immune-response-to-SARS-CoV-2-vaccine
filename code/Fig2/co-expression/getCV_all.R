# 计算差异系数(coefficient of variation)

library(data.table)

# celltypes <- c('Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
#                'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
#                 'CD8+ Tem', 'Eryth/Platelet', 'MAIT')

print('Loading all.csv file.')
t1 <- proc.time()
data <- fread('/database/huhuajie/findCoexpressionGenes/finalResults/all.csv')
t2 <- proc.time()
t <- t2 - t1
print(paste0('Down! Time taken: ', t[3][[1]] / 60, ' mins'))

# for (type in celltypes){
type <- 'Eryth/Platelet CD14+ Mono new'
print(paste0('Start processing ', type))
t3 <- proc.time()

# 选取特定类型(type)的细胞，删除不需要的列
df <- data[celltype %chin% c('Eryth/Platelet','CD14+ Mono'), 
           !c("V1", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res.0.3", "seurat_clusters",
              "Dlabel", "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", 
              "nucleosome_percentile", "TSS.enrichment", "TSS.percentile", "n200", "n200_k33", "annotation", "celltype")]
df1 <- df[Participants %chin% c('P1', 'P3', 'P4')]

# 对于P2，剔除数据质量差的Day1、Day3、Day171、Day185四个时间点的数据
df2 <- subset(df, Participants == 'P2' & Timepoints %chin% c('Day0', 'Day6','Day14','Day30','Day31','Day33','Day36','Day44'))
df <- rbind(df1, df2)
df <- df[, !c('Participants')]

# 求该细胞类型所有基因不同时间点的平均表达量及表达量标准差
meandf <- df[, lapply(.SD, mean, na.rm = TRUE), by = Timepoints]
sddf <- df[, lapply(.SD, sd, na.rm = TRUE), by = Timepoints]

# 按时间点升序排列
for (i in 1:nrow(meandf)){
  if (meandf[, 1][i] == 'Day0'){
    meandf[, 1][i] = 0
  }else if (meandf[, 1][i] == 'Day1'){
    meandf[, 1][i] = 1
  }else if (meandf[, 1][i] == 'Day3'){
    meandf[, 1][i] = 3
  }else if (meandf[, 1][i] == 'Day6'){
    meandf[, 1][i] = 6
  }else if (meandf[, 1][i] == 'Day14'){
    meandf[, 1][i] = 14
  }else if (meandf[, 1][i] == 'Day30'){
    meandf[, 1][i] = 30
  }else if (meandf[, 1][i] == 'Day31'){
    meandf[, 1][i] = 31
  }else if (meandf[, 1][i] == 'Day33'){
    meandf[, 1][i] = 33
  }else if (meandf[, 1][i] == 'Day36'){
    meandf[, 1][i] = 36
  }else if (meandf[, 1][i] == 'Day44'){
    meandf[, 1][i] = 44
  }
}
meandf$Timepoints <- as.numeric(meandf$Timepoints)
meandf <- meandf[order(meandf$Timepoints)]

for (i in 1:nrow(sddf)){
  if (sddf[, 1][i] == 'Day0'){
    sddf[, 1][i] = 0
  }else if (sddf[, 1][i] == 'Day1'){
    sddf[, 1][i] = 1
  }else if (sddf[, 1][i] == 'Day3'){
    sddf[, 1][i] = 3
  }else if (sddf[, 1][i] == 'Day6'){
    sddf[, 1][i] = 6
  }else if (sddf[, 1][i] == 'Day14'){
    sddf[, 1][i] = 14
  }else if (sddf[, 1][i] == 'Day30'){
    sddf[, 1][i] = 30
  }else if (sddf[, 1][i] == 'Day31'){
    sddf[, 1][i] = 31
  }else if (sddf[, 1][i] == 'Day33'){
    sddf[, 1][i] = 33
  }else if (sddf[, 1][i] == 'Day36'){
    sddf[, 1][i] = 36
  }else if (sddf[, 1][i] == 'Day44'){
    sddf[, 1][i] = 44
  }
}
sddf$Timepoints <- as.numeric(sddf$Timepoints)
sddf <- sddf[order(sddf$Timepoints)]

# 计算CV，CV保存原始差异系数值，CV_ED存放CV变化的等距斜率，CV_NED存放CV变化的不等距斜率
x <- meandf[, .(Timepoints)]
x <- as.numeric(unlist(x))
num <- 1
genes <- colnames(meandf)
CV <- data.table()
CV_ED <- data.table()
CV_NED <- data.table()
for (i in 2:length(genes)){
  mean <- meandf[[i]]
  mean <- as.numeric(unlist(mean))
  sd <- sddf[[i]]
  sd <- as.numeric(unlist(sd))
  for (j in 1:(length(mean)-1)){
    if (mean[j+1] == 0 | mean[j] == 0){
      tmp1 <- 0
    }
    else{
      tmp1 <- (sd[j+1]/mean[j+1] - sd[j]/mean[j]) / (x[j+1] - x[j])
    }
    if (j == 1){
      slope1 <- tmp1
    }else{
      slope1 <- c(slope1, tmp1)
    }
    if (mean[j+1] == 0 | mean[j] == 0){
      tmp2 <- 0
    }
    else{
      tmp2 <- (sd[j+1]/mean[j+1] - sd[j]/mean[j])
    }
    if (j == 1){
      slope2 <- tmp2
    }else{
      slope2 <- c(slope2, tmp2)
    }
  }
  CV_NED[, genes[i]:= slope1]
  CV_ED[, genes[i]:= slope2]
  
  for (k in 1:length(mean)){
    if (mean[k] == 0){
      tmp3 <- 0
    }
    else{
      tmp3 <- sd[k] / mean[k]
    }
    if (k == 1){
      slope3 <- tmp3
    }else{
      slope3 <- c(slope3, tmp3)
    }
  }
  CV[, genes[i]:= slope3]
  
  if (num %% 2000 == 0 | num == length(genes)){
    print(num)
  }
  num <- num + 1
}
type <- gsub(' ','_',type)
type <- gsub('/', '_', type)
path1 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getCV/','CV_',type,'_ED.csv')
path2 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getCV/','CV_',type,'_NED.csv')
path3 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getCV/','CV_',type,'.csv')
fwrite(CV_ED, path1)
fwrite(CV_NED, path2)
fwrite(CV, path3)
t4 <- proc.time()
t <- t4 - t3
print(paste0(type, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
rm(df, meandf, sddf, CV_ED, CV_NED, CV, path1, path2, path3)
gc()
# }


t5 <- proc.time()
t <- t5 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
