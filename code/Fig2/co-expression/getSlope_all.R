# 计算斜率数据，包括跨一个点的斜率

library(data.table)

# celltypes <- c('Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
#                  'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
#                  'CD8+ Tem', 'Eryth/Platelet', 'MAIT')

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
df <- data[celltype %chin% c('Eryth/Platelet', 'CD14+ Mono'), 
                 !c("V1", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res.0.3", "seurat_clusters",
                    "Dlabel", "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", 
                    "nucleosome_percentile", "TSS.enrichment", "TSS.percentile", "n200", "n200_k33", "annotation", "celltype")]
df1 <- df[Participants %chin% c('P1', 'P3', 'P4')]

# 对于P2，剔除数据质量差的Day1、Day3、Day171、Day185四个时间点的数据
df2 <- subset(df, Participants == 'P2' & Timepoints %chin% c('Day0', 'Day6','Day14','Day30','Day31','Day33','Day36','Day44'))
df <- rbind(df1, df2)
df <- df[, !c('Participants')]

# 求该细胞类型所有基因不同时间点的平均表达量
df <- df[, lapply(.SD, mean, na.rm = TRUE), by = Timepoints]

# 按时间点升序排列
for (i in 1:nrow(df)){
  if (df[, 1][i] == 'Day0'){
    df[, 1][i] = 0
  }else if (df[, 1][i] == 'Day1'){
    df[, 1][i] = 1
  }else if (df[, 1][i] == 'Day3'){
    df[, 1][i] = 3
  }else if (df[, 1][i] == 'Day6'){
    df[, 1][i] = 6
  }else if (df[, 1][i] == 'Day14'){
    df[, 1][i] = 14
  }else if (df[, 1][i] == 'Day30'){
    df[, 1][i] = 30
  }else if (df[, 1][i] == 'Day31'){
    df[, 1][i] = 31
  }else if (df[, 1][i] == 'Day33'){
    df[, 1][i] = 33
  }else if (df[, 1][i] == 'Day36'){
    df[, 1][i] = 36
  }else if (df[, 1][i] == 'Day44'){
    df[, 1][i] = 44
  }
}
df$Timepoints <- as.numeric(df$Timepoints)
df <- df[order(df$Timepoints),]

# 保存该细胞类型所有基因不同时间点的平均表达量为csv文件
type <- gsub(' ','_',type)
type <- gsub('/', '_', type)
path <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/gex_',type,'.csv')
fwrite(df, path)

# 计算斜率，slope_NED存放不等距斜率，slope_ED存放等距斜率
x <- df$Timepoints
num <- 1
genes <- colnames(df)
slope_NED <- data.table()
slope_ED <- data.table()
for (i in 2:length(genes)){
  y <- df[[i]]
  
  # 非跨点斜率，包含翻转
  for (j in 1:(length(y)-1)){
    tmp1 <- (y[j+1] - y[j]) / (x[j+1] - x[j])
    if (j == 1){
      slope1 <- tmp1
      slope1_flip <- -tmp1
    }else{
      slope1 <- c(slope1, tmp1)
      slope1_flip <- c(slope1_flip, -tmp1)
    }
    
    tmp2 <- (y[j+1] - y[j])
    if (j == 1){
      slope2 <- tmp2
      slope2_flip <- -tmp2
    }else{
      slope2 <- c(slope2, tmp2)
      slope2_flip <- c(slope2_flip, -tmp2)
    }
  }
  
  # 跨一个点的斜率，包括翻转
  for (j in 1:(length(y)-2)){
    tmp1 <- (y[j+2] - y[j]) / (x[j+2] - x[j])
    slope1 <- c(slope1, tmp1)
    slope1_flip <- c(slope1_flip, -tmp1)
    
    tmp2 <- (y[j+2] - y[j])
    slope2 <- c(slope2, tmp2)
    slope2_flip <- c(slope2_flip, -tmp2)
  }
  
  slope_NED[, genes[i]:= slope1]
  slope_NED[, paste0(genes[i],'_flip'):= slope1_flip]
  slope_ED[, genes[i]:= slope2]
  slope_ED[, paste0(genes[i],'_flip'):= slope2_flip]
  if (num %% 2000 == 0 | num == length(genes)){
    print(num)
  }
  num <- num + 1
}

path1 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getSlope/','Slope_',type,'_ED.csv')
path2 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getSlope/','Slope_',type,'_NED.csv')
fwrite(slope_ED, path1)
fwrite(slope_NED, path2)
t4 <- proc.time()
t <- t4 - t3
print(paste0(type, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
rm(df, df1, df2, slope_ED, slope_NED, path1, path2)
gc()
# }

t5 <- proc.time()
t <- t5 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
