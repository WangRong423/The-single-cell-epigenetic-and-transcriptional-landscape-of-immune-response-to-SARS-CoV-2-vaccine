# 计算零值百分比

library(data.table)

celltypes <- c('Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
               'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
               'CD8+ Tem', 'Eryth/Platelet', 'MAIT')


# 计算一组值的零值百分比，table函数统计每个值出现的频次
getZero <- function(x){
  if (names(table(x))[1] == '0'){
    nums_zero <- table(x)[['0']]
    nums_all <- length(x)
    res <- 1 - (nums_zero / nums_all)
  }else{
    res <- 0
  }
  return(res)
}

print('Loading all.csv file.')
t1 <- proc.time()
data <- fread('/database/huhuajie/findCoexpressionGenes/finalResults/all.csv')
t2 <- proc.time()
t <- t2 - t1
print(paste0('Down! Time taken: ', t[3][[1]] / 60, ' mins'))

for (type in celltypes){
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

  # 求该细胞类型所有基因不同时间点的零值百分比
  zerodf <- df[, lapply(.SD, getZero), by = Timepoints]

  # 按时间点升序排列 
  for (i in 1:nrow(zerodf)){
    if (zerodf[, 1][i] == 'Day0'){
      zerodf[, 1][i] = 0
    }else if (zerodf[, 1][i] == 'Day1'){
      zerodf[, 1][i] = 1
    }else if (zerodf[, 1][i] == 'Day3'){
      zerodf[, 1][i] = 3
    }else if (zerodf[, 1][i] == 'Day6'){
      zerodf[, 1][i] = 6
    }else if (zerodf[, 1][i] == 'Day14'){
      zerodf[, 1][i] = 14
    }else if (zerodf[, 1][i] == 'Day30'){
      zerodf[, 1][i] = 30
    }else if (zerodf[, 1][i] == 'Day31'){
      zerodf[, 1][i] = 31
    }else if (zerodf[, 1][i] == 'Day33'){
      zerodf[, 1][i] = 33
    }else if (zerodf[, 1][i] == 'Day36'){
      zerodf[, 1][i] = 36
    }else if (zerodf[, 1][i] == 'Day44'){
      zerodf[, 1][i] = 44
    }
  }
  zerodf$Timepoints <- as.numeric(zerodf$Timepoints)
  zerodf <- zerodf[order(zerodf$Timepoints)]

  # 计算零值百分比，
  # zeroPT保存原始零值百分比值，zeroPT_ED存放零值百分比变化的等距斜率，zeroPT_NED存放零值百分比变化的不等距斜率  
  x <- zerodf[, .(Timepoints)]
  x <- as.numeric(unlist(x))
  num <- 1
  genes <- colnames(zerodf)
  zeroPT <- data.table()
  zeroPT_ED <- data.table()
  zeroPT_NED <- data.table()

  for (i in 2:length(genes)){
    zero <- zerodf[[i]]
    zero <- as.numeric(unlist(zero))
    zeroPT[, genes[i]:= zero]

    for (j in 1:(length(zero)-1)){

      tmp1 <- (zero[j+1] - zero[j]) / (x[j+1] - x[j])
      if (j == 1){
        zero1 <- tmp1
      }else{
        zero1 <- c(zero1, tmp1)
      }

      tmp2 <- (zero[j+1] - zero[j])
      if (j == 1){
        zero2 <- tmp2
      }else{
        zero2 <- c(zero2, tmp2)
      }
    }
    zeroPT_NED[, genes[i]:= zero1]
    zeroPT_ED[, genes[i]:= zero2]

    if (num %% 2000 == 0 | num == length(genes)){
      print(num)
    }
    num <- num + 1
  }
  type <- gsub(' ','_',type)
  type <- gsub('/', '_', type)
  path1 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getZero/','ZeroPT_',type,'_ED.csv')
  path2 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getZero/','ZeroPT_',type,'_NED.csv')
  path3 <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getZero/','ZeroPT_',type,'.csv')
  fwrite(zeroPT_ED, path1)
  fwrite(zeroPT_NED, path2)
  fwrite(zeroPT, path3)
  t4 <- proc.time()
  t <- t4 - t3
  print(paste0(type, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
  rm(df, zerodf, zeroPT_ED, zeroPT_NED, zeroPT, path1, path2, path3)
  # rm(df, zerodf, zeroPT, path3)
  gc()
}
  
t5 <- proc.time()
t <- t5 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
