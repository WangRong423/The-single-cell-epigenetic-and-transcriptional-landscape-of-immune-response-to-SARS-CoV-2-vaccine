library(data.table)

people <- c('P1', 'P2', 'P3', 'P4')
# celltypes <- c('Naïve CD4+ T cells', 'Naïve CD8+ T cells', 'NK/NKT', 'Plasmablasts/Memrary B', 'Naïve B',
#                'Intermediate B', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'Plasmacytoid DC', 'Treg', 
#                'CD8+ Tem', 'Eryth/Platelet', 'MAIT')

print('Loading all.csv file.')
t1 <- proc.time()
data <- fread('/database/huhuajie/findCoexpressionGenes/finalResults/all.csv')
t2 <- proc.time()
t <- t2 - t1
print(paste0('Down! Time taken: ', t[3][[1]] / 60, ' mins'))
for (person in people){
  print(paste0('Start processing ', person))
  t3 <- proc.time()
  persondata <- data[Participants == person]
  # for (type in celltypes){
  type <- 'Eryth/Platelet CD14+ Mono new'
  print(paste0('Start processing ', type, ' of ', person))
  t4 <- proc.time()
  df <- persondata[celltype %chin% c('Eryth/Platelet','CD14+ Mono'), 
            !c("V1", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res.0.3", "seurat_clusters",
                "Dlabel", "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "Participants",
                "nucleosome_percentile", "TSS.enrichment", "TSS.percentile", "n200", "n200_k33", "annotation", "celltype")]
  df <- df[, lapply(.SD, mean, na.rm = TRUE), by = Timepoints]
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
    }else if (df[, 1][i] == 'Day171'){
      df[, 1][i] = 171
    }else if (df[, 1][i] == 'Day185'){
      df[, 1][i] = 185
    }
  }
  df$Timepoints <- as.numeric(df$Timepoints)
  df <- df[order(df$Timepoints)]
  type <- gsub(' ','_',type)
  type <- gsub('/', '_', type)
  path <- paste0('/database/huhuajie/findCoexpressionGenes/finalResults/getGEX/byPeople/',person,'_',type,'.csv')
  fwrite(df, path)
  t5 <- proc.time()
  t <- t5 - t4
  print(paste0(type, ' of ', person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
# }
  t6 <- proc.time()
  t <- t6 - t3
  print(paste0(person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
}

t7 <- proc.time()
t <- t7 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
