library(data.table)

people <- c('P2', 'P1', 'P3', 'P4')
annotations <- c('CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ Tc1', 'CD8+ Tc1/Tc17', 
                 'CD8+ Tem', 'GC B', 'Macrophage', 'Memery B', 'Naïve B', 'Naïve CD4+ T cells', 'Naïve CD8+ T cells', 
                 'Neutrophil', 'NK', 'Plasmablasts', 'Plasmacytoid DC', 'ProB/Immature B', 'Treg')

print('Loading all.csv file.')
t1 <- proc.time()
data <- fread('/database/huhuajie/findCoexpressionGenes/results/all.csv')
t2 <- proc.time()
t <- t2 - t1
print(paste0('Down! Time taken: ', t[3][[1]] / 60, ' mins'))

for (person in people){
  print(paste0('Start processing ', person))
  t3 <- proc.time()
  persondata <- data[Participants == person]
  for (type in annotations){
    print(paste0('Start processing ', type, ' of ', person))
    t4 <- proc.time()
    df <- persondata[annotation == type, 
                     !c('V1', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 
                        'Dlabel', 'Participants', 'n200', 'n200_k34', 'annotation')]
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
    
    x <- df$Timepoints
    num <- 1
    genes <- colnames(df)
    slope_ED <- data.table()
    slope_NED <- data.table()
    for (i in 2:length(genes)){
      y <- df[[i]]
      for (j in 1:(length(y)-1)){
        tmp1 <- (y[j+1] - y[j]) / (x[j+1] - x[j])
        if (j == 1){
          slope1 <- tmp1
        }else{
          slope1 <- c(slope1, tmp1)
        }
        
        tmp2 <- (y[j+1] - y[j])
        if (j == 1){
          slope2 <- tmp2
        }else{
          slope2 <- c(slope2, tmp2)
        }
      }
      
      for (j in 1:(length(y)-2)){
        tmp1 <- (y[j+2] - y[j]) / (x[j+2] - x[j])
        slope1 <- c(slope1, tmp1)
        
        tmp2 <- (y[j+2] - y[j])
        slope2 <- c(slope2, tmp2)
      }
      
      slope_NED[, genes[i]:= slope1]
      slope_ED[, genes[i]:= slope2]
      if (num %% 2000 == 0){
        print(num)
      }
      num <- num + 1
    }
    type <- gsub(' ','_',type)
    type <- gsub('/', '_', type)
    path1 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/',person,'_',type,'_ED.csv')
    path2 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/',person,'_',type,'_NED.csv')
    fwrite(slope_ED, path1)
    fwrite(slope_NED, path2)
    t5 <- proc.time()
    t <- t5 - t4
    print(paste0(type, ' of ', person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
    rm(df, slope_ED, slope_NED, path1, path2)
    gc()
  }
  t6 <- proc.time()
  t <- t6 - t3
  print(paste0(person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
  rm(persondata)
  gc()
}


t7 <- proc.time()
t <- t7 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
