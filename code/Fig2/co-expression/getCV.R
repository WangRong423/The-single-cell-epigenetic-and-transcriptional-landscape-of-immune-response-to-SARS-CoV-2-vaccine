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
    meandf <- df[, lapply(.SD, mean, na.rm = TRUE), by = Timepoints]
    sddf <- df[, lapply(.SD, sd, na.rm = TRUE), by = Timepoints]
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
      }else if (meandf[, 1][i] == 'Day171'){
        meandf[, 1][i] = 171
      }else if (meandf[, 1][i] == 'Day185'){
        meandf[, 1][i] = 185
      }
    }
    meandf$Timepoints <- as.numeric(meandf$Timepoints)
    meandf <- meandf[order(meandf$Timepoints),]
    
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
      }else if (sddf[, 1][i] == 'Day171'){
        sddf[, 1][i] = 171
      }else if (sddf[, 1][i] == 'Day185'){
        sddf[, 1][i] = 185
      }
    }
    sddf$Timepoints <- as.numeric(sddf$Timepoints)
    sddf <- sddf[order(sddf$Timepoints),]
    
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
        if (mean[j+1] == 0 || mean[j] == 0){
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
        if (mean[j+1] == 0 || mean[j] == 0){
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
      
      if (num %% 2000 == 0){
        print(num)
      }
      num <- num + 1
    }
    type <- gsub(' ','_',type)
    type <- gsub('/', '_', type)
    path1 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getCV/',person,'_',type,'_ED.csv')
    path2 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getCV/',person,'_',type,'_NED.csv')
    path3 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getCV/',person,'_',type,'.csv')
    fwrite(CV_ED, path1)
    fwrite(CV_NED, path2)
    fwrite(CV, path3)
    t5 <- proc.time()
    t <- t5 - t4
    print(paste0(type, ' of ', person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
    rm(df, meandf, sddf, CV_ED, CV_NED, CV, path1, path2, path3)
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
