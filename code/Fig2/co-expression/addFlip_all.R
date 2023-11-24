library(data.table)

t0 <- proc.time()
annotations <- c('NK', 'Naïve CD8+ T cells', 'Naïve CD4+ T cells', 'CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'CD4+ Tem',
                 'CD8+ Tc1', 'CD8+ Tc1/Tc17', 'CD8+ Tem', 'GC B', 'Macrophage', 'Memery B', 'Naïve B',
                 'Neutrophil', 'Plasmablasts', 'Plasmacytoid DC', 'ProB/Immature B', 'Treg')

for (type in annotations){
  print(paste0('Start processing ', type))
  t1 <- proc.time()
  type <- gsub(' ','_',type)
  type <- gsub('/', '_', type)
  EDfilepath <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/allPeople/ignore/','all_',type,'_ED.csv')
  NEDfilepath <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/allPeople/ignore/','all_',type,'_NED.csv')
  
  ED <- fread(EDfilepath)
  NED <- fread(NEDfilepath)
  
  num <- 1
  genes <- colnames(ED)
  
  for (i in 1:length(genes)){
    
    y1 <- ED[[i]]
    y1 <- as.numeric(unlist(y1))
    y1 <- -y1
    ED[, paste0(genes[i],'_flip'):= y1]
    
    y2 <- NED[[i]]
    y2 <- as.numeric(unlist(y2))
    y2 <- -y2
    NED[, paste0(genes[i],'_flip'):= y2]
    
    if (num %% 2000 == 0 | num == length(genes)){
      print(num)
    }
    num <- num + 1
  }
  path1 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/allPeople/ignore/addFlip/','all_',type,'_ED.csv')
  path2 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/allPeople/ignore/addFlip/','all_',type,'_NED.csv')
  fwrite(ED, path1)
  fwrite(NED, path2)
  t2 <- proc.time()
  t <- t2 - t1
  print(paste0(type, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
  rm(EDfilepath, NEDfilepath, genes, ED, NED, path1, path2)
  gc()
}

t3 <- proc.time()
t <- t3 - t0
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
