library(data.table)

t0 <- proc.time()
people <- c('P2', 'P1', 'P3', 'P4')
annotations <- c('CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ Tc1', 'CD8+ Tc1/Tc17', 
                 'CD8+ Tem', 'GC B', 'Macrophage', 'Memery B', 'Naïve B', 'Naïve CD4+ T cells', 'Naïve CD8+ T cells', 
                 'Neutrophil', 'NK', 'Plasmablasts', 'Plasmacytoid DC', 'ProB/Immature B', 'Treg')

for (person in people){
  print(paste0('Start processing ', person))
  t1 <- proc.time()
  
  for (type in annotations){
    print(paste0('Start processing ', type, ' of ', person))
    t2 <- proc.time()
    type <- gsub(' ','_',type)
    type <- gsub('/', '_', type)
    EDfilepath <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/',person,'_',type,'_ED.csv')
    NEDfilepath <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/',person,'_',type,'_NED.csv')
    
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
      
      if (num %% 2000 == 0){
        print(num)
      }
      num <- num + 1
    }
    path1 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/Slope/',person,'-',type,'_ED.csv')
    path2 <- paste0('/database/huhuajie/findCoexpressionGenes/results/getSlope/Slope/',person,'_',type,'_NED.csv')
    fwrite(ED, path1)
    fwrite(NED, path2)
    t3 <- proc.time()
    t <- t3 - t2
    print(paste0(type, ' of ', person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
    rm(EDfilepath, NEDfilepath, genes, ED, NED, path1, path2)
    gc()
  }
  t4 <- proc.time()
  t <- t4 - t1
  print(paste0(person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
}

t5 <- proc.time()
t <- t5 - t0
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
