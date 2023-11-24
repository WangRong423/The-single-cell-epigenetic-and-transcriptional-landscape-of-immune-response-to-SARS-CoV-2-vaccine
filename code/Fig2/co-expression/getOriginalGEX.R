library(data.table)

people <- c('P2', 'P1', 'P3', 'P4')
annotations <- c('CD14+ Mono', 'CD16+ Mono', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ Tc1', 'CD8+ Tc1/Tc17', 
'CD8+ Tem', 'GC B', 'Macrophage', 'Memery B', 'Naïve B', 'Naïve CD4+ T cells', 'Naïve CD8+ T cells', 
'Neutrophil', 'NK', 'Plasmablasts', 'Plasmacytoid DC', 'ProB/Immature B', 'Treg')

print('Loading all_new.csv file.')
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
                     !c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 
                        'Dlabel', 'n200', 'n200_k34')]
    
    type <- gsub(' ','_',type)
    type <- gsub('/', '_', type)
    path <- paste0('/database/huhuajie/findCoexpressionGenes/results/getGEX/original/',person,'_',type,'.csv')
    fwrite(df, path)
    t5 <- proc.time()
    t <- t5 - t4
    print(paste0(type, ' of ', person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
  }
  t6 <- proc.time()
  t <- t6 - t3
  print(paste0(person, ' was down! Time taken: ',  t[3][[1]] / 60, ' mins'))
}

t7 <- proc.time()
t <- t7 - t1
print(paste0('All time spend: ',  t[3][[1]] / 60, ' mins'))
