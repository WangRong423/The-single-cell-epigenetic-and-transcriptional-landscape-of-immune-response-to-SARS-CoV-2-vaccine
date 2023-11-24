library(log4r)
library(Seurat)
library(data.table)


## We divide the large matrix into smaller chunks to save it as a CSV file.
write_sparse_csv <- function(x, file, ..., chunk = 1000){
  passes <- ncol(x) %/% chunk
  remaining <- ncol(x) %% chunk
  if(passes > 0){
    idx <- seq_len(chunk)
    y <- x[ , idx, drop = FALSE]
    y <- as.matrix(y)
    passes <- passes - 1L
    for(i in seq_len(passes)){
      idx <- idx + chunk
      tmp <- x[ , idx, drop = FALSE]
      tmp <- as.matrix(tmp)
      y <- data.frame(y, tmp)
    }
    if(remaining > 0){
      p <- idx[chunk] + 1
      q <- idx[chunk] + remaining
      tmp <- x[ , p:q, drop = FALSE]
      tmp <- as.matrix(tmp)
      y <- data.frame(y, tmp)
      y <- t(y)
      y <- as.data.frame(y)
      fwrite(y, file, sep = ",", col.names = TRUE, row.names = TRUE)
    }
  } 
}

## extract the necessary data from the RDS file, including:
## 1. gene expression data
## 2. samples metadata
rnaData <- readRDS('/database/findCoexpressionGenes/data/5.18.annotation.signac.RDS')
geneExpression <- rnaData@assays$RNA
geneExpression <- geneExpression@data
metadata <- rnaData@meta.data
## delete the variable, free up memory.
rm(rnaData)
gc()

## create a log to record the outputs
path = 'database/findCoexpressionGenes/logs/01_dataPreprocess.txt'
create.logger(path, level = 2) 
appender = file_appender(path)

## save data as CSV files
appender(level = 'INFO', 'Start saving geneExpression.csv.')
write_sparse_csv(geneExpression, '/database/findCoexpressionGenes/outputs/geneExpression.csv')
appender(level = 'INFO', 'geneExpression.csv was saved.')

appender(level = 'INFO', 'Start saving metaData.csv.')
write.csv(metadata, '/database/findCoexpressionGenes/outputs/metaData.csv')
appender(level = 'INFO', 'metaData.csv was saved.')

rm(geneExpression, metadata)
gc()

## merge data by barcode information
appender(level = 'INFO', 'Start merging data.')
data <- fread('/database/findCoexpressionGenes/outputs/geneExpression.csv')
meta <- fread('/database/findCoexpressionGenes/outputs/metaData.csv')
meta$V1 <- gsub("-","\\.",meta$V1)
all <- merge.data.table(meta, data, by = 'V1')
fwrite(all, '/database/findCoexpressionGenes/outputs/all.csv', col.names = TRUE)
appender(level = 'INFO', 'Merged-data was saved!.')