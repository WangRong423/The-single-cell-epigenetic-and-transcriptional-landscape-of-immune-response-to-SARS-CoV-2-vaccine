library(clusterProfiler)
library(org.Hs.eg.db)


cell_type<-readRDS(file="/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.6.celltype.rds")
cell_type<-gsub('/','_',cell_type)
cell_type<-gsub(' ','_',cell_type)
KEGG_DEG<-list()
GO_DEG <- list()
clustering=list()
filterd=list()
filterd1=list()
trans_KEGG_DEG=list()
trans_GO_DEG=list()
temp <- list()


for( i in 1:length(cell_type)){
print(cell_type[i])
annotation<-read.table(paste0('/database/wangrong/Results/0712_ATAC+RNA/gene_clustering/filtered/filtered_',cell_type[i],'.csv'),header = TRUE,sep = ",")
colnames(annotation)<-gsub(' ','.',colnames(annotation))
del <- grep('_flip', annotation$gene, value = F)
annotation <- annotation[-del,]
print(dim(annotation))
clustering[[i]]<-levels(factor(annotation$clustering.results))
print(paste(cell_type[i],' is filtered'))


for( j in 1:length(clustering[[i]])){
	print(clustering[[i]][[j]])
	clustering_gene_table<-annotation[annotation$clustering.results==clustering[[i]][[j]],]
	clustering_gene_table$gene<-as.character(clustering_gene_table$gene)
tryCatch(expr = {
  df <-(bitr(clustering_gene_table$gene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db"))
KEGG_DEG[[j]] <- enrichKEGG(df$ENTREZID,organism = "hsa",pAdjustMethod = "none", pvalueCutoff = 0.05)
temp[[j]] <- setReadable(KEGG_DEG[[j]],OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
trans_KEGG_DEG[[j]]<-temp[[j]]
	trans_KEGG_DEG[[j]]@result$clustering<-rep(clustering[[i]][[j]], nrow(trans_KEGG_DEG[[j]]@result))
	print(paste0(cell_type[i],'_',clustering[[i]][[j]],'KEGG is done'))
},
error = function(e){          
    print('bug')
    }
  )}
write.csv(do.call(rbind,lapply(trans_KEGG_DEG,as.data.frame)),file=paste('/database/wangrong/Results/0712_ATAC+RNA/gene_clustering/enrichment_kegg/',cell_type[i],'_KEGG.csv'))

}  

KEGG_list <- list()
for( k in 1:length(cell_type)){
KEGG_list[[k]] <- read.table(paste('/database/wangrong/Results/0712_ATAC+RNA/gene_clustering/enrichment_kegg/6.23/6.7_',cell_type[k],'_KEGG.csv'),header = TRUE,sep = ",")
KEGG_list[[k]]$cell_type <- rep(cell_type[k], nrow(KEGG_list[[k]]))
}
ALL_KEGG <- do.call(rbind,KEGG_list)
filter_all_KEGG <- ALL_KEGG[ALL_KEGG$Count>2,]
sord_all_KEGG <- filter_all_KEGG[order(filter_all_KEGG$pvalue),]
 write.csv(sord_all_KEGG,file='/database/wangrong/Results/0712_ATAC+RNA/gene_clustering/enrichment_kegg/6.23/sord_all_KEGG_.csv')
 