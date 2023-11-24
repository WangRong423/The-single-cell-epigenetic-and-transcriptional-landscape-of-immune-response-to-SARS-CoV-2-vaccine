library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(limma)
library(edgeR)
library(mashr)


counts <- pbmc@assays$RNA@counts
metadata<-pbmc@meta.data
metadata<-metadata[,c("orig.ident","Participants","Timepoints","celltype")]
metadata$orig.ident <- gsub("-","",metadata$orig.ident)
metadata$Timepoints <- as.factor(metadata$Timepoints)
metadata$Participants <- as.factor(metadata$Participants)
metadata$orig.ident <- as.factor(metadata$orig.ident)

sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
#groups <- colData(sce)[, c("celltype", "orig.ident")]
Timepoints <- purrr::set_names(levels(sce$Timepoints))
nt <- length(Timepoints)
celltype <- purrr::set_names(levels(sce$celltype))
na<- length(celltype)
Participants <- purrr::set_names(levels(sce$Participants))
np<- length(Participants)
orig.ident <- purrr::set_names(levels(sce$orig.ident))
nr<- length(orig.ident)
table(sce$orig.ident)				
nr_cells <- as.numeric(table(sce$orig.ident))
mr <- match(orig.ident, sce$orig.ident)
eia <- data.frame(colData(sce)[mr, ], nr_cells, row.names = NULL) 
groups <- colData(sce)[, c("celltype","orig.ident")]
pb1 <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
class(pb1)
dim(pb1)
pb1[1:6, 1:6]
splitf <- sapply(stringr::str_split(rownames(pb1), pattern = "_",  n = 2), `[`, 1)
pb2 <- split.data.frame(pb1, factor(splitf)) %>%lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
str(pb2)
options(width = 100)
table(sce$celltype, sce$orig.ident)

get_sample_ids <- function(x){
        pb2[[x]] %>%colnames()
}
de_samples <- map(1:length(celltype ), get_sample_ids) %>%unlist()

samples_list <- map(1:length(celltype), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb2)[x], each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(celltype), get_cluster_ids) %>%unlist()

gg_df <- data.frame(celltype = de_cluster_ids,orig.ident = de_samples)
gg_df <- left_join(gg_df, eia[, c("orig.ident","Timepoints", "Participants")]) 
metadata <- gg_df %>%dplyr::select(celltype , orig.ident, Timepoints, Participants)         

metadata$celltype <- as.factor(metadata$celltype)
clusters <- levels(metadata$celltype)
save(pb2,metadata, clusters, file = "/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.21.psudocelltype.matrix.RData")



load( file = "/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.21.psudocelltype.matrix.RData")

clusters<-gsub("/","_",clusters)
metadata$celltype<-gsub("/","_",metadata$celltype)
names(pb2)<-gsub("/","_",names(pb2))

for( i in 1:length(clusters)){
print(clusters[i])

cluster_metadata <- metadata[which(metadata$celltype == clusters[i]), ]
cluster_metadata <- cluster_metadata[which(cluster_metadata$Timepoints != "Day171"),]
cluster_metadata <- cluster_metadata[which(cluster_metadata$Timepoints != "Day185"),]
print(head(cluster_metadata))
rownames(cluster_metadata) <- cluster_metadata$orig.ident
counts <- pb2[[clusters[i]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])   
dim(cluster_counts)

condition <- factor(cluster_metadata$Timepoints, levels = c(c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44") ) )
design <- model.matrix(~0 + condition ,cluster_metadata)

y <- DGEList(cluster_counts)
y <- calcNormFactors(y)

v = voom(y, design)

timepoints.gene <- lmFit( v , design)
timepoints.Bayes.gene <- eBayes(timepoints.gene)

timepoints.condition.gene.mean <-timepoints.Bayes.gene$coefficients[, 1:10]
timepoints.condition.gene.se <- ( timepoints.Bayes.gene$stdev.unscaled * sqrt(timepoints.Bayes.gene$s2.post) )[,1:10]

colnames(timepoints.condition.gene.mean) <- colnames(timepoints.condition.gene.se) <- gsub("condition","", colnames(timepoints.condition.gene.se))
#machr
data.gene = mash_set_data(timepoints.condition.gene.mean, timepoints.condition.gene.se)
data.L.gene = mash_update_data(data.gene, ref = 1)
U.c.gene = cov_canonical(data.L.gene)
print(head(names(U.c.gene)))
mashcontrast.model.gene = mash(data.L.gene, U.c.gene  )

pdf(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/7.15_timepoints_',clusters[[i]],'.pairwise_sharing.pdf'))
corrplot(get_pairwise_sharing(mashcontrast.model.gene, factor=0.5, lfsr_thresh = 0.1) , method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude\n(< Factor of 2)', mar = c(4,0,4,0))
dev.off()
print(head( get_pm(mashcontrast.model.gene) ))

est_pi <- data.frame( Type = factor( names(get_estimated_pi(mashcontrast.model.gene) ), levels = c(names(get_estimated_pi(mashcontrast.model.gene) )) ), 
                      estimates =  get_estimated_pi(mashcontrast.model.gene)  )

pdf(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/7.15_timepoints_',clusters[[i]],'.est_pi.pdf'))
ggplot( est_pi, aes( x = Type, y = estimates) )+
	geom_bar(stat = "identity")+
	ylab(expression( pi) )+
	theme_classic()+
	theme(text = element_text(color = "black", size = 18), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = -45, vjust = 0), plot.margin = margin(r = 10, b = 10) )+
	scale_y_continuous(expand = c(0,0))
dev.off()
sig.DEG_ID <- get_significant_results(mashcontrast.model.gene,thresh = 0.1) 

mashResult.gene <- get_pm(mashcontrast.model.gene)
colnames(mashResult.gene) <- gsub("-Day0",".log2FC",colnames(mashResult.gene))
mashResult.gene <- cbind(mashResult.gene, get_lfsr( mashcontrast.model.gene ) )
colnames(mashResult.gene) <- gsub("-Day0",".lfsr",colnames(mashResult.gene))
sig.mashResults.gene<- mashResult.gene[sig.DEG_ID,]
sig.mashResults.gene<-as.data.frame(sig.mashResults.gene)
sig.mashResults.gene$total.lfsr<-rowSums(sig.mashResults.gene[, c(10,11,12,13,14,15,16,17,18)])
sig.mashResults.gene<-sig.mashResults.gene[order(sig.mashResults.gene$total.lfsr),]
sig.mashResults.gene<-sig.mashResults.gene[1:25,]
write.csv(sig.mashResults.gene,file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/7.15_timepoints_',clusters[[i]],'.sig.mashResults.gene.top25.csv'))


sig.mashResults.gene_melt <- data.frame( geneName = rep( rownames(sig.mashResults.gene), 9 ),
                                                    posterior.log2FC = c(sig.mashResults.gene$Day1.log2FC,
                                                                        sig.mashResults.gene$Day3.log2FC,
                                                                        sig.mashResults.gene$Day6.log2FC,
                                                                        sig.mashResults.gene$Day14.log2FC,
																		sig.mashResults.gene$Day30.log2FC,
                                                                        sig.mashResults.gene$Day31.log2FC,
                                                                        sig.mashResults.gene$Day33.log2FC,
                                                                        sig.mashResults.gene$Day36.log2FC,
																		sig.mashResults.gene$Day44.log2FC),
                                                    lfsr = c( sig.mashResults.gene$Day1.lfsr,
                                                              sig.mashResults.gene$Day3.lfsr,
                                                              sig.mashResults.gene$Day6.lfsr,
                                                              sig.mashResults.gene$Day14.lfsr,
															  sig.mashResults.gene$Day30.lfsr,
                                                              sig.mashResults.gene$Day31.lfsr,
                                                              sig.mashResults.gene$Day33.lfsr,
                                                              sig.mashResults.gene$Day36.lfsr,
															  sig.mashResults.gene$Day44.lfsr),
                                                    condition = c( rep("Day1", nrow(sig.mashResults.gene) ),
                                                                   rep("Day3", nrow(sig.mashResults.gene) ),
                                                                   rep("Day6", nrow(sig.mashResults.gene) ),
                                                                   rep("Day14", nrow(sig.mashResults.gene) ), 
                                                                   rep("Day30", nrow(sig.mashResults.gene) ),
                                                                   rep("Day31", nrow(sig.mashResults.gene) ),
                                                                   rep("Day33", nrow(sig.mashResults.gene) ),
                                                                   rep("Day36", nrow(sig.mashResults.gene) ),
																   rep("Day44", nrow(sig.mashResults.gene) )))

sig.mashResults.gene_melt$condition <- factor(sig.mashResults.gene_melt$condition , levels = c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44") ) 
pdf(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/7.15_timepoints_',clusters[[i]],'.sig.mashResults.gene.pdf'))
ggplot(sig.mashResults.gene_melt, aes( x = condition, y = geneName, size = -log(lfsr), color = posterior.log2FC ) )+geom_point()+scale_color_gradientn( colors=c("cyan","red" ))+scale_size_continuous(breaks= -log(c(0.4,0.1,0.05,0.01,0.001)),labels=c(0.4,0.1,0.05,0.01,0.001), range = c(0.5,8),  name = "lfsr" )+geom_point(data =sig.mashResults.gene_melt[sig.mashResults.gene_melt$lfsr<0.1,], aes( x = condition, y = geneName, size = -log(lfsr) ), shape = 1, color = "black" )+ggtitle("Significant in at least one conditions")+theme_bw()+theme(
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1),
                   plot.title = element_text(size = 16, hjust = 0.5),
                   axis.title=element_blank(),
                   legend.title=element_text(size=15 ),legend.text = element_text(size = 20),
                   axis.text.x = element_text(size = 15,angle = -20, vjust = 0.3 , color = "black") ,axis.text.y = element_text(size = 15, color = "black" ) )
dev.off()
print(paste0(clusters[i],'is done!'))
}

marshr.gene<-list()
for (i in 1:length(clusters)){
print(clusters[i])
marshr.gene[[i]] <- read.table(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/6.24_timepoints_',clusters[[i]],'.sig.mashResults.gene.top25.csv'),sep = ",",header=TRUE)
marshr.gene[[i]]$celltype<-rep(clusters[[i]],25)
}
all_mashr<-do.call(rbind,marshr.gene)