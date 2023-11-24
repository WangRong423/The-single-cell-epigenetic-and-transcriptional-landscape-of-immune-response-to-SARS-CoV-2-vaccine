library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)

#rename_pbmc<-readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.18.call_peak.assay.RDS")
#rename_pbmc
#rename_pbmc <- SCTransform(rename_pbmc)
#saveRDS(rename_pbmc,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.18.call_peak.sct.RDS",compress=F)
#rename_pbmc <-readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.18.call_peak.sct.RDS")

#DefaultAssay(rename_pbmc) <- "SCT"
#rename_pbmc <- RunPCA(rename_pbmc)
#rename_pbmc <- RunUMAP(rename_pbmc, dims = 1:50, reduction.name = "umap.rna")
#rename_pbmc <- FindNeighbors(rename_pbmc, dims = 1:50)
#rename_pbmc <- FindClusters(rename_pbmc, algorithm = 3)
#print("umap.rna is done!")
#DefaultAssay(rename_pbmc) <- "peaks"
#rename_pbmc <- FindTopFeatures(rename_pbmc, min.cutoff = 10)
#rename_pbmc <- RunTFIDF(rename_pbmc)
#rename_pbmc <- RunSVD(rename_pbmc)
#rename_pbmc <- RunUMAP(rename_pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.atac")
#rename_pbmc <- FindNeighbors(rename_pbmc, reduction = "lsi", dims = 2:40)
#rename_pbmc <- FindClusters(rename_pbmc, algorithm = 3)
#print("umap.atac is done!")
rename_pbmc<-readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.19.umap.RDS")
DefaultAssay(rename_pbmc) <- "peaks"
rename_pbmc <- RegionStats(rename_pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)



b6<-brewer.pal(12,"Paired")
b7<-brewer.pal(8,"Set2")
c7<-c(b6[c(1,2,3,4,5,7,8,9,10,11,12)],b7)
b8<-brewer.pal(9,"Pastel1")
b9<-brewer.pal(11,"Spectral")
c33<-c(b6,b7,b8,b9)

DimPlot(rename_pbmc, reduction = "umap.rna", label = TRUE, repel = TRUE,group.by = 'celltype',cols=c7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.rna", label = TRUE, repel = TRUE,group.by = 'annotation',cols=c33)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.rna", label = TRUE, repel = TRUE,group.by = 'Timepoints',cols=c7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.rna", label = TRUE, repel = TRUE,group.by = 'Participants',cols=b7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE,group.by = 'celltype',cols=c7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE,group.by = 'annotation',cols=c33)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE,group.by = 'Timepoints',cols=c7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

DimPlot(rename_pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE,group.by = 'Participants',cols=b7)+
annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# link peaks to genes
rename_pbmc <- LinkPeaks(object = rename_pbmc,peak.assay = "peaks",expression.assay = "SCT")
#Testing 29845 genes and 269029 peaks,Found gene coordinates for 18233 genes


saveRDS(rename_pbmc,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.19.peak_to_gene.RDS",compress=F)
print("all done!")

Links_pbmc <- Links(rename_pbmc)
# ratio of positive to negative links
sum(Links_pbmc $score < 0) / length(Links_pbmc ) * 100
#2.908547
# total over 100 kb
lnk<-Links_pbmc
sum(width(lnk) > 100000) / length(lnk)
#0.4743392
link.df <- as.data.frame (Links_pbmc)
# for each gene, find the number of linked peaks
links_per_gene <- link.df %>% 
  mutate(pos_link = score > 0) %>% 
  group_by(gene) %>% 
  summarise(positive_links = sum(pos_link), negative_links = sum(!pos_link))
mean(links_per_gene$positive_links + links_per_gene$negative_links)
#7.034951
sd(links_per_gene$positive_links + links_per_gene$negative_links)
#8.606771

# total links per gene
link_per_gene_plot <- links_per_gene %>%
  group_by(positive_links, negative_links) %>%
  summarise(count = n()) %>% 
  ggplot(data = ., aes(x = positive_links, y = negative_links, fill = log10(count+1))) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c() +
  ylab("Total negative links") +
  xlab("Total positive links") +
  ggtitle("Number of linked peaks per gene")
  
# number of linked genes per peak
genes_per_link <- link.df %>% 
  mutate(pos_link = score > 0) %>% 
  group_by(peak) %>% 
  summarise(positive_links = sum(pos_link), negative_links = sum(!pos_link))
mean(genes_per_link$positive_links + genes_per_link$negative_links)
#1.292599 
sd(genes_per_link$positive_links + genes_per_link$negative_links)
#0.781051

# total links per gene
gene_per_link_plot <- genes_per_link %>%
  group_by(positive_links, negative_links) %>%
  summarise(count = n()) %>% 
  ggplot(data = ., aes(x = positive_links, y = negative_links, fill = log10(count+1))) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(breaks = 0:10) +
  scale_fill_viridis_c() +
  ylab("Total negative links") +
  xlab("Total positive links") +
  ggtitle("Number of linked genes per peak")
  
  # distance from peak to tss
distplot_positive <- ggplot(data = link.df[link.df$score > 0, ], aes(x = width)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("") +
  ylab("Count") +
  ggtitle("Positive gene associations")

distplot_negative <- ggplot(data = link.df[link.df$score < 0, ], aes(x = width)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("Distance to gene TSS (bp)") +
  ylab("Count") +
  ggtitle("Negative gene associations")

pval_dist <- ggplot(data = link.df, mapping = aes(x = pvalue)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("p-value") +
  ylab("Count") +
  ggtitle("p-value distribution")
  
lnkplot <- gene_per_link_plot / link_per_gene_plot
distances <- pval_dist / distplot_positive / distplot_negative
ggsave(filename = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.25.lnkplot.pdf", plot = lnkplot, width = 5, height = 6, units = 'in')
ggsave(filename = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.25.distances.pdf", plot = distances, width = 4, height = 6, units = 'in')

gene.activities <- GeneActivity(rename_pbmc)
rename_pbmc[['RNA_score']] <- CreateAssayObject(counts = gene.activities)
rename_pbmc <- NormalizeData(object = rename_pbmc,assay = 'RNA_score',normalization.method = 'LogNormalize',scale.factor = median(rename_pbmc$nCount_RNA))

DefaultAssay(rename_pbmc) <- "RNA"
pbmc.marker <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/hsPBMC_markers_v3--wangrong.xlsx")
 
 rename_pbmc$celltype<-gsub("Eryth/Platelet","CD14+ Mono",rename_pbmc$celltype)
 rename_pbmc$annotation<-gsub("Eryth/Platelet","CD14+ Mono_4",rename_pbmc$annotation)
 rename_pbmc$annotation<-gsub("NKNK","NK",rename_pbmc$annotation)
 rename_pbmc$annotation<-gsub("Memrary","Memory",rename_pbmc$annotation)
 
 saveRDS(rename_pbmc,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.13.signac_peak_to_gene.rds",compress = F)