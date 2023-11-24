---
title: "HIPPO"
author: "wangrong"
date: '2022-05-16'
output: workflowr::wflow_html
code_folding: hide
editor_options:
  
  chunk_output_type: console
---

library(knitr)
library(rmarkdown)
library(tidyverse)
library(Seurat)
library(DoubletDecon)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(Seurat)
library(limma)
library(DoubletDecon)
library(scater)
library(HIPPO)
library(lightHippo)
library(irlba)

co_all <-  readRDS("/database/wangrong/Results/0712_ATAC+RNA/Signac/co_joint/5.14_merge_cojoint.rds")
n200 <- readRDS("/database/jiangxinyu/Result/Data/0516_Total_cojoint_res_200cells_random5000_improve0.8.rds")
n400 <- readRDS("/database/jiangxinyu/Result/Data/0516_Total_cojoint_res_400cells_random5000_improve0.8.rds")
visualize_hippo_hierarchy(n200)
visualize_hippo_hierarchy(n400)
co_all@meta.data$n200 <- n200$next_round_IDs
table(co_all@meta.data$n200)
Idents(co_all) <- co_all$n200
pbmc.marker <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/hsPBMC_markers_v3--wangrong.xlsx")
DotPlot(co_all, features = unique(pbmc.marker$Gene) ,col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "pink",low ="white", midpoint = 0)


new.clusters_33 <- cut_hierarchy(n200, K = 33, cut_sequence = TRUE)
visualize_hippo_hierarchy(new.clusters_33)
co_all@meta.data$n200_k33 <- new.clusters_33$labels
table(co_all@meta.data$n200_k33)
Idents(co_all) <- co_all$n200_k33
DotPlot(co_all, features = unique(pbmc.marker$Gene) ,col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
DotPlot(co_all, features = unique(pbmc.marker[31:60,]$Gene),col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
DotPlot(co_all, features = c("CD3D","CD8A","CD8B","CD4","CD40LG","CCR7","CD27","MK167","CXCR5","FOXP3","GZMB","IFNG","TRDV2","TRGV9","TRAV1","NCAM1","NCR1","FCGR3A","IL2","IL1A","IL1B","IL10","IL22","TNF","LTA","IL4","IL5","IL13","IL17A","IL17F","IL21","IL12A","IL12B"),scale = FALSE,col.min=-2, col.max=2,idents=c('3','4','5','6','7','10','9','12','13','16','17','18','19','22','23','25','26','27','30','32','33'))+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)


reference <- LoadH5Seurat("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/pbmc_multimodal.h5seurat")
pbmc3k <- SCTransform(co_all, verbose = FALSE)
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc3k,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
pbmc3k <- MapQuery(
  anchorset = anchors,
  query = pbmc3k,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 <- DimPlot(pbmc3k, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()+ggsci::scale_color_igv()
p2 <- DimPlot(pbmc3k, reduction = "ref.umap", group.by = "n200_k33", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()+ggsci::scale_color_igv()
p1 + p2
table(pbmc3k$n200_k33,pbmc3k$predicted.celltype.l2)
b <- table(pbmc3k$n200_k33,pbmc3k$predicted.celltype.l2)
b <- as.data.frame(b)
colnames(b)<-c("ours", "azimuth", 'value')
b<-b %>%group_by(ours) %>%mutate(freq = value / sum(value))
ggplot(a,aes(x=ours,y=azimuth,fill=freq))+
  geom_tile(colour="white",size=0.2)+
  scale_fill_distiller(palette = "Spectral")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_fixed()

#c
c <- table(pbmc3k$n200,pbmc3k$predicted.celltype.l2)
c <- as.data.frame(c)
colnames(c)<-c("ours", "azimuth", 'value')
c<-c %>%group_by(ours) %>%mutate(freq = value / sum(value))
ggplot(c,aes(x=ours,y=azimuth,fill=freq))+
  geom_tile(colour="white",size=0.2)+
  scale_fill_distiller(palette = "Spectral")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_fixed()

head(Idents(co_all))
My_levels <- c('1','14','20','31','11','2','8','28','21','29','6','15','24','9','27','4','7','12','16','23','33','5','18','32','22','3','17','25','26','10','13','30','19')
Idents(co_all) <- factor(Idents(co_all), levels= My_levels)
new.cluster.ids <- c("CD16+ Mono", "CD14+ Mono_1","CD14+ Mono_2","CD14+ Mono_3","Plasmacytoid DC",   "Naïve B_1", "Naïve B_2","Naïve B_3","Plasmablasts/Memrary B"," Intermediate B", "NK/NKT_1", "NK/NKT_2", "NK/NKT_3", "Naïve CD8+ T cells_1", "Naïve CD8+ T cells_2", "CD8+ Tem_1", "CD8+ Tem_2", "CD8+ Tem_3", "CD8+ Tem_4","CD8+ Tem_5", "CD8+ Tem_6", "Naïve CD4+ T cells_1","Naïve CD4+ T cells_2","Naïve CD4+ T cells_3","Treg", "CD4+ Tcm_1","CD4+ Tcm_2", "CD4+ Tcm_3", "CD4+ Tcm_4","CD4+ Tcm_5","CD4+ Tcm_6", "MAIT","Eryth/Platelet")
names(new.cluster.ids) <- levels(co_all)
rename_pbmc <- RenameIdents(co_all, new.cluster.ids)

new.cluster.ids <- c("CD16+ Mono", "CD14+ Mono","CD14+ Mono","CD14+ Mono","Plasmacytoid DC",   "Naïve B", "Naïve B","Naïve B","Plasmablasts/Memrary B"," Intermediate B", "NK/NKT", "NK/NKT", "NK/NKT", "Naïve CD8+ T cells", "Naïve CD8+ T cells", "CD8+ Tem", "CD8+ Tem", "CD8+ Tem", "CD8+ Tem","CD8+ Tem", "CD8+ Tem_6", "Naïve CD4+ T cells","Naïve CD4+ T cells","Naïve CD4+ T cells","Treg", "CD4+ Tcm","CD4+ Tcm", "CD4+ Tcm", "CD4+ Tcm","CD4+ Tcm","CD4+ Tcm", "MAIT","Eryth/Platelet")
names(new.cluster.ids) <- levels(rename_pbmc)
rename_pbmc <- RenameIdents(rename_pbmc, new.cluster.ids)
saveRDS(rename_pbmc,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.17.annotation.signac.RDS")


#call peak
 
DefaultAssay(rename_pbmc) <- "ATAC"
peaks <- CallPeaks(rename_pbmc, macs2.path = "/home/wangrong/.local/bin/macs2",group.by="celltype")
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(fragments = Fragments(rename_pbmc),features = peaks,cells = colnames(rename_pbmc))

rename_pbmc[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,fragments = Fragments(rename_pbmc),annotation = annotation)
											
saveRDS(rename_pbmc,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.17.annotation.call_peak.RDS")

##
DefaultAssay(rename_pbmc) <- "RNA"
rename_pbmc <- SCTransform(rename_pbmc)
rename_pbmc <- RunPCA(rename_pbmc)
rename_pbmc <- RunUMAP(rename_pbmc, dims = 1:50, reduction.name = "umap.rna")
rename_pbmc <- FindNeighbors(rename_pbmc, dims = 1:50)
rename_pbmc <- FindClusters(rename_pbmc, algorithm = 3)

DefaultAssay(rename_pbmc) <- "peaks"
rename_pbmc <- FindTopFeatures(rename_pbmc, min.cutoff = 10)
rename_pbmc <- RunTFIDF(rename_pbmc)
rename_pbmc <- RunSVD(rename_pbmc)
rename_pbmc <- RunUMAP(rename_pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.atac")
rename_pbmc <- FindNeighbors(rename_pbmc, reduction = "lsi", dims = 2:40)
rename_pbmc <- FindClusters(rename_pbmc, algorithm = 3)

DefaultAssay(rename_pbmc) <- "peaks"
rename_pbmc <- RegionStats(rename_pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
rename_pbmc <- LinkPeaks(object = rename_pbmc,peak.assay = "peaks",expression.assay = "SCT")