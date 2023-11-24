#####
#Multiome scRNA-seq/scATAC-seq analysis with ArchR
#DATE:2022-6-14
######

library(ArchR)
library(pheatmap)
library(Seurat) 
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(magick,lib.loc = "/home/zhang_lab/R/x86_64-pc-linux-gnu-library/4.2")
addArchRThreads(threads = 96)
addArchRGenome("hg38")

batch_list=list("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
                "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
                "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
                "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")

name_list=c("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
            "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
            "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
            "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")

arrowlist<-c()
for(i in 1:length(batch_list)){
  arrowlist<-c(arrowlist,paste0('/database/wangrong/Results/0712_ATAC+RNA/atac_output/ArrowFiles/',batch_list[[i]],'.arrow'))
}										
archrproj <- ArchRProject(ArrowFiles = arrowlist,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/atac_output",copyArrows = FALSE)
archrproj

atac_input<-c()
for(i in 1:length(batch_list)){
  atac_input<-c(atac_input,paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/atac_fragments.tsv.gz'))
}
ArrowFiles<- createArrowFiles(inputFiles = atac_input,sampleNames = name_list)

rna_input<-c()
for(i in 1:length(batch_list)){
  rna_input<-c(rna_input,paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/filtered_feature_bc_matrix.h5'))
}
seRNA <- import10xFeatureMatrix(input = rna_input,names = name_list,featureType = "Gene Expression")

projpbmc <- ArchRProject(ArrowFiles = ArrowFiles,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/atac_output")
new_ArrowFiles<-addGeneExpressionMatrix(input=projpbmc, seRNA=seRNA)
proj <- new_ArrowFiles
saveRDS(proj,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.24.arch.RDS",compress =  F)
metadata<-read.table('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/sample_metadata.csv', sep = ",", header = TRUE,row.names = 1)
rownames(metadata) <- gsub("_","#",rownames(metadata))
length(rownames(metadata))
a <- rownames(metadata)
b <- a[!is.na(match(a,getCellNames(proj)))]
proj <- subsetArchRProject(proj,cells = b,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/5.24_ArchRSubset")
dim(metadata[!is.na(match(rownames(metadata),b)),])
c=metadata[!is.na(match(rownames(metadata),b)),]
proj <- addCellColData(ArchRProj = proj,data = c$annotation,name = "annotation",cells = getCellNames(proj),force = FALSE)
proj <- addCellColData(ArchRProj = proj,data = c$celltype,name = "celltype",cells = getCellNames(proj),force = FALSE)

proj <- addIterativeLSI(ArchRProj = proj, clusterParams = list(resolution = 0.2, sampleCells = 10000,n.start = 10),saveIterations = FALSE,useMatrix = "TileMatrix", depthCol = "nFrags",name = "LSI_ATAC")

proj <- addIterativeLSI(ArchRProj = proj, clusterParams = list(resolution = 0.2, sampleCells = 10000,n.start = 10),saveIterations = FALSE,useMatrix = "GeneExpressionMatrix", depthCol = "Gex_nUMI",varFeatures = 2500,firstSelection = "variable",binarize = FALSE,name = "LSI_RNA")
#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)

p1 <- plotEmbedding(archrproj, name = "celltype", embedding = "UMAP_ATAC", pal=c13_match,colorBy = "cellColData",size = 0.8, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(archrproj, name = "celltype", embedding = "UMAP_RNA",pal=c13_match,colorBy = "cellColData",size = 0.8, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(archrproj, name = "celltype", embedding = "UMAP_Combined",pal=c13_match,colorBy = "cellColData",size = 0.8, labelAsFactors=F, labelMeans=F)
plotPDF(p1, p2, p3, name = "6.30_UMAP-scATAC-scRNA-Combined_celltype", addDOC = FALSE)

#callpeak
archrproj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype",force = TRUE)
pathToMacs2 <- "/home/wangrong/.local/bin/macs2"
archrproj <- addReproduciblePeakSet(ArchRProj = archrproj,groupBy = "celltype",pathToMacs2 = pathToMacs2)

archrproj <- addPeakMatrix(archrproj)
saveRDS(archrproj,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.14.Archr_call_peak.RDS",compress = F)



p4 <- plotEmbedding(proj, name = "annotation", embedding = "UMAP_ATAC", pal=c33,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p5 <- plotEmbedding(proj, name = "annotation", embedding = "UMAP_RNA", pal=c33,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p6 <- plotEmbedding(proj, name = "annotation", embedding = "UMAP_Combined",pal=c33, colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
plotPDF(p4, p5, p6, name = "6.14_UMAP-scATAC-scRNA-Combined_annotation", addDOC = FALSE)

p2g_archr$Participants <- p2g_archr$Sample
p2g_archr$Participants <- substr(p2g_archr$Participants,1,2)
p2g_archr$Participants <- gsub("M1","P1",p2g_archr$Participants)
p2g_archr$Participants <- gsub("M2","P2",p2g_archr$Participants)
p2g_archr$Participants <- gsub("M3","P3",p2g_archr$Participants)
p2g_archr$Participants <- gsub("M5","P4",p2g_archr$Participants)

p2g_archr$Timepoints <- p2g_archr$Sample
p2g_archr$Timepoints <- substring(p2g_archr$Timepoints,4)
p2g_archr$Timepoints <- gsub("11","Day171",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("12","Day185",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("10","Day44",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^1$","Day0",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^2$","Day1",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^3$","Day3",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^4$","Day6",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^5$","Day14",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^6$","Day30",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^7$","Day31",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^8$","Day33",p2g_archr$Timepoints)
p2g_archr$Timepoints <- gsub("^9$","Day36",p2g_archr$Timepoints)

p7 <- plotEmbedding(p2g_archr, name = "Participants", embedding = "UMAP_ATAC", pal=b7,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p8 <- plotEmbedding(p2g_archr, name = "Participants", embedding = "UMAP_RNA", pal=b7,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p9 <- plotEmbedding(p2g_archr, name = "Participants", embedding = "UMAP_Combined",pal=b7, colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
plotPDF(p7, p8, p9, name = "6.14_UMAP-scATAC-scRNA-Combined_Participants", addDOC = FALSE)

p10 <- plotEmbedding(p2g_archr, name = "Timepoints", embedding = "UMAP_ATAC", pal=c7,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p11 <- plotEmbedding(p2g_archr, name = "Timepoints", embedding = "UMAP_RNA", pal=c7,colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p12 <- plotEmbedding(p2g_archr, name = "Timepoints", embedding = "UMAP_Combined",pal=c7, colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
plotPDF(p10, p11, p12, name = "6.14_UMAP-scATAC-scRNA-Combined_Timepoints", addDOC = FALSE)

p18 <- plotEmbedding(archrproj, name = "Sample", embedding = "UMAP_ATAC",colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p19 <- plotEmbedding(archrproj, name = "Sample", embedding = "UMAP_RNA",colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
p20 <- plotEmbedding(archrproj, name = "Sample", embedding = "UMAP_Combined", colorBy = "cellColData",size = 1, labelAsFactors=F, labelMeans=F)
plotPDF(p18, p19, p20, name = "6.30_UMAP-scATAC-scRNA-Combined_sample", addDOC = FALSE)
##TSSEnrichment
p13 <- plotGroups(ArchRProj = p2g_archr, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "ridges")
p14 <- plotGroups(ArchRProj = p2g_archr, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
##nFrags
p15 <- plotGroups(ArchRProj = p2g_archr, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
##Fragment size distributions 
p16 <- plotFragmentSizes(ArchRProj = p2g_archr)
p17 <- plotTSSEnrichment(ArchRProj = p2g_archr)


table(p2g_archr$celltype)
##finding and Visualizing Marker Genes on an Embedding
markersGE <- getMarkerFeatures(ArchRProj = archrproj, useMatrix = "GeneExpressionMatrix", groupBy = "celltype",testMethod = "wilcoxon")
markerList <- getMarkers(markersGE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
pbmc.marker <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/hsPBMC_markers_v3--wangrong.xlsx")
heatmapGE <- markerHeatmap(seMarker = markersGE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = unique(pbmc.marker$Gene,transpose = TRUE))
ComplexHeatmap::draw(heatmapGE, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGE, name = "6.14_Expression-Marker-Heatmap", width = 8, height = 10, ArchRProj = archrproj, addDOC = TRUE)

markersPeaks <- getMarkerFeatures(ArchRProj = archrproj, useMatrix = "PeakMatrix", groupBy = "celltype",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList_peak <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
plotPDF(heatmapPeaks, name = "6.14_Peak-Marker-Heatmap", width = 8, height = 10, ArchRProj = archrproj, addDOC = TRUE)

##Motif and Feature Enrichment 
motif_archr<- addMotifAnnotations(ArchRProj = archrproj, motifSet = "cisbp", name = "Motif",force = TRUE)
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = motif_archr,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "6.14_Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = motif_archr, addDOC = FALSE)
##ChromVAR Deviatons Enrichment 
motif_archr <- addBgdPeaks(motif_archr)
motif_archr <- addDeviationsMatrix(ArchRProj = motif_archr, peakAnnotation = "Motif",force = TRUE)
plotVarDev <- getVarDeviations(motif_archr, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = motif_archr, addDOC = FALSE)


##Co-accessibility
motif_archr <- addCoAccessibility(ArchRProj = motif_archr,reducedDims = "LSI_Combined")
cA <- getCoAccessibility(ArchRProj = motif_archr ,corCutOff = 0.5, resolution = 1,returnLoops = TRUE)
saveRDS(motif_archr,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.30.Archr_down.RDS",compress = F)
save(motif_archr,markersGE , markersPeaks,archrproj,motifPositions,cA, file = "/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.14.archr_motif.RData")
archrproj <- addPeak2GeneLinks(ArchRProj = archrproj,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
archrproj$celltype <- factor(archrproj$celltype,levels = level)
p2g <- getPeak2GeneLinks(ArchRProj = archrproj,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
metadata(p2g)[[1]]

saveRDS(archrproj,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.25.Peak2Gene.RDS",compress = F)
pd <- paletteDiscrete(archrproj$celltype)
c13 <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A", "#33A02C" ,"#FB9A99" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A", "#FFFF99", "#B15928", "#66C2A5" ,"#FC8D62")
names(c13) <- names(pd)
plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "celltype",palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2gmatrix <- plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "celltype",returnMatrices = TRUE,nPlot = 70000,palRNA= paletteer_c("grDevices::TealRose", 30))
p2gmatrix_byp <- plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "Participants",returnMatrices = TRUE,nPlot = 50000,palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "Participants",nPlot = 50000,palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))

p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
p2g.df.null <- as.data.frame(p2g)
hist(p2g.df.null$Correlation,col = "lightblue",main = "Histogram of peak-to-gene correlations",xlab = "Correlation")
hist(p2g.df.null$FDR,col="lightblue",main = "Histogram peak-to-gene FDR",xlab = "FDR")
hist(table(p2g.df.null$idxRNA),main="Distribution of peaks per gene",xlab = "Number of peaks per gene")
hist(table(p2g.df.null$idxATAC),main="Distribution of gene per peaks",xlab = "Number of gene per peaks")
mean(table(p2g.df.null$idxRNA))
mean(table(p2g.df.null$idxATAC))

p2g.df.null$kmeans <- p2gmatrix$RNA$kmeansId

##celltype specific p2g
p2g.df.null$group <- p2g.df.null$kmeans
myeliod_B <- C(2,3,4,5)
myeliod <-c(1,6,7,8,9,10)
B <- c(11,12,14,16)
B_CD4 <- c(13)
B_CD8 <- c(15)
CD4T <- c(21,22,23)
CD8T <- c(17,18,19,20)
T <- c(24,25)
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^24$","^25$"),to = "T")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^17$","^18$","^19$","^20$"),to = "CD8T")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^21$","^22$","^23$"),to = "CD4T")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^15$"),to = "B_CD8")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^13$"),to = "B_CD4")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^11$","^12$","^14$","^16$"),to = "B")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^1$","^6$","^7$","^8$","^9$","^10$"),to = "myeliod")
p2g.df.null$group <- Replace(data=p2g.df.null$group,from = c("^2$","^3$","^4$","^5$"),to = "myeliod_B")

sig.mashResults.gene <- list()
for( i in 1:length(clusters)){
  print(clusters[i])
  sig.mashResults.gene[[i]] <- read.table(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/marsh/6.24_timepoints_',clusters[[i]],'.sig.mashResults.gene.top25.csv'),header = T,sep = ",")
  sig.mashResults.gene[[i]]$celltype <- rep(clusters[i],nrow(sig.mashResults.gene[[i]]))
}
ALL_sig.gene <- do.call(rbind,sig.mashResults.gene)
rownames(ALL_sig.gene) <- ALL_sig.gene$X
sig.gene <- ALL_sig.gene[which(ALL_sig.gene$X %in% intersect(unique(p2g.df.null$geneName),unique(ALL_sig.gene$X))),]
sig.gene.p2g <- p2g.df.null[which(p2g.df.null$geneName %in% intersect(unique(p2g.df.null$geneName),unique(ALL_sig.gene$X))),]
length(unique(sig.gene.p2g$geneName))
sig.gene.p2g<-arrange(sig.gene.p2g,sig.gene.p2g$kmeans)


p <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = unique(sig.gene.p2g[1:1087,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_myeliod.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_1 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1:1087,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_1, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_myeliod.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_1 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1:1087,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_1, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_myeliod.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

p_c_2 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = unique(sig.gene.p2g[1088:1332,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_c_2, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_B.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_2 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1088:1332,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_2, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_B.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

p_c_3 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = unique(sig.gene.p2g[1333:1590,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_c_3, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_CD8.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_3 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1333:1590,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_3, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_CD8.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_c_4 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = unique(sig.gene.p2g[1591:1744,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_c_4, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_CD4.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_4 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1591:1744,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_4, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_CD4.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

p_c_5 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = unique(sig.gene.p2g[1745:1811,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_c_5, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_T.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
p_t_5 <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = unique(sig.gene.p2g[1745:1811,]$geneName), upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = p_t_5, 
        name = "6.28.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-timepoints_T.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
save(p,p_t_1,p_c_2,p_t_2,p_c_3,p_t_3,p_c_4,p_t_4,p_c_5,p_t_5,file="/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.29.plotBrowserTrack.RData")

#Identification of Positive TF-Regulators
seGroupMotif <- getGroupSE(ArchRProj = motif_archr, useMatrix = "MotifMatrix", groupBy = "celltype")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(ArchRProj = motif_archr,useMatrix1 = "GeneExpressionMatrix",useMatrix2 = "MotifMatrix",reducedDims = "LSI_Combined")
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05))+
ggrepel::geom_text_repel(aes(label=GeneExpressionMatrix_name,color=TFRegulator),data.frame(corGSM_MM),
  size = 4, #注释文本的字体大小
  box.padding = 0.5, #字到点的距离
  point.padding = 0.8, #字到点的距离，点周围的空白宽度
  min.segment.length = 0.5, #短线段可以省略
  segment.color = "black", #显示线段
  show.legend = F)
time_markersPeaks <- getMarkerFeatures(
  ArchRProj = motif_archr, useMatrix = "PeakMatrix", groupBy = "Timepoints",bias = c("TSSEnrichment", "log10(nFrags)","Sample"),testMethod = "wilcoxon")
save(corGSM_MM,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/7.12.Positive_TF-Regulators.RData")

positive_tf <-plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = 'SPIB', upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = positive_tf, 
        name = "7.19.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_ptf.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
positive_tf_1 <-plotBrowserTrack(ArchRProj = archrproj, groupBy = "Timepoints", geneSymbol = 'SPIB', upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(archrproj))
plotPDF(plotList = positive_tf_1, 
        name = "7.19.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-celltype_ptf_1.pdf", 
        ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)
##Motif Footprinting
motifPositions <- getPositions(archrproj)
motifs <- c('SPIB', 'FOSL2', 'FOSL1', 'API1', 'BCL11AN', 'TCF4', 'JDP2', 'ZEB1', 'POU2F2', 'EBF1', 'ETS1' , 'LWF1','NR4A1','EOMES','LTS1','LEF1','TCF7')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
seFoot <- getFootprints(ArchRProj = archrproj, positions = motifPositions[markerMotifs], groupBy = "celltype")
plotFootprints(seFoot = seFoot,ArchRProj = archrproj, 
  normMethod = "Subtract",plotName = "Footprints-Subtract-Bias",addDOC = FALSE,smoothWindow = 5)
archrproj_t <- addGroupCoverages(ArchRProj = archrproj, groupBy = "Timepoints")
seFoot_t <- getFootprints(ArchRProj = archrproj_t, positions = motifPositions[markerMotifs], groupBy = "Timepoints")
plotFootprints(seFoot = seFoot_t,ArchRProj = archrproj_t, 
               normMethod = "Subtract",plotName = "Footprints-Subtract-Bias",addDOC = FALSE,smoothWindow = 5)

##participants specific p2g
library(ChIPpeakAnno)
library(stats)

P1_cells<-archrproj@cellColData[archrproj$Participants=="P1",]
P1_cells <- as.data.frame(P1_cells)
P1_cells <- rownames(P1_cells)
proj_P1 <- subsetArchRProject(motif_archr,cells = P1_cells,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/ArchR/P1")
proj_P1 <- addPeak2GeneLinks(ArchRProj = proj_P1,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
plotPeak2GeneHeatmap(ArchRProj = proj_P1, groupBy = "celltype",palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p1 <- getPeak2GeneLinks(ArchRProj = proj_P1,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_p1$geneName <- mcols(metadata(p2g_p1)$geneSet)$name[p2g_p1$idxRNA]
p2g_p1$peakName <- (metadata(p2g_p1)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_p1$idxATAC]
p2g_p1.df <- as.data.frame(p2g_p1)
p2g.df.null$paired_id <- paste0(p2g.df.null$peakName,p2g.df.null$geneName)
p2g_p1.df$paired_id <- paste0(p2g_p1.df$peakName,p2g_p1.df$geneName)
length(intersect(p2g_p1.df$paired_id,p2g.df.null$paired_id))
#31230
P2_cells<-archrproj@cellColData[archrproj$Participants=="P2",]
P2_cells <- as.data.frame(P2_cells)
P2_cells <- rownames(P2_cells)
proj_P2 <- subsetArchRProject(motif_archr,cells = P2_cells,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/ArchR/P2")
proj_P2 <- addPeak2GeneLinks(ArchRProj = proj_P2,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
plotPeak2GeneHeatmap(ArchRProj = proj_P2, groupBy = "celltype",palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p2 <- getPeak2GeneLinks(ArchRProj = proj_P2,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_p2$geneName <- mcols(metadata(p2g_p2)$geneSet)$name[p2g_p2$idxRNA]
p2g_p2$peakName <- (metadata(p2g_p2)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_p2$idxATAC]
p2g_p2.df <- as.data.frame(p2g_p2)
p2g_p2.df$paired_id <- paste0(p2g_p2.df$peakName,p2g_p2.df$geneName)
length(intersect(p2g_p2.df$paired_id,p2g.df.null$paired_id))
#33255
P3_cells<-archrproj@cellColData[archrproj$Participants=="P3",]
P3_cells <- as.data.frame(P3_cells)
P3_cells <- rownames(P3_cells)
proj_P3 <- subsetArchRProject(motif_archr,cells = P3_cells,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/ArchR/P3",force = TRUE)
proj_P3 <- addPeak2GeneLinks(ArchRProj = proj_P3,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
plotPeak2GeneHeatmap(ArchRProj = proj_P3, groupBy = "celltype",palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p3 <- getPeak2GeneLinks(ArchRProj = proj_P3,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_p3$geneName <- mcols(metadata(p2g_p3)$geneSet)$name[p2g_p3$idxRNA]
p2g_p3$peakName <- (metadata(p2g_p3)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_p3$idxATAC]
p2g_p3.df <- as.data.frame(p2g_p3)
p2g_p3.df$paired_id <- paste0(p2g_p3.df$peakName,p2g_p3.df$geneName)
length(intersect(p2g_p3.df$paired_id,p2g.df.null$paired_id))
#31461
P4_cells<-archrproj@cellColData[archrproj$Participants=="P4",]
P4_cells <- as.data.frame(P4_cells)
P4_cells <- rownames(P4_cells)
proj_P4 <- subsetArchRProject(motif_archr,cells = P4_cells,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/ArchR/P4",force = TRUE)
proj_P4 <- addPeak2GeneLinks(ArchRProj = proj_P4,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
plotPeak2GeneHeatmap(ArchRProj = proj_P4, groupBy = "celltype",palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p4 <- getPeak2GeneLinks(ArchRProj = proj_P4,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_p4$geneName <- mcols(metadata(p2g_p4)$geneSet)$name[p2g_p4$idxRNA]
p2g_p4$peakName <- (metadata(p2g_p4)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_p4$idxATAC]
p2g_p4.df <- as.data.frame(p2g_p4)
p2g_p4.df$paired_id <- paste0(p2g_p4.df$peakName,p2g_p4.df$geneName)
length(intersect(p2g_p4.df$paired_id,p2g.df.null$paired_id))
#38068
#韦恩图
library(VennDiagram)
library(RColorBrewer)
venn_ploy <-venn.diagram(x = list(all_Participants = p2g.df.null$paired_id,P1 = p2g_p1.df$paired_id,P2 = p2g_p2.df$paired_id,P3 = p2g_p3.df$paired_id,P4 = p2g_p4.df$paired_id),filename = "/home/wangrong/results/6.20.P2g.venn.tiff",fill = brewer.pal(5, "Pastel2"))
save(proj_P1,proj_P2,proj_P3,proj_P4,archrproj,file = "/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.20.archr_p2g.RData")

##top p2g in intersect group
co_p2g<- Reduce(intersect, list(p2g.df.null$paired_id,p2g_p1.df$paired_id,p2g_p2.df$paired_id,p2g_p3.df$paired_id,p2g_p4.df$paired_id))
co_p2g_df <- p2g.df.null[!is.na(match(p2g.df.null$paired_id,co_p2g)),]
number_p <- as.data.frame(table(co_p2g_df$geneName))
mean(table(co_p2g_df$idxRNA))#6.584461
number_p <- number_p[number_p$Freq>=7,]
co_p2g_df_filter <- co_p2g_df[which(co_p2g_df$geneName %in% number_p$Var1),]
co_p2g_df_arrange <- arrange(co_p2g_df_filter,desc(co_p2g_df_filter$Correlation))
write_csv(co_p2g_df_arrange,file="/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.20_co_p2g_df_arrange.csv")
co_GO <- bitr(unique(co_p2g_df_arrange$geneName[1:186]), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db")
GO_df <- enrichGO(co_GO$ENTREZID,OrgDb = org.Hs.eg.db,ont = "all",pAdjustMethod = "none", pvalueCutoff = 0.05,readable =T)
-log(p.adjust)
dotplot(GO_df, showCategory = 15)+scale_color_gradient(low = '#008B8B', high ='#FFF8DC' )
GO_data <- GO_df@result[1:15,]

ggplot(GO_data,aes(x = Count, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_colour_gradient(low='#FFF8DC',high='#008B8B')+
  labs(color=expression(-log10(p.adjust)),
    size=" Count Number",x="Gene Count")+theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))

##散点图 number of enhance 和correlation
mc <- tapply(co_p2g_df_filter$Correlation,co_p2g_df_filter$geneName, median)
np <- as.data.frame(table(co_p2g_df_filter$geneName))
np$Correlation <- mc
mean(np$Correlation)
#[1] 0.6709494
np1 <- np[np$Correlation>0.8,]
np2 <- np1[np1$Freq>10,]
np2$group <- "sig"
np2$gene <- np2$Var1
np_g <- merge(np,np2,all.x = TRUE)
np_g[is.na(np_g)] <-  "non_sig"


ggplot(np_g, aes(x = Freq, y = Correlation)) +geom_point(aes(colour = group,))+
       scale_color_manual(values=c("non_sig" = "#80B1D3","sig" = "#FB8072"))+
       geom_vline(xintercept=10,lty=2,col="black",lwd=1) +
       theme(legend.background=element_blank(), legend.key=element_blank(),
                         legend.title = element_blank(),
                         panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme_bw()+
       ggrepel::geom_text_repel(
             aes(label=gene,color=group),np_g,
             size = 4, #注释文本的字体大小
             box.padding = 0.5, #字到点的距离
             point.padding = 0.8, #字到点的距离，点周围的空白宽度
             min.segment.length = 0.5, #短线段可以省略
             segment.color = "black", #显示线段
            show.legend = F)
top_gene <- c("PID1","VCAN","IL1B","IL1A","RBM47","FCAR","TLR2","CLEC7A","SLC8A1")
p <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "celltype", geneSymbol = top_gene, upstream = 50000,downstream = 50000,
  loops = getPeak2GeneLinks(archrproj))
grid::grid.draw(p$IL1B)

p_p <- plotBrowserTrack(ArchRProj = archrproj, groupBy = "Participants", geneSymbol = top_gene, upstream = 50000,downstream = 50000,
                        loops = getPeak2GeneLinks(archrproj))

grid::grid.draw(p_p$IL1B)                                                
##P2_specific p2g

P2_specific_p2g <- Reduce(setdiff, list(p2g_p2.df$paired_id,p2g.df.null$paired_id,p2g_p1.df$paired_id,p2g_p3.df$paired_id,p2g_p4.df$paired_id))
p2gmatrix_p2 <- plotPeak2GeneHeatmap(ArchRProj = proj_P2, groupBy = "celltype",returnMatrices = TRUE,nPlot = 70000,palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p2.df$kmeansId <- p2gmatrix_p2$RNA$kmeansId
p2_p2g_df <- p2g_p2.df[!is.na(match(p2g_p2.df$paired_id,P2_specific_p2g )),]
#number_p <- as.data.frame(table(p2_p2g_df$geneName))
mean(table(p2_p2g_df$idxRNA))#2.272436
mean(table(p2_p2g_df$idxATAC))
#number_p <- number_p[number_p$Freq>=7,]
#co_p2g_df_filter <- co_p2g_df[which(co_p2g_df$geneName %in% number_p$Var1),]
p2_p2g_df_arrange <- arrange(p2_p2g_df,desc(p2_p2g_df$Correlation))
#write_csv(co_p2g_df_arrange,file="/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.20_co_p2g_df_arrange.csv")
p2_GO <- bitr(unique(p2_p2g_df_arrange$geneName[1:127]), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db")
p2_GO_df <- enrichGO(p2_GO$ENTREZID,OrgDb = org.Hs.eg.db,ont = "all",pAdjustMethod = "none", pvalueCutoff = 0.05,readable =T)
p2_GO_data <- p2_GO_df@result[1:15,]

ggplot(p2_GO_data,aes(x = Count, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_colour_gradient(low='#FFF8DC',high='#008B8B')+
  labs(color=expression(-log10(p.adjust)),
       size=" Count Number",x="Gene Count")+theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))
mc_p2 <- tapply(p2_p2g_df$Correlation,p2_p2g_df$geneName, median)
np_p2 <- as.data.frame(table(p2_p2g_df$geneName))
np_p2$Correlation <- mc_p2
mean(np_p2$Correlation)
#[1] 0.5035295
np1_p2 <- np_p2[np_p2$Correlation>0.5,]
np2_p2 <- np1_p2[np1_p2$Freq>3,]
np2_p2$group <- "sig"
np2_p2$gene <- np2_p2$Var1
np_g_p2 <- merge(np_p2,np2_p2,all.x = TRUE)
np_g_p2[is.na(np_g_p2)] <-  "non_sig"


ggplot(np_g_p2, aes(x = Freq, y = Correlation)) +geom_point(aes(colour = group,))+
  scale_color_manual(values=c("non_sig" = "#80B1D3","sig" = "#FB8072"))+
  geom_vline(xintercept=10,lty=2,col="black",lwd=1) +
  theme(legend.background=element_blank(), legend.key=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme_bw()+
  ggrepel::geom_text_repel(
    aes(label=gene,color=group),np_g_p2,
    size = 4, #注释文本的字体大小
    box.padding = 0.5, #字到点的距离
    point.padding = 0.8, #字到点的距离，点周围的空白宽度
    min.segment.length = 0.5, #短线段可以省略
    segment.color = "black", #显示线段
    show.legend = F)



##P3_specific p2g

P3_specific_p2g <- Reduce(setdiff, list(p2g_p3.df$paired_id,p2g.df.null$paired_id,p2g_p1.df$paired_id,p2g_p2.df$paired_id,p2g_p4.df$paired_id))
p2gmatrix_p3 <- plotPeak2GeneHeatmap(ArchRProj = proj_P3, groupBy = "celltype",returnMatrices = TRUE,nPlot = 70000,palGroup = c13,palRNA= paletteer_c("grDevices::TealRose", 30))
p2g_p3.df$kmeansId <- p2gmatrix_p3$RNA$kmeansId
p3_p2g_df <- p2g_p3.df[!is.na(match(p2g_p3.df$paired_id,P3_specific_p2g )),]
#number_p <- as.data.frame(table(p2_p2g_df$geneName))
mean(table(p3_p2g_df$idxRNA))
mean(table(p3_p2g_df$idxATAC))
#number_p <- number_p[number_p$Freq>=7,]
#co_p2g_df_filter <- co_p2g_df[which(co_p2g_df$geneName %in% number_p$Var1),]
p3_p2g_df_arrange <- arrange(p3_p2g_df,desc(p3_p2g_df$Correlation))
#write_csv(co_p2g_df_arrange,file="/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.20_co_p2g_df_arrange.csv")
p3_GO <- bitr(unique(p3_p2g_df_arrange$geneName[1:122]), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db")
p3_GO_df <- enrichGO(p3_GO$ENTREZID,OrgDb = org.Hs.eg.db,ont = "all",pAdjustMethod = "none", pvalueCutoff = 0.05,readable =T)
p3_GO_data <- p3_GO_df@result[1:15,]

ggplot(p3_GO_data,aes(x = Count, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_colour_gradient(low='#FFF8DC',high='#008B8B')+
  labs(color=expression(-log10(p.adjust)),
       size=" Count Number",x="Gene Count")+theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))

mc_p3 <- tapply(p3_p2g_df$Correlation,p3_p2g_df$geneName, median)
np_p3 <- as.data.frame(table(p3_p2g_df$geneName))
np_p3$Correlation <- mc_p3
mean(np_p3$Correlation)
#0.5008309
np1_p3 <- np_p3[np_p3$Correlation>0.5,]
np2_p3 <- np1_p3[np1_p3$Freq>3,]
np2_p3$group <- "sig"
np2_p3$gene <- np2_p3$Var1
np_g_p3 <- merge(np_p3,np2_p3,all.x = TRUE)
np_g_p3[is.na(np_g_p3)] <-  "non_sig"

ggplot(np_g_p3, aes(x = Freq, y = Correlation)) +geom_point(aes(colour = group,))+
  scale_color_manual(values=c("non_sig" = "#80B1D3","sig" = "#FB8072"))+
  geom_vline(xintercept=10,lty=2,col="black",lwd=1) +
  theme(legend.background=element_blank(), legend.key=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme_bw()+
  ggrepel::geom_text_repel(
    aes(label=gene,color=group),np_g_p3,
    size = 4, #注释文本的字体大小
    box.padding = 0.5, #字到点的距离
    point.padding = 0.8, #字到点的距离，点周围的空白宽度
    min.segment.length = 0.5, #短线段可以省略
    segment.color = "black", #显示线段
    show.legend = F)


signac_pbmc <- readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.13.signac_peak_to_gene.rds")
VlnPlot(signac_pbmc, features = "IL1B",pt.size = 0 , group.by = "Participants")
VlnPlot(signac_pbmc, features = "IL1B",pt.size = 0 , group.by = "celltype")

cd14mono_sub <- subset(signac_pbmc,idents=c("CD14+ Mono"))
VlnPlot(cd14mono_sub, features = "IL1B",pt.size = 0 , group.by = "Participants")

NK_sub<- subset(signac_pbmc,idents=c("NK/NKT"))
Treg_sub<- subset(signac_pbmc,idents=c("Treg"))
Treg_p_marker <- FindMarkers(Treg_sub, ident.1 = "P2",ident.2 = "P3" ,group.by = "Participants")
Treg_p_GO <- bitr(rownames(Treg_p_marker), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db")
Treg_p_GO_df <- enrichGO(Treg_p_GO$ENTREZID,OrgDb = org.Hs.eg.db,ont = "all",pAdjustMethod = "none", pvalueCutoff = 0.05,readable =T)
Treg_p_GO_data <- Treg_p_GO_df@result[1:15,]

ggplot(Treg_p_GO_data,aes(x = Count, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_colour_gradient(low='#FFF8DC',high='#008B8B')+
  labs(color=expression(-log10(p.adjust)),
       size=" Count Number",x="Gene Count")+theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))
Treg_p_KEGG_df <- enrichKEGG(Treg_p_GO$ENTREZID,organism = "hsa",keyType = "kegg", pAdjustMethod = "none", pvalueCutoff = 0.05)
Treg_p_KEGG_data <- Treg_p_KEGG_df@result[1:15,]
ggplot(Treg_p_KEGG_data,aes(x = Count, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_colour_gradient(low='#FFF8DC',high='#008B8B')+
  labs(color=expression(-log10(p.adjust)),
       size=" Count Number",x="Gene Count")+theme_bw()+
  theme(axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_blank()
  )+ scale_size(range=c(5, 10))

markersPeaks_p <- getMarkerFeatures(ArchRProj = archrproj, useMatrix = "PeakMatrix", groupBy = "celltype_participants",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList_peak_p <- getMarkers(markersPeaks_p, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heat_mapPeaks_p <- markerHeatmap(seMarker = markersPeaks_p, cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE,labelRows = FALSE)
heatmapPeaks_p <- markerHeatmap(seMarker = markersPeaks_p, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",returnMatrix = TRUE,labelRows = FALSE)
plotPDF(heat_mapPeaks_p, name = "7.11_Peak-Marker-Heatmap", width = 10, height = 10, ArchRProj = archrproj, addDOC = TRUE)

markersPeaks_p_Treg <- getMarkerFeatures(ArchRProj = archrproj, useMatrix = "PeakMatrix", groupBy = "celltype_participants",useGroups = c("Treg_P1","Treg_P2","Treg_P3","Treg_P4"),bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
#markerList_peak_p <- getMarkers(markersPeaks_p, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heatmapPeaks_p_Treg <- markerHeatmap(seMarker = markersPeaks_p_Treg, cutOff = "FDR <= 0.01 & Log2FC >= 0.5",labelRows = FALSE)
plotPDF(heatmapPeaks_p_Treg, name = "7.11_Peak-Marker-Heatmap-Treg", width = 10, height = 10, ArchRProj = archrproj, addDOC = TRUE)

m_cells<-archrproj@cellColData[archrproj$celltype==c("CD14+ Mono","CD16+ Mono"),]
m_cells <- as.data.frame(m_cells)
m_cells <- rownames(m_cells)
proj_m <- subsetArchRProject(archrproj,cells = m_cells,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/ArchR/m_cells",force = TRUE)
proj_m_p <- plotBrowserTrack(ArchRProj = proj_m, groupBy = "Participants", geneSymbol = 'IL1B', upstream = 50000,downstream = 50000,loops = getPeak2GeneLinks(proj_m))
plotPDF(plotList = proj_m_p, 
        name = "7.20.Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-m-p-IL1B-3.pdf", 
        ArchRProj = proj_p2g_m, addDOC = FALSE, width = 5, height = 5)

p2g_m <- getPeak2GeneLinks(ArchRProj = proj_m,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_m$geneName <- mcols(metadata(p2g_m)$geneSet)$name[p2g_m$idxRNA]
p2g_m$peakName <- (metadata(p2g_m)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_m$idxATAC]
p2g_m.df <- as.data.frame(p2g_m)
markersPeaks_m <- getMarkerFeatures(ArchRProj = proj_m, useMatrix = "PeakMatrix", groupBy = "Participants",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList_peak_m <- getMarkers(markersPeaks_m, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heatmapPeaks_m <- markerHeatmap(seMarker = markersPeaks_m, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
plotPDF(heatmapPeaks_m, name = "7.21_Peak-Marker-Heatmap", width = 8, height = 10, ArchRProj = proj_m, addDOC = TRUE)
