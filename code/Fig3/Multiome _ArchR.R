#####
#Multiome scRNA-seq/scATAC-seq analysis with ArchR
#DATE:2022-5-24
######

library(ArchR)
library(pheatmap)
library(Seurat) 
library(ggplot2)
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
metadata<-read.table('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/sample_metadata.csv', sep = ",", header = TRUE)
rownames(metadata) <- gsub("_","#",rownames(metadata))
length(rownames(metadata))
a <- rownames(metadata)
b <- a[!is.na(match(a,getCellNames(proj)))]
proj <- subsetArchRProject(proj,cells = b,outputDirectory = "/database/wangrong/Results/0712_ATAC+RNA/atac_output/ArchRSubset")
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

p1 <- plotEmbedding(proj, name = "celltype", embedding = "UMAP_ATAC", colorBy = "cellColData",size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "celltype", embedding = "UMAP_RNA", colorBy = "cellColData",size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "celltype", embedding = "UMAP_Combined", colorBy = "cellColData",size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)

archrproj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype")
pathToMacs2 <- "/home/wangrong/.local/bin/macs2"
archrproj <- addReproduciblePeakSet(ArchRProj = archrproj,groupBy = "celltype",pathToMacs2 = pathToMacs2)

archrproj <- addPeakMatrix(archrproj)
saveRDS(archrproj,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.24.Archr_call_peak.RDS",compress = F)
archrproj <- addPeak2GeneLinks(ArchRProj = archrproj,reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = archrproj,corCutOff = 0.45,resolution = 1,returnLoops = FALSE)
p2g_list<-as.data.frame(p2g )
write.csv(p2g_list,'/database/wangrong/Results/0712_ATAC+RNA/atac_output/p2g_list.csv')

plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "celltype")
saveRDS(archrproj,"/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.24.Peak2Gene.RDS",compress = F)