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

batch_list=list("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
                "M2-1","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
                "M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
                "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")
Timepoints=list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44",
                "Day0","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185",
                "Day6","Day14","Day30","Day31","Day33","Day36","Day44",
                "Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44")
Participants=list("P1","P1","P1","P1","P1","P1","P1","P1","P1","P1",
                  "P2","P2","P2","P2","P2","P2","P2","P2","P2","P2","P2",
                  "P3","P3","P3","P3","P3","P3","P3",
                  "P4","P4","P4","P4","P4","P4","P4","P4","P4","P4")
for( i in 1:length(batch_list))
{
print(batch_list[[i]])
dir=paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/filtered_feature_bc_matrix')
print(dir)
s_object_data=Read10X(data.dir = dir)

s_object=CreateSeuratObject(counts =s_object_data$`Gene Expression`, project = batch_list[[i]])
s_object[["percent.mt"]] <- PercentageFeatureSet(s_object, pattern = "^MT-")
print(dim(s_object))
a=summary(s_object$nFeature_RNA)
b=summary(s_object$nCount_RNA)
c=summary(s_object$percent.mt)
print(a)
print(b)
print(c)
#QC percent.mt
s_object <- subset(s_object, subset =  percent.mt <10)
print(dim(s_object))
a1=summary(s_object$nFeature_RNA)
b1=summary(s_object$nCount_RNA)
c1=summary(s_object$percent.mt)
print(a1)
print(b1)
print(c1)
#Normalize
s_object <- NormalizeData(s_object, normalization.method = "LogNormalize", scale.factor = 10000)
#FindVariableFeatures
s_object <- FindVariableFeatures(s_object, selection.method = "vst", nfeatures = 2000)
#Scale
s_object <- ScaleData(s_object, features = rownames(s_object))
#PCA
s_object <- RunPCA(s_object, features = VariableFeatures(s_object),npcs = 20)
#cluster
name = paste0('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/',batch_list[[i]], '_bef_double_umap.pdf')
pdf(file = name)
s_object <- FindNeighbors(s_object, dims = 1:20)
s_object <- FindClusters(s_object, resolution = 0.3)
s_object <- RunUMAP(s_object, dims = 1:20)
p1<-DimPlot(s_object,reduction = "umap",pt.size = 0.5,label = TRUE)

location="/database/wangrong/Results/0712_ATAC+RNA/HIPPO/"
filename=paste0(batch_list[[i]],'_PBMC_example')
newFiles=Improved_Seurat_Pre_Process(s_object, num_genes=50, write_files=FALSE)

results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename=filename, 
                           location=location,
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="hsa", 
                           rhop=1, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=4,
                           nCores=6)
dou=length(row.names(results$Final_doublets_groups))
sing=length(row.names(results$Final_nondoublets_groups))
p=dou/(sing+dou)
print(p)

DRS_doublet_table=results$DRS_doublet_table
addmargins(table(DRS_doublet_table$isADoublet))
DRS_doublet_table$Dlabel=factor(DRS_doublet_table$isADoublet,c(FALSE,TRUE),ordered = T, labels=c("Singlet","Doublet"))
seuratObject=AddMetaData(s_object,DRS_doublet_table$Dlabel, col.name = "Dlabel")
Idents(seuratObject)=seuratObject@meta.data$Dlabel
dp1=DimPlot(seuratObject,reduction="umap")
dp1=dp1 +
  theme(legend.position="top") +
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.text = element_text(size=9))
print(dp1)
print(p1)
dev.off()
saveRDS(seuratObject,file=paste0('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/',batch_list[[i]],'_BeforeDoublet.rds'))
seuratObject@meta.data$Timepoints <- Timepoints[[i]]
seuratObject@meta.data$Participants <- Participants[[i]]
seuratObject_RemoveDoublet<-subset(x = seuratObject, Dlabel=="Singlet")
print(dim(seuratObject_RemoveDoublet))
a2=summary(seuratObject_RemoveDoublet$nFeature_RNA)
b2=summary(seuratObject_RemoveDoublet$nCount_RNA)
c2=summary(seuratObject_RemoveDoublet$percent.mt)
print(a2)
print(b2)
print(c2)
saveRDS(seuratObject_RemoveDoublet,file=paste0('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/',batch_list[[i]],'_RemoveDoublet.rds'))

}