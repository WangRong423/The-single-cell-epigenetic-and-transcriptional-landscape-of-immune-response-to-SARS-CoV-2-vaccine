library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

batch_list=list("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
"M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
"M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
"M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")
gex<-list()
for( i in 1:length(batch_list))
{
print(batch_list[[i]])
dir=paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/filtered_feature_bc_matrix')
counts=Read10X(data.dir = dir)
metadata <- read.csv(file = paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/per_barcode_metrics.csv'),header = TRUE,row.names = 1)
fragment.path <- paste0('/database/wangrong/Results/0712_ATAC+RNA/',batch_list[[i]],'/outs/atac_fragments.tsv.gz')
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'

pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragment.path,annotation = annotation)
DefaultAssay(pbmc) <- "ATAC"

gex[[i]] = readRDS(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/HIPPO/hippo/',batch_list[[i]],'_RemoveDoublet.rds'))
print("gex is down!")
pbmc <- subset(pbmc,cells=rownames(gex[[i]]@meta.data))
print("RNA subset is doen!")
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

pbmc<-subset(x = pbmc,subset = nCount_ATAC < 100000 &nCount_ATAC > 1000 &nucleosome_signal < 5 &TSS.enrichment > 1)
sub_gex <- subset(gex[[i]],cells=rownames(pbmc@meta.data))
pbmc@meta.data <- cbind(sub_gex@meta.data,pbmc@meta.data)
saveRDS(pbmc,paste0('/database/wangrong/Results/0712_ATAC+RNA/Signac/co_joint/',batch_list[[i]],'_cojoint.rds'))
}


#2022.5.14
cell_name_ids=C("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
                "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
                "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
                "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")
batch_list=list("M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
                "M2-1","M2-2","M2-3","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10","M2-11","M2-12",
                "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
                "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10")
all_p=list()           
for( i in 1:length(batch_list))
{
  all_p[[i]] = readRDS(file=paste0('/database/wangrong/Results/0712_ATAC+RNA/Signac/co_joint/',batch_list[[i]],'_cojoint.rds'))
  DefaultAssay(all_p[[i]]) <- "ATAC"
}
PBMCmerge <- merge(all_p[[1]],all_p[2:length(batch_list)],add.cell.ids = cell_name_ids)
saveRDS(PBMCmerge,file="/database/wangrong/Results/0712_ATAC+RNA/Signac/co_joint/5.14_merge_cojoint.rds")