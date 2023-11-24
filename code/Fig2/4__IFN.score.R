pbmc<-readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/5.18.annotation.signac.RDS")

pbmc.marker <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/hsPBMC_markers_v3--wangrong.xlsx")
                                                                                                                                   
#View(pbmc.marker)
pbmc.marker<-pbmc.marker[-54,]
pbmc.marker<-pbmc.marker[-53,]
#DotPlot(pbmc, features = unique(pbmc.marker$Gene) ,col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
levels=c("Naïve B"," Intermediate B","Plasmablasts/Memrary B",
          "Naïve CD4+ T cells","CD4+ Tcm","Treg",
          "Naïve CD8+ T cells","CD8+ Tem","MAIT","NK/NKT",
          "CD14+ Mono","CD16+ Mono","Plasmacytoid DC")
Idents(pbmc)<- factor(pbmc$celltype, levels= levels)
RNA_dotplot <- DotPlot(pbmc, features = unique(pbmc.marker$Gene) ,col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

RNA_dotplot_df <- RNA_dotplot$data
A<-RNA_dotplot_df 
ave <- AverageExpression(pbmc,assays ="RNA",slot = "counts",group.by = "ident",features = unique(pbmc.marker$Gene))
ave <- melt(ave)
df<-cbind(A,ave)
df$features.plot<- as.factor(A$features.plot)
#A$features.plot <- fct_inorder(A$features.plot)

#A$features.plot <- fct_inorder(A$features.plot)#处理因子，使排序不变

p <- ggplot(df,aes(x=id,y= features.plot)) + 
  geom_point(aes(size=pct.exp, color=value)) +#散点图
  geom_point(shape=21,aes(size= pct.exp), color="black")+#叠加一个图层，使气泡有外边框
  scale_size(rang = c(1.5,6)) +#点大小范围调整
  labs(x=NULL,y=NULL,size="pct.exp")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),#标题居中显示
        axis.text =element_text(size = 10, color = "black"),#坐标轴字体设置
        axis.text.y = element_text(face="italic"),#y轴字体设置为斜体
        axis.text.x=element_text(angle=90,hjust = 0.5,vjust=0.5))+#x轴字体设置
  scale_color_gradient2(low='navy',high='firebrick3', mid="white", midpoint = 0)+#修改填充颜色
  #geom_hline(yintercept=c(5.5, 13.5))#添加横格线


GWAS_gene<-c("OAS1","OAS2","CCR2","CXCR6","XCR1","FYCO1","IFNAR2","SLC6A20","TYK2","DPP9")
GWAS_dotplot <- DotPlot(pbmc, features = GWAS_gene ,col.min=-2, col.max=2)+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
GWAS_dotplot_df <- GWAS_dotplot$data
AA<-GWAS_dotplot_df 
ave_GWAS <- AverageExpression(pbmc,assays ="RNA",slot = "counts",group.by = "ident",features = GWAS_gene)
ave_GWAS <- melt(ave_GWAS)
df_GWAS<-cbind(AA,ave_GWAS)
df_GWAS$features.plot<- as.factor(df_GWAS$features.plot)
ggplot(df_GWAS,aes(x=id,y= features.plot)) + 
  geom_point(aes(size=pct.exp, color=value)) +#散点图
  geom_point(shape=21,aes(size= pct.exp), color="black")+#叠加一个图层，使气泡有外边框
  scale_size(rang = c(1.5,6)) +#点大小范围调整
  labs(x=NULL,y=NULL,size="pct.exp")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),#标题居中显示
        axis.text =element_text(size = 10, color = "black"),#坐标轴字体设置
        axis.text.y = element_text(face="italic"),#y轴字体设置为斜体
        axis.text.x=element_text(angle=90,hjust = 0.5,vjust=0.5))+#x轴字体设置
  scale_color_gradient2(low='navy',high='firebrick3', mid="white", midpoint = 0)

GWAS_dotplot_t <- DotPlot(pbmc, features = GWAS_gene ,col.min=-2, col.max=2,group.by="Timepoints")+ RotatedAxis()+ scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
GWAS_dotplot_df_t <- GWAS_dotplot_t$data
AAA<-GWAS_dotplot_df_t 
ave_GWAS_t <- AverageExpression(pbmc,assays ="RNA",slot = "counts",group.by = "Timepoints",features = GWAS_gene)
ave_GWAS_t <- melt(ave_GWAS_t)
df_GWAS_t<-cbind(AAA,ave_GWAS_t)
df_GWAS_t$features.plot<- as.factor(df_GWAS_t$features.plot)
df_GWAS_t$id<- factor(df_GWAS_t$id,levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
ggplot(df_GWAS_t,aes(x=id,y= features.plot)) + 
  geom_point(aes(size=pct.exp, color=value)) +#散点图
  geom_point(shape=21,aes(size= pct.exp), color="black")+#叠加一个图层，使气泡有外边框
  scale_size(rang = c(1.5,6)) +#点大小范围调整
  labs(x=NULL,y=NULL,size="pct.exp")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),#标题居中显示
        axis.text =element_text(size = 10, color = "black"),#坐标轴字体设置
        axis.text.y = element_text(face="italic"),#y轴字体设置为斜体
        axis.text.x=element_text(angle=45,hjust = 0.5,vjust=0.5))+#x轴字体设置
  scale_color_gradient2(low='navy',high='firebrick3', mid="white", midpoint = 0)




  # B <- unique(A$features.plot)
# B <- as.data.frame(B)#提取gene并转化为数据框
# gene_group1 <- B$B %>% as.data.frame() %>%
  # mutate(group=rep(c("Immune supression","T cell activation","Immune reponse"),
                   # times=c(5,8,7))) %>% #分组，注释基因
  # mutate(p="")%>%
  # ggplot(aes(p,.,fill=group))+
  # geom_tile() + 
  # scale_y_discrete(position="left") +
  # theme_minimal()+xlab(NULL) + ylab(NULL) +
  # theme(panel.grid = element_blank(),
    # axis.text.y = element_blank(),
        # axis.text.x =element_text(
          # angle =90,hjust =0.5,vjust = 0.5),
        # legend.position = 'none')+
  # scale_fill_manual(values = alpha(c('#852f88',
                               # '#eb990c',
                               # '#0f8096'), 0.5))#替换分组填充颜色，并设置透明度为0.5
  


# gene_group2 <- B$B %>% as.data.frame() %>%
  # mutate(group=rep(c("Immune supression","T cell activation","Immune reponse"),
                   # times=c(5,8,7))) %>%
  # mutate(p="")%>%
  # ggplot(aes(p,.,fill=group))+
  # geom_tile() + 
  # scale_y_discrete(position="left") +
  # theme_minimal()+xlab(NULL) + ylab(NULL) +
  # theme(panel.grid = element_blank(),
        # axis.text.y = element_blank(),
        # axis.text.x =element_text(
          # angle =90,hjust =0.5,vjust = 0.5),
        # legend.position = 'none')+
  # scale_fill_manual(values = alpha(c('#852f88',
                                     # '#eb990c',
                                     # '#0f8096'), 0))+#设置透密度为0
  # annotate("text", label="Immune supression", x=1.2, y=18)+
  # annotate("text", label="T cell activation", x=1.2, y=10)+
  # annotate("text", label="Immune reponse", x=1.2, y=5)#添加分组标题至合适的位置




# left <- ggplotGrob(gene_group1)
# p1 <- p+annotation_custom(left,xmin=-5,xmax=1,ymin=-0.5,ymax=21)

# library(aplot)
# p1%>%insert_left(gene_group2, width = 0.5)

#https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_RESPONSE_TO_TYPE_I_INTERFERON
IFN_genes = ["ABCE1", 	"ADAR", 	"BST2", 	"CACTIN", 	"CDC37", 	"CNOT7", 	"DCST1", 	"EGR1", 	"FADD", 	"GBP2", 

#IFN response（across 时间点和细胞类型）

x = readLines("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -c(1:2))	

timepoints <- list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
subset <-list()
day_aggregate<-list()
for (i in 1:length(timepoints)){
print(timepoints[[i]])
subset[[i]]<- subset(pbmc, Timepoints==timepoints[[i]])
subset[[i]] <- NormalizeData(subset[[i]])
subset[[i]] <- FindVariableFeatures(subset[[i]], selection.method = "vst")
subset[[i]] <- ScaleData(subset[[i]], features = rownames(subset[[i]]))
gene.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)

day_s <- AddModuleScore(subset[[i]],assay="RNA",features=gene.list,nbin = 24,min.cells=1,ctrl = 100,name = "INF_score")

day_aggregate[[i]]<-aggregate(day_s$INF_score1, by=list(type=day_s$celltype),mean)
print(paste0(timepoints[i], " is down"))
}

names(day_aggregate) <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
day_aggregate<-day_aggregate[-c(11,12)]
df <- as.data.frame(day_aggregate)


df1 <- df[,c(2,4,6,8,10,12,14,16,18,20)]/df[,2]
df1$celltype<-df$Day0.type
df1 <- df1[,-1]
df1<-as.data.frame(t(df1))
colnames(df1) <- df1[10,]
df1 <- df1[-10,]
rownames(df1) <- gsub(".x","",rownames(df1))
df2<-melt(df1,value.name = "celltype")

df3 <- df1
df1 <- as.data.frame(lapply(df1,as.numeric))
rownames(df1) <- rownames(df3)
colnames(df1) <- colnames(df3)
#ggplot(df2,aes(x=celltype,y=timepoints))+ xlab("celltype") + ylab("timepoints") + theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank()) + theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + theme(legend.position="top") +  geom_tile(aes(fill=value)) + scale_fill_gradient(low = "white", high = "red")

				
pheatmap(df1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 40,cellheight = 20,legend_labels ="Fold change",annotation_legend=T)	
save(df1,df2,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/ifn_score.Rdata")				

#IFN response（across 时间点和人）
 x = readLines("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -c(1:2))	

timepoints <- list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
subset <-list()
day_aggregate<-list()
for (i in 1:length(timepoints)){
print(timepoints[[i]])
subset[[i]]<- subset(signac_pbmc, Timepoints==timepoints[[i]])
subset[[i]] <- NormalizeData(subset[[i]])
subset[[i]] <- FindVariableFeatures(subset[[i]], selection.method = "vst")
subset[[i]] <- ScaleData(subset[[i]], features = rownames(subset[[i]]))
gene.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)

day_s <- AddModuleScore(subset[[i]],assay="RNA",features=gene.list,nbin = 24,min.cells=1,ctrl = 100,name = "INF_score")

day_aggregate[[i]]<-aggregate(day_s$INF_score1, by=list(type=day_s$Participants),mean)
print(paste0(timepoints[i], " is down"))
}

names(day_aggregate) <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
day_aggregate<-day_aggregate[-c(11,12)]
df <- as.data.frame(day_aggregate)


df1 <- df[,c(2,4,6,8,10,12,14,16,18,20)]/df[,2]
df1$Participants<-df$Day0.type
df1 <- df1[,-1]
df1<-as.data.frame(t(df1))
colnames(df1) <- df1[10,]
df1 <- df1[-10,]
rownames(df1) <- gsub(".x","",rownames(df1))
df2<-melt(df1,value.name = "Participants")

df3 <- df1
df1 <- as.data.frame(lapply(df1,as.numeric))
rownames(df1) <- rownames(df3)
colnames(df1) <- colnames(df3)
			
pheatmap(df1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 40,cellheight = 20,legend_labels ="Fold change",annotation_legend=T)	

save(df1,df2,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/p_t_ifn_score.Rdata")				

#IFN response m（across 时间点和人）
x = readLines("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -c(1:2))	

timepoints <- list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
subset <-list()
day_aggregate<-list()
for (i in 1:length(timepoints)){
print(timepoints[[i]])
subset[[i]]<- subset(atac_m_sub, Timepoints==timepoints[[i]])
subset[[i]] <- NormalizeData(subset[[i]])
subset[[i]] <- FindVariableFeatures(subset[[i]], selection.method = "vst")
subset[[i]] <- ScaleData(subset[[i]], features = rownames(subset[[i]]))
gene.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)

day_s <- AddModuleScore(subset[[i]],assay="RNA",features=gene.list,nbin = 24,min.cells=1,ctrl = 100,name = "INF_score")

day_aggregate[[i]]<-aggregate(day_s$INF_score1, by=list(type=day_s$Participants),mean)
print(paste0(timepoints[i], " is down"))
}

names(day_aggregate) <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
day_aggregate<-day_aggregate[-c(11,12)]
df <- as.data.frame(day_aggregate)


df1 <- df[,c(2,4,6,8,10,12,14,16,18,20)]/df[,2]
df1$Participants<-df$Day0.type
df1 <- df1[,-1]
df1<-as.data.frame(t(df1))
colnames(df1) <- df1[10,]
df1 <- df1[-10,]
rownames(df1) <- gsub(".x","",rownames(df1))
df2<-melt(df1,value.name = "Participants")

df3 <- df1
df1 <- as.data.frame(lapply(df1,as.numeric))
rownames(df1) <- rownames(df3)
colnames(df1) <- colnames(df3)
			
pheatmap(df1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 40,cellheight = 20)	

save(df1,df2,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/7.22_m_ifn_score.Rdata")

#cytokine, chemokine and growth factors（across 时间点和细胞类型）
library(readxl,lib.loc = "/home/shisusanna/R/x86_64-pc-linux-gnu-library/4.2")
cytokine <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/48因子与网络节点信息.xlsx")
colnames(cytokine) <- cytokine[1,]
cytokine <- cytokine[-1,]
length(unique(cytokine$`Gene Name`))

ave_cytokine_c <- AverageExpression(signac_pbmc,assays ="RNA",slot = "counts",group.by = "celltype",features = unique(cytokine$`Gene Name`))
ave_cytokine_c$RNA <- t(ave_cytokine_c$RNA)
pheatmap(ave_cytokine_c$RNA,border="black", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c( "aliceblue", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10,legend_labels ="Fold change",annotation_legend=T)

ave_cytokine <- AverageExpression(signac_pbmc,assays ="RNA",slot = "counts",group.by = "Participants",features = unique(cytokine$`Gene Name`))
ave_cytokine$RNA <- t(ave_cytokine$RNA)
pheatmap(ave_cytokine$RNA,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10,legend_labels ="Fold change",annotation_legend=T)	

ave_cytokine_t <- AverageExpression(rename_pbmc,assays ="RNA",slot = "counts",group.by = "Timepoints",features = unique(cytokine$`Gene Name`))
ave_cytokine_t$RNA <- t(ave_cytokine_t$RNA)
pheatmap(ave_cytokine_t$RNA,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)

rename_pbmc <- readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.13.signac_peak_to_gene.rds")
m_sub <- subset(signac_pbmc,idents=c("CD14+ Mono","CD16+ Mono"))	
ave_cytokine_m <- AverageExpression(m_sub,assays ="RNA",slot = "counts",group.by = "Participants",features = unique(cytokine$`Gene Name`))
ave_cytokine_m$RNA <- t(ave_cytokine_m$RNA)
pheatmap(ave_cytokine_m$RNA,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10,legend_labels ="Fold change",annotation_legend=T)	

ave_cytokine_m_t <- AverageExpression(m_sub,assays ="RNA",slot = "counts",group.by = "Timepoints",features = unique(cytokine$`Gene Name`))
ave_cytokine_m_t$RNA <- t(ave_cytokine_m_t$RNA)
rownames(ave_cytokine_m_t$RNA)<- factor(rownames(ave_cytokine_m_t$RNA),levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
pheatmap(ave_cytokine_m_t$RNA,border="white", cluster_cols = T, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)



#cytokine, chemokine and growth factors（across 时间点和细胞类型）

library(readxl,lib.loc = "/home/shisusanna/R/x86_64-pc-linux-gnu-library/4.2")
cytokine <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/48因子与网络节点信息.xlsx")
colnames(cytokine) <- cytokine[1,]
cytokine <- cytokine[-1,]
length(unique(cytokine$`Gene Name`))

ave_cytokine_c <- AverageExpression(rename_pbmc,assays ="RNA",slot = "counts",group.by = "celltype",features = unique(cytokine$`Gene Name`))
ave_cytokine_c$RNA <- t(ave_cytokine_c$RNA)
pheatmap(ave_cytokine_c$RNA,border="black", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c( "aliceblue", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10,legend_labels ="Fold change",annotation_legend=T)

ave_cytokine <- AverageExpression(rename_pbmc,assays ="RNA",slot = "counts",group.by = "Participants",features = unique(cytokine$`Gene Name`))
ave_cytokine$RNA <- t(ave_cytokine$RNA)
pheatmap(ave_cytokine$RNA,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10,legend_labels ="Fold change",annotation_legend=T)	

ave_cytokine_t <- AverageExpression(rename_pbmc,assays ="RNA",slot = "counts",group.by = "Timepoints",features = unique(cytokine$`Gene Name`))
ave_cytokine_t$RNA <- t(ave_cytokine_t$RNA)
x <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
#ave_cytokine_t1 <- ave_cytokine_t[c(2,3,4,5,6,7,8,9,10,11,12),]/ave_cytokine_t[2,]
ave_cytokine=list()
time <- c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
for (i in 1:length(time)) {
ave_cytokine[[i]]<-FoldChange(rename_pbmc, group.by = 'Timepoints', ident.1="Day0", ident.2=time[i], features = unique(cytokine$`Gene Name`))
}
all_ave_cytokine<-do.call(rbind,ave_cytokine)
all_log2FC_cytokine <- all_ave_cytokine[,c(1,4,7,10,13,16,19,22,25,28,31)]
colnames(all_log2FC_cytokine)<-c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
pheatmap(all_log2FC_cytokine,border="white", cluster_cols = F, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)

rename_pbmc <- readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.13.signac_peak_to_gene.rds")
t_sub <- subset(rename_pbmc,idents=c("CD4+ Tcm","CD8+ Tem","NK/NKT"))	
ave_cytokine_t <- AverageExpression(t_sub,assays ="RNA",slot = "counts",group.by = "Participants",features = unique(cytokine$`Gene Name`))
ave_cytokine_t$RNA <- t(ave_cytokine_t$RNA)
pheatmap(ave_cytokine_t$RNA,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)	

ave_cytokine_m_t <- AverageExpression(m_sub,assays ="RNA",slot = "counts",group.by = "Timepoints",features = unique(cytokine$`Gene Name`))
ave_cytokine_m_t$RNA <- t(ave_cytokine_m_t$RNA)
rownames(ave_cytokine_m_t$RNA)<- factor(rownames(ave_cytokine_m_t$RNA),levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
pheatmap(ave_cytokine_m_t$RNA,border="white", cluster_cols = T, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)


my_comparisons <- list(c("P1", "P2"), 
c("P2", "P3"), 
c("P3", "P4"),
c("P1","P3"),
c("P1","P4"),
c("P2","P4"))

plots_violins <- VlnPlot(m_sub_nor, cols = c("#A6CEE3", "#B2DF8A","Bisque","#FB9A99"),pt.size = 0,group.by = "Participants",features = c('FOSL2','HMBOX1','SREBF2','ZNFR652','JUND','SPI1','ATF4','LIN54','LEF1','MAFF','FOSL1'), log = FALSE)
> plist <- list()
plist <- list()
for(i in 1:length(plots_violins)) {
    data  <- plots_violins[[i]]$data
    colnames(data) <- c('RNA expression', 'Participants')
    p <- ggviolin(data,x="Participants",y="RNA expression",fill="Participants",palette = c("#A6CEE3", "#B2DF8A","Bisque","#FB9A99"),add = "boxplot", add.params = list(fill="white"), legend = NULL)+ stat_compare_means(comparisons = my_comparisons)
    p2 <- ggpar(p,legend =  "none")
	plist[[i]] <- p2
}

CombinePlots(plist, ncol = 4)

plotGroups(ArchRProj = proj_m,groupBy = "Participants",colorBy = "GeneExpressionMatrix",name = c("IL1B","BACH1","FOS","JUN","KLF4","REL",'FOSL2','HMBOX1','SREBF2','ZNFR652','JUND','SPI1','ATF4','LIN54','LEF1','MAFF','FOSL1'),plotAs = "violin",addBoxPlot = TRUE)

for(i in 1:length(plots_violins)) {
    data  <- plots_violins[[i]]$data
    colnames(data) <- c('RNA expression', 'Participants')
    pdf(paste0("/home/wangrong/Plots/",features[i],"_4.5.vlnplot.pdf"),width=5,height=5)
	p <- ggviolin(data,x="Participants",y="RNA expression",fill="Participants",palette = c("#A6CEE3", "#B2DF8A","Bisque","#FB9A99"),add = "boxplot", add.params = list(fill="white"), legend = NULL)+ stat_compare_means(comparisons = my_comparisons)
    p2 <- ggpar(p,legend =  "none")
	plist[[i]] <- p2
	dev.off()
}



geneintegration= getMatrixFromProject(archrproj, useMatrix="GeneExpressionMatrix")
geneintedata=assay(geneintegration)
genedata=rowData(geneintegration)
geneinteinfo=cbind(genedata, geneintedata)
geneinteinfo2=as.data.frame(t(as.data.frame(geneinteinfo[,7:length(colnames(geneinteinfo))])))
colnames(geneinteinfo2)=genedata$name
rownames(geneinteinfo2)=colnames(geneinteinfo)[7:length(colnames(geneinteinfo))]
for(i in rownames(geneinteinfo2)){
  geneinteinfo2[i,"celltype"]=archrproj$celltype[which(archrproj$cellNames==i)]
  geneinteinfo2[i,"Samples"]=archrproj$Sample[which(archrproj$cellNames==i)]
}

bubble_plot_info=data.frame()
for(i in gene_list){
  for(k in 1:length(unique(geneinteinfo2$celltype))){
    a=nrow(bubble_plot_info)
    l=unique(geneinteinfo2$celltype)[k]
    bubble_plot_info[a+1,"gene_name"]=i
    bubble_plot_info[a+1,"celltype"]=l
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'pct_exp']=(length(geneinteinfo2[(geneinteinfo2$",i,">0 & geneinteinfo2$celltype=='",l,"'),'",i,"'])/nrow(geneinteinfo2[geneinteinfo2$celltype=='",l,"',]))*100", sep=""))))
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'avg_exp']=mean(geneinteinfo2[geneinteinfo2$",i,">0 & geneinteinfo2$celltype=='",l,"','",i,"'])", sep=""))))
  }
}

pdf("/home/wangrong/Plots/4.6.Bubble-plot.pdf", width=7, height=7)
ggplot(data = bubble_plot_info, mapping = aes_string(x = 'gene_name', y = 'celltype')) +
  geom_point(mapping = aes_string(size = 'pct_exp', color = "avg_exp")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  labs(
    x = 'gene_name',
    y = 'celltype'
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0))
dev.off()



for(i in 1:length(plots_violins)) {
data  <- plots_violins[[i]]$data
colnames(data) <- c('RNA expression', 'Participants')
pdf(paste0("/home/wangrong/Plots/",features[i],"_4.14.vlnplot.pdf"),width=5,height=5)
p <- ggviolin(data,x="Participants",y="RNA expression",fill="Participants",palette = c("#A6CEE3", "#B2DF8A","Bisque","#FB9A99"),add = "mean", add.params = list(fill="white"), legend = NULL)+ stat_compare_means(comparisons = my_comparisons)
p2 <- ggpar(p,legend =  "none")
print(p2)
dev.off()}

#2023.6.1 
innate_immune_cells <- subset(rename_pbmc,idents=c("Plasmacytoid DC","CD14+ Mono","CD16+ Mono"))
library(readxl,lib.loc = "/home/shisusanna/R/x86_64-pc-linux-gnu-library/4.2")
cytokine <- read_excel("/database/wangrong/Results/0712_ATAC+RNA/48因子与网络节点信息.xlsx")
colnames(cytokine) <- cytokine[1,]
cytokine <- cytokine[-1,]
length(unique(cytokine$`Gene Name`))

inna_cytokine=list()
time <- c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
for (i in 1:length(time)) {
inna_cytokine[[i]]<-FoldChange(innate_immune_cells, group.by = 'Timepoints', ident.1="Day0", ident.2=time[i], features = unique(cytokine$`Gene Name`))
}
all_inna_cytokine<-do.call(cbind,inna_cytokine)
inna_log2FC_cytokine <- all_inna_cytokine[,c(1,4,7,10,13,16,19,22,25,28,31)]
colnames(inna_log2FC_cytokine)<-c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
pheatmap(inna_log2FC_cytokine,border="white", cluster_cols = F, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)

inna_cytokine=list()
time <- c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
for (i in 1:length(time)) {
inna_cytokine[[i]]<-FoldChange(innate_immune_cells, group.by = 'Timepoints', ident.1="Day0", ident.2=time[i], features = c("CXCL8","IL1B","IL1A","IL15","TNF","CCL4","CCL5","CCL3","CSF2RA"))
}
all_inna_cytokine<-do.call(cbind,inna_cytokine)
inna_log2FC_cytokine <- all_inna_cytokine[,c(1,4,7,10,13,16,19,22,25,28,31)]
colnames(inna_log2FC_cytokine)<-c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
pheatmap(inna_log2FC_cytokine,border="white", cluster_cols = F, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)

inna_cytokine=list()
time <- c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
for (i in 1:length(time)) {
inna_cytokine[[i]]<-FoldChange(innate_immune_cells, group.by = 'Timepoints', ident.1="Day0", ident.2=time[i], features = c("CXCL8","IL1B","IL1A","IL15","TNF","CCL4","CCL5","CCL3","CSF2RA"))
}
all_inna_cytokine<-do.call(cbind,inna_cytokine)
inna_log2FC_cytokine <- all_inna_cytokine[,c(1,4,7,10,13,16,19,22,25,28,31)]
colnames(inna_log2FC_cytokine)<-c("Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
pheatmap(inna_log2FC_cytokine,border="white", cluster_cols = F, cluster_rows = T,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 8,fontsize_col = 8,
cellwidth = 10,cellheight = 10)


x = readLines("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -c(1:2))	

timepoints <- list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44")
subset <-list()
day_aggregate<-list()
for (i in 1:length(timepoints)){
print(timepoints[[i]])
subset[[i]]<- subset(rename_pbmc, Timepoints==timepoints[[i]])
subset[[i]] <- NormalizeData(subset[[i]])
subset[[i]] <- FindVariableFeatures(subset[[i]], selection.method = "vst")
subset[[i]] <- ScaleData(subset[[i]], features = rownames(subset[[i]]))
gene.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)

day_s <- AddModuleScore(subset[[i]],assay="RNA",features=gene.list,nbin = 24,min.cells=1,ctrl = 100,name = "INF_score")

day_aggregate[[i]]<-aggregate(day_s$INF_score1, by=list(type=day_s$celltype),mean)
print(paste0(timepoints[i], " is down"))
}

names(day_aggregate) <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44")
df <- as.data.frame(day_aggregate)

df1 <- df[,c(2,4,6,8,10,12,14,16,18,20)]/df[,2]
df1$celltype<-df$Day0.type
df1 <- df1[,-1]
df1<-as.data.frame(t(df1))
colnames(df1) <- df1[10,]
df1 <- df1[-10,]
rownames(df1) <- gsub(".x","",rownames(df1))
df2<-melt(df1,value.name = "celltype")

df3 <- df1
df1 <- as.data.frame(lapply(df1,as.numeric))
rownames(df1) <- rownames(df3)
colnames(df1) <- colnames(df3)

pheatmap(df1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 40,cellheight = 20,legend_labels ="Fold change",annotation_legend=T)	
save(df1,df2,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.3.ifn_score.Rdata")	



write.table(genelist, "/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/6.21.mTOR_signaling.csv",    row.names=FALSE,col.names=TRUE,sep=",") 
timepoints <- list("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44")
subset <-list()
day_aggregate<-list()
for (i in 1:length(timepoints)){
print(timepoints[[i]])
subset[[i]]<- subset(rename_pbmc, Timepoints==timepoints[[i]])
subset[[i]] <- NormalizeData(subset[[i]])
subset[[i]] <- FindVariableFeatures(subset[[i]], selection.method = "vst")
subset[[i]] <- ScaleData(subset[[i]], features = rownames(subset[[i]]))
gene.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)

day_s <- AddModuleScore(subset[[i]],assay="RNA",features=gene.list,nbin = 24,min.cells=1,ctrl = 100,name = "INF_score")

day_aggregate[[i]]<-aggregate(day_s$INF_score1, by=list(type=day_s$Participants),mean)
print(paste0(timepoints[i], " is down"))
}

names(day_aggregate) <- c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185")
day_aggregate<-day_aggregate[-c(11,12)]
df <- as.data.frame(day_aggregate)


df1 <- df[,c(2,4,6,8,10,12,14,16,18,20)]/df[,2]
df1$Participants<-df$Day0.type
df1 <- df1[,-1]
df1<-as.data.frame(t(df1))
colnames(df1) <- df1[10,]
df1 <- df1[-10,]
rownames(df1) <- gsub(".x","",rownames(df1))
df2<-melt(df1,value.name = "Participants")

df3 <- df1
df1 <- as.data.frame(lapply(df1,as.numeric))
rownames(df1) <- rownames(df3)
colnames(df1) <- colnames(df3)
			
pheatmap(df1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 40,cellheight = 20)	

save(df1,df2,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/7.22_m_ifn_score.Rdata")