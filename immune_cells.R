library(CytoTRACE)
library(Seurat)
library(dplyr)
library(harmony)
library(patchwork)
library(ggplot2)
library(cowplot)
load("E:/Bowen-files/out/v5_addT90_new/scRNA_harmony_0.2_named.Rdata")
Idents(sce)<- "celltype"
sce <-subset(x = sce, idents = "Immune cells")
setwd("E:/Bowen-files/out/v5_addT90_new/Immune cell/26_3_11")
table(sce@meta.data$celltype)
Epi_sce =CreateSeuratObject(counts = GetAssayData(sce, assay="RNA",layer='counts'),  # 使用提取的细胞构建新的Seurat对象
                            meta.data =sce@meta.data) # 保留meta.data所有信息
#标注化、归一化、高变基因、pca
Epi_sce <- NormalizeData(Epi_sce)%>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
Epi_sce <-RunHarmony(Epi_sce, group.by.var = "orig.ident")#降维聚类
Epi_sce
p <-ElbowPlot(Epi_sce, ndims=50, reduction="pca")
pdf("ElbowPlot.pdf",width= 7,height = 7)
p
dev.off()
Epi_sce <-FindNeighbors(Epi_sce, reduction = "pca", dims =1:20)#reduction="ha
Epi_sce =FindClusters(Epi_sce,resolution = 0.05)
table(Epi_sce@meta.data$seurat_clusters) 
Epi_sce <-RunUMAP(Epi_sce, reduction = "harmony", dims =1:15)##reduction="harmony
#EPI_umap图_未注释
Epi_sce$group
# 将Epi_sce$group转换为因子，并指定水平顺序
Epi_sce$group<- factor(Epi_sce$group, levels = c("E90", "E110","P7")) # 验证水平顺序
levels(Epi_sce$group)

plot1 =DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE)
plot2 = DimPlot(Epi_sce, reduction = "umap", group.by='orig.ident',raster=FALSE)
plot3 = DimPlot(Epi_sce, reduction = "umap",split.by = "group",label = T,raster=FALSE)
plot4 = DimPlot(Epi_sce, reduction = "umap",group.by = "group",shuffle = T,raster=FALSE)
pdf("Germ_cells_3group.pdf",width = 21,height = 7)
plot3
dev.off()

pdf("Immune_cell_group.pdf",width = 7,height = 7)
plot4
dev.off()

pdf("Immune_cell_umap.pdf",width = 7,height = 7)
plot1
dev.off()
save.image("Immune_cell.Rdata")
