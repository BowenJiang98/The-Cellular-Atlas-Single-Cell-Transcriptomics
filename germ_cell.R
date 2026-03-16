library(CytoTRACE)
library(Seurat)
library(dplyr)
library(harmony)
library(patchwork)
library(ggplot2)
library(cowplot)
load("E:/Bowen-files/out/v5_addT90_new/scRNA_harmony_0.2_named.Rdata")
Idents(sce)<- "celltype"
sce <-subset(x = sce, idents = "Germ cells")
setwd("E:/Bowen-files/out/v5_addT90_new/Germ_cell/26_3_11")
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

pdf("Germ_cells_group.pdf",width = 7,height = 7)
plot4
dev.off()

pdf("Germ_cells_umap.pdf",width = 7,height = 7)
plot1
dev.off()
save.image("Germ_cell.Rdata")
saveRDS(Epi_sce,file="Tcell.RDS")
rm(list=ls())
library(monocle)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
Tcell <- readRDS("Tcell.RDS")
DefaultAssay(Tcell) <- "RNA"
Tcell <- JoinLayers(Tcell)
# 从亚群seurat对象中提取拟时分析所需的数据
exp <- Tcell[["RNA"]]$counts#提取表达量，counts就是一个稀疏矩阵
##提取feature data基因信息
fdata <- data.frame(
  gene_short_name = row.names(Tcell), 
  row.names = row.names(Tcell)
)


pdata <- Tcell@meta.data 

## 1.2 CDS对象的创建
# 将基因特征文件和细胞表型文件重新写一个对象
fd <- new("AnnotatedDataFrame", data = fdata)
pd <- new("AnnotatedDataFrame", data = pdata)

# 构建monocle分析的CDS对象
CDS <- newCellDataSet(cellData = exp, phenoData = pd, featureData = fd)

#计算size factors 和 dispersions（离差），用于后期分析；
#结果：在phenoData表格添加1列Size_Factor；
CDS <- estimateSizeFactors(CDS)
memory.limit(100000)
CDS <- estimateDispersions(CDS)
gc()

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS, file = "CDS_germ.Rdata")
rm(list = ls())
gc()
load("CDS_germ.Rdata")
dim(CDS)

CDS <- detectGenes(CDS, min_expr = 0.1)
#过滤低表达的基因，以降低噪音和减少差异分析计算量，fdata中细胞数小于20去除，大于保留
expressed_genes <- fData(CDS) %>%
  subset(num_cells_expressed >= 50) %>%
  rownames()
length(expressed_genes)

#CDS <- clusterCells(CDS, verbose = F)

#不同细胞类型的差异分析,只提供p和Qvalue，p小于0.05保留
clustering_DEG_genes <- differentialGeneTest(
  CDS[expressed_genes, ], 
  fullModelFormulaStr = '~seurat_clusters',
  cores=10
)
gc()


## 3. 时间轨迹及其差异分析
## 3.1 构建细胞轨迹
# 第一步：选择用于构建细胞轨迹的基因集；
# 选择Top800差异基因作为排序基因；
ordering_genes <- rownames(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]


#将差异基因储存到CDS对象中；
CDS <- setOrderingFilter(CDS, ordering_genes = ordering_genes)
plot_ordering_genes(CDS)

#第二步: 数据降维
#降维函数与上文的t-SNE一致，但降维算法这里用的是“DDRTree”，
CDS <- reduceDimension(CDS, method = 'DDRTree')

#第三步: 构建细胞分化轨迹
#按照分化轨迹排序细胞；
CDS <- orderCells(CDS)

#绘制细胞分化轨迹：
#按“Pseudotime”分组；
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")
ggsave(Pseudotime_clusters,filename = "germ_ProPseudotime.pdf", width = 7, height = 7,units = "in",dpi = 300)
#按“State”分组；
State_clusters <-plot_cell_trajectory(CDS, color_by = "State")
ggsave(State_clusters,filename = "germ_State.pdf", width = 7, height = 7,units = "in",dpi = 300)

#按seurat分群结果分组
clusters <-plot_cell_trajectory(CDS, color_by = "seurat_clusters")
ggsave(clusters,filename = "germ_clusters.pdf", width = 7, height = 7,units = "in",dpi = 300)
#按细胞类型分组
#plot_cell_trajectory(CDS, color_by = "celltype")
#按照ori_ident分类
clusters <-plot_cell_trajectory(CDS, color_by = "group")
ggsave(clusters,filename = "germ_orig.ident.pdf", width = 7, height = 7,units = "in",dpi = 300)

save(CDS, file = "CDS_germ_pseudotime.Rda")

sig_gene_names2 <- c('LGALS1','ART3','YBX1','HSF2BP','EIF4A2','PABPC1','SKP1')

for (gene in sig_gene_names2) {
  
  # 计算表达量
  pData(CDS)[[gene]] <- log2(exprs(CDS)[gene, ] + 1)
  
  # 画图
  p <- plot_cell_trajectory(
    CDS,
    color_by = gene,
    size = 1,
    show_backbone = TRUE
  ) +
    theme(panel.border = element_rect(fill = NA,
                                      color = "black",
                                      size = 1,
                                      linetype = "solid")) +
    scale_color_gradient2(
      low = 'grey',
      mid = 'orange',
      high = 'blue'
    )
  
  # 保存
  ggsave(
    filename = paste0("pseudotime_", gene, ".png"),
    plot = p,
    width = 7,
    height = 7,
    units = "in",
    dpi = 300
  )
}
