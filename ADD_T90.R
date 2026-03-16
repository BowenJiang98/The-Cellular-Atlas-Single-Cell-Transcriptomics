#加载包
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90")

E90.1  <- Read10X(data.dir = "E:/Bowen-files/out/T90_2/filtered_feature_bc_matrix")
E90.1=CreateSeuratObject(counts = E90.1,project = "E90.1", min.cells = 3, min.features = 200)
E90.2  <- Read10X(data.dir = "E:/Bowen-files/out/T90_3/filtered_feature_bc_matrix")
E90.2=CreateSeuratObject(counts = E90.2,project = "E90.2", min.cells = 3, min.features = 200)
E90.3  <- Read10X(data.dir = "E:/Bowen-files/out/T90_4/filtered_feature_bc_matrix")
E90.3=CreateSeuratObject(counts = E90.3,project = "E90.3", min.cells = 3, min.features = 200)
#E90.4  <- Read10X(data.dir = "E:/Bowen-files/out/T90_1/filtered_feature_bc_matrix")
#E90.4=CreateSeuratObject(counts = E90.4,project = "E90.4", min.cells = 3, min.features = 200)

E110.1 <- Read10X(data.dir = "E:/Bowen-files/out/P1101/filtered_feature_bc_matrix")
E110.1=CreateSeuratObject(counts = E110.1,project = "E110.1", min.cells = 3, min.features = 200)
E110.2 <- Read10X(data.dir = "E:/Bowen-files/out/P1102/filtered_feature_bc_matrix")
E110.2=CreateSeuratObject(counts = E110.2,project = "E110.2", min.cells = 3, min.features = 200)
E110.3 <- Read10X(data.dir = "E:/Bowen-files/out/220925A_T110_T110")
E110.3=CreateSeuratObject(counts = E110.3,project = "E110.3", min.cells = 3, min.features = 200)

P7.1 <- Read10X(data.dir = "E:/Bowen-files/out/P71/filtered_feature_bc_matrix")
P7.1=CreateSeuratObject(counts = P7.1,project = "P7.1", min.cells = 3, min.features = 200)
P7.2 <- Read10X(data.dir = "E:/Bowen-files/out/P72/filtered_feature_bc_matrix")
P7.2=CreateSeuratObject(counts = P7.2,project = "P7.2", min.cells = 3, min.features = 200)

#合并多个样本策略
#seurat_object <- merge(x,y=c(a,b,c), add.cell.ids = c("x","a","b","c"), project = "four_merged")
seurat_object <- merge(E90.1,y=c(E90.2,E90.3,E110.1,E110.2,E110.3,P7.1,P7.2), add.cell.ids = c("E90.1","E90.2","E90.3","E110.1","E110.2","E110.3","P7.1","P7.2"), project = "10_merged")
#计算每个细胞的线粒体基因转录本数的百分比（%），使用[[ ]] 操作符存放到 metadata 中； 
seurat_object [["percent.mt"]] <- PercentageFeatureSet(seurat_object, features = c("ND1","COX1","ND2","COX2","ATP8","ATP6","ND3","COX3","ND5","ND4","ND4L","ND6","CYTB"))
#seurat_object <- subset(seurat_object, subset = classifications=="Singlet")
#nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例。
pdf("nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例.pdf",height = 7,width = 21)
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#根据4分位数的1.5倍确定nCount_RNA和nFeature_RNA
res_nCount_RNA <- seurat_object@meta.data[["nCount_RNA"]]
res_nCount_RNA<-quantile(res_nCount_RNA, probs = c(0,0.25,0.5,0.75,1))
res_nCount_RNA

res_nFeature_RNA <- seurat_object@meta.data[["nFeature_RNA"]]
res_nFeature_RNA<-quantile(res_nFeature_RNA, probs = c(0,0.25,0.5,0.75,1))
res_nFeature_RNA
#过滤细胞：保留 gene 数大于200 小于 2500 的细胞；目的是去掉空 GEMs 和 1 个 GEMs 包 含 2 个以上细胞的数据；
#而保留线粒体基因的转录本数低于20%的细胞，为了过滤掉死细胞 等低质量的细胞数据。 
seurat_object <- subset(seurat_object, subset = 5013 > nFeature_RNA & nFeature_RNA > 200 & nCount_RNA < 14321 & percent.mt < 20)

## 对过滤后的 QC metrics 进行可视化（绘制散点图）；
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("plot2.pdf", width = 28, height = 25, units = "cm")

#表达量数据标准化：LogNormalize 的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ) 
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因（feature selection），用于下游分析，PCA
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# 提取表达量变变化最高的 10 个基因； 
top10 <- head(VariableFeatures(seurat_object), 10)
top10
# 绘制带有和不带有标签的变量特征的散点图
plot3 <- VariableFeaturePlot(seurat_object)+NoLegend()
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge=0, ynudge=0)
plot3+plot4
ggsave("plot3.pdf", width = 28, height = 25, units = "cm")
#排除细胞周期基因对分群和降维的影响
#细胞周期marker基因加载
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
##You can try to merge the layer together with data.filt <- JoinLayers(data.filt)
seurat_object <- JoinLayers(seurat_object)
#计算细胞周期分数
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(seurat_object[[]])
#去除细胞周期对分群和降维的影响
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(seurat_object))
#而对所有基因进行标准化的方法如下： 
#seurat_object <- ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = row.names(seurat_object)) ##耗时2min
#线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定； 
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) 
# 检查 PCA 分群结果， 这里只展示前 12 个 PC,每个 PC 只显示 3 个基因； 
print(seurat_object[["pca"]], dims = 1:20, nfeatures = 3) 
#绘制 pca 散点图； 去除图例
DimPlot(seurat_object, reduction = "pca")
#画前 2 个主成分的热图； 
DimHeatmap(seurat_object, dims = 1:2, cells = 500, balanced = TRUE) 
#ggsave("plot2.png", width = 15, height = 13, units = "cm")


#确定数据集的分群个数 
##方法 1：Jackstraw 置换检验算法；重复取样（原数据的 1%），重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因 scores; 
seurat_object <- JackStraw(seurat_object, num.replicate = 100,dims = 30)  ##耗时3min
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:30) #dims必须小于等于上一步
JackStrawPlot(seurat_object, dims = 1:30) #dims必须在上一步的范围内
#方法 2：肘部图（碎石图），基于每个主成分对方差解释率的排名； 
ElbowPlot(seurat_object,ndims = 30)

###方法二：Harmony矫正
#harmony矫正
seurat_object = RunHarmony(seurat_object,"orig.ident", plot_convergence = TRUE)#耗时1min
#细胞聚类
seurat_object <- seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.4)
#降维可视化
#seurat_object = RunTSNE(seurat_object,dims = 1:20)
p1 <- DimPlot(seurat_object, reduction = "umap", shuffle = T,group.by = "orig.ident", pt.size = 0.5)
p2 <- DimPlot(seurat_object, reduction = "umap",label=TRUE,pt.size = 0.4)
p1+p2
pdf("umap_0.4.pdf") 
DimPlot(seurat_object, reduction = "umap",pt.size = 0.5,label = T)
DimPlot(seurat_object, reduction = "umap",label=TRUE,pt.size = 0.5,group.by = "orig.ident")
dev.off()

#p3 <- DimPlot(seurat_object, reduction = "tsne", shuffle = T,group.by = "orig.ident", pt.size = 0.5)
#p4 <- DimPlot(seurat_object, reduction = "tsne",label=TRUE,pt.size = 0.5)
#p3+p4

table(seurat_object$orig.ident)
prop.table(table(Idents(seurat_object)))
table(Idents(seurat_object), seurat_object$orig.ident)

saveRDS(seurat_object,file="0.4_unnamed_6sample.rds")
#生殖细胞
UCHL1 <- FeaturePlot(seurat_object, features = c('UCHL1'),label=TRUE,order = TRUE)
ggsave(UCHL1,filename = "UCHL1.png",width = 7, height = 7,units = "in",dpi = 300) 

#膜细胞和基质细胞
ICFBP5 <- FeaturePlot(seurat_object, features = c('HSD3B'),label=TRUE,order = TRUE)
ggsave(ICFBP5,filename = "FATE1.png",width = 7, height = 7,units = "in",dpi = 300)


DCN <- FeaturePlot(seurat_object, features = c('DCN'),label=TRUE,order = TRUE)
ggsave(DCN,filename = "DCN.png",width = 7, height = 7,units = "in",dpi = 300) 
COL1A1 <- FeaturePlot(seurat_object, features = c('COL1A1'),label=TRUE,order = TRUE)
ggsave(COL1A1,filename = "COL1A1.png",width = 7, height = 7,units = "in",dpi = 300) 
COL3A1 <- FeaturePlot(seurat_object, features = c('COL3A1'),label=TRUE,order = TRUE)
ggsave(COL3A1,filename = "COL3A1.png",width = 7, height = 7,units = "in",dpi = 300) 

#平滑肌细胞：
MUSTN1 <- FeaturePlot(seurat_object, features = c('MUSTN1'),label=TRUE,order = TRUE)
ggsave(MUSTN1,filename = "MUSTN1.png",width = 7, height = 7,units = "in",dpi = 300) 

ACTA2 <- FeaturePlot(seurat_object, features = c('ACTA2'),label=TRUE,order = TRUE)
ggsave(ACTA2,filename = "ACTA2.png",width = 7, height = 7,units = "in",dpi = 300) 

#内皮细胞:
TM4SF1 <- FeaturePlot(seurat_object, features = c('TM4SF1'),label=TRUE,order = TRUE)
ggsave(TM4SF1,filename = "TM4SF1.png",width = 7, height = 7,units = "in",dpi = 300) 
VWF <- FeaturePlot(seurat_object, features = c('VWF'),label=TRUE,order = TRUE)
ggsave(VWF,filename = "VWF.png",width = 7, height = 7,units = "in",dpi = 300) 
#单核细胞：
TYROBP <- FeaturePlot(seurat_object, features = c('TYROBP'),label=TRUE,order = TRUE)
ggsave(TYROBP,filename = "TYROBP'.png",width = 7, height = 7,units = "in",dpi = 300)
IFI30 <- FeaturePlot(seurat_object, features = c('IFI30'),label=TRUE,order = TRUE)
ggsave(IFI30,filename = "IFI30.png",width = 7, height = 7,units = "in",dpi = 300)
CD163 <- FeaturePlot(seurat_object, features = c('CD163'),label=TRUE,order = TRUE)
ggsave(CD163,filename = "CD163.png",width = 7, height = 7,units = "in",dpi = 300)
C1QB <- FeaturePlot(seurat_object, features = c('C1QB'),label=TRUE,order = TRUE)
ggsave(C1QB,filename = "C1QB.png",width = 7, height = 7,units = "in",dpi = 300)


#自然杀伤细胞：
CCL5 <- FeaturePlot(seurat_object, features = c('CCL5'),label=TRUE,order = TRUE)
ggsave(CCL5,filename = "CCL5.png",width = 7, height = 7,units = "in",dpi = 300)
NKG7 <- FeaturePlot(seurat_object, features = c('NKG7'),label=TRUE,order = TRUE)
ggsave(NKG7,filename = "NKG7.png",width = 7, height = 7,units = "in",dpi = 300)
#T淋巴细胞：
KLRB1 <- FeaturePlot(seurat_object, features = c('KLRB1'),label=TRUE,order = TRUE)
ggsave(KLRB1,filename = "KLRB1.png",width = 7, height = 7,units = "in",dpi = 300)

#支持细胞
INHA<- FeaturePlot(seurat_object, features = c('INHA'),label=TRUE,order = TRUE)
ggsave(INHA,filename = "INHA.png",width = 7, height = 7,units = "in",dpi = 300)

PRND <- FeaturePlot(seurat_object, features = c('PRND'),label=TRUE,order = TRUE)
ggsave(PRND,filename = "PRND.png",width = 7, height = 7,units = "in",dpi = 300)
#
WNT5A <- FeaturePlot(seurat_object, features = c('WNT5A'),label=TRUE,order = TRUE)
ggsave(WNT5A,filename = "WNT5A.png",width = 7, height = 7,units = "in",dpi = 300)
#间质细胞干细胞 ITGAV
CD51 <- FeaturePlot(seurat_object, features = c('CD51'),label=TRUE,order = TRUE)
ggsave(CD51,filename = "CD51.png",width = 7, height = 7,units = "in",dpi = 300)

ITGAV <- FeaturePlot(seurat_object, features = c('ITGAV'),label=TRUE,order = TRUE)
ggsave(ITGAV,filename = "ITGAV.png",width = 7, height = 7,units = "in",dpi = 300)

HSD3B <- FeaturePlot(seurat_object, features = c('HSD3B7'),label=TRUE,order = TRUE)
ggsave(HSD3B,filename = "HSD3B.png",width = 7, height = 7,units = "in",dpi = 300)

PDGFRA <- FeaturePlot(seurat_object, features = c('PDGFRA'),label=TRUE,order = TRUE)
ggsave(PDGFRA,filename = "PDGFRA.png",width = 7, height = 7,units = "in",dpi = 300)

NES <-FeaturePlot(seurat_object, features = c('NES'),label=TRUE,order = TRUE)
ggsave(NES,filename = "NES.png",width = 7, height = 7,units = "in",dpi = 300)

NR2F2 <-FeaturePlot(seurat_object, features = c('NR2F2'),label=TRUE,order = TRUE)
ggsave(NR2F2,filename = "NR2F2.png",width = 7, height = 7,units = "in",dpi = 300)

P75NTR <-FeaturePlot(seurat_object, features = c('P75NTR'),label=TRUE,order = TRUE)
ggsave(P75NTR,filename = "P75NTR.png",width = 7, height = 7,units = "in",dpi = 300)

THY1 <-FeaturePlot(seurat_object, features = c('THY1'),label=TRUE,order = TRUE)
ggsave(THY1,filename = "THY1.png",width = 7, height = 7,units = "in",dpi = 300)

new.cluster.ids <- c("Stroma_cell",#0
                     "Sertoli_cell",  #1
                     "Leydig_cell", #2
                     "Stroma_cell", #3
                     "Sertoli_cell", #4
                     "Sertoli_cell", #5
                     "Stroma_cell",#6
                     "Macrophage", #7
                     "Vascular_cell", #8
                     "Lymphocyte",#9
                     "Sertoli_cell", #10
                     "Myoid_cell", #11
                     "Macrophage",#12
                     "Stroma_cell",#13
                     "S'gonia", #14
                     "Vascular_cell",#15
                     "Sertoli_cell", #16
                     "Macrophage", #17
                     "Leydig_cell",#18
                     "Macrophage",#19
                     "Leydig_cell",#20
                     "Macrophage",#21
                     "Macrophage"#22
)#25
names(new.cluster.ids) 
levels(seurat_object)
#将seurat_object的水平属性赋值给new.cluster.ids的names属性； 
names(new.cluster.ids) <- levels(seurat_object)
names(new.cluster.ids) 
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 

#重画umap
namedmap <- DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()
#自定义颜色
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
p3<-namedmap+scale_color_manual(values = color)
p3
#自定义主题
mytheme<-theme_classic()+
  theme(panel.background = element_rect(fill = 'white',color = 'white'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(color = 'black',size = 10),
        axis.title = element_text(color = 'black',size = 12),
        plot.title = element_text(size = 18,hjust = 0.5)
  )
p4<-p3+ggtitle('Cell Type')+mytheme
pdf("umap_0.4_named.pdf") 
p4
dev.off()
ggsave(filename = "namedumap.png",plot = p4,width = 9, height = 7, units = "in",dpi = 300)


gene_to_check3=c('COL1A1','COL3A1','DCN','LUM','COL3A1','PDGFRA','CLU','AMH','GATA4','SOX9','CYP17A1','CYP11A1','INSL3','STAR','PTPRC','CST3','C1QA','TYROBP','PECAM1','VWF','RAMP2','CDH5','IL7R','CD52','MYH11','ACTA2','NOTCH3','RGS5','DDX4','UCHL1','DMRT1','SYCP3')

gene_to_check3 <- as.data.frame(gene_to_check3)
markerdata<-ScaleData(seurat_object,features = as.character(unique(gene_to_check3$gene_to_check3)))
pdf("气泡图_named.pdf") 
DotPlot(markerdata, features = as.character(unique(gene_to_check3$gene_to_check3)))+coord_flip()+ theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.x = element_text(angle = 45, vjust = 1))
dev.off()

saveRDS(seurat_object,file='seurat_object_named_0.4_harmony.Rds')
#差异基因
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(markers,"markers_allclusters.csv",row.names=FALSE,col.names=TRUE,sep=",")
heatmap10 <- DoHeatmap(seurat_object, features = top10$gene,size = 0.5,slot = "data")+NoLegend()
ggsave(heatmap10,filename = "heatmap.png",width = 12, height = 20,units = "in",dpi = 300)
#画出来是空白的；很奇怪
#seurat_object[["RNA4"]] <- as(object = seurat_object[["RNA"]], Class = "Assay4")

#top20heat <- DoHeatmap(seurat_object,features = top20$gene,group.by = "seurat_clusters",assay = 'RNA',label=F,slot = "data",
#                       group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
#  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

#ggsave(top20heat,filename = "heatmaptop20.png",width = 12, height = 20,units = "in",dpi = 300)
seurat_object <- seurat_object_named_0.4_harmony
seurat_object[["celltype"]]<-Idents(seurat_object)
#seurat_object <- seurat_object_named_0.4_harmony
rm(seurat_object_named_0.4_harmony)
DefaultAssay(seurat_object)='RNA'
seurat_object<- NormalizeData(seurat_object)
seurat_object<- FindVariableFeatures(seurat_object)
seurat_object<-ScaleData(seurat_object,verbose = FALSE,features = rownames(seurat_object),vars.to.regress = c("S.Score", "G2M.Score"))
seurat_object@assays[["RNA"]]@layers[["scale.data"]] <- scale(seurat_object@assays[["RNA"]]@layers[["data"]], scale = TRUE)

markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

heatmap20 <- DoHeatmap(seurat_object,features = c(top20$gene),group.by = "celltype",assay = 'RNA',label=F, group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
#尝试使用pdf保存
pdf("heatmap20.pdf")
heatmap20
dev.off()
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap10 <- DoHeatmap(seurat_object,features = c(top10$gene),group.by = "celltype",assay = 'RNA',label=F,group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
#这个可以画出来，但是会omit一些基因
DoHeatmap(seurat_object, label = F , # 不加label
          features = as.character(unique(top10$gene)),   
          group.by = "ident",  
          assay = "RNA",  
          group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B","#FFA500"))+ #设置组别颜色
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))#设置热图颜色

#画出来是空白
ggsave(heatmap10,filename = "heatmaptop10.png",width = 12, height = 20,units = "in",dpi = 300)
#尝试使用pdf保存
pdf("heatmap10.pdf")
heatmap10
dev.off()

#接着把间质细胞提取出来做一个拟时序分析
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/monocle_leydigcell")
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/seurat_object_named_0.4_harmony.Rds")
Tcell <- subset(x = seurat_object, idents = "Leydig_cell")

TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG

rm(seurat_object)
testA <- subset(Tcell,orig.ident=="E90")
testA.seu <- testA[["RNA"]]$counts
testA.seu=CreateSeuratObject(counts = testA.seu,project = "E90", min.cells = 3, min.features = 200)
testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)

testB <- subset(Tcell,orig.ident=="E110")
testB.seu <- testB[["RNA"]]$counts
testB.seu=CreateSeuratObject(counts = testB.seu,project = "E110", min.cells = 3, min.features = 200)
testB.seu <- NormalizeData(testB.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testB.seu <- FindVariableFeatures(testB.seu, selection.method = "vst", nfeatures = 2000)


testC <- subset(Tcell,orig.ident=="P7")
testC.seu <- testC[["RNA"]]$counts
testC.seu=CreateSeuratObject(counts = testC.seu,project = "P7", min.cells = 3, min.features = 200)
testC.seu <- NormalizeData(testC.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testC.seu <- FindVariableFeatures(testC.seu, selection.method = "vst", nfeatures = 2000)


testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu,testC.seu), dims = 1:20)#可设置不同的k.anchor = 5,k.filter = 30参数，这里使用默认的k.filter =200
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:20)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.1)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:20)
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_Leydigs_1.png",width = 14, height = 7,units = "in",dpi = 300)

Tcell <- testAB.integrated
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
testAB.integrated <- Tcell  
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_Leydigs.png",width = 14, height = 7,units = "in",dpi = 300)

CYP17A1 <- FeaturePlot(Tcell, features = c('CYP17A1'),label=TRUE,order = TRUE)
ggsave(CYP17A1,filename = "CYP17A1.png",width = 7, height = 7,units = "in",dpi = 300)

save(Tcell, file="unnamed_leydig.Rda")
rm(list=ls())
load(file="unnamed_leydig.Rda")
library(monocle)
library(Seurat)
library(AUCell)
library(patchwork)
library(ggplot2)
library(DOSE)#人疾病
library(clusterProfiler)#富集分析
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Rn.eg.db)
#library(org.At.tair.db)
library(dplyr)
library(GO.db)
library(rjson)
library(stringr)
library("reshape2")
#dim(Tcell)
#Tcell <- subset(Tcell, orig.ident == "E110.1")
#Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(0,2,20,5,1,14,9,10,3,4)]
#sertoli cell

# 为了减少内存占用，可以保存提取数据的seurat对象，并且删除原来的seurat对象
#save(Tcell, file = "sertolinamed.Rda")
#rm(Tcell)
gc()
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

CDS <- estimateDispersions(CDS)
gc()

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS, file = "CDS_Leydig.Rdata")
rm(list = ls())
gc()
load("CDS_Leydig.Rdata")
dim(CDS)

################################################################################
## 2. 差异分析寻找高变基因
# detectGenes()函数：同时统计表达当前基因的细胞数量和细胞表达的基因数；
# min_expr参数用于设置检测阈值，比如min_expr = 0.1表示当前基因的表达量超过0.1才会
# 纳入统计；
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

################################################################################
## 3. 时间轨迹及其差异分析
## 3.1 构建细胞轨迹
# 第一步：选择用于构建细胞轨迹的基因集；
# 选择Top200差异基因作为排序基因；
ordering_genes <- rownames(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]
flc_candidates <- c("PDGFRA","CYP11A1", "STAR","CYP17A1","WNT5A","NES","NR5A1", "DLK1")
sig_genes <- union(flc_candidates,ordering_genes)
#fData(CDS)$use_for_ordering <-fData(CDS)$num_cells_expressed > 0.1 * ncol(CDS)

#将差异基因储存到CDS对象中；
CDS <- setOrderingFilter(CDS, ordering_genes = sig_genes)
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
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.png", width = 7, height = 7,units = "in",dpi = 300)
#按“State”分组；
State_clusters <-plot_cell_trajectory(CDS, color_by = "State")
ggsave(State_clusters,filename = "leydig_State.png", width = 7, height = 7,units = "in",dpi = 300)

#按seurat分群结果分组
clusters <-plot_cell_trajectory(CDS, color_by = "seurat_clusters")
ggsave(clusters,filename = "leydig_clusters.png", width = 7, height = 7,units = "in",dpi = 300)
#按细胞类型分组
plot_cell_trajectory(CDS, color_by = "celltype")
#按照ori_ident分类
clusters <-plot_cell_trajectory(CDS, color_by = "orig.ident")
ggsave(clusters,filename = "leydig_orig.ident.png", width = 7, height = 7,units = "in",dpi = 300)

save(CDS, file = "CDS_leydig_pseudotime.Rda")

##3.2 比较细胞分化轨迹进程中功能基因的表达差异
#主要用到sm.ns()函数根据表达量拟合曲线；用拟时序做基因差异表达分析
diff_test_res <- differentialGeneTest(
  CDS[expressed_genes, ], 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
head(diff_test_res[, c("gene_short_name", "pval", "qval")])
save.image(file = "monocle_leydig.RData")
# 按q值从小到大排序后，查看最显著的前6个基因的拟时间趋势
sig_gene_names1 <- rownames(diff_test_res[order(diff_test_res$qval)[1:4], ])
seurat_cluster <- plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'Pseudotime')
ggsave(seurat_cluster,filename = "Leydig_seurat_cluster.png", width = 14, height = 14,units = "in",dpi = 300)
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

# 也可以查看比较关注基因的拟时间趋势
sig_gene_names1 <- c("DDX4", "UCHL1", "DNMT3L","PCNA","PIWIL4","ETV4","ETV5","KDM1B","GFRA1")
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

sig_gene_names2 <- c('CYP17A1','CYP11A1','INSL3','STAR')
Leydigmarker <- plot_genes_in_pseudotime(CDS[sig_gene_names2, ], color_by = 'State')
ggsave(Leydigmarker,filename = "pseudotime_Leydigmarker.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$PDGFRA=log2(exprs(CDS)['PDGFRA',]+1)
pseudotime_PDGFRA <- plot_cell_trajectory(CDS,color_by="PDGFRA", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_PDGFRA,filename = "pseudotime_PDGFRA.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$CYP11A1=log2(exprs(CDS)['CYP11A1',]+1)
pseudotime_CYP11A1 <- plot_cell_trajectory(CDS,color_by="CYP11A1", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_CYP11A1,filename = "pseudotime_CYP11A1.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$STAR=log2(exprs(CDS)['STAR',]+1)
pseudotime_STAR <- plot_cell_trajectory(CDS,color_by="STAR", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_STAR,filename = "pseudotime_STAR.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$INSL3=log2(exprs(CDS)['INSL3',]+1)
pseudotime_INSL3 <- plot_cell_trajectory(CDS,color_by="INSL3", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_INSL3,filename = "pseudotime_INSL3.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$KIT=log2(exprs(CDS)['KIT',]+1)
pseudotime_KIT <- plot_cell_trajectory(CDS,color_by="KIT", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_KIT,filename = "pseudotime_KIT.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$DNMT3L=log2(exprs(CDS)['DNMT3L',]+1)
pseudotime_DNMT3L <- plot_cell_trajectory(CDS,color_by="DNMT3L", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_DNMT3L,filename = "pseudotime_DNMT3L.png",width = 7, height = 7,units = "in",dpi = 300)


#更换根节点：#按照分化轨迹排序细胞；
CDS <- orderCells(CDS,root_state = 2)
#按“Pseudotime”分组；
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")+scale_color_gsea()+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+scale_color_gradient2(low="grey",mid="orange",high="blue")
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.png", width = 7, height = 7,units = "in",dpi = 300)


##3.4 差异基因的拟时表达模式聚类分析
#提取差异基因；
sig_gene_names2 <- row.names(subset(diff_test_res, qval < 0.001))
#绘制拟时间差异基因表达谱热图；#重新聚类
plot_pseudotime_heatmap(
  CDS[sig_gene_names2, ], 
  num_clusters = 3, #cluster数量
  show_rownames = F,
  cores = 8
)

Time_diff <- differentialGeneTest(CDS[sig_gene_names2,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
write.csv(Time_diff,"Tiff_diff_all.csv",row.names = F)

Time_genes <- top_n(Time_diff, n = 2000, desc(qval)) %>% pull(gene_short_name) %>% as.character()  

#基因太多了画不出来就选择前2000个基因，所有基因Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p<-plot_pseudotime_heatmap(CDS[sig_gene_names2,], num_clusters=3, show_rownames=F, return_heatmap=T)  #字体太大，看不清，可以导出
ggsave("Time_heatmap.png", p, width = 5, height = 8)
ggsave("Time_heatmap.pdf", p, width = 5, height = 10)

clusters<-cutree(p$tree_row,k=3)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
write.csv(clustering,"间质细胞2000个拟时序基因热图分3个cluster.csv")
save.image(file = "间质细胞-monocle2.Rdata") 
################################################################################
## 4.单细胞轨迹分支分析 
# 当细胞分化轨迹出现分支的时候，意味着细胞将面临不同的分化命运“fate”，接下来主要
# 分析分支事件，比如沿着分化轨迹，基因的表达量如何变化？
# 不同分支之间的差异基因有哪些？

# Monocle 提供一个特定的统计检验方法: branched expression analysis modeling（BEAM）.
## 4.1 BEAM检验
# 使用BEAM()函数对基因进行注释；
# BEAM函数的输入对象： 完成拟时间分析的CellDataSet且轨迹中包含1个分支点；
# 返回一个包含每个基因significance scores 的表格,若得分是显著的则表明该基因的表
# 达是与分支相关的（branch-dependent）。

BEAM_res <- BEAM(CDS[expressed_genes, ], branch_point = 1,cores=10,progenitor_method = "duplicate")
saveRDS(BEAM_res, file = "BEAM_res.RDS")

#按照qval升序排列；
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
head(BEAM_res)

##4.2轨迹分支表达分析
#使用pheatmap包绘制分支表达量热图；
sig_gene_names3 <- row.names(subset(BEAM_res, qval < 1e-5))
plot_genes_branched_heatmap(
  CDS[sig_gene_names3, ],
  branch_point = 1,#设置分支节点
  num_clusters = 3,
  use_gene_short_name = FALSE,
  show_rownames = FALSE)

#使用 plot_genes_branched_pseudotime() 函数绘制拟合曲线，不同线表示不同分支
plot_genes_branched_pseudotime(
  CDS[rownames(BEAM_res)[1:4], ], 
  branch_point = 1, 
  color_by = "State",
  ncol = 1
)

##细胞通讯分析画图
setwd("E:/Bowen-files/out/v5_addT90_new/cellchat")
#setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/cellchat")
rm(list=ls())
library(dplyr)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(monocle)
library(tidyverse)
library(Rcpp)
library(ggplot2)
library(clustree)
library(cowplot)
library(ggsci)
library(CellChat)
library(ggalluvial)
library(svglite)
library(NMF)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
#pbmc3k.final <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/seurat_object_named_0.4_harmony.Rds")
load("E:/Bowen-files/out/v5_addT90_new/scRNA_harmony_0.2_named.Rdata")
Idents(sce)
pbmc3k.final <- sce
pbmc3k.final@commands$FindClusters 
pbmc3k.final$cell_annotations <- Idents(pbmc3k.final) 
pbmc3k.final
pbmc3k.final [["RNA3"]] <- as(object = pbmc3k.final[["RNA"]], Class = "Assay")
DefaultAssay(pbmc3k.final) <- "RNA3"
pbmc3k.final[["RNA"]] <- NULL
pbmc3k.final <- RenameAssays(object = pbmc3k.final, RNA3 = 'RNA')


#cellchat <- createCellChat(object=pbmc3k.final,group.by = "seurat_clusters")
cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "celltype")
summary(cellchat)
str(cellchat)
levels(cellchat@idents)                   
groupSize <- as.numeric(table(cellchat@idents))#
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB) 
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
p <- showDatabaseCategory(CellChatDB)
ggsave(filename = "showDatabaseCategory.pdf",plot = p,width = 21, height = 7, units = "in",dpi = 300)


unique(CellChatDB$interaction$annotation)#
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat)
###这一步报错了
future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used cellchatDB

cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 20)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")


cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2),xpd =TRUE)

pdf(file="number_weight_plot.pdf")
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Number of interactions")  
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name ="Interaction weights")
dev.off()

Number_of_interactions<- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T ,
                 label.edge= F, title.name = "Number of interactions")  
ggsave(Number_of_interactions,filename = "Number_of_interactions.pdf",width = 7, height = 7,units = "in",dpi = 300)


mat <- cellchat@net$count

par (mfrow = c(3,3),xpd = TRUE)
pdf(file="number_plot.pdf")
for (i in 1:nrow(mat)) {
  # i = 1 
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()


mat <- cellchat@net$weight
par (mfrow = c(3,3), xpd=T)
pdf(file="weight_plot.pdf")
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0,nrow = nrow(mat), ncol= ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
  
}
dev.off()

cellchat@netP$pathways ##you can check the all the passway 
pathways.show <- c("CD40")  ##choose one of the pass way 
pathways.show.all <- cellchat@netP$pathways##
levels(cellchat@idents)
##
vertex.receiver = c(1,2,3,4,5,6) 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy.pdf

##
levels(cellchat@idents)    
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
##
pathways.show <- c("IGF")
pdf("netVisual_heatmap_IGF.pdf",height = 14,width = 14)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("netVisual_chord_PTN.pdf",height = 14,width = 14)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

##bubble map
levels(cellchat@idents)    
#netVisual_bubble(cellchat,sources.use = c(1),targets.use = c(2,4,6,9),remove.isolate = FALSE,font.size =20,line.size = 0.8,font.size.title = 20)
pdf("netVisual_bubbblemap.pdf",height = 14,width = 14)
netVisual_bubble(cellchat,sources.use = c(1,2,3,4,5,6,7),targets.use = c(4),signaling = pathways.show.all,remove.isolate = FALSE)
dev.off()

pdf("netVisual_bubbblemap_2.pdf",height = 14,width = 14)
netVisual_bubble(cellchat,sources.use = c(4),targets.use = c(1,2,3,4,5,6,7),signaling = pathways.show.all,remove.isolate = FALSE)
dev.off()


pdf("netVisual_heatmap.pdf",height = 14,width = 14)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

pdf("Diff_number_strength_heatmap.pdf",height = 7,width = 10)
par(mfrow = c(1,2)) 
netVisual_heatmap(cellchat,color.heatmap = "Reds")
netVisual_heatmap(cellchat,measure = "weight",color.heatmap = "Reds")
dev.off()

##
netAnalysis_contribution(cellchat, signaling =pathways.show)
plotGeneExpression(cellchat,signaling = c("FSH","LHB"))
plotGeneExpression(cellchat,signaling = c("KIT","PTN"))
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling =pathways.show[1],  
                    vertex.receiver = vertex.receiver,layout = 'hierarchy')
save.image(file = "cellchat_named.RData")
par(mfrow = c(1,2)) 
netVisual_diffInteraction(cellchat, weight.scale = T) 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#做一个比较的出生前后的
named_seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/4sample_0.4_harmony/named_seurat_object.rds")
named_seurat_object <-sce
TCELL_ORIG<-named_seurat_object@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
named_seurat_object@meta.data[["orig.ident"]]<-TCELL_ORIG

named_seurat_object [["RNA3"]] <- as(object = named_seurat_object[["RNA"]], Class = "Assay")
DefaultAssay(named_seurat_object) <- "RNA3"
named_seurat_object[["RNA"]] <- NULL
named_seurat_object <- RenameAssays(object = named_seurat_object, RNA3 = 'RNA')

cellchat.E110 <- subset(named_seurat_object , orig.ident == "E110")
cellchat.E110$cell_annotations <- Idents(cellchat.E110)
cellchat.E110 <- createCellChat(cellchat.E110@assays$RNA@data, meta = cellchat.E110@meta.data, group.by = "celltype")
summary(cellchat.E110)
str(cellchat.E110)
levels(cellchat.E110@idents)                   
groupSize <- as.numeric(table(cellchat.E110@idents))#
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB) 
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() 
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#"Secreted Signaling"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat.E110@DB <- CellChatDB.use # set the used database in the object

cellchat.E110 <- subsetData(cellchat.E110)
future::plan("multiprocess", workers = 4)
cellchat.E110 <- identifyOverExpressedGenes(cellchat.E110)
cellchat.E110 <- identifyOverExpressedInteractions(cellchat.E110) #Identify over-expressed ligand-receptor interactions (pairs) within the used cellchatDB
cellchat.E110 <- projectData(cellchat.E110, PPI.human)

cellchat.E110 <- computeCommunProb(cellchat.E110, raw.use = FALSE, population.size = TRUE)#????????????????????????PPI??????????????????raw.use = TRUE?????????
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.E110 <- filterCommunication(cellchat.E110, min.cells = 10)
df.net <- subsetCommunication(cellchat.E110)
write.csv(df.net, "cellchat.E110_net_lr.csv")

cellchat.E110 <- computeCommunProbPathway(cellchat.E110)
df.netp <- subsetCommunication(cellchat.E110,slot.name = "netP")
write.csv(df.netp, "net_pathway_cellchat.E110.csv")

cellchat.E110 <- aggregateNet(cellchat.E110)


cellchat.P7 <- subset(named_seurat_object , orig.ident == "P7")
cellchat.P7$cell_annotations <- Idents(cellchat.P7)
cellchat.P7 <- createCellChat(cellchat.P7@assays$RNA@data, meta = cellchat.P7@meta.data, group.by = "celltype")
summary(cellchat.P7)
str(cellchat.P7)
levels(cellchat.P7@idents)                   
groupSize <- as.numeric(table(cellchat.P7@idents))#
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB) #
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new()
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#"Secreted Signaling"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat.P7@DB <- CellChatDB.use # set the used database in the object



cellchat.P7 <- subsetData(cellchat.P7)
future::plan("multiprocess", workers = 1)
cellchat.P7 <- identifyOverExpressedGenes(cellchat.P7)
cellchat.P7 <- identifyOverExpressedInteractions(cellchat.P7) #Identify over-expressed ligand-receptor interactions (pairs) within the used cellchatDB
cellchat.P7 <- projectData(cellchat.P7, PPI.human)

cellchat.P7 <- computeCommunProb(cellchat.P7, raw.use = FALSE, population.size = TRUE)#????????????????????????PPI??????????????????raw.use = TRUE?????????
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.P7 <- filterCommunication(cellchat.P7, min.cells = 10)
df.net <- subsetCommunication(cellchat.P7)
write.csv(df.net, "cellchat.p7_net_lr.csv")


cellchat.P7 <- computeCommunProbPathway(cellchat.P7)
df.netp <- subsetCommunication(cellchat.P7,slot.name = "netP")
write.csv(df.netp, "net_pathway_cellchat.P7.csv")
cellchat.P7 <- aggregateNet(cellchat.P7)

cellchat.E90 <- subset(named_seurat_object , orig.ident == "E90")
cellchat.E90$cell_annotations <- Idents(cellchat.E90)
cellchat.E90 <- createCellChat(cellchat.E90@assays$RNA@data, meta = cellchat.E90@meta.data, group.by = "celltype")
summary(cellchat.E90)
str(cellchat.E90)
levels(cellchat.E90@idents)                   
groupSize <- as.numeric(table(cellchat.E90@idents))#
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB) #
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new()
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#"Secreted Signaling"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat.E90@DB <- CellChatDB.use # set the used database in the object



cellchat.E90 <- subsetData(cellchat.E90)
#future::plan("multiprocess", workers = 1)
cellchat.E90 <- identifyOverExpressedGenes(cellchat.E90)
cellchat.E90 <- identifyOverExpressedInteractions(cellchat.E90) #Identify over-expressed ligand-receptor interactions (pairs) within the used cellchatDB
cellchat.E90 <- projectData(cellchat.E90, PPI.human)

cellchat.E90 <- computeCommunProb(cellchat.E90, raw.use = FALSE, population.size = TRUE)#????????????????????????PPI??????????????????raw.use = TRUE?????????
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.E90 <- filterCommunication(cellchat.E90, min.cells = 10)
df.net <- subsetCommunication(cellchat.E90)
write.csv(df.net, "cellchat.E90_net_lr.csv")


cellchat.E90 <- computeCommunProbPathway(cellchat.E90)
df.netp <- subsetCommunication(cellchat.E90,slot.name = "netP")
write.csv(df.netp, "net_pathway_cellchat.E90.csv")
cellchat.E90 <- aggregateNet(cellchat.E90)

object.list <- list(E90=cellchat.E90,E110 = cellchat.E110, P7 = cellchat.P7)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.cellchat#> An object of class CellChat created from a merged object with multiple datasets #>  555 signaling genes.#>  7563 cells.
#比较交互总数和交互强度
pdf(file="test_plot_counts_weight_2.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2
dev.off()

pdf(file="test_plot_counts_weight_LINE_2.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf(file="test_plot_heatmap_2.pdf",width = 16,height = 8)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file="test_plot_counts_LINE_E90_E110_P7_2.pdf")
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
netVisual_circle(cellchat, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name ="Interaction weights")

pdf(file="number_weight_plot.pdf")
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Number of interactions")  
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name ="Interaction weights")
dev.off()

pdf(file="test_plot_counts_LINE_E90_E110_P7_3.pdf")
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

pdf(file="test_PTPRM_P7.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.P7, signaling = pathways.show, color.heatmap = "Reds")
# save as TIL/CXCL_heatmap.pdf
dev.off()


# 分别绘制各时间点的热图
pdf(file = "timepoint_heatmaps.pdf", width = 18, height = 12)
for (i in 1:length(object.list)) {
  timepoint_name <- names(object.list)[i]
  cat("Plotting", timepoint_name, "\n")
  
  # 数量热图
  p_count <- netVisual_heatmap(object.list[[i]], 
                               title = paste(timepoint_name, " - Count"),
                               color.heatmap = "Reds")
  
  # 强度热图
  p_weight <- netVisual_heatmap(object.list[[i]], 
                                measure = "weight",
                                title = paste(timepoint_name, " - Weight"),
                                color.heatmap = "Reds")
  
  # 组合并绘制
  print(p_count + p_weight)
}
dev.off()



pdf(file="test_plot_leydigcellsource_2.pdf",width = 10,height = 35)
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = c(4),  comparison = c(1, 2,3), angle.x = 45)
#> Comparing communications on a merged object
dev.off()


pdf(file="test_plot_leydigcelltargets_2.pdf",width = 10,height = 35)
netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(1,2,3,4,5,6,7),  comparison = c(1, 2,3), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

pdf(file="test_plot_stromal_celltargets_1.pdf",width = 10,height = 20)
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = c(7),  comparison = c(1, 2,3), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

#cellchat气泡图优化
cellchat@netP$pathways

netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(7), 
                 signaling = c("PTN","MK"), remove.isolate = FALSE)

pdf("bubble_cellchat_stromal.pdf",width = 14,height = 20) 
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = c(7), 
                 signaling = c("FGF","TGFb","NOTCH","PDGF", "COLLAGEN","LAMININ"), remove.isolate = FALSE, comparison = c(1, 2,3), angle.x = 45)
dev.off()
pdf("bubble_cellchat_PTN_MK_5.pdf") 
netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(4), 
                 signaling = c("MK","PTN","PTPRM"), remove.isolate = FALSE, comparison = c(1, 2,3), angle.x = 45)
dev.off()
pdf("bubble_cellchat_PTN_MK_6.pdf") 
netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(1,2,3,4,5,6,7), 
                 signaling = c("PTPRM"), remove.isolate = FALSE, comparison = c(1, 2,3), angle.x = 45)
dev.off()
pathways.show <- c("PTPRM")
pdf("netVisual_chord_PTPRM_E90.pdf",height = 14,width = 14)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.E90, signaling = pathways.show, layout = "chord")
dev.off()

pdf("netVisual_chord_PTPRM_E110.pdf",height = 14,width = 14)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.E110, signaling = pathways.show, layout = "chord")
dev.off()

pdf("netVisual_chord_PTPRM_P7.pdf",height = 14,width = 14)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.P7, signaling = pathways.show, layout = "chord")
dev.off()

levels(cellchat@idents[["E90"]])
cellchat@netP[["E90"]][["pathways"]]

save.image(file="cellchat.Rdata")




#重新整一下气泡图，查看双细胞的情况
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/0.4_unnamed_8sample.rds")
gene_to_check3=c('UCHL1','DMRT1','FANCD2','MSH2','FANCI','STRA8','SPO11','MEIOC','SYCP2','SYCE3','SPAG6','ACRV1','TNP1','PRM2','PRM3','SOX9','FATE1','PRND','INHA','COL1A1','COL3A1','DCN','ACTA2','MYH11','CYP26B1','CYP11A1','CYP17A1','STAR','INSL3','VWF','PECAM1','CD93','CD34','NOTCH3','RGS5','ITGB1','PTPRC','CD163','C1QB','TYROBP','CTSB','CCL14','CCL5','GZMB','CTSW')
gene_to_check3 <- as.data.frame(gene_to_check3)
markerdata<-ScaleData(seurat_object,features = as.character(unique(gene_to_check3$gene_to_check3)))
pdf("气泡图_unnamed.pdf") 
DotPlot(seurat_object, features = as.character(unique(gene_to_check3$gene_to_check3)))+coord_flip()+ theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+theme(axis.text.x = element_text(angle = 45, vjust = 1))
dev.off()

table(seurat_object$orig.ident)
prop.table(table(Idents(seurat_object)))
table(Idents(seurat_object), seurat_object$orig.ident)

Tcell <- seurat_object
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
seurat_object <- Tcell  
p1<- DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(seurat_object, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "named_umap.png",width = 14, height = 7,units = "in",dpi = 300)

P7 <- subset(seurat_object, subset = orig.ident=="P7")
E90 <- subset(seurat_object, subset = orig.ident=="E90")
E110 <- subset(seurat_object, subset = orig.ident=="E110")

pdf("umap_0.4_3_time_8sample.pdf") 
DimPlot(P7, reduction = "umap",pt.size = 0.5,label = T)+ggtitle('P7')
DimPlot(E90, reduction = "umap",pt.size = 0.5,label = T)+ggtitle('E90')
DimPlot(E110, reduction = "umap",pt.size = 0.5,label = T)+ggtitle('E110')
dev.off()

table(seurat_object$orig.ident)
prop.table(table(Idents(seurat_object)))
table(Idents(seurat_object), seurat_object$orig.ident)

SOX9 <- FeaturePlot(E90, features = c('SOX9'),label=TRUE,order = TRUE)
ggsave(SOX9,filename = "SOX9_E90.png",width = 7, height = 7,units = "in",dpi = 300)
SOX9 <- FeaturePlot(E110, features = c('SOX9'),label=TRUE,order = TRUE)
ggsave(SOX9,filename = "SOX9_E110.png",width = 7, height = 7,units = "in",dpi = 300)
SOX9 <- FeaturePlot(P7, features = c('SOX9'),label=TRUE,order = TRUE)
ggsave(SOX9,filename = "SOX9_P7.png",width = 7, height = 7,units = "in",dpi = 300)

CYP17A1 <- FeaturePlot(E90, features = c('CYP17A1'),label=TRUE,order = TRUE)
ggsave(CYP17A1,filename = "CYP17A1_E90.png",width = 7, height = 7,units = "in",dpi = 300)
CYP17A1 <- FeaturePlot(E110, features = c('CYP17A1'),label=TRUE,order = TRUE)
ggsave(CYP17A1,filename = "CYP17A1_E110.png",width = 7, height = 7,units = "in",dpi = 300)
CYP17A1 <- FeaturePlot(P7, features = c('CYP17A1'),label=TRUE,order = TRUE)
ggsave(CYP17A1,filename = "CYP17A1_P7.png",width = 7, height = 7,units = "in",dpi = 300)

CYP11A1 <- FeaturePlot(E90, features = c('CYP11A1'),label=TRUE,order = TRUE)
ggsave(CYP11A1,filename = "CYP11A1_E90.png",width = 7, height = 7,units = "in",dpi = 300)
CYP11A1 <- FeaturePlot(E110, features = c('CYP11A1'),label=TRUE,order = TRUE)
ggsave(CYP11A1,filename = "CYP11A1_E110.png",width = 7, height = 7,units = "in",dpi = 300)
CYP11A1 <- FeaturePlot(P7, features = c('CYP11A1'),label=TRUE,order = TRUE)
ggsave(CYP11A1,filename = "CYP11A1_P7.png",width = 7, height = 7,units = "in",dpi = 300)

STAR <- FeaturePlot(E90, features = c('STAR'),label=TRUE,order = TRUE)
ggsave(STAR,filename = "STAR_E90.png",width = 7, height = 7,units = "in",dpi = 300)
STAR <- FeaturePlot(E110, features = c('STAR'),label=TRUE,order = TRUE)
ggsave(STAR,filename = "STAR_E110.png",width = 7, height = 7,units = "in",dpi = 300)
STAR <- FeaturePlot(P7, features = c('STAR'),label=TRUE,order = TRUE)
ggsave(STAR,filename = "STAR_P7.png",width = 7, height = 7,units = "in",dpi = 300)

#修改细胞占比图
install.packages("RImagePalette")

#提取细胞在细胞亚群和样本的分布数据,细胞堆叠图(整三个时期的)
seurat_object<-sce 
Idents(seurat_object) <- "celltype"
TCELL_ORIG<-seurat_object@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
seurat_object@meta.data[["orig.ident"]]<-TCELL_ORIG
seurat_object[["celltype"]]<-Idents(seurat_object)
df<-table(seurat_object$orig.ident,seurat_object$celltype)
df<-as.data.frame(df)
colnames(df)<-c('sample','celltype','Freq')
table(df$sample)

#如果你使用的是tidyverse包，也可以使用以下代码：
A <- df %>%
  group_by(sample) %>%
  mutate(total = sum(Freq)) %>%
  mutate(ratio = Freq / total) %>%
  ungroup() %>%
  print()
A <- A[,-3:-4]
colnames(A)[colnames(A) == "ratio"] <- "Freq"
colnames(A)[colnames(A) == "sample"] <- "stage"
tem <- A$stage
A[,1] <-A[,2]
A[,2]<-tem
colnames(A)[1] <- "celltype"
colnames(A)[2] <- "stage"
#更换行从E90  E110  P7的顺序来
tem_E110 <-A[c(1,4,7,10,13,16,19),]
tem_E90 <-A[c(2,5,8,11,14,17,20),]
tem_P7 <-A[c(3,6,9,12,15,18,21),]
rm(B)
B <- rbind(tem_E90,tem_E110,tem_P7)
B <- as.data.frame(B)
B$stage <- factor(B$stage, levels = c("E90", "E110", "P7"))
colnames(B) <- c("celltype", "stage", "Freq")
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

p1<-ggplot(B,aes(x = stage,y =Freq,
             group=celltype))+
  stat_summary(geom = 'line',fun='mean',cex=1,col='white')+
  geom_area(data =B,aes(fill=celltype))+
  scale_fill_manual(values=cluster_cols)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text = element_text(color = "black",size = 10))+
  geom_vline(aes(xintercept ="E90"),linetype="dashed", size=1.2, colour="white")+
  geom_vline(aes(xintercept ="E110"),linetype="dashed", size=1.2, colour="white")+
  geom_vline(aes(xintercept ="P7"),linetype="dashed", size=1.2, colour="white")
p1
ggsave(p1,filename = "xibaozhanbi.pdf",width = 7, height = 7,units = "in",dpi = 300)

#cellchat气泡图优化
cellchat@netP$pathways

netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(3), 
                 signaling = c("PTN","MK"), remove.isolate = FALSE)

pdf("bubble_cellchat_PTN_MK.pdf") 
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(3), 
                 signaling = c("PTN","MK"), remove.isolate = FALSE)
dev.off()


#检测双细胞
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90")
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/seurat_object_named_0.4_harmony.Rds")
#seurat_object[["RNA3"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
#DefaultAssay(seurat_object) <- "RNA3"
scRNA_harmony <- seurat_object
rm(seurat_object)
options(future.globals.maxSize = 600 * 1024^2)
sweep.res.list <- paramSweep_v3(scRNA_harmony, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(scRNA_harmony$seurat_clusters) 
nExp_poi <- round(0.075 *nrow(scRNA_harmony@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量
scRNA_harmony <- doubletFinder_v3(scRNA_harmony, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                                  nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
DimPlot(scRNA_harmony, reduction = "UMAP", group.by = "DF.classifications_0.25_33_9717", raster = FALSE)
##太大了，爆内存 使用服务器上跑
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/Stroma_cell")
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/0.4_unnamed_8sample.rds")
##重新聚类基质细胞：提取基质细胞
Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(0,3,6,13)]
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG

rm(seurat_object)
testA <- subset(Tcell,orig.ident=="E90")
testA.seu <- testA[["RNA"]]$counts
testA.seu=CreateSeuratObject(counts = testA.seu,project = "E90", min.cells = 3, min.features = 200)
testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)

testB <- subset(Tcell,orig.ident=="E110")
testB.seu <- testB[["RNA"]]$counts
testB.seu=CreateSeuratObject(counts = testB.seu,project = "E110", min.cells = 3, min.features = 200)
testB.seu <- NormalizeData(testB.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testB.seu <- FindVariableFeatures(testB.seu, selection.method = "vst", nfeatures = 2000)


testC <- subset(Tcell,orig.ident=="P7")
testC.seu <- testC[["RNA"]]$counts
testC.seu=CreateSeuratObject(counts = testC.seu,project = "P7", min.cells = 3, min.features = 200)
testC.seu <- NormalizeData(testC.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testC.seu <- FindVariableFeatures(testC.seu, selection.method = "vst", nfeatures = 2000)


testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu,testC.seu), dims = 1:20)#可设置不同的k.anchor = 5,k.filter = 30参数，这里使用默认的k.filter =200
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:20)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.1)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:20)
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_stromal.png",width = 14, height = 7,units = "in",dpi = 300)

Tcell <- testAB.integrated
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
testAB.integrated <- Tcell  
p1<- DimPlot(Tcell, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(Tcell, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_Leydigs_2.png",width = 14, height = 7,units = "in",dpi = 300)
save.image(file="stromal_cell.Rdata")
pdf("stromal_0.1.pdf") 
DimPlot(Tcell, reduction = "umap",pt.size = 0.5,label = T)
DimPlot(Tcell, reduction = "umap",label=TRUE,pt.size = 0.5,group.by = "orig.ident")
dev.off()
#看一下文献中间质细胞干细胞的marker 
#间质细胞干细胞 ITGAV
CD51 <- FeaturePlot(Tcell, features = c('CD51'),label=TRUE,order = TRUE)
ggsave(CD51,filename = "CD51.png",width = 7, height = 7,units = "in",dpi = 300)

ITGAV <- FeaturePlot(Tcell, features = c('ITGAV'),label=TRUE,order = TRUE)
ggsave(ITGAV,filename = "ITGAV.png",width = 7, height = 7,units = "in",dpi = 300)

HSD3B <- FeaturePlot(Tcell, features = c('HSD3B7'),label=TRUE,order = TRUE)
ggsave(HSD3B,filename = "HSD3B.png",width = 7, height = 7,units = "in",dpi = 300)

PDGFRA <- FeaturePlot(Tcell, features = c('PDGFRA'),label=TRUE,order = TRUE)
ggsave(PDGFRA,filename = "PDGFRA.png",width = 7, height = 7,units = "in",dpi = 300)

NES <-FeaturePlot(Tcell, features = c('NES'),label=TRUE,order = TRUE)
ggsave(NES,filename = "NES.png",width = 7, height = 7,units = "in",dpi = 300)

NR2F2 <-FeaturePlot(Tcell, features = c('NR2F2'),label=TRUE,order = TRUE)
ggsave(NR2F2,filename = "NR2F2.png",width = 7, height = 7,units = "in",dpi = 300)

P75NTR <-FeaturePlot(Tcell, features = c('P75NTR'),label=TRUE,order = TRUE)
ggsave(P75NTR,filename = "P75NTR.png",width = 7, height = 7,units = "in",dpi = 300)

THY1 <-FeaturePlot(Tcell, features = c('THY1'),label=TRUE,order = TRUE)
ggsave(THY1,filename = "THY1.png",width = 7, height = 7,units = "in",dpi = 300)

##查看IGF1的表达情况
IGF1 <-FeaturePlot(seurat_object, features = c('IGF1'),label=TRUE,order = TRUE)
ggsave(IGF1,filename = "IGF1.png",width = 7, height = 7,units = "in",dpi = 300)

#支持细胞亚群分析
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/sertoli_cell")
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/0.4_unnamed_8sample.rds")
##重新聚类基质细胞：提取基质细胞
Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(1,4,5,10,16)]
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG

rm(seurat_object)
testA <- subset(Tcell,orig.ident=="E90")
testA.seu <- testA[["RNA"]]$counts
testA.seu=CreateSeuratObject(counts = testA.seu,project = "E90", min.cells = 3, min.features = 200)
testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)

testB <- subset(Tcell,orig.ident=="E110")
testB.seu <- testB[["RNA"]]$counts
testB.seu=CreateSeuratObject(counts = testB.seu,project = "E110", min.cells = 3, min.features = 200)
testB.seu <- NormalizeData(testB.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testB.seu <- FindVariableFeatures(testB.seu, selection.method = "vst", nfeatures = 2000)


testC <- subset(Tcell,orig.ident=="P7")
testC.seu <- testC[["RNA"]]$counts
testC.seu=CreateSeuratObject(counts = testC.seu,project = "P7", min.cells = 3, min.features = 200)
testC.seu <- NormalizeData(testC.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testC.seu <- FindVariableFeatures(testC.seu, selection.method = "vst", nfeatures = 2000)


testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu,testC.seu), dims = 1:20)#可设置不同的k.anchor = 5,k.filter = 30参数，这里使用默认的k.filter =200
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:20)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.08)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:20)
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_sertoli.png",width = 14, height = 7,units = "in",dpi = 300)

Tcell <- testAB.integrated
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
testAB.integrated <- Tcell  
p1<- DimPlot(Tcell, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(Tcell, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_sertoli_2.png",width = 14, height = 7,units = "in",dpi = 300)
save.image(file="sertoli_cell.Rdata")
pdf("sertoli_0.08.pdf") 
DimPlot(Tcell, reduction = "umap",pt.size = 0.5,label = T)
DimPlot(Tcell, reduction = "umap",label=TRUE,pt.size = 0.5,group.by = "orig.ident")
dev.off()

#支持细胞亚群分析
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/sertoli_cell")
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/0.4_unnamed_8sample.rds")
Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(1,4,5,10,16)]
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG

rm(seurat_object)
testA <- subset(Tcell,orig.ident=="E90")
testA.seu <- testA[["RNA"]]$counts
testA.seu=CreateSeuratObject(counts = testA.seu,project = "E90", min.cells = 3, min.features = 200)
testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)

testB <- subset(Tcell,orig.ident=="E110")
testB.seu <- testB[["RNA"]]$counts
testB.seu=CreateSeuratObject(counts = testB.seu,project = "E110", min.cells = 3, min.features = 200)
testB.seu <- NormalizeData(testB.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testB.seu <- FindVariableFeatures(testB.seu, selection.method = "vst", nfeatures = 2000)


testC <- subset(Tcell,orig.ident=="P7")
testC.seu <- testC[["RNA"]]$counts
testC.seu=CreateSeuratObject(counts = testC.seu,project = "P7", min.cells = 3, min.features = 200)
testC.seu <- NormalizeData(testC.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testC.seu <- FindVariableFeatures(testC.seu, selection.method = "vst", nfeatures = 2000)


testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu,testC.seu), dims = 1:20)#可设置不同的k.anchor = 5,k.filter = 30参数，这里使用默认的k.filter =200
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:20)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.08)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:20)
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_sertoli.png",width = 14, height = 7,units = "in",dpi = 300)

Tcell <- testAB.integrated
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
testAB.integrated <- Tcell  
p1<- DimPlot(Tcell, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(Tcell, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_sertoli_2.png",width = 14, height = 7,units = "in",dpi = 300)
save.image(file="sertoli_cell.Rdata")
pdf("sertoli_0.08.pdf") 
DimPlot(Tcell, reduction = "umap",pt.size = 0.5,label = T)
DimPlot(Tcell, reduction = "umap",label=TRUE,pt.size = 0.5,group.by = "orig.ident")
dev.off()

library(monocle)
library(Seurat)
library(AUCell)
library(patchwork)
library(ggplot2)
library(DOSE)#人疾病
library(clusterProfiler)#富集分析
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Rn.eg.db)
#library(org.At.tair.db)
library(dplyr)
library(GO.db)
library(rjson)
library(stringr)
library("reshape2")
#dim(Tcell)
#Tcell <- subset(Tcell, orig.ident == "E110.1")
#Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(0,2,20,5,1,14,9,10,3,4)]
#sertoli cell

# 为了减少内存占用，可以保存提取数据的seurat对象，并且删除原来的seurat对象
#save(Tcell, file = "sertolinamed.Rda")
#rm(Tcell)
gc()
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

CDS <- estimateDispersions(CDS)
gc()

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS, file = "CDS_sertoli.Rdata")
rm(list = ls())
gc()
load("CDS_sertoli.Rdata")
dim(CDS)

################################################################################
## 2. 差异分析寻找高变基因
# detectGenes()函数：同时统计表达当前基因的细胞数量和细胞表达的基因数；
# min_expr参数用于设置检测阈值，比如min_expr = 0.1表示当前基因的表达量超过0.1才会
# 纳入统计；
CDS <- detectGenes(CDS, min_expr = 0.1)

#过滤低表达的基因，以降低噪音和减少差异分析计算量，fdata中细胞数小于20去除，大于保留
expressed_genes <- fData(CDS) %>%
  subset(num_cells_expressed >= 20) %>%
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

################################################################################
## 3. 时间轨迹及其差异分析
## 3.1 构建细胞轨迹
# 第一步：选择用于构建细胞轨迹的基因集；
# 选择Top2000差异基因作为排序基因；
ordering_genes <- rownames(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]

#fData(CDS)$use_for_ordering <-fData(CDS)$num_cells_expressed > 0.1 * ncol(CDS)

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
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.png", width = 7, height = 7,units = "in",dpi = 300)
#按“State”分组；
State_clusters <-plot_cell_trajectory(CDS, color_by = "State")
ggsave(State_clusters,filename = "leydig_State.png", width = 7, height = 7,units = "in",dpi = 300)

#按seurat分群结果分组
clusters <-plot_cell_trajectory(CDS, color_by = "seurat_clusters")
ggsave(clusters,filename = "leydig_clusters.png", width = 7, height = 7,units = "in",dpi = 300)
#按细胞类型分组
save.image(file="sertoli_cell_monocle.Rdata")

#间质细胞提取
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("E:/Bowen-files/pig-process/merged/peixun/addT90/leydig_cell")
seurat_object <- readRDS("E:/Bowen-files/pig-process/merged/peixun/addT90/0.4_unnamed_8sample.rds")
Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(2,18,20)]
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG

rm(seurat_object)
testA <- subset(Tcell,orig.ident=="E90")
testA.seu <- testA[["RNA"]]$counts
testA.seu=CreateSeuratObject(counts = testA.seu,project = "E90", min.cells = 3, min.features = 200)
testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)

testB <- subset(Tcell,orig.ident=="E110")
testB.seu <- testB[["RNA"]]$counts
testB.seu=CreateSeuratObject(counts = testB.seu,project = "E110", min.cells = 3, min.features = 200)
testB.seu <- NormalizeData(testB.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testB.seu <- FindVariableFeatures(testB.seu, selection.method = "vst", nfeatures = 2000)


testC <- subset(Tcell,orig.ident=="P7")
testC.seu <- testC[["RNA"]]$counts
testC.seu=CreateSeuratObject(counts = testC.seu,project = "P7", min.cells = 3, min.features = 200)
testC.seu <- NormalizeData(testC.seu, normalization.method = "LogNormalize", scale.factor = 10000)
testC.seu <- FindVariableFeatures(testC.seu, selection.method = "vst", nfeatures = 2000)


testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu,testC.seu), dims = 1:20)#可设置不同的k.anchor = 5,k.filter = 30参数，这里使用默认的k.filter =200
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:20)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.08)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:20)
p1<- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(testAB.integrated, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_leydig.png",width = 14, height = 7,units = "in",dpi = 300)

Tcell <- testAB.integrated
TCELL_ORIG<-Tcell@meta.data[["orig.ident"]]
TCELL_ORIG <- substr(TCELL_ORIG,1,nchar(TCELL_ORIG)-2)
Tcell@meta.data[["orig.ident"]]<-TCELL_ORIG
testAB.integrated <- Tcell  
p1<- DimPlot(Tcell, reduction = "umap", group.by = "orig.ident")
p2<- DimPlot(Tcell, reduction = "umap", label=TRUE,group.by = "ident")
unnamedumap <-  p1+p2
ggsave(unnamedumap,filename = "unnamed_umap_leydig_2.png",width = 14, height = 7,units = "in",dpi = 300)
save.image(file="leydig_cell.Rdata")
pdf("leydig_亚群.pdf") 
DimPlot(Tcell, reduction = "umap",pt.size = 0.5,label = T)
DimPlot(Tcell, reduction = "umap",label=TRUE,pt.size = 0.5,group.by = "orig.ident")
dev.off()
#间质细胞的monocle
library(monocle)
library(Seurat)
library(AUCell)
library(patchwork)
library(ggplot2)
library(DOSE)#人疾病
library(clusterProfiler)#富集分析
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Rn.eg.db)
#library(org.At.tair.db)
library(dplyr)
library(GO.db)
library(rjson)
library(stringr)
library("reshape2")
#dim(Tcell)
#Tcell <- subset(Tcell, orig.ident == "E110.1")
#Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(0,2,20,5,1,14,9,10,3,4)]
#sertoli cell

# 为了减少内存占用，可以保存提取数据的seurat对象，并且删除原来的seurat对象
#save(Tcell, file = "sertolinamed.Rda")
#rm(Tcell)
gc()
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

CDS <- estimateDispersions(CDS)
gc()

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS, file = "CDS_Leydig_2.Rdata")
rm(list = ls())
gc()
load("CDS_Leydig_2.Rdata")
dim(CDS)

################################################################################
## 2. 差异分析寻找高变基因
# detectGenes()函数：同时统计表达当前基因的细胞数量和细胞表达的基因数；
# min_expr参数用于设置检测阈值，比如min_expr = 0.1表示当前基因的表达量超过0.1才会
# 纳入统计；
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

################################################################################
## 3. 时间轨迹及其差异分析
## 3.1 构建细胞轨迹
# 第一步：选择用于构建细胞轨迹的基因集；
# 选择Top2000差异基因作为排序基因；
ordering_genes <- rownames(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]

#fData(CDS)$use_for_ordering <-fData(CDS)$num_cells_expressed > 0.1 * ncol(CDS)
genes <- c("BTG2", "CEBPB", "DNAJA1", "KLF4", "TSC22D1", "ALDH1A1", "APOE", "CARD19", "CLU", "DUSP1", "EGR1", "INHBA", "JUN", "NR4A1", "RHOB", "TIMP3", "VNN1", "ZFP36","PCNA", "TOP2A", "MCM6", "MKI67","CYP17A1","CYP11A1","STAR","NR5A1","PDGFRA","NR2F2","ITGAV")
#将差异基因储存到CDS对象中；
sig_genes <- union(ordering_genes, genes)
CDS <- setOrderingFilter(CDS, ordering_genes = sig_genes)

#第二步: 数据降维
#降维函数与上文的t-SNE一致，但降维算法这里用的是“DDRTree”，
CDS <- reduceDimension(CDS, method = 'DDRTree')

#第三步: 构建细胞分化轨迹
#按照分化轨迹排序细胞；
CDS <- orderCells(CDS)

#绘制细胞分化轨迹：
#按“Pseudotime”分组；
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")
ggsave(Pseudotime_clusters,filename = "leydig_2_ProPseudotime_200gene.png", width = 7, height = 7,units = "in",dpi = 300)
#按“State”分组；
State_clusters <-plot_cell_trajectory(CDS, color_by = "State")
ggsave(State_clusters,filename = "leydig_2_State_200gene.png", width = 7, height = 7,units = "in",dpi = 300)

#按seurat分群结果分组
clusters <-plot_cell_trajectory(CDS, color_by = "seurat_clusters")
ggsave(clusters,filename = "leydig_2_clusters.png", width = 7, height = 7,units = "in",dpi = 300)
save.image(file="leydig_2_200genes_monocle.Rdata")
#按细胞类型分组
plot_cell_trajectory(CDS, color_by = "celltype")
#按照ori_ident分类
clusters <-plot_cell_trajectory(CDS, color_by = "orig.ident")
ggsave(clusters,filename = "leydig_2_200genesorig.ident.png", width = 7, height = 7,units = "in",dpi = 300)

save(CDS, file = "CDS_leydig_pseudotime.Rda")

##3.2 比较细胞分化轨迹进程中功能基因的表达差异
#主要用到sm.ns()函数根据表达量拟合曲线；用拟时序做基因差异表达分析
diff_test_res <- differentialGeneTest(
  CDS[expressed_genes, ], 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
head(diff_test_res[, c("gene_short_name", "pval", "qval")])
save.image(file = "monocle_leydig.RData")
# 按q值从小到大排序后，查看最显著的前6个基因的拟时间趋势
sig_gene_names1 <- rownames(diff_test_res[order(diff_test_res$qval)[1:4], ])
seurat_cluster <- plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'Pseudotime')
ggsave(seurat_cluster,filename = "Leydig_seurat_cluster.png", width = 14, height = 14,units = "in",dpi = 300)
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

# 也可以查看比较关注基因的拟时间趋势
sig_gene_names1 <- c("DDX4", "UCHL1", "DNMT3L","PCNA","PIWIL4","ETV4","ETV5","KDM1B","GFRA1")
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

sig_gene_names2 <- c('CYP17A1','CYP11A1','INSL3','STAR')
Leydigmarker <- plot_genes_in_pseudotime(CDS[sig_gene_names2, ], color_by = 'State')
ggsave(Leydigmarker,filename = "pseudotime_Leydigmarker.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$PDGFRA=log2(exprs(CDS)['PDGFRA',]+1)
pseudotime_PDGFRA <- plot_cell_trajectory(CDS,color_by="PDGFRA", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_PDGFRA,filename = "pseudotime_PDGFRA.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$CYP11A1=log2(exprs(CDS)['CYP11A1',]+1)
pseudotime_CYP11A1 <- plot_cell_trajectory(CDS,color_by="CYP11A1", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_CYP11A1,filename = "pseudotime_CYP11A1.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$STAR=log2(exprs(CDS)['STAR',]+1)
pseudotime_STAR <- plot_cell_trajectory(CDS,color_by="STAR", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_STAR,filename = "pseudotime_STAR.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$INSL3=log2(exprs(CDS)['INSL3',]+1)
pseudotime_INSL3 <- plot_cell_trajectory(CDS,color_by="INSL3", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_INSL3,filename = "pseudotime_INSL3.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$KIT=log2(exprs(CDS)['KIT',]+1)
pseudotime_KIT <- plot_cell_trajectory(CDS,color_by="KIT", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_KIT,filename = "pseudotime_KIT.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$DNMT3L=log2(exprs(CDS)['DNMT3L',]+1)
pseudotime_DNMT3L <- plot_cell_trajectory(CDS,color_by="DNMT3L", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_DNMT3L,filename = "pseudotime_DNMT3L.png",width = 7, height = 7,units = "in",dpi = 300)


#更换根节点：#按照分化轨迹排序细胞；
CDS <- orderCells(CDS,root_state = 2)
#按“Pseudotime”分组；
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")+scale_color_gsea()+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+scale_color_gradient2(low="grey",mid="orange",high="blue")
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.png", width = 7, height = 7,units = "in",dpi = 300)


##3.4 差异基因的拟时表达模式聚类分析
#提取差异基因；
sig_gene_names2 <- row.names(subset(diff_test_res, qval < 0.001))
#绘制拟时间差异基因表达谱热图；#重新聚类
plot_pseudotime_heatmap(
  CDS[sig_gene_names2, ], 
  num_clusters = 3, #cluster数量
  show_rownames = F
)

Time_diff <- differentialGeneTest(CDS[sig_gene_names2,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
write.csv(Time_diff,"Tiff_diff_all.csv",row.names = F)

Time_genes <- top_n(Time_diff, n = 2000, desc(qval)) %>% pull(gene_short_name) %>% as.character()  

#基因太多了画不出来就选择前2000个基因，所有基因Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p<-plot_pseudotime_heatmap(CDS[sig_gene_names2,], num_clusters=3, show_rownames=F, return_heatmap=T)  #字体太大，看不清，可以导出
ggsave("Time_heatmap.png", p, width = 5, height = 8)
ggsave("Time_heatmap.pdf", p, width = 5, height = 10)

clusters<-cutree(p$tree_row,k=3)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
write.csv(clustering,"间质细胞2000个拟时序基因热图分3个cluster.csv")
save.image(file = "间质细胞-monocle2.Rdata") 
################################################################################
## 4.单细胞轨迹分支分析 
# 当细胞分化轨迹出现分支的时候，意味着细胞将面临不同的分化命运“fate”，接下来主要
# 分析分支事件，比如沿着分化轨迹，基因的表达量如何变化？
# 不同分支之间的差异基因有哪些？

# Monocle 提供一个特定的统计检验方法: branched expression analysis modeling（BEAM）.
## 4.1 BEAM检验
# 使用BEAM()函数对基因进行注释；
# BEAM函数的输入对象： 完成拟时间分析的CellDataSet且轨迹中包含1个分支点；
# 返回一个包含每个基因significance scores 的表格,若得分是显著的则表明该基因的表
# 达是与分支相关的（branch-dependent）。

BEAM_res <- BEAM(CDS[expressed_genes[1:200], ], branch_point = 1,cores=10,progenitor_method = "duplicate")
saveRDS(BEAM_res, file = "BEAM_res.RDS")

#按照qval升序排列；
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
head(BEAM_res)

##4.2轨迹分支表达分析
#使用pheatmap包绘制分支表达量热图；
sig_gene_names3 <- row.names(subset(BEAM_res, qval < 1e-5))
plot_genes_branched_heatmap(
  CDS[sig_gene_names3, ],
  branch_point = 1,#设置分支节点
  num_clusters = 3,
  use_gene_short_name = FALSE,
  show_rownames = FALSE)

#使用 plot_genes_branched_pseudotime() 函数绘制拟合曲线，不同线表示不同分支
plot_genes_branched_pseudotime(
  CDS[rownames(BEAM_res)[1:5], ], 
  branch_point = 1, 
  color_by = "State",
  ncol = 1
)
write.csv(sig_gene_names3,file="BEAM_res.csv")
-#间质细胞的周期性
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# 对细胞周期基因做PCA，揭示了细胞完全按阶段分离
Tcell <- RunPCA(Tcell, features = c(s.genes, g2m.genes))
DimPlot(Tcell)
Tcell <- CellCycleScoring(Tcell, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 查看细胞周期分数和阶段分配
head(Tcell[[]])

# 可视化整个细胞周期标记的分布
RidgePlot(Tcell, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# 对细胞周期基因做PCA，揭示了细胞完全按阶段分离
DimPlot(Tcell)
p1 <- DimPlot(Tcell, reduction = "umap", shuffle = T, pt.size = 0.5)
ggsave(p1,filename = "周期.png",width = 7, height = 7,units = "in",dpi = 300) 

RidgePlot(Tcell, 
          features = c("PCNA", "TOP2A", "C", "MKI67"), 
          cols = pal_npg("nrc", alpha = 0.7)(3),
          ncol = 2)
PCNA <- FeaturePlot(Tcell, features = c('PCNA'),label=TRUE,order = TRUE)
ggsave(PCNA,filename = "PCNA.png",width = 7, height = 7,units = "in",dpi = 300)

TOP2A <- FeaturePlot(Tcell, features = c('TOP2A'),label=TRUE,order = TRUE)
ggsave(TOP2A,filename = "TOP2A.png",width = 7, height = 7,units = "in",dpi = 300)

MCM6 <- FeaturePlot(Tcell, features = c('MCM6'),label=TRUE,order = TRUE)
ggsave(MCM6,filename = "MCM6.png",width = 7, height = 7,units = "in",dpi = 300)

save.image(file="leydig_cell.Rdata")

pData(CDS)$PCNA=log2(exprs(CDS)['PCNA',]+1)
pseudotime_PCNA <- plot_cell_trajectory(CDS,color_by="PCNA", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_PCNA,filename = "pseudotime_PCNA.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$TOP2A=log2(exprs(CDS)['TOP2A',]+1)
pseudotime_TOP2A <- plot_cell_trajectory(CDS,color_by="TOP2A", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_TOP2A,filename = "pseudotime_TOP2A.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$MCM6=log2(exprs(CDS)['MCM6',]+1)
pseudotime_MCM6 <- plot_cell_trajectory(CDS,color_by="MCM6", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_MCM6,filename = "pseudotime_MCM6.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$WNT5A=log2(exprs(CDS)['WNT5A',]+1)
pseudotime_WNT5A <- plot_cell_trajectory(CDS,color_by="WNT5A", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_WNT5A,filename = "pseudotime_WNT5A.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$NG2=log2(exprs(CDS)['NG2',]+1)
pseudotime_NG2 <- plot_cell_trajectory(CDS,color_by="NG2", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_NG2,filename = "pseudotime_NG2.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$EGR1=log2(exprs(CDS)['EGR1',]+1)
pseudotime_EGR1 <- plot_cell_trajectory(CDS,color_by="EGR1", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_EGR1,filename = "pseudotime_EGR1.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$DLK1=log2(exprs(CDS)['DLK1',]+1)
pseudotime_DLK1 <- plot_cell_trajectory(CDS,color_by="DLK1", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_DLK1,filename = "pseudotime_DLK1.png",width = 7, height = 7,units = "in",dpi = 300)
gene_names <- c("PCNA", "TOP2A", "MCM6", "MKI67")
# 遍历每个基因
for (gene in features) {
  # 计算基因的 log2 表达值
  pData(CDS)[[gene]] <- log2(exprs(CDS)[gene, ] + 1)
  # 绘制细胞轨迹图
  pseudotime_plot <- plot_cell_trajectory(CDS, color_by = gene, size = 1, show_backbone = TRUE) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
    scale_color_gradient2(low = 'grey', mid = 'orange', high = 'blue')
  # 保存图片
  filename <- paste0("pseudotime_", gene, ".png")
  ggsave(pseudotime_plot, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}
for (gene in features) {
  # 计算基因的 log2 表达值
  pData(CDS)[[gene]] <- exprs(CDS)[gene, ]
  # 绘制细胞轨迹图
  pseudotime_plot <- plot_cell_trajectory(CDS, color_by = gene, size = 1, show_backbone = TRUE) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
    scale_color_gradient2(low = 'grey', mid = 'orange', high = 'blue')
  # 保存图片
  filename <- paste0("pseudotime_200genes", gene, ".png")
  ggsave(pseudotime_plot, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}

# 定义基因列表
genes <- c("PECAM1","VWF","CD34","NOS3", "TEK","KDR","ICAM2","FLT1","TIE1","PTPRB",#血管内皮 Endothelial Cells
           "ACTA2","CNN1","MYH11","MYLK","TAGLN","SNCG","MYL9","PLN", # Smooth Muscle Cells, SMCs
           "PDGFRB","CSPG4","DES","CD248","ABCC9","ANPEP","VTN","RGS5","S1PR3", # 周细胞（Pericytes）
           "VIM","PDGFRA","COL5A1","LUM","COL8A2","DCN","COL12A1","COL3A1","MMP2","CYP11A1", "STAR","CYP17A1","EGR1","WNT5A","NES","NR5A1", "DLK1", "HES1"# Fibroblasts
)
features = c("PCNA", "TOP2A", "MCM6", "MKI67")
# 遍历基因列表
for (gene in genes) {
  # 绘制特征图
  p <- FeaturePlot(Tcell, features = gene, label = TRUE, order = TRUE)
  
  # 生成包含基因名的文件名
  filename <- paste0(gene, ".png")
  
  # 保存图片
  ggsave(p, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}


CSF1 
states_de <- differentialGeneTest(CDS[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]
write.csv(states_de, 
          file = "Leydig_cell_State_diff_genes.csv")

genes <- c("PCNA", "TOP2A", "MCM6")
genes <- c("BTG2", "CEBPB", "DNAJA1", "KLF4", "TSC22D1", "ALDH1A1", "APOE", "CARD19", "CLU", "DUSP1", "EGR1", "INHBA", "JUN", "NR4A1", "RHOB", "TIMP3", "VNN1", "ZFP36")
for (gene in genes) {
  # 绘制特征图
  p <- FeaturePlot(Tcell, features = gene, label = TRUE, order = TRUE)
  
  # 生成包含基因名的文件名
  filename <- paste0(gene, "_200genes.png")
  
  # 保存图片
  ggsave(p, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}
p <- FeaturePlot(Tcell, features = "INSL3", label = TRUE, order = TRUE)
ggsave(p, filename = "INSL3.png", width = 7, height = 7, units = "in", dpi = 300)

markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(markers,"markers.csv",row.names=FALSE,col.names=TRUE,sep=",")
heatmap10 <- DoHeatmap(seurat_object, features = top10$gene,size = 0.5,slot = "data")+NoLegend()
ggsave(heatmap10,filename = "heatmap.png",width = 12, height = 20,units = "in",dpi = 300)
#monocle3的安装：
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))


marker_genes <- c(
  'ARX','PDGFRA','WNT5A','TCF21','NR2F2','ADAM23','INHBA','MYLK',
  'MYH11','ACTA2','ALDH1A2',
  'PDGFC','PODXL', 'ITGA8', 'NTRK3','RXFP2','LGI1','DGKB','NDST3',
  'MCAM','STEAP4','CSPG4','FSHR','LHCGR','NR5A1')

DotPlot(sce, features = marker_genes) + 
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend(order = 3)) + 
  scale_color_gradientn(
    values = seq(0, 1, 0.2), 
    colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
  )

for (gene in marker_genes) {
  if (gene %in% rownames(sce)) {
    p <- FeaturePlot(sce, features = c(gene), label = TRUE, order = TRUE)
    print(p) # 在屏幕上显示
    tryCatch({
      ggsave(plot = p, filename = paste0(gene, ".png"), width = 7, height = 7, units = "in", dpi = 300)
      print(paste("已保存", gene, ".png"))
    }, error = function(e){
      warning(paste("无法保存", gene, ".png:", e$message))
    })
  } else {
    warning(paste("Marker 基因", gene, "不在 Seurat 对象的基因列表中，跳过绘图。"))
  }
}


##某个基因的三个趋势
library(Seurat)
library(ggplot2)
library(dplyr)

gene_of_interest <- "AMH"
celltype_of_interest <- "Sertoli cells"

# 过滤细胞
sce_sub <- subset(sce, subset = celltype == celltype_of_interest)

# 获取表达和group信息
expr_df <- FetchData(sce_sub, vars = c(gene_of_interest, "group"))
expr_df$group <- factor(expr_df$group, levels = c("E90", "E110", "P7"))

#使用数据绘制小提琴图
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
ggplot(data = expr_df,aes(x = group,y = get(gene_of_interest)))+
  geom_violin(aes(fill = group),alpha = 0.8,scale = 'count')+
  scale_fill_manual(values = color)+
  geom_jitter(width=0.2, size=0.5, alpha=0.3) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill="red") +
  labs(title=paste0(gene_of_interest, " expression in ", celltype_of_interest),
       y=paste(gene_of_interest, "expression"),
       x="Group (Timepoint)") +
  theme_classic()

# 统计检验
print(kruskal.test(expr_df[[gene_of_interest]] ~ expr_df$group))

# 两两检验
print(pairwise.wilcox.test(expr_df[[gene_of_interest]], expr_df$group, p.adjust.method = "BH"))

# 线性趋势
lm_res <- lm(expr_df[[gene_of_interest]] ~ as.numeric(expr_df$group))
summary(lm_res)

##箱线图
library(ggplot2)
pairwise_res <- pairwise.wilcox.test(expr_df[[gene_of_interest]], expr_df$group,
                                     p.adjust.method = "BH")
print(pairwise_res)



library(ggpubr)
ggplot(data = expr_df, aes(x = group, y = get(gene_of_interest), fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
#  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = color) +
  labs(title = paste0(gene_of_interest, " expression in ", celltype_of_interest),
       y = paste(gene_of_interest, "expression"),
       x = "Group (Timepoint)") +
  theme_classic() +
  stat_compare_means(comparisons = list(c("E90", "E110"), c("E110", "P7"), c("E90", "P7")),
                     method = "wilcox.test",
                     label = "p.signif")

#####################################################################################
P150 <- Read10X(data.dir = "E:/Bowen-files/out/GSM5328101_testis_P150_P150")
P150 =CreateSeuratObject(counts = P150,project = "P150", min.cells = 3, min.features = 200)
P150[["mt_percent"]] <- PercentageFeatureSet(P150, features = c("ND1","COX1","ND2","COX2","ATP8","ATP6","ND3","COX3","ND5","ND4","ND4L","ND6","CYTB"))
P150 <- NormalizeData(P150)
P150 <- FindVariableFeatures(P150)
P150 <- ScaleData(P150, vars.to.regress = c("mt_percent"))#消除线粒体影响
P150 <- RunPCA(P150, verbose=F)
ElbowPlot(P150, ndims = 50)
P150 <- FindNeighbors(P150, reduction = "pca", dims = 1:20)

P150 <- FindClusters(P150, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
P150 <- RunUMAP(P150, dims = 1:20, reduction = "pca")
save(P150,file = "P150.Rdata")

library(clustree)
pdf("clustree_不同分辨率.pdf",width = 9,height = 7)
clustree(P150)
dev.off()

pdf("P150_umap.pdf",width = 7,height = 7)
DimPlot(P150, reduction = "umap")+
  ggtitle("P150")
dev.off()
Idents(P150)
P150$RNA_snn_res.0.1#resolution 0.1分成16个cluster
Idents(P150) <- "RNA_snn_res.0.1"
markers=c('COL1A1','COL3A1','DCN','PDGFRA',#基质细胞
          'CLU','AMH','SOX9','HSD17B3',#支持细胞
          'CYP17A1','CYP11A1','INSL3','STAR',#间质细胞
          'PTPRC','CST3','C1QA','CD68','TYROBP','IL7R','CD52','CD3D',#免疫细胞
          'PECAM1','EPCAM','VWF','RAMP2','CDH5',#内皮外皮细胞
          'MYH11','ACTA2','NOTCH3','RGS5',#肌样细胞
          'DDX4','UCHL1','KIT','SYCP3','SYCP1','PIWIL1','ACRV1','TNP1','PRM1')#生殖细胞

DotPlot(P150, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))
table(Idents(P150))
leidig_P150 <- subset(P150,subset = seurat_clusters %in% c(11))
leidig_P150 <- subset(P150, idents = 11)
leidig_P150
##接着把leidig_胚胎的数据导入
Idents(sce) <- "celltype"
leydig_we <- subset(sce, subset = celltype %in% c("Leydig cells"))

head(leydig_we@meta.data)
library(Seurat)
library(harmony)

leidig_P150$group <- "P150"

merged <- merge(leydig_we, y = leidig_P150)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))

merged <- RunHarmony(merged, group.by.vars = "group")

merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.1)

DimPlot(merged, reduction = "umap", group.by = "group")
DimPlot(merged, reduction = "umap", label = TRUE)

markers=c('COL1A1','COL3A1','DCN','PDGFRA',#基质细胞
          'CLU','AMH','SOX9','HSD17B3',#支持细胞
          'CYP17A1','CYP11A1','INSL3','STAR',#间质细胞
          'PTPRC','CST3','C1QA','CD68','TYROBP','IL7R','CD52','CD3D',#免疫细胞
          'PECAM1','EPCAM','VWF','RAMP2','CDH5',#内皮外皮细胞
          'MYH11','ACTA2','NOTCH3','RGS5',#肌样细胞
          'DDX4','UCHL1','KIT','SYCP3','SYCP1','PIWIL1','ACRV1','TNP1','PRM1')#生殖细胞

DotPlot(merged, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

table(Idents(merged), merged$group)

save.image(file="P150.Rdata")
for (gene in markers) {
  if (gene %in% rownames(merged)) {
    p <- FeaturePlot(merged, features = c(gene), label = TRUE, order = TRUE)
    print(p) # 在屏幕上显示
    tryCatch({
      ggsave(plot = p, filename = paste0(gene, ".png"), width = 7, height = 7, units = "in", dpi = 300)
      print(paste("已保存", gene, ".png"))
    }, error = function(e){
      warning(paste("无法保存", gene, ".png:", e$message))
    })
  } else {
    warning(paste("Marker 基因", gene, "不在 Seurat 对象的基因列表中，跳过绘图。"))
  }
}
FeaturePlot(merged, features = c("CYP17A1","CYP11A1"), label = TRUE, order = TRUE)
#######################################################################################
P150 <- Read10X(data.dir = "E:/Bowen-files/out/GSM5328101_testis_P150_P150")
P150 =CreateSeuratObject(counts = P150,project = "P150", min.cells = 3, min.features = 200)
P150[["mt_percent"]] <- PercentageFeatureSet(P150, features = c("ND1","COX1","ND2","COX2","ATP8","ATP6","ND3","COX3","ND5","ND4","ND4L","ND6","CYTB"))

P150 <- NormalizeData(P150)
P150 <- FindVariableFeatures(P150)
P150 <- ScaleData(P150, vars.to.regress = c("mt_percent"))#消除线粒体影响
P150 <- RunPCA(P150, verbose=F)
ElbowPlot(P150, ndims = 50)
P150 <- FindNeighbors(P150, reduction = "pca", dims = 1:20)

P150 <- FindClusters(P150, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
P150 <- RunUMAP(P150, dims = 1:20, reduction = "pca")
save(P150,file = "P90.Rdata")

library(clustree)
pdf("clustree_不同分辨率.pdf",width = 9,height = 7)
clustree(P150)
dev.off()

pdf("P90_umap.pdf",width = 7,height = 7)
DimPlot(P150, reduction = "umap")+
  ggtitle("P150")
dev.off()
Idents(P150)
P150$RNA_snn_res.0.1#resolution 0.1分成16个cluster
Idents(P150) <- "RNA_snn_res.0.1"
markers=c('COL1A1','COL3A1','DCN','PDGFRA',#基质细胞
          'CLU','AMH','SOX9','HSD17B3',#支持细胞
          'CYP17A1','CYP11A1','INSL3','STAR',#间质细胞
          'PTPRC','CST3','C1QA','CD68','TYROBP','IL7R','CD52','CD3D',#免疫细胞
          'PECAM1','EPCAM','VWF','RAMP2','CDH5',#内皮外皮细胞
          'MYH11','ACTA2','NOTCH3','RGS5',#肌样细胞
          'DDX4','UCHL1','KIT','SYCP3','SYCP1','PIWIL1','ACRV1','TNP1','PRM1')#生殖细胞

DotPlot(P150, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))
table(Idents(P150))
#leidig_P150 <- subset(P150,subset = seurat_clusters %in% c(11))
leidig_P150 <- subset(P150, idents = 9)
leidig_P150
load("E:/Bowen-files/out/v5_addT90_new/scRNA_harmony_0.2_named.Rdata")
##接着把leidig_胚胎的数据导入
Idents(sce) <- "celltype"
leydig_we <- subset(sce, subset = celltype %in% c("Leydig cells"))

head(leydig_we@meta.data)
library(Seurat)
library(harmony)

leidig_P150$group <- "P90"

merged <- merge(leydig_we, y = leidig_P150_90)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))

merged <- RunHarmony(merged, group.by.vars = "group")
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.1)

DimPlot(merged, reduction = "umap", group.by = "group")
DimPlot(merged, reduction = "umap", label = TRUE)

markers=c('COL1A1','COL3A1','DCN','PDGFRA',#基质细胞
          'CLU','AMH','SOX9','HSD17B3',#支持细胞
          'CYP17A1','CYP11A1','INSL3','STAR',#间质细胞
          'PTPRC','CST3','C1QA','CD68','TYROBP','IL7R','CD52','CD3D',#免疫细胞
          'PECAM1','EPCAM','VWF','RAMP2','CDH5',#内皮外皮细胞
          'MYH11','ACTA2','NOTCH3','RGS5',#肌样细胞
          'DDX4','UCHL1','KIT','SYCP3','SYCP1','PIWIL1','ACRV1','TNP1','PRM1')#生殖细胞

DotPlot(merged, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

table(Idents(merged), merged$group)

save.image(file="P90.Rdata")
for (gene in markers) {
  if (gene %in% rownames(merged)) {
    p <- FeaturePlot(merged, features = c(gene), label = TRUE, order = TRUE)
    print(p) # 在屏幕上显示
    tryCatch({
      ggsave(plot = p, filename = paste0(gene, ".png"), width = 7, height = 7, units = "in", dpi = 300)
      print(paste("已保存", gene, ".png"))
    }, error = function(e){
      warning(paste("无法保存", gene, ".png:", e$message))
    })
  } else {
    warning(paste("Marker 基因", gene, "不在 Seurat 对象的基因列表中，跳过绘图。"))
  }
}
FeaturePlot(merged, features = c("CYP17A1","CYP11A1"), label = TRUE, order = TRUE)
##
library(Seurat)
library(dplyr)

# 1. 确认Idents设为聚类信息
Idents(merged) <- "seurat_clusters"

# 2. 找到cluster 7的细胞ID
cluster7_cells <- WhichCells(merged, idents = 7)

# 3. 提取cluster 7细胞的表达矩阵（这里用RNA表达）
expr_mat <- GetAssayData(merged, assay = "RNA", slot = "data")[c("PRM1", "TNP1"), cluster7_cells]

# 4. 判断表达的细胞（表达量 > 0）
cells_expressing_PRM1_or_TNP1 <- colnames(expr_mat)[colSums(expr_mat > 0) > 0]

# 5. 从整个Seurat对象中去除这些细胞
cells_to_keep <- setdiff(colnames(merged), cells_expressing_PRM1_or_TNP1)

# 6. 构建新的Seurat对象，保留筛选后的细胞
seurat_obj_filtered <- subset(merged, cells = cells_to_keep)
table(Idents(seurat_obj_filtered), seurat_obj_filtered$group)


markers_2_vs_7 <- FindMarkers(seurat_obj_filtered, ident.1 = 2, ident.2 = 7)
write.table(markers_2_vs_7,"markers_2_vs_7_filted.csv",row.names=T,col.names=TRUE,sep=",")

library(dplyr)

# 假设markers已经读取
# 过滤掉基因名包含 "ENSSS" 的行
markers_filtered <- markers_2_vs_7[!grepl("ENSSS", rownames(markers_2_vs_7)), ]
markers <- markers_filtered 
# 计算 -log10(p_val_adj)
markers <- markers %>%
  mutate(Gene = rownames(markers)) %>%
  mutate(log10_pval = -log10(p_val_adj + 1e-300),  # 避免log(0)错误
         log2FC = avg_log2FC)                       # 你数据里的差异倍数列名是这个

# 自定义上调/下调/非显著基因
markers <- markers %>%
  mutate(Significance = case_when(
    p_val_adj < 0.05 & log2FC > 1  ~ "Up",
    p_val_adj < 0.05 & log2FC < 1 ~ "Down",
    TRUE ~ "Not Significant"
  ))

# 画火山图
ggplot(markers, aes(x = log2FC, y = log10_pval, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = expression(log[2]*" Fold Change"),
       y = expression(-log[10]*" Adjusted P-value")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme(legend.title = element_blank())


##Leydig 间质细胞在不同时期高表达的marker################################
load("E:/Bowen-files/out/zhao_out/zhao_out/P_150_90_e90_110_7_leydigcell.Rdata")
setwd("E:/Bowen-files/out/mfuzz/leydigcell")
library(Mfuzz)
# 重新设置因子的levels顺序
merged@meta.data[["group"]] <- factor(
  merged@meta.data[["group"]],
  levels = c("E90", "E110", "P7", "P90", "P150")
)

# 验证排序结果
table(merged@meta.data[["group"]])

# 加载必要的包
library(Biobase)
library(Seurat)

plot_mfuzz <- function(obj,assay='RNA',mode='knn',nclusters=NULL,prefix='out',gene.use=NULL,nfeatures=3000,group.by='group'){
  library(Mfuzz)
  if (is.null(gene.use)) {
    obj <- FindVariableFeatures(obj,nfeatures=nfeatures)
    gene.use <- VariableFeatures(obj)
  }
  ave <- AverageExpression(obj,group.by = group.by,assays=c(assay),features=gene.use)
  ingene <- ave$RNA
  gene_tpm <- data.matrix(ingene)
  eset <- new("ExpressionSet",exprs = gene_tpm)
  gene.r <- filter.NA(eset, thres=0.25)
  #gene.f <- fill.NA(gene.r,mode="mean")
  gene.f <- fill.NA(gene.r,mode=mode)
  #gene.f <- fill.NA(gene.r,mode="wknn")
  tmp <- filter.std(gene.f,min.std=0)
  gene.s <- standardise(tmp)
  
  if (is.null(nclusters)) {
    nclusters=length(colnames(gene_tpm))
  }
  print(paste0('nclusters is:',nclusters))
  c <- nclusters
  m <- mestimate(gene.s)
  cl <- mfuzz(gene.s, c = c, m = m)
  saveRDS(cl,paste0(prefix,'.rds'))
  write.table(cl$cluster,"output.txt",quote=F,row.names=T,col.names=F,sep="\t")
  
  nrow <- ceiling(nclusters/3)
  
  pdf(paste0(prefix,'.pdf'),18,ceiling(nclusters/6)*10)
  mfuzz.plot(gene.s,cl,mfrow=c(nrow,3),new.window= FALSE)
  #mfuzz.plot(gene.s,cl,mfrow=c(2,3),new.window= FALSE,time.labels=as.vector(colnames(gene_tpm)))
  mfuzz.plot2(gene.s,cl,mfrow=c(nrow,3),time.labels=colnames(gene_tpm),x11 = FALSE)
  dev.off()
}


library(Seurat)
obj <- merged

plot_mfuzz(obj,assay='RNA',mode='knn',nclusters=NULL,prefix='out',gene.use=NULL,nfeatures=3000,group.by='group')
#简化版
plot_mfuzz(obj)

sce@meta.data[["group"]] <- factor(
  sce@meta.data[["group"]],
  levels = c("E90", "E110", "P7")
)
pdf("harmony_named_umap_3group.pdf",width = 21,height = 7)
DimPlot(sce, reduction = "umap",split.by = "group")
dev.off()

pdf("harmony_named_umap_3group_no_split.pdf",width = 7,height = 7)
DimPlot(sce, reduction = "umap", group.by = "group")
dev.off()

for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                layer = "counts",
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent"),
                                pt.size = 0.01,
                                ncol = 3)
}
pdf("violin_after_nFeature_nCount_mt_HB.pdf")
violin_before# 输出质控前的小提琴图，有多少样本就会输出多少张图片
dev.off()

scRNA_harmony@meta.data$group <- factor(
  scRNA_harmony@meta.data$group,
  levels = c("E90", "E110", "P7")
)
Idents(scRNA_harmony) <- "group"
pdf("violin_after_nFeature_nCount_mt_group.pdf",width = 16,height = 4)
VlnPlot(scRNA_harmony,
        features = c("nFeature_RNA", "nCount_RNA", "mt_percent"),
        split.by = "group",
        layer = "counts",
        pt.size = 0.01,
        ncol = 4)
dev.off()

scRNA_harmony@meta.data$orig.ident <- factor(
  scRNA_harmony@meta.data$orig.ident,
  levels = c("E90_1", "E90_2", "E90_3", "E110_1", "E110_2", "E110_3", "P7_1", "P7_2")
)

Idents(scRNA_harmony) <- "orig.ident"
pdf("umap_sample_8sample.pdf",width = 7,height = 7)
DimPlot(scRNA_harmony, reduction = "umap", shuffle = T,group.by = "orig.ident", pt.size = 0.5)
dev.off()

#immune_cell#######################################################################
library(CytoTRACE)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
load("E:/Bowen-files/out/v5_addT90_new/scRNA_harmony_0.2_named.Rdata")
Idents(sce) <- "celltype"
sce <- subset(x = sce, idents = "Endothelial cells")

sce$group
# 将Epi_sce$group转换为因子，并指定水平顺序
sce$group <- factor(sce$group, levels = c("E90", "E110", "P7"))
# 验证水平顺序
levels(sce$group)

table(sce@meta.data$celltype)
#亚群常规方法是提取counts来走作标准流程，这样做肯定没问题；
# 创建SeuratObject对象
Epi_sce = CreateSeuratObject(counts = GetAssayData(sce, assay="RNA", layer='counts'),  # 使用提取的细胞构建新的Seurat对象
                             meta.data = sce@meta.data) # 保留meta.data所有信息
#标注化、归一化、高变基因、pca
Epi_sce <- NormalizeData(Epi_sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
#对数据进行归一化、找高变基因、均一化、进行PCA降维
#重点！！！！！！！！！！
#这里要注意是否选择去批次，可选可不选都能解释通（建议都跑一遍）
# 比如我这里是上皮细胞里有癌细胞，本身存在很强的样本异质性，去了批次反而去掉了样本间的差异
# 免疫细胞也可以不去批次存在较强异质性
# 我这里就不运行
#Epi_sce <- RunHarmony(Epi_sce, group.by.var = "orig.ident")
#降维聚类
Epi_sce
ElbowPlot(Epi_sce, ndims=50, reduction="pca")
Epi_sce <- FindNeighbors(Epi_sce, reduction = "pca", dims = 1:20)#reduction="ha
Epi_sce = FindClusters(Epi_sce,resolution = 0.05)
table(Epi_sce@meta.data$seurat_clusters)

Epi_sce <- RunUMAP(Epi_sce, reduction = "pca", dims = 1:20)##reduction="harmony
#EPI_umap图_未注释
levels(Epi_sce$group)
plot1 =DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE)
plot1

plot2 = DimPlot(Epi_sce, reduction = "umap", group.by='orig.ident',raster=FALSE)
plot2
plot3 = DimPlot(Epi_sce, reduction = "umap",split.by = "group",label = T,raster=FALSE)
plot3
plot4 = DimPlot(Epi_sce, reduction = "umap",group.by = "group",shuffle = T,raster=FALSE)
plot4
pdf("sertoli_cell_3group.pdf",width = 21,height = 7)
plot3
dev.off()
pdf("sertoli_cell_group.pdf",width = 7,height = 7)
plot4
dev.off()

pdf("sertoli_cell_cluster.pdf",width = 7,height = 7)
plot1
dev.off()
save.image("sertoli.Rdata")
#######################区分一下淋巴系和髓系
markers=c('CD19',#B细胞
          'CD4',#T细胞
          'CD14','CXCR3',#单核细胞
          'CST3','C1QA','TYROBP','IL7R','CD52','CD3D',#免疫细胞
          'CD68','CD163',#巨噬细胞
           "CSF3R", "MMP9", "S100A8",#粒细胞
          'CD9')#巨核细胞/血小板

DotPlot(Epi_sce, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

for (gene in markers) {
  # 绘制特征图
  p <- FeaturePlot(Epi_sce, features = gene, label = TRUE, order = TRUE)
  # 生成包含基因名的文件名
  filename <- paste0(gene, "_IMMUNE.pdf")
  # 保存图片
  ggsave(p, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}
table(Epi_sce@meta.data[["RNA_snn_res.0.05"]])
table(Epi_sce$group)
prop.table(table(Idents(Epi_sce)))
table(Idents(Epi_sce), Epi_sce$group)
markers <- FindAllMarkers(Epi_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(markers,"markers.csv",row.names=FALSE,col.names=TRUE,sep=",")

############cytotrace###################
#提取髓系和淋巴系细胞：
#sce <- subset(Epi_sce,idents = c(1,6,9))
sce <- Epi_sce
Idents(sce) <- "group"
exp1 <- as.matrix(GetAssayData(sce, assay = "RNA", layer = "counts"))
exp1 <- exp1[apply(exp1 > 0, 1, sum) >= 5,]
# 使用 CytoTRACE 进行分析，设置 ncores = 1 表示使用一个 CPU 核心进行计算
results <- CytoTRACE(exp1, ncores = 1)
# 提取细胞类型注释信息，并将其转换为字符向量
sce$group <- factor(sce$group, levels = c("E90", "E110", "P7"))
Idents(sce) <- "group"
phenot <- sce$group
phenot <- as.character(phenot)
# 将细胞类型注释的名字设置为元数据中的行名
names(phenot) <- rownames(sce@meta.data)
# 提取 UMAP 降维后的细胞嵌入信息
emb <- sce@reductions[["umap"]]@cell.embeddings
# 使用 CytoTRACE 结果、细胞类型注释和 UMAP 嵌入信息绘制 CytoTRACE 图，并将结果保存到指定目录
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './CytoTRACE结果/')
# 绘制 CytoTRACE 基因表达图，显示前 30 个基因，并将结果保存到指定目录
plotCytoGenes(results, numOfGenes = 30, outputDir = './CytoTRACE结果/')
##### monocle2
library(monocle)
library(Seurat)
library(AUCell)
library(patchwork)
library(ggplot2)
library(DOSE)#人疾病
library(clusterProfiler)#富集分析
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(org.Rn.eg.db)
#library(org.At.tair.db)
library(dplyr)
library(GO.db)
library(rjson)
library(stringr)
library("reshape2")
#dim(Tcell)
#Tcell <- subset(Tcell, orig.ident == "E110.1")
#Tcell = seurat_object[,seurat_object@meta.data$seurat_clusters%in% c(0,2,20,5,1,14,9,10,3,4)]
#sertoli cell

# 为了减少内存占用，可以保存提取数据的seurat对象，并且删除原来的seurat对象
#save(Tcell, file = "sertolinamed.Rda")
#rm(Tcell)
gc()
Tcell<- Epi_sce

rm(sce)
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

CDS <- estimateDispersions(CDS)
gc()

#fData()函数用于提取CDS对象中的基因注释表格，得到的结果为数据框；
#pData()函数作用类似，提取CDS对象中的细胞表型表格；
head(pData(CDS))
head(fData(CDS))
head(dispersionTable(CDS))

#保存创建的CDS对象
save(CDS, file = "CDS_sertolicells.Rdata")
rm(list = ls())
gc()
load("CDS_sertolicells.Rdata")
dim(CDS)

################################################################################
## 2. 差异分析寻找高变基因
# detectGenes()函数：同时统计表达当前基因的细胞数量和细胞表达的基因数；
# min_expr参数用于设置检测阈值，比如min_expr = 0.1表示当前基因的表达量超过0.1才会
# 纳入统计；
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

################################################################################
## 3. 时间轨迹及其差异分析
## 3.1 构建细胞轨迹
# 第一步：选择用于构建细胞轨迹的基因集；
# 选择Top200差异基因作为排序基因；
ordering_genes <- rownames(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
#flc_candidates <- c("PDGFRA","CYP11A1", "STAR","CYP17A1","WNT5A","NES","NR5A1", "DLK1")
#sig_genes <- union(flc_candidates,ordering_genes)
#fData(CDS)$use_for_ordering <-fData(CDS)$num_cells_expressed > 0.1 * ncol(CDS)

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
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.pdf", width = 7, height = 7,units = "in",dpi = 300)
#按“State”分组；
State_clusters <-plot_cell_trajectory(CDS, color_by = "State")
ggsave(State_clusters,filename = "leydig_State.pdf", width = 7, height = 7,units = "in",dpi = 300)

#按seurat分群结果分组
clusters <-plot_cell_trajectory(CDS, color_by = "seurat_clusters")
ggsave(clusters,filename = "leydig_clusters.pdf", width = 7, height = 7,units = "in",dpi = 300)
#按细胞类型分组
plot_cell_trajectory(CDS, color_by = "celltype")
#按照ori_ident分类
clusters <-plot_cell_trajectory(CDS, color_by = "group")
ggsave(clusters,filename = "leydig_orig.ident.pdf", width = 7, height = 7,units = "in",dpi = 300)

save(CDS, file = "CDS_leydig_pseudotime.Rda")

##3.2 比较细胞分化轨迹进程中功能基因的表达差异
#主要用到sm.ns()函数根据表达量拟合曲线；用拟时序做基因差异表达分析
diff_test_res <- differentialGeneTest(
  CDS[expressed_genes, ], 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
head(diff_test_res[, c("gene_short_name", "pval", "qval")])
save.image(file = "monocle_leydig.RData")
# 按q值从小到大排序后，查看最显著的前6个基因的拟时间趋势
sig_gene_names1 <- rownames(diff_test_res[order(diff_test_res$qval)[1:4], ])
seurat_cluster <- plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'Pseudotime')
ggsave(seurat_cluster,filename = "Leydig_seurat_cluster.pdf", width = 14, height = 14,units = "in",dpi = 300)
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

# 也可以查看比较关注基因的拟时间趋势
sig_gene_names1 <- c("DDX4", "UCHL1", "DNMT3L","PCNA","PIWIL4","ETV4","ETV5","KDM1B","GFRA1")
plot_genes_in_pseudotime(CDS[sig_gene_names1, ], color_by = 'State')

sig_gene_names2 <- c('CYP17A1','CYP11A1','INSL3','STAR')
Leydigmarker <- plot_genes_in_pseudotime(CDS[sig_gene_names2, ], color_by = 'State')
ggsave(Leydigmarker,filename = "pseudotime_Leydigmarker.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$PDGFRA=log2(exprs(CDS)['PDGFRA',]+1)
pseudotime_PDGFRA <- plot_cell_trajectory(CDS,color_by="PDGFRA", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_PDGFRA,filename = "pseudotime_PDGFRA.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$CYP11A1=log2(exprs(CDS)['CYP11A1',]+1)
pseudotime_CYP11A1 <- plot_cell_trajectory(CDS,color_by="CYP11A1", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_CYP11A1,filename = "pseudotime_CYP11A1.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$STAR=log2(exprs(CDS)['STAR',]+1)
pseudotime_STAR <- plot_cell_trajectory(CDS,color_by="STAR", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_STAR,filename = "pseudotime_STAR.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$INSL3=log2(exprs(CDS)['INSL3',]+1)
pseudotime_INSL3 <- plot_cell_trajectory(CDS,color_by="INSL3", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_INSL3,filename = "pseudotime_INSL3.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$KIT=log2(exprs(CDS)['KIT',]+1)
pseudotime_KIT <- plot_cell_trajectory(CDS,color_by="KIT", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_KIT,filename = "pseudotime_KIT.png",width = 7, height = 7,units = "in",dpi = 300)

pData(CDS)$DNMT3L=log2(exprs(CDS)['DNMT3L',]+1)
pseudotime_DNMT3L <- plot_cell_trajectory(CDS,color_by="DNMT3L", size=1,show_backbone=TRUE)+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+ scale_color_gradient2(low='grey', mid = 'orange', high = 'blue')
ggsave(pseudotime_DNMT3L,filename = "pseudotime_DNMT3L.png",width = 7, height = 7,units = "in",dpi = 300)


#更换根节点：#按照分化轨迹排序细胞；
CDS <- orderCells(CDS,root_state = 2)
#按“Pseudotime”分组；
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")+scale_color_gsea()+theme(panel.border = element_rect(fill = NA,color = "black",size = 1,linetype = "solid"))+scale_color_gradient2(low="grey",mid="orange",high="blue")
Pseudotime_clusters <- plot_cell_trajectory(CDS, color_by = "Pseudotime")
ggsave(Pseudotime_clusters,filename = "leydig_ProPseudotime.pdf", width = 7, height = 7,units = "in",dpi = 300)


##3.4 差异基因的拟时表达模式聚类分析
#提取差异基因；
sig_gene_names2 <- row.names(subset(diff_test_res, qval < 0.001))
#绘制拟时间差异基因表达谱热图；#重新聚类
plot_pseudotime_heatmap(
  CDS[sig_gene_names2, ], 
  num_clusters = 2, #cluster数量
  show_rownames = F,
  cores = 1
)

Time_diff <- differentialGeneTest(CDS[sig_gene_names2,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
write.csv(Time_diff,"Tiff_diff_all.csv",row.names = F)

Time_genes <- top_n(Time_diff, n = 2000, desc(qval)) %>% pull(gene_short_name) %>% as.character()  

#基因太多了画不出来就选择前2000个基因，所有基因Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p<-plot_pseudotime_heatmap(CDS[sig_gene_names2,], num_clusters=3, show_rownames=F, return_heatmap=T)  #字体太大，看不清，可以导出
ggsave("Time_heatmap.png", p, width = 5, height = 8)
ggsave("Time_heatmap.pdf", p, width = 5, height = 10)

clusters<-cutree(p$tree_row,k=3)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
write.csv(clustering,"间质细胞2000个拟时序基因热图分3个cluster.csv")
save.image(file = "间质细胞-monocle2.Rdata") 
################################################################################
## 4.单细胞轨迹分支分析 
# 当细胞分化轨迹出现分支的时候，意味着细胞将面临不同的分化命运“fate”，接下来主要
# 分析分支事件，比如沿着分化轨迹，基因的表达量如何变化？
# 不同分支之间的差异基因有哪些？

# Monocle 提供一个特定的统计检验方法: branched expression analysis modeling（BEAM）.
## 4.1 BEAM检验
# 使用BEAM()函数对基因进行注释；
# BEAM函数的输入对象： 完成拟时间分析的CellDataSet且轨迹中包含1个分支点；
# 返回一个包含每个基因significance scores 的表格,若得分是显著的则表明该基因的表
# 达是与分支相关的（branch-dependent）。

BEAM_res <- BEAM(CDS[expressed_genes, ], branch_point = 1,cores=1,progenitor_method = "duplicate")
saveRDS(BEAM_res, file = "BEAM_res.RDS")

#按照qval升序排列；
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
head(BEAM_res)

##4.2轨迹分支表达分析
#使用pheatmap包绘制分支表达量热图；
sig_gene_names3 <- row.names(subset(BEAM_res, qval < 1e-5))
plot_genes_branched_heatmap(
  CDS[sig_gene_names3, ],
  branch_point = 1,#设置分支节点
  num_clusters = 3,
  use_gene_short_name = FALSE,
  show_rownames = FALSE)

#使用 plot_genes_branched_pseudotime() 函数绘制拟合曲线，不同线表示不同分支
plot_genes_branched_pseudotime(
  CDS[rownames(BEAM_res)[1:4], ], 
  branch_point = 1, 
  color_by = "State",
  ncol = 1
)

#######################区分一下淋巴系和髓系
markers=c("SOX9","AMH","EGR3","TOMM7","ENO1","BEX1","CITED1","CST9L","DEFB119","FATE1")

DotPlot(Epi_sce, features = markers) + coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

for (gene in markers) {
  # 绘制特征图
  p <- FeaturePlot(Epi_sce, features = gene, label = TRUE, order = TRUE)
  # 生成包含基因名的文件名
  filename <- paste0(gene, "_endothelial_CELL.pdf")
  # 保存图片
  ggsave(p, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}

#细胞周期marker基因加载
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
##You can try to merge the layer together with data.filt <- JoinLayers(data.filt)
Epi_sce <- JoinLayers(Epi_sce)
#计算细胞周期分数
Epi_sce <- CellCycleScoring(Epi_sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(seurat_object[[]])

Idents(Epi_sce) <- "RNA_snn_res.0.1"
plot2 =DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE)
plot2
ggsave("plot2_sertoli_Phase.pdf", width = 7, height = 7, units = "in",dpi=300)


markers <- c("CCL21","PROX1","HEY1","IGFBP3","CD36","CA4","ACKR1")

for (gene in markers) {
  # 绘制特征图
  p <- FeaturePlot(Epi_sce, features = gene, label = TRUE, order = TRUE)
  # 生成包含基因名的文件名
  filename <- paste0(gene, "_PTM_CELL.pdf")
  # 保存图片
  ggsave(p, filename = filename, width = 7, height = 7, units = "in", dpi = 300)
}

FeaturePlot(Epi_sce, features = "CD93", label = TRUE, order = TRUE)

library(ClusterGVis)
library(Seurat)

levels(cellchat@idents)



library(dplyr)
library(ggplot2)
library(scales)

df <- as.data.frame(table(sce$group, sce$celltype))
colnames(df) <- c("stage","celltype","count")

df_ratio <- df %>%
  group_by(stage) %>%
  mutate(freq = count/sum(count))

cell_cols <- c(
  "Endothelial cells"="#E64B35",
  "Germ cells"="#FFB000",
  "Immune cells"="#4DBBD5",
  "Leydig cells"="#00A087",
  "PTM cells"="#3C5488",
  "Sertoli cells"="#F39B7F",
  "Stromal cells"="#8491B4"
)

p <- ggplot(df_ratio,
            aes(x=stage,
                y=freq,
                fill=celltype)) +
  
  geom_bar(stat="identity",
           width=0.75,
           color="black",
           size=0.2) +
  
  scale_y_continuous(labels=percent) +
  
  scale_fill_manual(values=cell_cols) +
  
  labs(x=NULL,
       y="Cell proportion") +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size=11,color="black"),
    axis.title = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    legend.position="right"
  )
p +
  geom_text(aes(label=ifelse(freq>0.05,
                             percent(freq,accuracy=0.1),
                             "")),
            position=position_stack(vjust=0.5),
            size=3,
            color="white")


