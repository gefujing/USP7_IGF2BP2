
# 准备环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(harmony)
load(file = "./4-scRNA/1-sce.all.filt.Rdata")
load(file = "./0-RawData/my_color.Rdata")

# 标准化数据
sce.all.filt = NormalizeData(sce.all.filt, normalization.method = "LogNormalize", scale.factor = 1e4) 

# 高变异基因的选择
sce.all.filt = FindVariableFeatures(sce.all.filt)

# 归一化数据
sce.all.filt = ScaleData(sce.all.filt)

# PCA分析
sce.all.filt = RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))

# 确定数据集的"主成分个数"
ElbowPlot(sce.all.filt)

# 整合批次效应
sce.all.filt = RunHarmony(sce.all.filt, group.by.vars = "orig.ident")
sce.all.filt = Runtsne(sce.all.filt, dims = 1:20, reduction = "harmony")
sce.all.filt = RunTSNE(sce.all.filt, dims = 1:20, reduction = "harmony")
DimPlot(object = sce.all.filt, group.by = "orig.ident", reduction = "umap")
DimPlot(object = sce.all.filt, group.by = "orig.ident", reduction = "tsne")

# 细胞聚类
sce.all.filt = FindNeighbors(sce.all.filt, reduction = "harmony", dims = 1:20)
sce.all.filt = FindClusters(sce.all.filt, resolution = 1, algorithm = 1)

# 取子集
table(sce.all.filt@meta.data$seurat_clusters)
DimPlot(object = sce.all.filt, group.by = "seurat_clusters", reduction = "tsne")
sce.all.filt = subset(x = sce.all.filt, seurat_clusters %in% c(0:19))
sce.all.filt@meta.data$seurat_clusters = factor(sce.all.filt@meta.data$seurat_clusters, levels = 0:19)

DimPlot(sce.all.filt, reduction = "tsne", group.by = "orig.ident", label = F) + scale_color_manual(values = col_vector)
ggsave(filename = "./4-scRNA/3-harmony_orig.pdf", width = 4.7, height = 3.2)

DimPlot(sce.all.filt, reduction = "tsne", group.by = "seurat_clusters", label = F) + scale_color_manual(values = col_vector)
ggsave(filename = "./4-scRNA/4-tsne_cluster.pdf", width = 4, height = 3.3)

# 保存数据
sce.all.filt = SetIdent(sce.all.filt, value = "seurat_clusters")
table(sce.all.filt@active.ident)
sce.all.filt.int = sce.all.filt
save(sce.all.filt.int, file = "./4-scRNA/2-sce.all.filt.int.Rdata")
