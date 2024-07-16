
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(harmony)
library(mMCPcounter)
library(CARD)
library(patchwork)
library(hdf5r)
library(SeuratDisk)

load(file = "./1-data-prepare/st.all.Rdata")
load(file = "./0-RawData/mycol.Rdata")

## 标记基因
SpatialFeaturePlot(st, features = c("Igf2bp2", "Acta2"), slot = "counts", ncol = 2)
ggsave(filename = "./2-celltype/0-st.markers2.pdf", width = 6.5, height = 6.5)


## 标记基因
SpatialFeaturePlot(st, features = c("Usp7", "Igf2bp2", "Pdgfa"), slot = "counts", ncol = 3)
ggsave(filename = "./2-celltype/0-st.markers.pdf", width = 9.5, height = 6.5)

# 数据归一化
st = SCTransform(st, assay = "Spatial", vst.flavor = "v1", verbose = FALSE)
st = RunPCA(st, assay = "SCT", verbose = FALSE)
st = RunUMAP(st, reduction = "pca", dims = 1:30)

st = FindNeighbors(st, reduction = "pca", dims = 1:30)
st = FindClusters(st, verbose = FALSE)

# 可视化
DimPlot(st, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + scale_color_manual(values = col_vec)
p1 = SpatialDimPlot(st, label = TRUE, images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)
p2 = SpatialDimPlot(st, label = TRUE, images = "slice1.1", label.size = 3) + scale_fill_manual(values = col_vec)
(p1 | p2) + plot_layout(guides = 'collect')
ggsave(filename = "./2-celltype/1-cluster.pdf", width = 7.4, height = 3.7)

# 细胞注释
## 构建markerlist
markerlist = read.csv(file = "./2-celltype/nature.markers.csv")
markerlist$dp = markerlist$pct.1 - markerlist$pct.2
markerlist = markerlist[markerlist$dp > 0.5,]
table(markerlist$cluster)

markerlist2 = list(tumor = markerlist$gene[markerlist$cluster %in% c("Ductal", "Tumor", "Acinar")],
                   stroma = markerlist$gene[markerlist$cluster %in% c("CAF","B_cells", "DCs", "MonoMacro", "T_NK")])

for (i in names(markerlist2)) {
  intersect_sign = intersect(markerlist2[[i]], rownames(st))
  cell_types = make.names(paste(i, "GeneMean", sep = "_"))
  cell_types_Seurat = make.names(paste(i, "Seurat", sep = "_"))
  st = AddMetaData(st, apply(as.matrix(st@assays[["SCT"]]@data[intersect_sign,]), 2, mean), col.name = cell_types) 
  st = AddModuleScore(object = st, features = list(intersect_sign), name = cell_types_Seurat)
}

head(st)


# 检查总体基因情况
# 点图
genes_to_check = colnames(st@meta.data)[c(9,11)]
genes_to_check
p = DotPlot(st, features = unique(genes_to_check), group.by = "seurat_clusters", assay='SCT') + coord_flip() + theme_bw() + scale_color_gradient(low = "grey", high = "red")
p
ggsave(filename = "./2-celltype/1-cluster.markers.pdf", width = 6.8, height = 1.2)

# 小提琴图
p1 = VlnPlot(st, features = "tumor_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)
p2 = VlnPlot(st, features = "stroma_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)

(p1 | p2) + plot_layout(guides = 'collect')
ggsave(filename = "./2-celltype/1-cluster.markers1.pdf", width = 8, height = 3.3)

# 特征图
SpatialPlot(st, features = "tumor_Seurat1")
ggsave(filename = "./2-celltype/2-tumor-feature.pdf", width = 8, height = 4.7)

SpatialPlot(st, features = "stroma_Seurat1")
ggsave(filename = "./2-celltype/2-tumor-feature2.pdf", width = 8, height = 4.7)

# 重命名
# 细胞注释
celltype = data.frame(ClusterID = 0:15,
                      celltype = 0:15) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,1,2,4,6,8,9,11,14,15), 2] = 'Tumor'
celltype[celltype$ClusterID %in% c(3,5,7,10,12,13), 2] = 'Stroma'

## 写入细胞亚群
table(celltype$celltype)
st@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  st@meta.data[which(st@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(st@meta.data$celltype)
st@meta.data$celltype = factor(st@meta.data$celltype, levels = c('Tumor', 'Stroma'))

p1 = SpatialDimPlot(st, label = TRUE, group.by = "celltype", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)
p2 = SpatialDimPlot(st, label = TRUE, group.by = "celltype", images = "slice1.1", label.size = 3) + scale_fill_manual(values = col_vec)

(p1 | p2) + plot_layout(guides = 'collect')
ggsave(filename = "./2-celltype/3-celltype.pdf", width = 7.4, height = 3.7)

# 保存为HD5
stmatrix = st@assays$SCT@counts
rownames(stmatrix) = rownames(st)
colnames(stmatrix) = colnames(st)

st@meta.data$celltype2 = ifelse((st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] > mean(stmatrix["Igf2bp2",])), "Iptumor",
                                ifelse((st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] < mean(stmatrix["Igf2bp2",])), "Intumor",
                                       ifelse((st@meta.data$celltype == "Stroma") & (stmatrix["Acta2",] > mean(stmatrix["Acta2",])), "Apstroma", "Anstroma")))

SpatialDimPlot(st, label = TRUE, group.by = "celltype2", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)
save(st, file = "./2-celltype/st.all.Rdata")

# 保存为HD5
SaveH5Seurat(st, filename = "./2-celltype/st.all.h5Seurat")
Convert("./2-celltype/st.all.h5Seurat", dest = "h5ad")




