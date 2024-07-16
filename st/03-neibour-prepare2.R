
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(hdf5r)
library(SeuratDisk)
load(file = "./0-RawData/mycol.Rdata")

# 数据读取1
p2.sce = Read10X(data.dir = "./0-RawData/P2/filtered_feature_bc_matrix/")
p2.image = Read10X_Image(image.dir = file.path("./0-RawData/P2/", "spatial"), filter.matrix = TRUE)
p2.st = CreateSeuratObject(counts = p2.sce, assay = "Spatial")

p2.image = p2.image[Cells(x = p2.st)]
DefaultAssay(p2.st = p2.image) = "Spatial"
p2.st[["slice1"]] = p2.image
p2.st$orig.ident = "P2"

# 数据归一化
p2.st = SCTransform(p2.st, assay = "Spatial", vst.flavor = "v1", verbose = FALSE)
p2.st = NormalizeData(p2.st, verbose = FALSE, assay = "Spatial")
p2.st = RunPCA(p2.st, assay = "SCT", verbose = FALSE)
p2.st = RunUMAP(p2.st, reduction = "pca", dims = 1:30)

p2.st = FindNeighbors(p2.st, reduction = "pca", dims = 1:30)
p2.st = FindClusters(p2.st, verbose = FALSE)

# 可视化
DimPlot(p2.st, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + scale_color_manual(values = col_vec)

# 细胞注释
## 构建markerlip2.st
markerlist = read.csv(file = "./2-celltype/nature.markers.csv")
markerlist$dp = markerlist$pct.1 - markerlist$pct.2
markerlist = markerlist[markerlist$dp > 0.5,]
table(markerlist$cluster)

markerlist2 = list(tumor = markerlist$gene[markerlist$cluster %in% c("Ductal", "Tumor", "Acinar")],
                   stroma = markerlist$gene[markerlist$cluster %in% c("CAF","B_cells", "DCs", "MonoMacro", "T_NK")])

for (i in names(markerlist2)) {
  intersect_sign = intersect(markerlist2[[i]], rownames(p2.st))
  cell_types = make.names(paste(i, "GeneMean", sep = "_"))
  cell_types_Seurat = make.names(paste(i, "Seurat", sep = "_"))
  p2.st = AddMetaData(p2.st, apply(as.matrix(p2.st@assays[["SCT"]]@data[intersect_sign,]), 2, mean), col.name = cell_types) 
  p2.st = AddModuleScore(object = p2.st, features = list(intersect_sign), name = cell_types_Seurat)
}

head(p2.st)


# 检查总体基因情况
# 点图
genes_to_check = colnames(p2.st@meta.data)[c(9,11)]
genes_to_check
DotPlot(p2.st, features = unique(genes_to_check), group.by = "seurat_clusters", assay='SCT') + coord_flip() + theme_bw() + scale_color_gradient(low = "grey", high = "red")

# 小提琴图
VlnPlot(p2.st, features = "tumor_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)
VlnPlot(p2.st, features = "stroma_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)

# 特征图
SpatialPlot(p2.st, features = "tumor_Seurat1")
SpatialPlot(p2.st, features = "stroma_Seurat1")
SpatialDimPlot(p2.st, group.by = "seurat_clusters", label = T)

# 重命名
# 细胞注释
celltype = data.frame(ClusterID = 0:13,
                      celltype = 0:13) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(1,2,3,4,5,6,7,11,13), 2] = 'Tumor'
celltype[celltype$ClusterID %in% c(0,8,9,10,12), 2] = 'stroma'

## 写入细胞亚群
table(celltype$celltype)
p2.st@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  p2.st@meta.data[which(p2.st@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(p2.st@meta.data$celltype)
p2.st@meta.data$celltype = factor(p2.st@meta.data$celltype, levels = c('Tumor', 'stroma'))

SpatialDimPlot(p2.st, label = TRUE, group.by = "celltype", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)

# 保存为HD5
stmatrix = p2.st@assays$SCT@counts
rownames(stmatrix) = rownames(p2.st)
colnames(stmatrix) = colnames(p2.st)

p2.st@meta.data$celltype2 = ifelse((p2.st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] > mean(stmatrix["Igf2bp2",])), "Iptumor",
                                   ifelse((p2.st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] < mean(stmatrix["Igf2bp2",])), "Intumor",
                                          ifelse((p2.st@meta.data$celltype == "stroma") & (stmatrix["Acta2",] > mean(stmatrix["Acta2",])), "Apstroma", "Anstroma")))

SpatialDimPlot(p2.st, label = TRUE, group.by = "celltype2", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)

p2.st@meta.data$celltype3 = ifelse(p2.st@meta.data$celltype2 %in% c("Iptumor", "Intumor"), "Tumor", p2.st@meta.data$celltype2)
SpatialDimPlot(p2.st, label = TRUE, group.by = "celltype3", images = "slice1", label.size = 3)

save(p2.st, file = "./3-neibour/p2.st.all.Rdata")

# 导出数据
anno = p2.st@meta.data
metadata = read.csv('./3-neibour/metadata.csv', row.names = 1)
metadata$celltype1 = anno$celltype[match(metadata$barcode, rownames(anno))]
metadata$celltype2 = anno$celltype2[match(metadata$barcode, rownames(anno))]
metadata$celltype3 = anno$celltype3[match(metadata$barcode, rownames(anno))]

table(metadata$celltype1)
table(metadata$celltype2)
table(metadata$celltype3)

write.csv(metadata, './3-neibour/annotated_metadata.csv')

