
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(hdf5r)
library(SeuratDisk)
load(file = "./0-RawData/mycol.Rdata")

# 数据读取1
p1.sce = Read10X(data.dir = "./0-RawData/P1/filtered_feature_bc_matrix/")
p1.image = Read10X_Image(image.dir = file.path("./0-RawData/P1/", "spatial"), filter.matrix = TRUE)
p1.st = CreateSeuratObject(counts = p1.sce, assay = "Spatial")

p1.image = p1.image[Cells(x = p1.st)]
DefaultAssay(p1.st = p1.image) = "Spatial"
p1.st[["slice1"]] = p1.image
p1.st$orig.ident = "P1"

# 数据归一化
p1.st = SCTransform(p1.st, assay = "Spatial", vst.flavor = "v1", verbose = FALSE)
p1.st = NormalizeData(p1.st, verbose = FALSE, assay = "Spatial")
p1.st = RunPCA(p1.st, assay = "SCT", verbose = FALSE)
p1.st = RunUMAP(p1.st, reduction = "pca", dims = 1:30)

p1.st = FindNeighbors(p1.st, reduction = "pca", dims = 1:30)
p1.st = FindClusters(p1.st, verbose = FALSE)

# 可视化
DimPlot(p1.st, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + scale_color_manual(values = col_vec)

# 细胞注释
## 构建markerlip1.st
markerlist = read.csv(file = "./2-celltype/nature.markers.csv")
markerlist$dp = markerlist$pct.1 - markerlist$pct.2
markerlist = markerlist[markerlist$dp > 0.5,]
table(markerlist$cluster)

markerlist2 = list(tumor = markerlist$gene[markerlist$cluster %in% c("Ductal", "Tumor", "Acinar")],
                   stroma = markerlist$gene[markerlist$cluster %in% c("CAF","B_cells", "DCs", "MonoMacro", "T_NK")])

for (i in names(markerlist2)) {
  intersect_sign = intersect(markerlist2[[i]], rownames(p1.st))
  cell_types = make.names(paste(i, "GeneMean", sep = "_"))
  cell_types_Seurat = make.names(paste(i, "Seurat", sep = "_"))
  p1.st = AddMetaData(p1.st, apply(as.matrix(p1.st@assays[["SCT"]]@data[intersect_sign,]), 2, mean), col.name = cell_types) 
  p1.st = AddModuleScore(object = p1.st, features = list(intersect_sign), name = cell_types_Seurat)
}

head(p1.st)


# 检查总体基因情况
# 点图
genes_to_check = colnames(p1.st@meta.data)[c(9,11)]
genes_to_check
DotPlot(p1.st, features = unique(genes_to_check), group.by = "seurat_clusters", assay='SCT') + coord_flip() + theme_bw() + scale_color_gradient(low = "grey", high = "red")

# 小提琴图
VlnPlot(p1.st, features = "tumor_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)
VlnPlot(p1.st, features = "stroma_Seurat1", group.by = "seurat_clusters", pt.size = 0) + scale_fill_manual(values = col_vec)

# 特征图
SpatialPlot(p1.st, features = "tumor_Seurat1")
SpatialPlot(p1.st, features = "stroma_Seurat1")
SpatialDimPlot(p1.st, group.by = "seurat_clusters", label = T)

# 重命名
# 细胞注释
celltype = data.frame(ClusterID = 0:12,
                      celltype = 0:12) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,1,2,7,8,9,12), 2] = 'Tumor'
celltype[celltype$ClusterID %in% c(3,4,5,6,10,11), 2] = 'stroma'

## 写入细胞亚群
table(celltype$celltype)
p1.st@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  p1.st@meta.data[which(p1.st@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(p1.st@meta.data$celltype)
p1.st@meta.data$celltype = factor(p1.st@meta.data$celltype, levels = c('Tumor', 'stroma'))

SpatialDimPlot(p1.st, label = TRUE, group.by = "celltype", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)

# 保存为HD5
stmatrix = p1.st@assays$SCT@counts
rownames(stmatrix) = rownames(p1.st)
colnames(stmatrix) = colnames(p1.st)

p1.st@meta.data$celltype2 = ifelse((p1.st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] > mean(stmatrix["Igf2bp2",])), "Iptumor",
                                ifelse((p1.st@meta.data$celltype == "Tumor") & (stmatrix["Igf2bp2",] < mean(stmatrix["Igf2bp2",])), "Intumor",
                                       ifelse((p1.st@meta.data$celltype == "stroma") & (stmatrix["Acta2",] > mean(stmatrix["Acta2",])), "Apstroma", "Anstroma")))

SpatialDimPlot(p1.st, label = TRUE, group.by = "celltype2", images = "slice1", label.size = 3) + scale_fill_manual(values = col_vec)

p1.st@meta.data$celltype3 = ifelse(p1.st@meta.data$celltype2 %in% c("Iptumor", "Intumor"), "Tumor", p1.st@meta.data$celltype2)
SpatialDimPlot(p1.st, label = TRUE, group.by = "celltype3", images = "slice1", label.size = 3)

save(p1.st, file = "./3-neibour/p1.st.all.Rdata")

# 导出数据
anno = p1.st@meta.data
metadata = read.csv('./3-neibour/metadata.csv', row.names = 1)
metadata$celltype1 = anno$celltype[match(metadata$barcode, rownames(anno))]
metadata$celltype2 = anno$celltype2[match(metadata$barcode, rownames(anno))]
metadata$celltype3 = anno$celltype3[match(metadata$barcode, rownames(anno))]

table(metadata$celltype1)
table(metadata$celltype2)
table(metadata$celltype3)

write.csv(metadata, './3-neibour/annotated_metadata.csv')

