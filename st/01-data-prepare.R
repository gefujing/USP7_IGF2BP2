
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(hdf5r)
library(SeuratDisk)

# 数据读取1
p1.sce = Read10X(data.dir = "./0-RawData/P1/filtered_feature_bc_matrix/")
p1.image = Read10X_Image(image.dir = file.path("./0-RawData/P1/", "spatial"), filter.matrix = TRUE)
p1.st = CreateSeuratObject(counts = p1.sce, assay = "Spatial")

p1.image = p1.image[Cells(x = p1.st)]
DefaultAssay(p1.st = p1.image) = "Spatial"
p1.st[["slice1"]] = p1.image
p1.st$orig.ident = "P1"

# 数据读取2
p2.sce = Read10X(data.dir = "./0-RawData/P2/filtered_feature_bc_matrix/")
p2.image = Read10X_Image(image.dir = file.path("./0-RawData/P2/", "spatial"), filter.matrix = TRUE)
p2.st = CreateSeuratObject(counts = p2.sce, assay = "Spatial")

p2.image = p2.image[Cells(x = p2.st)]
DefaultAssay(p2.st = p2.image) = "Spatial"
p2.st[["slice1"]] = p2.image
p2.st$orig.ident = "P2"

# 数据整合
st = merge(x = p1.st, y = p2.st)
save(st, file = "./1-data-prepare/st.all.Rdata")

SpatialFeaturePlot(st, features = c("Krt19", "Igf2bp2", "Pdgfa", "Acta2", "Usp7"), ncol = 5)
ggsave(filename = "./1-data-prepare/st.markers.pdf", width = 11, height = 6)


















