
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(ggheatmap)
library(ggsci)
library(Nebulosa)
library(RColorBrewer)

# 设置数据
load(file = "./4-scRNA/2-sce.all.filt.int.Rdata")
load(file = "./0-RawData/my_color.Rdata")

# 检测分群相关性
table(sce.all.filt.int$seurat_clusters)
av = AverageExpression(sce.all.filt.int, group.by = "seurat_clusters", assays = "RNA")
av = av[[1]]
cg = names(tail(sort(apply(av, 1, sd)), 2000))
kk = cor(as.matrix(av[cg,]), method = 'spearman')

ggheatmap(kk,
          scale = "none",
          cluster_rows = T,
          cluster_cols = T,
          border = "black",
          show_cluster_cols = F,
          show_cluster_rows = F,
          color = c("blue", "white", "red"))

ggsave(filename = "./4-scRNA/5-cluster.cor.heatmap.pdf", width = 4, height = 3.3)

DimPlot(object = sce.all.filt.int, group.by = "seurat_clusters", reduction = "tsne", label = T)  + scale_color_manual(values = col_vector)
ggsave(file = "./4-scRNA/4-tsne_cluster.pdf", width = 4, height = 3.3)

# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
                   'CD19', 'CD79A', 'MS4A1', # B cells
                   'IGHG1', 'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'TPSAB1', 'TPSB2', # mast cells,
                   'RGS5', 'CD73', 'CD105', 'CD44', # perivascular
                   'CD14', 'S100A9', 'S100A8', 'MMP19', # monocyte
                   'FCGR3A', 'FGFBP2', 'CX3CR1', 'KLRB1', 'NCR1', # NK cells
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'FGF7','MME', 'ACTA2', ## human Fibroblasts 
                   'DCN', 'LUM', 'GSN' , ## mouse PDAC Fibroblasts 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  'CDH5', ## Endothelial cells
                   'AMY1', 'AMY2A2', 'PRSS1',  ## Acinar cells
                   'EPCAM' , 'KRT19', 'KRT7', 'PROM1', 'ALDH1A1', 'CD24', # epithelial or tumor
                   'CHGB' ## Endocrine cells
)

genes_to_check
p <- DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()
p
ggsave(filename = "./4-scRNA/6-cluster.markers.pdf", width = 8, height = 10)


# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', # T Cells 
                   'CD79A', 'MS4A1', # B cells
                   'MZB1', # Plasma cells
                   'CD68', 'CD163', 'AIF1',# macrophages
                   'RGS5', # perivascular
                   'S100A9', 'S100A8', # monocyte
                   'FCGR3A', 'KLRB1', # NK cells
                   'CD1E','CD1C', # DC2
                   'FGF7', 'ACTA2', ## human Fibroblasts 
                   'DCN', 'LUM', ## mouse PDAC Fibroblasts 
                   'VWF',  'CDH5', ## Endothelial cells
                   'AMY2A', 'PRSS1',  ## Acinar cells
                   'AMBP', 'EPCAM' , 'KRT19', 'KRT7', # epithelial or tumor
                   'TOP2A', 'MKI67'
)

genes_to_check

p = DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA', cols = c("#41b6e6", "#f6003c")) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = "./4-scRNA/7-cluster.markers4.pdf", width = 8, height = 3.5)

# 细胞注释
celltype = data.frame(ClusterID = 0:19, celltype = 0:19) 

## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,1,8,13,16), 2] = 'PCC'
celltype[celltype$ClusterID %in% c(2,12), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(3,5,17), 2] = 'Fibroblast'
celltype[celltype$ClusterID %in% c(4,6), 2] = 'Stellate'
celltype[celltype$ClusterID %in% c(7,11), 2] = 'Macrophage'
celltype[celltype$ClusterID %in% c(14), 2] = 'DCs'
celltype[celltype$ClusterID %in% c(9), 2] = 'Bcell' 
celltype[celltype$ClusterID %in% c(19), 2] = 'Plasma'
celltype[celltype$ClusterID %in% c(10,15), 2] = 'Tcell'
celltype[celltype$ClusterID %in% c(18), 2] = 'Acinar'

## 写入细胞亚群
table(celltype$celltype)
sce.all.filt.int@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.filt.int@meta.data[which(sce.all.filt.int@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(sce.all.filt.int@meta.data$celltype)
sce.all.filt.int$celltype = factor(sce.all.filt.int$celltype, levels = c('PCC', 'Acinar', 'Fibroblast', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))
save(sce.all.filt.int, file = "./4-scRNA/3-sce.all.filt.int.celltype.Rdata")

# 查看细胞亚群
DimPlot(sce.all.filt.int, reduction = "tsne", group.by = "celltype", label = T, cols = col_vector) 
ggsave(filename = "./4-scRNA/8-tsne.celltype.pdf", width = 5.5, height = 4)

plot_density(object = sce.all.filt.int, features = "IGF2BP2", reduction = "tsne")
ggsave(filename = "./4-scRNA/9-IGF2BP2.cell.pdf", width = 4.2, height = 3.3)


plot_density(object = sce.all.filt.int, features = c("KRT19", "IGF2BP2"), joint = TRUE)[[3]]
ggsave(filename = "./4-scRNA/9-IGF2BP2.cell2.pdf", width = 4.5, height = 3.3)

