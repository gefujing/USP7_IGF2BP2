
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(ggheatmap)
library(ggsci)
library(org.Hs.eg.db)
library(clusterProfiler)
library(CellChat)
library(scRNAtoolVis)

# 设置数据
load(file = "./4-scRNA/3-sce.all.filt.int.celltype.Rdata")
load(file = "./0-RawData/my_color.Rdata")

# 分组
Idents(sce.all.filt.int) = sce.all.filt.int@meta.data$celltype
sce.all.filt.int@meta.data$celltype = as.character(sce.all.filt.int@meta.data$celltype)
episub = sce.all.filt.int[ ,Idents(sce.all.filt.int) %in% c("PCC")]
table(sce.all.filt.int$celltype)

epimatrix = episub@assays$RNA@layers$counts
rownames(epimatrix) = rownames(episub)
colnames(epimatrix) = colnames(episub)

episub$celltype = ifelse(epimatrix["IGF2BP2",] > 0, "IGF2BP2.pos", "IGF2BP2.neg")
table(episub$celltype)

sce.all.filt.int$celltype[match(colnames(episub), colnames(sce.all.filt.int))] = episub$celltype 
table(sce.all.filt.int$celltype)

sce.all.filt.int@meta.data$celltype = factor(sce.all.filt.int@meta.data$celltype, levels = c('IGF2BP2.pos', 'IGF2BP2.neg', 'Acinar', 'Fibroblast', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))
DimPlot(sce.all.filt.int, reduction = "tsne", group.by = 'celltype', label = TRUE, repel = T, pt.size = 1, cols = col_vector)
ggsave(filename = "./4-scRNA/10-igf2.tsne.pdf", width = 5.2, height = 3.9)

save(sce.all.filt.int, file = "./4-scRNA/4-sce.all.filt.int.IGF2BP2.Rdata")

VlnPlot(object = sce.all.filt.int, features = "IGF2BP2", group.by = "celltype", split.by = "celltype", pt.size = 0, cols = col_vector)
ggsave(filename = "./4-scRNA/11-igf2.vln.pdf", width = 5, height = 3.5)


# 热图
genes_to_check = c('IGHG1', 'MZB1', 'SDC1', # Plasma cells
                   'CD79A', 'MS4A1', 'CD79A', # B cells
                   'CD3D', 'CD3E', 'CD2',  # T Cells
                   'CD1C', 'CD1E', #DC
                   'CD68', 'CD163', 'CD14',# macrophages
                   'VWF',  'CDH5', 'PLVAP', ## Endothelial cells
                   'RGS5', 'ADIRF', 'ACTA2', # Stellate cell
                   'LUM', 'DCN', 'SFRP2', ## human Fibroblasts
                   'PRSS1', 'CTRB2', 'CELA3A', ## Acinar cell                  
                   'EPCAM' , 'KRT19', 'KRT7', 'IGF2BP2' # epithelial or tumor
)

genes_to_check = unique(genes_to_check)
Idents(sce.all.filt.int) = sce.all.filt.int@meta.data$celltype
sce.all.filt.int@meta.data$celltype = factor(sce.all.filt.int@meta.data$celltype, levels = c('IGF2BP2.pos', 'IGF2BP2.neg', 'Acinar', 'Fibroblast', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))

p = DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "celltype", assay = 'RNA', cols = c("#41b6e6", "#f6003c")) + 
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p
ggsave(filename = "./4-scRNA/12-cluster.markers4.pdf", width = 4.5, height = 5.4)


# 细胞通讯分析
rm(list = ls()) 
options(stringsAsFactors = F)
library(CellChat)
load(file = "./4-scRNA/4-sce.all.filt.int.IGF2BP2.Rdata")
table(sce.all.filt.int@meta.data$celltype)

## 构建CELLCHAT对象
data.input = sce.all.filt.int@assays$RNA@layers$data
rownames(data.input) = rownames(sce.all.filt.int)
colnames(data.input) = colnames(sce.all.filt.int)

meta.data = sce.all.filt.int@meta.data

meta.data$celltype = factor(meta.data$celltype, levels = c('IGF2BP2.pos', 'IGF2BP2.neg', 'Acinar', 'Fibroblast', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))
cellchat = createCellChat(object = data.input, meta = meta.data, group.by = "celltype")

rm(list = c("sce.all.filt.int", "data.input", "meta.data"))

## 加载CellChatDB数据库
CellChatDB = CellChatDB.human
CellChatDB.use = CellChatDB
cellchat@DB = CellChatDB.use

## 对表达数据进行预处理
cellchat = subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) # do parallel

cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)

## 计算通讯概率，推断细胞通讯网络
cellchat = computeCommunProb(cellchat, population.size = F)
cellchat = filterCommunication(cellchat, min.cells = 10)

## 提取预测的细胞通讯网络为data frame
df.net = subsetCommunication(cellchat)
write.csv(x = df.net, file = "./4-scRNA/1-cellchat.df.csv")
df.pathway = subsetCommunication(cellchat, slot.name = "netP")

levels(cellchat@idents)
# df.net = subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(3,4))
# df.net = subsetCommunication(cellchat, signaling = c("PDGF", "SPP1"))

## 在信号通路水平推断细胞通讯
cellchat = computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)

## 计算加和的cell-cell通讯网络
cellchat = aggregateNet(cellchat)
groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize,
                 weight.scale = T, 
                 label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,
                 weight.scale = T, 
                 label.edge= F,
                 title.name = "Interaction weights/strength")


## 保存数据
save(cellchat, file = "./4-scRNA/5-sce.all.filt.int.cellchat.Rdata")


