
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(dplyr)
library(CellChat)
library(Nebulosa)

# 设置数据
load(file = "./4-scRNA/5-sce.all.filt.int.cellchat.Rdata")
load(file = "./0-RawData/my_color.Rdata")


# CELLCHAT可视化

## 使用层次图（Hierarchical plot），圆圈图（Circle plot）或和弦图（Chord diagram）可视化每个信号通路
pathways.show = c("PDGF") 
vertex.receiver = c(1,2,3,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", color.use = col_vector)

par(mfrow = c(1,3))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", label.edge= T)

pairLR.PDGF = extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show = pairLR.PDGF[1,]
LR.show
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = "PDGFA_PDGFRA", layout = "circle")

## 热图
par(mfrow = c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

## 计算每个配体-受体对L-R pairs对整个信号通路的贡献
netAnalysis_contribution(cellchat, signaling = pathways.show)

## 观察多种配体受体或信号通路介导的细胞-细胞通讯
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(4:5), remove.isolate = T)


# 可视化
load(file = "./4-scRNA/3-sce.all.filt.int.celltype.Rdata")
plot_density(object = sce.all.filt.int, features = c("PDGFRA", "PDGFRB"), reduction = "tsne")
ggsave(filename = "./4-scRNA/19-PDGFR.cell.pdf", width = 8.4, height = 3.3)

plot_density(object = sce.all.filt.int, features = c("PDGFA", "IGF2BP2"), joint = TRUE)[[3]]
ggsave(filename = "./4-scRNA/20-PDGFA.cell2.pdf", width = 4.5, height = 3.3)

plot_density(object = sce.all.filt.int, features = c("IGF2BP2", "LAMA4"), joint = TRUE)[[3]]
ggsave(filename = "./4-scRNA/21-LAMA4.cell.pdf", width = 4.5, height = 3.3)
