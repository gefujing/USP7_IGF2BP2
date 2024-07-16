
## 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(scRNAstat) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
load(file = "./1-Data-prepare/2-PDAC.singlecell.Rdata")
load(file = "./0-RawData/my_color.Rdata")

## 计算线粒体基因比例
mito_genes = rownames(sce.all)[grep("^MT-", rownames(sce.all))] # 人和鼠的基因名字稍微不一样 
mito_genes #13个线粒体基因
sce.all = PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")

## 计算核糖体基因比例
ribo_genes = rownames(sce.all)[grep("^RP[SL]", rownames(sce.all))]
ribo_genes
sce.all = PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")

## 计算红血细胞基因比例
rownames(sce.all)[grep("^HB[^(P)]", rownames(sce.all))]
sce.all = PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")

# 根据上述指标，过滤低质量细胞/基因
sce.all = JoinLayers(sce.all)
selected_c = WhichCells(sce.all, expression = (nFeature_RNA>300))
selected_f = rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@layers$counts > 0) > 3]
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)

dim(sce.all) 
dim(sce.all.filt) 

table(sce.all@meta.data$orig.ident) 
table(sce.all.filt@meta.data$orig.ident) 
## 可以看到，主要是过滤了基因，其次才是细胞

## 过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 10)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1)
length(selected_hb)
length(selected_ribo)
length(selected_mito)

sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 

## 过滤指标3:过滤特定基因
## 过滤线粒体基因
dim(sce.all.filt)
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt), ignore.case = T), ]
sce.all.filt <- sce.all.filt[!grepl("^RP[SL]", rownames(sce.all.filt), ignore.case = T), ]
sce.all.filt <- sce.all.filt[!grepl("HB[^(P)]", rownames(sce.all.filt), ignore.case = T), ]

dim(sce.all.filt) 
## 当然，还可以过滤更多

## 可视化过滤后的情况
sce.all.filt@meta.data$orig.ident = factor(sce.all.filt@meta.data$orig.ident, levels = c(paste0("T", 1:24), 
                                                                                         "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10",
                                                                                         "P12", "P13", "P14", "P15", "P19", "P20", "P22", "P23", "P26"))

feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered = VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, cols = col_vector) + NoLegend()
p1_filtered
ggsave(p1_filtered, filename="./4-scRNA/1-Filter_Vlnplot1.pdf", width = 10, height = 3.5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered = VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, cols = col_vector) + NoLegend()
p2_filtered
ggsave(p2_filtered, filename="./4-scRNA/2-Filter_Vlnplot2.pdf", width = 15, height = 3.5)

## 保存数据
save(sce.all.filt, file = "./4-scRNA/1-sce.all.filt.Rdata")


