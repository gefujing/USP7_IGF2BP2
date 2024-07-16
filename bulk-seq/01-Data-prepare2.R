
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(dplyr)
library(stringr)
library(data.table)

# 导入CRA001160数据
rawcount = fread(input = "./0-RawData/3-singlecell/CRA001160/count-matrix.txt", data.table = F)
rownames(rawcount) = rawcount[ ,1]
rawcount = rawcount[,-1]
sce.all.CRA001160 = CreateSeuratObject(counts = rawcount)

## 检查表达矩阵
as.data.frame(sce.all.CRA001160@assays$RNA@layers$counts[1:10, 1:2])
head(sce.all.CRA001160@meta.data, 20)
tail(sce.all.CRA001160@meta.data, 20)
table(sce.all.CRA001160@meta.data$orig.ident)
sce.all.CRA001160 = subset(x = sce.all.CRA001160, orig.ident %in% c(paste0("T", 1:24)))
sce.all.CRA001160@meta.data$orig.ident = factor(sce.all.CRA001160@meta.data$orig.ident, levels = c(paste0("T", 1:24)))
sce.all.CRA001160@meta.data$dataset = "CRA001160"
sce.all.CRA001160

## 删除无用数据
rm(list = c("rawcount"))

# # 导入GSE205013数据
# assays = dir("./0-RawData/3-singlecell/GSE205013/")
# dir = paste0("./0-RawData/3-singlecell/GSE205013/", assays)
# 
# ## 按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
# samples_name = assays
# 
# ## 批量读取
# scRNAlist = list()
# for(i in 1:length(dir)){
#   counts = Read10X(data.dir = dir[i])
#   scRNAlist[[i]] = CreateSeuratObject(counts, project = samples_name[i], min.cells = 10, min.features = 500)
#   scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])  #给细胞barcode加个前缀，防止合并后barcode重名
# }
# 
# names(scRNAlist) = samples_name
# 
# ## 合并
# sce.all.GSE205013 = merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
# sce.all.GSE205013@meta.data$dataset = "GSE205013"
# sce.all.GSE205013
# 
# 
# ## 合并两个数据集
# sce.all = merge(x = sce.all.CRA001160, y = sce.all.GSE205013)

sce.all = sce.all.CRA001160
## 保存数据
save(sce.all, file = "./1-Data-prepare/2-PDAC.singlecell.Rdata")

