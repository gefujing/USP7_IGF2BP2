
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggvenn)
library(ggVennDiagram)

# 读取数据
kd.list = read.csv(file = "./2-KD-Degs/1-IGF2BP2-KD-DEGS.csv")
rip.list = read.csv(file = "./0-RawData/2-RIP/GSE90639_IGF2BPs_RIP_bingding_gene_and_exp_info.csv")
cellchat.list = read.csv(file = "./4-scRNA/1-cellchat.df-genes.csv")

# 数据列表
kdlist = kd.list[kd.list$change == "DOWN", "X"]
riplist = rip.list[rip.list$log2_fold_change > 1, "geneName"]
cellchatlist = c(cellchat.list$ligand, cellchat.list$receptor, cellchat.list$receptor2)
cellchatlist = unique(cellchatlist)


# 取交集
symbol = intersect(kdlist, intersect(riplist, cellchatlist))
write.csv(symbol, file = "./5-merge/1-merge.csv")


# 韦恩图

x <- list(
  IGF2BP2_KO = kdlist, 
  CellChat = cellchatlist, 
  IGF2BP2_RIP = riplist
)

ggvenn(data = x,
       show_elements = F,
       show_percentage = F,
       fill_color = c("#1ee3cf", "#6b48ff", "#0d3f67"),
       stroke_size = 0.5,
       set_name_color = c("#1ee3cf", "#6b48ff", "#0d3f67"),
       set_name_size = 5
) 
ggsave(filename = "./5-merge/1-IGF2BP2.cellchat.venn.plot.pdf", width = 4.3, height = 3.6)
