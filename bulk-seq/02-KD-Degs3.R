

## 热图
rm(list = ls())
options(stringsAsFactors = F)
library(ggheatmap)
library(ggsci)
library(pheatmap)
library(tinyarray)
library(ggplot2)
load(file = "./2-KD-Degs/IGF2BP2-KD-fpkm.degs.Rdata")

table(grouplist)
grouplist = factor(grouplist, levels = c("SH", "NC"))
dgenes = c("IGF2BP2", "PDGFA", "SORT1", "LAMA4", "ITGB1")
heatdata = exp.fpkm[rownames(exp.fpkm) %in% dgenes, ]
rownames(heatdata) = factor(rownames(heatdata), 
                            levels = dgenes)

pdf(file = "./2-KD-Degs/4-shIGF2BP2_CCC_heatmap2.pdf", width = 3.4, height = 2.25)
draw_heatmap(n = heatdata,
                  group_list = grouplist,
                  cluster_cols = F,
                  n_cutoff = 2, 
                  legend = T, 
                  split_column = T,
                  annotation_legend = T,
                  show_rownames = T)
dev.off()

ggsave(hp, filename = "./2-KD-Degs/4-shIGF2BP2_CCC_heatmap2.pdf", width = 3.4, height = 2.25)
