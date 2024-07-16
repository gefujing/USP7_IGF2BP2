
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(tinyarray)
library(AnnoProbe)

# 读取数据
exp.fpkm = read.csv(file = "./0-RawData/1-KD/genes_expression.fpkm.csv")
exp.fpkm = na.omit(exp.fpkm)

# ID转换
exp.fpkm = exp.fpkm[!duplicated(exp.fpkm$gene_id), ]
rownames(exp.fpkm) = exp.fpkm$gene_id
exp.fpkm = exp.fpkm[,-1]

exp.fpkm = trans_exp(exp.fpkm, mrna_only = F)
exp.fpkm = exp.fpkm[,-1]
exp.fpkm = exp.fpkm[!str_detect(string = rownames(exp.fpkm), pattern = "-"),]
exp.fpkm = exp.fpkm[!str_detect(string = rownames(exp.fpkm), pattern = "[.]"),]

# 数据保存
save(exp.fpkm, file = "./1-Data-prepare/1-IGF2BP2.KD.exp.fpkm.Rdata")

