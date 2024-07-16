
# 数据读取
library(ggplot2)
library(stringr)
library(dplyr)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

rm(list = ls())
options(stringsAsFactors = F)

RIP = read.csv(file = "./0-RawData/2-RIP/GSE90639_IGF2BPs_RIP_bingding_gene_and_exp_info.csv")
colnames(RIP) = c("ID", "LogFC")
RIP = RIP[order(RIP$LogFC, decreasing = T) ,]
RIP$rank = 1:5486

# 数据预处理
## RIP排序
labgenes = c("ITGB1", "SORT1", "LAMA4", "PDGFA")
labdata  = RIP[RIP$ID %in% labgenes, ]

# 绘图
labdata$foldchange = 2^(labdata$LogFC)
labdata$ID = factor(labdata$ID, levels = c("LAMA4", "SORT1", "PDGFA", "ITGB1"))

ggplot(labdata, aes(x = foldchange, y = ID, fill = foldchange)) + 
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "red",high = "yellow" ) +
  labs(title = "IGF2BP2-RIP",
       x = "Fold Change", 
       y = "Gene Name")+
  theme(axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))+
  theme_bw()
ggsave(filename = "./3-RIP/2-RIP.genes2.pdf", width = 3.1, height = 2.1)




