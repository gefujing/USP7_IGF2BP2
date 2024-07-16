
# 导入数据
rm(list = ls())
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(msigdbr)
library(GSEABase)
library(GseaVis)
library(aPEAR)
library(dplyr)
load(file = "./2-KD-Degs/IGF2BP2-KD-fpkm.degs.Rdata")

# ID转换
s2e = bitr(DEG$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
DEG = inner_join(x = DEG, y = s2e, by = c("symbol" = "SYMBOL"))

# GO富集
cg = DEG[DEG$change == "DOWN", "ENTREZID"]

ego = enrichGO(gene = cg, OrgDb = org.Hs.eg.db, ont = "ALL",pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
ego.kk = ego@result
write.csv(ego.kk, file = "./2-KD-Degs/2-IGF2BP2-KD-DOWN-GO.csv")
rownames(ego.kk) = 1:7746

GO = ego[c(6,8,12,16,26,51,174,561),c(2,3,1,6,9)]
GO$geneID = str_replace_all(GO$geneID,"/",",") ### 修改geneID这一列的分隔符号
names(GO) = c("ID","Term","category","adj_pval","genes")

genedata = data.frame(ID = DEG$symbol, logFC = DEG$logFC)
circ = circle_dat(GO, genedata)

# 弦图
gene = circ %>% group_by(term) %>% top_n(5, logFC)
chord = chord_dat(circ, genes = gene$genes)
GOChord(chord, gene.order = "logFC", ribbon.col = brewer.pal(8, "Set3"))

ggsave(filename = "./2-KD-Degs/3-DEGs-go-circ2.pdf", width = 3.7, height = 4.7)
