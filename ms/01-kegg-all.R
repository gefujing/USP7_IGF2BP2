
# 设置环境
rm(list = ls())
library(tidyverse)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

# ID转换
gene = read.csv(file = "./Rawdata/picture.csv", header = T)
colnames(gene)[1] = "Gene"
enid = bitr(geneID = gene$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gene = inner_join(gene, enid, by = c("Gene" = "SYMBOL"))

# GO富集

ego <- enrichGO(gene = gene$ENTREZID,
                OrgDb= org.Hs.eg.db,
                ont = "ALL",
                readable = TRUE)

kk = ego@result

write.csv(kk, file = "./GO.result.csv")

dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

# kegg富集
kk = enrichKEGG(gene = gene$ENTREZID, pvalueCutoff = 1, qvalueCutoff = 1)
kk.result = kk@result
write.csv(kk.result, file = "./Table/1-KEGG-enrichment-all1.csv")

# 可视化2
colnames(kk.result)
kk.result = kk.result[c(1:2, 4:9, 11, 12, 15, 17, 18, 20:22, 30:32, 35),]
kk.result$order = factor(rev(as.integer(1:20)),
                         labels = rev(kk.result$Description))
fz = str_split(kk.result$GeneRatio, pattern = "/", simplify = T)[,1]
fm = str_split(kk.result$GeneRatio, pattern = "/", simplify = T)[,2]
kk.result$GeneRatio = as.numeric(fz)/as.numeric(fm)

fz = str_split(kk.result$BgRatio, pattern = "/", simplify = T)[,1]
fm = str_split(kk.result$BgRatio, pattern = "/", simplify = T)[,2]
kk.result$BgRatio = as.numeric(fz)/as.numeric(fm)

kk.result$EnrichFactor = kk.result$GeneRatio/kk.result$BgRatio

p1 = ggplot(kk.result, aes(x = Count, y = order, fill = pvalue)) + 
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low = "red", high = "yellow")+
  theme_bw()+
  theme(text = element_text(colour = "black"))
p1
ggsave(p1, filename = "./Rplot/1-KEGG-all2.pdf", width = 5.2, height = 4.3)

