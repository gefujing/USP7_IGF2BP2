
# 导入数据
rm(list = ls())
library(DESeq2)
library(limma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggrepel)
library(tidyverse)
library(pheatmap)


load("./1-Data-prepare/1-IGF2BP2.KD.exp.fpkm.Rdata")
grouplist = c("NC", "NC", "NC", "SH", "SH", "SH")
table(grouplist)
grouplist = factor(grouplist, levels = c("NC", "SH"))

# PCA
dat = as.data.frame(t(exp.fpkm))
dat.pca = PCA(dat, graph = FALSE)
pca_plot = fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = grouplist, # color by groups
                        palette = c("#293462", "#fe5f55"),
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Groups")

pca_plot
ggsave(pca_plot, filename = "./2-KD-Degs/1-PCA.pdf", width = 4.5, height = 3.5)


# limma分析差异基因
design = model.matrix(~0 + grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp.fpkm)

dge = DGEList(counts = exp.fpkm)
dge = calcNormFactors(dge)

v = voom(dge, design, normalize = "quantile")
fit = lmFit(v, design)

constrasts = paste(rev(levels(grouplist)), collapse = "-")
cont.matrix = makeContrasts(contrasts = constrasts, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, coef = constrasts, n = Inf)
DEG = na.omit(DEG)

DEG$symbol = rownames(DEG)

# 差异基因统计分析
pvalue_t = 0.05
logFC_t = 0

k1 = (DEG$P.Value < pvalue_t)&(DEG$logFC < -logFC_t);table(k1)
k2 = (DEG$P.Value < pvalue_t)&(DEG$logFC > logFC_t);table(k2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

table(DEG$change)

## 保存数据
write.csv(DEG, file = "./2-KD-Degs/1-IGF2BP2-KD-DEGS.csv")
save(exp.fpkm, grouplist, DEG, file = "./2-KD-Degs/IGF2BP2-KD-fpkm.degs.Rdata")

# 可视化
# 导入数据
rm(list = ls())
load(file = "./2-KD-Degs/IGF2BP2-KD-fpkm.degs.Rdata")

# 火山图
pdata = DEG
pdata$symbol = rownames(pdata)
for_label = pdata[c("IGF2BP2", "PDGFA", "APLP2", "GPT2", "SLC1A5", "TAB3"),]

ggplot(data = pdata) + 
  geom_point(aes(x = logFC, y = -log10(P.Value), color = logFC, size = -log10(P.Value))) + 
  geom_text_repel(data =  for_label, aes(x = logFC, y = -log10(P.Value), label = symbol),
                  nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, direction = "y", hjust = "left", max.overlaps = 200 ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(log2(1)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  theme_bw() + 
  ggtitle(label = "Volcano Plot")
ggsave(filename = "./2-KD-Degs/2-vocanal.pdf", width = 4, height = 3)

# 热图
pdata = pdata[order(pdata$logFC),]
cg = pdata$symbol[pdata$change != "NOT"]
cg = cg[!str_detect(string = cg, pattern = "-")]
cg = cg[!str_detect(string = cg, pattern = "[.]")]
cg = cg[c(1:10, 987:996)]
cg = c(cg, "IL6", "CXCL1", "CXCL5")
cg = unique(cg)

exp = exp[cg,]

colanno = data.frame(row.names = colnames(exp),
                     group = factor(grouplist, levels = c("Parent", "Crizoleres")))
anncolor = list(group = c(Parent = "#293462", Crizoleres = "#a64942")) 

pheatmap(mat = exp,
         scale = "row",
         show_colnames = F,
         cluster_cols = F,
         annotation_col = colanno,
         annotation_colors = anncolor,
         gaps_col = c(3),
         breaks = seq(-1.5,1.5,length.out = 100),
         cellwidth = 20, cellheight = 20)



















