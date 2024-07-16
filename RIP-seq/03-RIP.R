
# 数据读取
rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

RIP = read.csv("./0-RawData/2-RIP/GSE90639_IGF2BPs_RIP_bingding_gene_and_exp_info.csv")
colnames(RIP) = c("ID", "LogFC")
RIP = RIP[order(RIP$LogFC, decreasing = T) ,]
RIP$rank = 1:5486

# 聚类
s2e = bitr(RIP$ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
DEG = inner_join(RIP, s2e, by = c("ID" = "SYMBOL"))
deg = DEG
cg = DEG$ENTREZID
ego_BP = enrichGO(gene = cg, OrgDb= org.Hs.eg.db, ont = "BP", readable = TRUE)
gosort = ego_BP@result
sortid = gosort[c("GO:0030199", "GO:0010715", "GO:0048144", "GO:0006260"),]
pdata = deg

# 标记
lab_gene1 = data.frame(ID = t(str_split(string = sortid[1,8], pattern = "/", simplify = T)),
                       Discription = sortid[1,2])
lab_gene2 = data.frame(ID = t(str_split(string = sortid[2,8], pattern = "/", simplify = T)),
                       Discription = sortid[2,2])
lab_gene3 = data.frame(ID = t(str_split(string = sortid[3,8], pattern = "/", simplify = T)),
                       Discription = sortid[3,2])
lab_gene4 = data.frame(ID = t(str_split(string = sortid[4,8], pattern = "/", simplify = T)),
                       Discription = sortid[4,2])  

lab_gene = rbind(lab_gene1, lab_gene2, lab_gene3, lab_gene4)
lab_gene = inner_join(lab_gene, pdata)
lab_gene = lab_gene[!duplicated(lab_gene$ID),]
kk = c("PDGFA", "CST3", "DDR1", "FSCN1", "TGFB1", "PDPN", "COL11A2", "MMP11", "SERPINH1",
       "CDT1", "MCRS1", "RECQL4", "MYC", "PTMS", "CDK9", "MCM3", "RPA2", "E2F7", "WDR18")

lab_gene = lab_gene[lab_gene$ID %in% kk,]
lab_gene$Discription = c(rep("Signal release and ECM regulation",10), rep("DNA replication", 9))


p = ggplot(data = pdata, aes(x = rank, y = LogFC)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p


p1 = p + 
  geom_point(data = lab_gene, size = 2, aes(color = Discription)) +
  geom_text_repel(aes(label = ID), data = lab_gene, 
                  max.overlaps = Inf, color = "black",
                  nudge_x = 2000, direction="y", hjust = 0,
                  fontface = "italic", size = 2)+
  scale_color_manual(values = c("#3b9a9c", "#dd0a35"))

p1
ggsave(p1, filename = "./3-RIP/1-RIP-genes.pdf", width = 4.5, height = 3.4)



