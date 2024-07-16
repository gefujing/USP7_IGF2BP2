
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(grid)
library(forestploter)
library(ggsignif)
library(ggrepel)
library(ggplot2)
library(ggsci)
library(dplyr)

# 读取数据
load(file = "./Rdata/01-E-MTAB-6134.UNICOX.Rdata")
load(file = "./Rdata/06-GSE71729.UNICOX.Rdata")
load(file = "./Rdata/08-GSE183795.UNICOX.Rdata")
load(file = "./Rdata/17-TCGA-PAAD-seq.UNICOX.Rdata")
load(file = "./Rdata/18-CPTAC-PDAC-pro.UNICOX.Rdata")

# 数据合并

kk = c("METTL3", "METTL14", "METTL16", "PCIF1", "DIMT1", "METTL4", "FTO", "ALKBH5", "YTHDC1",
       "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", "ELAVL1", "EIF3A", "HNRNPC", "RBMX", "HNRNPA2B1", 
       "IGF2BP1", "IGF2BP2", "IGF2BP3", "FMR1")

CPTAC = read.csv(file = "./Rtable/01-unicox/18-CPTAC-PDAC-pro.csv")
EMTAB6134 = read.csv(file = "./Rtable/01-unicox/01-E-MTAB-6134.csv")
GSE183795 = read.csv(file = "./Rtable/01-unicox/08-GSE183795.csv")
GSE71729 = read.csv(file = "./Rtable/01-unicox/06-GSE71729.csv")
TCGA = read.csv(file = "./Rtable/01-unicox/17-TCGA-PAAD-seq.csv")

CPTAC = CPTAC[CPTAC$characteristics %in% kk, c(2,3,5,4)]
EMTAB6134 = EMTAB6134[EMTAB6134$characteristics %in% kk, c(2,3,5,4)]
GSE183795 = GSE183795[GSE183795$characteristics %in% kk, c(2,3,5,4)]
GSE71729 = GSE71729[GSE71729$characteristics %in% kk, c(2,3,5,4)]
TCGA = TCGA[TCGA$characteristics %in% kk, c(2,3,5,4)]

CPTAC$dataset = "CPTAC"
EMTAB6134$dataset = "EMTAB6134"
GSE183795$dataset  = "GSE183795"
GSE71729$dataset  = "GSE71729"
TCGA$dataset  = "TCGA"

CPTAC$n = 135
EMTAB6134$n = 288
GSE183795$n = 189
GSE71729$n = 125
TCGA$n = 177

pdata = rbind(CPTAC, EMTAB6134, GSE183795, GSE71729, TCGA)

pdata$category = ifelse(pdata$characteristics %in% c("METTL3", "METTL14", "METTL16", "PCIF1", "DIMT1", "METTL4"), "Writter", 
                        ifelse(pdata$characteristics %in% c("FTO", "ALKBH5"), "Eraser", "Reader"))

pdata$logp = -log10(pdata$P.Value)
pdata2 = pdata %>%
  group_by(characteristics) %>%
  summarize(mean = mean(Hazard.Ratio))

pdata2 = as.data.frame(pdata2)
pdata2 = pdata2[order(pdata2$mean, decreasing = T),]
pdata$characteristics = factor(pdata$characteristics, levels = pdata2$characteristics)
pdata$dataset = factor(pdata$dataset, levels = c("GSE71729", "CPTAC", "TCGA", "GSE183795", "EMTAB6134"))

# 可视化
ggplot(data = pdata, aes(x = characteristics, y = Hazard.Ratio))+
  geom_bar(stat = "summary", fun = mean, position = "dodge", fill = "grey", alpha = 0.8)+
  geom_point(aes(size = dataset, fill = logp), shape = 21, alpha = 0.8)+
  geom_hline(yintercept = 1, color = "red", lty = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_viridis_c()
ggsave(filename = "./Rplot/3-scatter2.pdf", width = 7.7, height = 4.1)

# 可视化2
ggplot(data = pdata, aes(x = characteristics, y = Hazard.Ratio))+
  geom_bar(stat = "summary", fun = mean, position = "dodge", fill = "grey", alpha = 0.8)+
  geom_point(aes(size = logp, fill = dataset), shape = 21, alpha = 0.8)+
  geom_hline(yintercept = 1, color = "red", lty = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_nejm()
ggsave(filename = "./Rplot/3-scatter3.pdf", width = 6.4, height = 3.5)



