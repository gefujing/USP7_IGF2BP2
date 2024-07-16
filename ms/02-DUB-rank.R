
# 读取单基因表达矩阵信息
rm(list = ls())
options(stringsAsFactors = F) 
library(tidyverse)
library(ggrepel)

rawdata = read.csv(file = "./RawData/picture.csv")
rownames(rawdata) = rawdata$Gene
dublist = read.csv(file = "./RawData/DubList.csv", header = F)
dublist = dublist$V1

kk1 = intersect(rawdata$Gene, dublist)
write.csv(kk1, file = "./Dub-merge.csv")

kk2 = c("MATR3", "PABPC1", "DCP1A", "IGF2BP2")


# 可视化
pdata = rawdata

kk = c(kk1, kk2)


lab_gene = pdata[kk,]
lab_gene$Discription = c(rep("DUBs",10), rep("IBPs", 4))


p = ggplot(data = pdata, aes(x = Rank, y = LOG2)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 1, lty = 4, col = "black", lwd = 0.8) +
  ylim(c(10, 35))+
  theme_bw()
p


p1 = p + 
  geom_point(data = lab_gene, size = 2, aes(color = Discription)) +
  geom_text_repel(aes(label = Gene), data = lab_gene, 
                  max.overlaps = Inf, color = "black",
                  nudge_x = 700, direction = "y", hjust = 0,
                  fontface = "italic", size = 2)+
  scale_color_manual(values = c("#3b9a9c", "#dd0a35"))+
  coord_polar()

p1
ggsave(p1, filename = "./IG2-MS.pdf", width = 5.1, height = 2.7)



