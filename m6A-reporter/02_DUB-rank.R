
# 读取单基因表达矩阵信息
rm(list = ls())
options(stringsAsFactors = F) 
library(tidyverse)
library(ggrepel)
lname = load(file = "./data/Step01_data_group")
lname  

# fc计算
DubCount$fc = DubCount$LowGFP /   DubCount$HighGFP
filter_count$fc = filter_count$LowGFP /   filter_count$HighGFP

# 可视化预处理
filter_count = filter_count[!str_detect(string = rownames(filter_count), pattern = "mir"),]
filter_count = filter_count[order(filter_count$fc, decreasing = T),]

filter_count$ID = rownames(filter_count)
filter_count$rank = 1:18869


# 可视化
pdata = filter_count

kk = c("USP5", "UBL5", "USP15", "USP43", "USP7", "UFD1L", "UCHL3", "USP46", "USP30", "STAMBP",
       "METTL3", "METTL14", "WTAP", "ALKBH5", "FTO", "YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "IGF2BP1", "IGF2BP2", "IGF2BP3")

lab_gene = pdata[kk,]
lab_gene$Discription = c(rep("DUBs",10), rep("m6A", 12))


p = ggplot(data = pdata, aes(x = rank, y = fc)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 1, lty = 4, col = "black", lwd = 0.8) +
  ylim(c(0,11))+
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
ggsave(p1, filename = "./Rplot/DUB.m6a.pdf", width = 3.6, height = 4.4)



