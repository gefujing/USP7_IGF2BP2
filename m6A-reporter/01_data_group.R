
# 1.数据读取和分组

rm(list = ls()) # 魔幻操作，一键清空
options(stringsAsFactors = F)

## 查看表达谱
dir.create("./data")
rawcount = read.csv(file = "./GSE138776_data_crisprScreen.csv")
rawcount[1:4,1:3]
colnames(rawcount) = c("GeneID", "HighGFP", "LowGFP")
dim(rawcount)
rawcount = aggregate(.~GeneID, mean, data = rawcount)
rawcount[rawcount$GeneID == "USP7", ]
table(!duplicated(rawcount$GeneID))
rownames(rawcount) = rawcount$GeneID
rawcount = rawcount[ ,-1]
dim(rawcount)

## 分组
group_list = factor(c("High_GFP", "Low_GFP"), levels = c("High_GFP", "Low_GFP"))
keep = rowSums(rawcount>0) >= floor(ncol(rawcount))
table(keep)

## 数据过滤
filter_count = rawcount[keep, ]
filter_count[1:4,1:2]
dim(filter_count)

## DUB相关数据
DubList = read.csv(file = "./DubList.csv")
DubList = as.vector(DubList[, 1])
DubCount = filter_count[rownames(filter_count) %in% DubList, ]

## 保存表达矩阵和分组结果
save(filter_count,
     DubCount,
     group_list,
     file = "./data/Step01_data_group")
