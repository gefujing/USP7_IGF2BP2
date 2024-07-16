
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(data.table)
library(tinyarray)
library(AnnoProbe)
library(survival)
library(survminer)
library(tidyr)
library(plyr)
library(GEOquery)

# 数据读取
## 读取表达量
exp = fread(input = "./RawData/8-GSE183795/GSE183795.matrix.csv", data.table = F)
exp = exp[!duplicated(exp$ID),]
rownames(exp) = exp$ID
exp = exp[,-1]

phe = fread(input = "./RawData/8-GSE183795/GSE183795.phenotype.csv", data.table = F)
rownames(phe) = phe$sample
phe = phe[,c(1,7,8)]
colnames(phe) = c("Sample", "OS.time", "OS")

# 整理样本名
## 取交集
phe = phe[phe$OS %in% c("0", "1"),]
colnames(exp) = str_remove_all(string = colnames(exp), pattern = ".CEL")
id = intersect(colnames(exp), rownames(phe))

## 取子集
exp = exp[,id]
phe = phe[id,]

# 数据删除
rm(list = c("ids", "id"))

# 数据预处理
## 转置表达量矩阵
exp = as.data.frame(t(exp))
identical(rownames(exp), rownames(phe))

## 联合表型
df = cbind(phe, exp)
df = df[,!duplicated(colnames(df))]
rm(list = c("exp", "phe"))

## 预处理
table(df$OS)
df$OS.time = as.numeric(df$OS.time)
df$OS = as.integer(df$OS)

# 批量COX
BaSurv = Surv(time = df$OS.time, event = df$OS)

UniCox = function(x){
  FML = as.formula(paste0("BaSurv~", x))
  GCox = coxph(FML, data = df)
  GSum = summary(GCox)
  HR = round(GSum$coefficients[,2],2)
  PValue = round(GSum$coefficients[,5],10)
  CI = paste0(round(GSum$conf.int[,3:4],2), collapse = "-")
  Unicox = data.frame("characteristics" = x,
                      "Hazard Ratio" = HR,
                      "CI95" = CI,
                      "P Value" = PValue)
  return(Unicox)
}

## 挑一个基因试一下
UniCox(colnames(df)[9])

# 单因素分析
VarNames = colnames(df)[4:23307]
VarNames = VarNames[!str_detect(string = VarNames, pattern = " ")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "-")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "@")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "/")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "_")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "[.]")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "[|]")]
VarNames = VarNames[!str_detect(string = VarNames, pattern = "[,]")]

UniVar = lapply(VarNames, UniCox)
UniVar = ldply(UniVar, data.frame)
GetFactors = UniVar$characteristics[which(UniVar$Hazard.Ratio > 1 & UniVar$P.Value < 0.05)] %>% as.character()

# 保存数据
write.csv(UniVar, file = "./Rtable/01-unicox/08-GSE183795.csv")
write.csv(GetFactors, file = "./Rtable/02-getfactors/08-GSE183795-get.csv")

GSE183795.GetFactors = GetFactors
save(GSE183795.GetFactors, file = "./Rdata/08-GSE183795.UNICOX.Rdata")




