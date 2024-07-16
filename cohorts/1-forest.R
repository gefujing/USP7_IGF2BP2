
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(grid)
library(forestploter)

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

rownames(CPTAC) = CPTAC$characteristics
rownames(EMTAB6134) = EMTAB6134$characteristics
rownames(GSE183795) = GSE183795$characteristics
rownames(GSE71729) = GSE71729$characteristics
rownames(TCGA) = TCGA$characteristics

CPTAC = CPTAC[kk, c(2,3,5,4)]
EMTAB6134 = EMTAB6134[kk, c(2,3,5,4)]
GSE183795 = GSE183795[kk, c(2,3,5,4)]
GSE71729 = GSE71729[kk, c(2,3,5,4)]
TCGA = TCGA[kk, c(2,3,5,4)]

rownames(CPTAC) = kk
rownames(EMTAB6134) = kk
rownames(GSE183795) = kk
rownames(GSE71729) = kk
rownames(TCGA) = kk

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
write.csv(pdata, file = "./Rtable/0-mergeForest.csv")

# 数据预处理
pdata = read.csv(file = "./Rtable/0-mergeForest2.csv")
pdata$characteristics = ifelse(is.na(pdata$n), pdata$characteristics, paste0("    ", pdata$characteristics))
pdata$se = (log(pdata$CI95.high) - log(pdata$Hazard.Ratio))/1.96
pdata$` ` = paste(rep(" ", 12), collapse = " ")
pdata$HR = paste0(pdata$Hazard.Ratio, " [", pdata$CI95.low, "-", pdata$CI95.high, "] ")
pdata$P.Value = round(pdata$P.Value, 4)

# 可视化
tm = forest_theme(base_size = 10,           # 设置基础字体大小
                  refline_col = "red4",     # 设置参考线颜色为红色
                  arrow_type = "closed",    # 设置箭头类型为闭合箭头
                  footnote_col = "blue4")   # 设置脚注文字颜色为蓝色


# 绘制森林图
p <- forest(pdata[,c(1, 2, 8, 9, 4)],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
            est = pdata$Hazard.Ratio,          # 效应值，也就是HR列
            lower = pdata$CI95.low,  # 置信区间下限
            upper = pdata$CI95.high,  # 置信区间上限
            ci_column = 3,             # 在第3列（可信区间列）绘制森林图
            ref_line = 1,              # 添加参考线
            arrow_lab = c("Low risk", "High Risk"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
            xlim = c(-1, 5),          # 设置x轴范围
            ticks_at = c(-0.5, 1, 3, 5),  # 在指定位置添加刻度
            theme = tm                )  # 添加自定义主题
ggsave(p, filename = "./Rplot/4-forest.pdf", heigh = 50, width = 10, limitsize = FALSE)

# 数据预处理
pdata = read.csv(file = "./Rtable/0-mergeForest3.csv")
pdata$characteristics = ifelse(is.na(pdata$Hazard.Ratio), pdata$characteristics, paste0("    ", pdata$characteristics))
pdata$se = (log(pdata$CI95.high) - log(pdata$Hazard.Ratio))/1.96
pdata$` ` = paste(rep(" ", 12), collapse = " ")
pdata$HR = paste0(pdata$Hazard.Ratio, " [", pdata$CI95.low, "-", pdata$CI95.high, "] ")
pdata$P.Value = round(pdata$P.Value, 4)

# 可视化
# 定义一个简单的主题，大家可以随意发挥自己的审美！
tm = forest_theme(base_size = 10,           # 设置基础字体大小
                  refline_col = "red4",     # 设置参考线颜色为红色
                  arrow_type = "closed",    # 设置箭头类型为闭合箭头
                  footnote_col = "blue4")   # 设置脚注文字颜色为蓝色


# 绘制森林图
p <- forest(pdata[,c(1, 2, 8, 9, 4)],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
            est = pdata$Hazard.Ratio,          # 效应值，也就是HR列
            lower = pdata$CI95.low,  # 置信区间下限
            upper = pdata$CI95.high,  # 置信区间上限
            ci_column = 3,             # 在第3列（可信区间列）绘制森林图
            ref_line = 1,              # 添加参考线
            arrow_lab = c("Low risk", "High Risk"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
            xlim = c(-1, 4),          # 设置x轴范围
            ticks_at = c(-0.5, 1, 2, 4),  # 在指定位置添加刻度
            theme = tm                )  # 添加自定义主题
ggsave(p, filename = "./Rplot/4-forest2.pdf", heigh = 5, width = 5, limitsize = FALSE)
