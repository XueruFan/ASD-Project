# this script is used to analysis ABIDE subject composition
# Xue-Ru Fan 24 Oct 2023 @BNU
# 提取用于后续数据分析的被试做人口学统计（性别、年龄和智商）
# 统计做高斯混合聚类分析使用的被试（发育阶段年龄6.0~17.9岁的男性ASD，且脑数据完整）
###################################################
# csv
###################################################
rm(list=ls())
packages <- c("openxlsx")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result", "pred", "515")
plotDir <- file.path(resuDir, "Plot")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_515_", resDate, ".csv")
cluster <- read.csv(file.path(resuDir, name))

cabic_all <- merge(cluster, pheno, by = "participant")

asd_male <- subset(cabic_all, SEX == "M" & GROUP == "ASD")
asd_male <- subset(cabic_all, AGE >= 5 & AGE < 15)

# 统计分析
sample <- data.frame()
# number  FIQ VIQ PIQ  

sample[1,1] <- sum(asd_male$GROUP == "ASD") # asd个数
row.names(sample) <- "ASD组男性"

# 年龄组成
sample[1,2] <- paste0(round(mean(asd_male$AGE), 2), "±", round(sd(asd_male$AGE), 2))
sample[1,3] <- paste0(round(min(asd_male$AGE), 2), "~", round(max(asd_male$AGE), 2))

colnames(sample)[2] <- "年龄(M±SD)"
colnames(sample)[3] <- "年龄范围"

# FIQ组成
asd_male$FIQ[asd_male$FIQ <= 0] <- NA
sample[1,4] <- paste0(round(mean(na.omit(asd_male$FIQ)), 2), "±", round(sd(na.omit(asd_male$FIQ)), 
                                                                        2))
colnames(sample)[4] <- "FIQ(M±SD)"


# save result
name <- paste0("pheno_summary_", resDate, ".xlsx")
write.xlsx(sample, file.path(resuDir, name), sheetName = "Sheet1", rowNames = T)
