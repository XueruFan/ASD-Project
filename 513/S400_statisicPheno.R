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
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
resDate <- "240315"
newDate <- "240610"

abide_pheno <- read.csv(file.path(dataDir, paste0("abide_A_all_", resDate, ".csv")))


asd_male <- subset(abide_pheno, SEX == "1" & DX_GROUP == "1")
# 获取用于谱聚类分析的被试编号
name <- paste0("Cluster_", newDate, ".csv")
cluster_id <- read.csv(file.path(abideDir, "Analysis/Cluster/Spect513", name))[, 2]

asd_male <- asd_male[asd_male$Participant %in% cluster_id, ]

# 统计分析
sample <- data.frame()
# number  FIQ VIQ PIQ  

sample[1,1] <- sum(asd_male$DX_GROUP == "1") # asd个数
row.names(sample) <- "ASD组男性"

# 年龄组成
sample[1,2] <- paste0(round(mean(asd_male$AGE_AT_SCAN), 2), "±", round(sd(asd_male$AGE_AT_SCAN), 2))
sample[1,3] <- paste0(round(min(asd_male$AGE_AT_SCAN), 2), "~", round(max(asd_male$AGE_AT_SCAN), 2))

colnames(sample)[2] <- "年龄(M±SD)"
colnames(sample)[3] <- "年龄范围"

# FIQ组成
asd_male$FIQ[asd_male$FIQ <= 0] <- NA
sample[1,4] <- paste0(round(mean(na.omit(asd_male$FIQ)), 2), "±", round(sd(na.omit(asd_male$FIQ)), 
                                                                        2))
colnames(sample)[4] <- "FIQ(M±SD)"

# VIQ组成
asd_male$VIQ[asd_male$VIQ <= 0] <- NA
sample[1,5] <- paste0(round(mean(na.omit(asd_male$VIQ)), 2), "±", round(sd(na.omit(asd_male$VIQ)), 
                                                                        2))
colnames(sample)[5] <- "VIQ(M±SD)"

# PIQ组成
asd_male$PIQ[asd_male$PIQ <= 0] <- NA
sample[1,6] <- paste0(round(mean(na.omit(asd_male$PIQ)), 2), "±", round(sd(na.omit(asd_male$PIQ)),
                                                                        2))
colnames(sample)[6] <- "PIQ(M±SD)"

# save result
name <- paste0("pheno_summary_", newDate, ".xlsx")
write.xlsx(sample, file.path(resDir, name), sheetName = "Sheet1", rowNames = T)
