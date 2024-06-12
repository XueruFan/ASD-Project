# this script is used to analysis ABIDE subject composition
# Xue-Ru Fan 24 Oct 2023 @BNU
# 提取用于后续数据分析的被试做人口学统计（性别、年龄和智商）

rm(list=ls())
packages <- c("openxlsx")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Statistic")
resDate <- "240315"
newDate <- "240610"

abide_pheno <- read.csv(file.path(dataDir, paste0("abide_A_all_", resDate, ".csv")))

####################### Part 1: Calculate how many #################################################

# 统计分析
sample <- data.frame()
# asd (female/asd)  FIQ VIQ PIQ  
# nc  (female/asd)

sample[1,1] <- sum(abide_pheno$DX_GROUP == "1") # asd个数
sample[2,1] <- sum(abide_pheno$DX_GROUP == "2") # nc个数

# 性别组成
sample[1,2] <- sum(abide_pheno$DX_GROUP == "1" & abide_pheno$SEX == "2") # asd里的女性个数
sample[1,3] <- paste0(round(sample[1,2]/sample[1,1]*100, 2), "%") # asd里的女性占比
sample[1,4] <- paste0(sample[1,1], " (", sample[1,3], ")")

sample[2,2] <- sum(abide_pheno$DX_GROUP == "2" & abide_pheno$SEX == "2") # asd里的女性个数
sample[2,3] <- paste0(round(sample[2,2]/sample[2,1]*100, 2), "%") # nc里的女性占比
sample[2,4] <- paste0(sample[2,1], " (", sample[2,3], ")")

row.names(sample) <- c("ASD组", "对照组")
colnames(sample)[4] <- "人数(女性%)"

asd <- subset(abide_pheno, DX_GROUP == "1")
nc <- subset(abide_pheno, DX_GROUP == "2")

# 年龄组成
sample[1,5] <- paste0(round(mean(asd$AGE_AT_SCAN), 2), "±", round(sd(asd$AGE_AT_SCAN), 2))
sample[1,6] <- paste0(round(min(asd$AGE_AT_SCAN), 2), "~", round(max(asd$AGE_AT_SCAN), 2))

sample[2,5] <- paste0(round(mean(nc$AGE_AT_SCAN), 2), "±", round(sd(nc$AGE_AT_SCAN), 2))
sample[2,6] <- paste0(round(min(nc$AGE_AT_SCAN), 2), "~", round(max(nc$AGE_AT_SCAN), 2))

sample <- sample[, c(-1:-3)]
colnames(sample)[2] <- "年龄(M±SD)"
colnames(sample)[3] <- "年龄范围"

# FIQ组成
asd$FIQ[asd$FIQ <= 0] <- NA
nc$FIQ[nc$FIQ <= 0] <- NA
sample[1,4] <- paste0(round(mean(na.omit(asd$FIQ)), 2), "±", round(sd(na.omit(asd$FIQ)), 2))
sample[2,4] <- paste0(round(mean(na.omit(nc$FIQ)), 2), "±", round(sd(na.omit(nc$FIQ)), 2))
colnames(sample)[4] <- "FIQ(M±SD)"

# VIQ组成
asd$VIQ[asd$VIQ <= 0] <- NA
nc$VIQ[nc$VIQ <= 0] <- NA
sample[1,5] <- paste0(round(mean(na.omit(asd$VIQ)), 2), "±", round(sd(na.omit(asd$VIQ)), 2))
sample[2,5] <- paste0(round(mean(na.omit(nc$VIQ)), 2), "±", round(sd(na.omit(nc$VIQ)), 2))
colnames(sample)[5] <- "VIQ(M±SD)"

# PIQ组成
asd$PIQ[asd$PIQ <= 0] <- NA
nc$PIQ[nc$PIQ <= 0] <- NA
sample[1,6] <- paste0(round(mean(na.omit(asd$PIQ)), 2), "±", round(sd(na.omit(asd$PIQ)), 2))
sample[2,6] <- paste0(round(mean(na.omit(nc$PIQ)), 2), "±", round(sd(na.omit(nc$PIQ)), 2))
colnames(sample)[6] <- "PIQ(M±SD)"

# save result
name <- paste0("abide_A_all_summary_", newDate, ".xlsx")
write.xlsx(sample, file.path(resDir, name), sheetName = "Sheet1", rowNames = T)


####################### Part 2: Calculate dev ######################################################
# 仅统计发育阶段年龄6.0~17.9岁
abide_pheno <- subset(abide_pheno, AGE_AT_SCAN < 18 & AGE_AT_SCAN >= 6)

# 统计分析
sample <- data.frame()
# asd (female/asd)  FIQ VIQ PIQ  
# nc  (female/asd)

sample[1,1] <- sum(abide_pheno$DX_GROUP == "1") # asd个数
sample[2,1] <- sum(abide_pheno$DX_GROUP == "2") # nc个数

# 性别组成
sample[1,2] <- sum(abide_pheno$DX_GROUP == "1" & abide_pheno$SEX == "2") # asd里的女性个数
sample[1,3] <- paste0(round(sample[1,2]/sample[1,1]*100, 2), "%") # asd里的女性占比
sample[1,4] <- paste0(sample[1,1], " (", sample[1,3], ")")

sample[2,2] <- sum(abide_pheno$DX_GROUP == "2" & abide_pheno$SEX == "2") # asd里的女性个数
sample[2,3] <- paste0(round(sample[2,2]/sample[2,1]*100, 2), "%") # nc里的女性占比
sample[2,4] <- paste0(sample[2,1], " (", sample[2,3], ")")

row.names(sample) <- c("ASD组", "对照组")
colnames(sample)[4] <- "人数(女性%)"

asd <- subset(abide_pheno, DX_GROUP == "1")
nc <- subset(abide_pheno, DX_GROUP == "2")

# 年龄组成
sample[1,5] <- paste0(round(mean(asd$AGE_AT_SCAN), 2), "±", round(sd(asd$AGE_AT_SCAN), 2))
sample[1,6] <- paste0(round(min(asd$AGE_AT_SCAN), 2), "~", round(max(asd$AGE_AT_SCAN), 2))

sample[2,5] <- paste0(round(mean(nc$AGE_AT_SCAN), 2), "±", round(sd(nc$AGE_AT_SCAN), 2))
sample[2,6] <- paste0(round(min(nc$AGE_AT_SCAN), 2), "~", round(max(nc$AGE_AT_SCAN), 2))

sample <- sample[, c(-1:-3)]
colnames(sample)[2] <- "年龄(M±SD)"
colnames(sample)[3] <- "年龄范围"

# FIQ组成
asd$FIQ[asd$FIQ <= 0] <- NA
nc$FIQ[nc$FIQ <= 0] <- NA
sample[1,4] <- paste0(round(mean(na.omit(asd$FIQ)), 2), "±", round(sd(na.omit(asd$FIQ)), 2))
sample[2,4] <- paste0(round(mean(na.omit(nc$FIQ)), 2), "±", round(sd(na.omit(nc$FIQ)), 2))
colnames(sample)[4] <- "FIQ(M±SD)"

# VIQ组成
asd$VIQ[asd$VIQ <= 0] <- NA
nc$VIQ[nc$VIQ <= 0] <- NA
sample[1,5] <- paste0(round(mean(na.omit(asd$VIQ)), 2), "±", round(sd(na.omit(asd$VIQ)), 2))
sample[2,5] <- paste0(round(mean(na.omit(nc$VIQ)), 2), "±", round(sd(na.omit(nc$VIQ)), 2))
colnames(sample)[5] <- "VIQ(M±SD)"

# PIQ组成
asd$PIQ[asd$PIQ <= 0] <- NA
nc$PIQ[nc$PIQ <= 0] <- NA
sample[1,6] <- paste0(round(mean(na.omit(asd$PIQ)), 2), "±", round(sd(na.omit(asd$PIQ)), 2))
sample[2,6] <- paste0(round(mean(na.omit(nc$PIQ)), 2), "±", round(sd(na.omit(nc$PIQ)), 2))
colnames(sample)[6] <- "PIQ(M±SD)"

# save result
name <- paste0("abide_A_dev_summary_", newDate, ".xlsx")
write.xlsx(sample, file.path(resDir, name), sheetName = "Sheet1", rowNames = T)


####################### Part 3: Calculate asd male dev #############################################
# 仅统计发育阶段年龄6.0~17.9岁的男性ASD
asd_male <- subset(abide_pheno, SEX == "1" & DX_GROUP == "1")

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
name <- paste0("abide_A_asd_male_dev_summary_", newDate, ".xlsx")
write.xlsx(sample, file.path(resDir, name), sheetName = "Sheet1", rowNames = T)