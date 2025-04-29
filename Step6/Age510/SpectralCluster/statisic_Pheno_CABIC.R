# Analyze CABIC (individuals used for clustering analysis) composition
# Male, ASD, 5~9.9 years old
# Xue-Ru Fan 24 Oct 2023 @BNU
###################################################
rm(list=ls())
packages <- c("openxlsx", "readxl")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result", "pred", "510")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_510_", resDate, ".csv")
cluster <- read.csv(file.path(resuDir, name))

cabic_all <- merge(cluster, pheno, by = "participant")

asd_male <- subset(cabic_all, SEX == "M" & GROUP == "ASD")
asd_male <- subset(cabic_all, AGE < 13)

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
