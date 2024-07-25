# 把两组人群的异常流行率合并为一个表
# 聚类算法：谱聚类算法(34个脑区作为分类指标)
# Xue-Ru Fan 25 July 2024 @BNU
###################################################
# 按照L组的异常率从大到小排列脑区，csv
###################################################
rm(list=ls())
packages <- c("dplyr")
# packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
datDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster34DK")
datDate <- "240610"
newDate <- "240725"


# load data
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_1_Abno5_Rank_", datDate, ".csv")
AbP_L <- read.csv(file.path(datDir, name))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_2_Abno95_Rank_", datDate, ".csv")
AbP_H <- read.csv(file.path(datDir, name))
# merge
AbP <- full_join(AbP_L, AbP_H, by = "label", suffix = c(".L", ".H"))

# 保留小数点后四位
AbP$perc.L <- round(AbP$perc.L, digits = 4)
AbP$perc.H <- round(AbP$perc.H, digits = 4)

# 保存文件
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_Abno_Perc_", newDate, ".csv")
write.csv(AbP, file.path(datDir, name), row.names = F)
