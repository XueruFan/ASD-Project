rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
cabicDir <- 'E:/PhDproject/CABIC/result'
resuDir <- file.path(abideDir, "Analysis/Statistic/Spect513")


cabic_h_p <- read.csv(file.path(cabicDir, "cabic_P_value_Cluster_2_vs_3.csv"))
rownames(cabic_h_p) <- cabic_h_p[,1]
cabic_h_p <- cabic_h_p[,-1]

cabic_l_p <- read.csv(file.path(cabicDir, "cabic_P_value_Cluster_1_vs_3.csv"))
rownames(cabic_l_p) <- cabic_l_p[,1]
cabic_l_p <- cabic_l_p[,-1]

abide_h_p <- read.csv(file.path(resuDir, "abide_P_value_Cluster_2_vs_3.csv"))
rownames(abide_h_p) <- abide_h_p[,1]
abide_h_p <- abide_h_p[,-1]

abide_l_p <- read.csv(file.path(resuDir, "abide_P_value_Cluster_1_vs_3.csv"))
rownames(abide_l_p) <- abide_l_p[,1]
abide_l_p <- abide_l_p[,-1]

############ H vs CN
# 1. 找到两个数据框中都不为 NA 的位置
common_non_na <- !is.na(abide_h_p) & !is.na(cabic_h_p)

# 2. 提取这些共同的非 NA 元素位置
indices <- which(common_non_na, arr.ind = TRUE)

# 3. 根据索引获取行名和列名
row_names <- rownames(abide_h_p)[indices[, 1]]
col_names <- colnames(abide_h_p)[indices[, 2]]

HC <- cbind(row_names ,col_names)

############## L VS CN
# 1. 找到两个数据框中都不为 NA 的位置
common_non_na <- !is.na(abide_l_p) & !is.na(cabic_l_p)

# 2. 提取这些共同的非 NA 元素位置
indices <- which(common_non_na, arr.ind = TRUE)

# 3. 根据索引获取行名和列名
row_names <- rownames(abide_l_p)[indices[, 1]]
col_names <- colnames(abide_l_p)[indices[, 2]]

LC <- cbind(row_names ,col_names)


##################

cabic_h_z <- read.csv(file.path(cabicDir, "cabic_Z_diff_Cluster_2_vs_3.csv"))
rownames(cabic_h_z) <- cabic_h_z[,1]
cabic_h_z <- cabic_h_z[,-1]

cabic_l_z <- read.csv(file.path(cabicDir, "cabic_Z_diff_Cluster_1_vs_3.csv"))
rownames(cabic_l_z) <- cabic_l_z[,1]
cabic_l_z <- cabic_l_z[,-1]

abide_h_z <- read.csv(file.path(resuDir, "abide_Z_diff_Cluster_2_vs_3.csv"))
rownames(abide_h_z) <- abide_h_z[,1]
abide_h_z <- abide_h_z[,-1]

abide_l_z <- read.csv(file.path(resuDir, "abide_Z_diff_Cluster_1_vs_3.csv"))
rownames(abide_l_z) <- abide_l_z[,1]
abide_l_z <- abide_l_z[,-1]



# 初始化空数据框
result <- data.frame(cluster = character(), 
                     ROI1 = character(), 
                     ROI2 = character(), 
                     Z_abide = numeric(), 
                     p_abide = numeric(), 
                     Z_cabic = numeric(), 
                     p_cabic = numeric(),
                     stringsAsFactors = FALSE)

# 遍历 H vs CN 的配对
for (i in 1:nrow(HC)) {
  roi1 <- HC[i, 1]
  roi2 <- HC[i, 2]
  
  # 检查 ROI 是否在所有数据框中存在
  if (!is.na(abide_h_p[roi1, roi2]) && !is.na(cabic_h_p[roi1, roi2])) {
    # 提取 Z 和 P 值
    Z_abide <- abide_h_z[roi1, roi2]
    p_abide <- abide_h_p[roi1, roi2]
    Z_cabic <- cabic_h_z[roi1, roi2]
    p_cabic <- cabic_h_p[roi1, roi2]
    
    # 将结果添加到数据框中
    result <- rbind(result, data.frame(
      cluster = "H",
      ROI1 = roi1,
      ROI2 = roi2,
      Z_abide = Z_abide,
      p_abide = p_abide,
      Z_cabic = Z_cabic,
      p_cabic = p_cabic
    ))
  }
}

# 遍历 L vs CN 的配对
for (i in 1:nrow(LC)) {
  roi1 <- LC[i, 1]
  roi2 <- LC[i, 2]
  
  # 检查 ROI 是否在所有数据框中存在
  if (!is.na(abide_l_p[roi1, roi2]) && !is.na(cabic_l_p[roi1, roi2])) {
    # 提取 Z 和 P 值
    Z_abide <- abide_l_z[roi1, roi2]
    p_abide <- abide_l_p[roi1, roi2]
    Z_cabic <- cabic_l_z[roi1, roi2]
    p_cabic <- cabic_l_p[roi1, roi2]
    
    # 将结果添加到数据框中
    result <- rbind(result, data.frame(
      cluster = "L",
      ROI1 = roi1,
      ROI2 = roi2,
      Z_abide = Z_abide,
      p_abide = p_abide,
      Z_cabic = Z_cabic,
      p_cabic = p_cabic
    ))
  }
}

result <- subset(result, sign(Z_abide) == sign(Z_cabic))

write.csv(result, file.path(resuDir, "StruCova_Both.csv"), row.names = F)
