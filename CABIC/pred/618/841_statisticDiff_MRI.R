# 本代码用来分析高斯混合聚类的两组ASD男性的脑形态（百分位数）之间的差异
# 雪如 2024年2月27日于北师大办公室
################################
# Part 01: 站点
################
# 以上每部分都会保存统计的csv文件和绘图的png文件
################
# Part Z：保存P值文件csv，之后需要手动excel，筛选出P值显著（＜0.05）的位置，保存一个xlsx文件
################################


rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result", "Diff")
resDate <- "240928"

name <- paste0("cabic_cluster_predictions_", resDate, ".csv")
cluster <- read.csv(file.path(resuDir, name))

colnames(cluster)[4:(ncol(cluster)-1)] <- paste0(colnames(cluster)[4:(ncol(cluster)-1)], "_centile")

# evalu <- c("Median", "Mean", "SD") # 计算哪些统计值

# 新建一个空数据框用来存p值
Pvalue <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(Pvalue) <- c("t-test", "w-test")

for (i in 1:41) {
  
  var <- cluster[, c(ncol(cluster), i+3)]
  var <- na.omit(var)
  # var[, 2] <- scale(var[, 2])
  
  # 统计检验
  L <- subset(var, predicted_cluster == 1)
  H <- subset(var, predicted_cluster == 2)
  
  Pvalue[i, "t-test"] <- t.test(L[,2], H[,2])[["p.value"]]
  Pvalue[i, "w-test"] <- wilcox.test(L[,2], H[,2])[["p.value"]]

  
}

rownames(Pvalue) <- names(cluster)[4:(ncol(cluster)-1)]
Pvalue$`t-test` <- as.numeric(Pvalue$`t-test`)
Pvalue$`w-test` <- as.numeric(Pvalue$`w-test`)

name <- file.path(resuDir, paste0("statis_MRI_", resDate, ".csv"))
write.csv(Pvalue, name)
