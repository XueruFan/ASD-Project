# Analyze differences of OoS centiles of MRI measures between the clusters
# Male, ASD, <13 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################


rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect510")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect510")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect510")
resDate <- "240315"
newDate <- "240610"

name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
colnames(cluster)[3:ncol(cluster)] <- paste0(colnames(cluster)[3:ncol(cluster)], "_centile")

# evalu <- c("Median", "Mean", "SD") # 计算哪些统计值

# 新建一个空数据框用来存p值
Pvalue <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(Pvalue) <- c("t-test", "w-test")

for (i in 1:41) {
  
  var <- cluster[, c(1, i+2)]
  var <- na.omit(var)
  # var[, 2] <- scale(var[, 2])
  
  # 统计检验
  L <- subset(var, clusterID == 1)
  H <- subset(var, clusterID == 2)
  
  Pvalue[i, "t-test"] <- t.test(L[,2], H[,2])[["p.value"]]
  Pvalue[i, "w-test"] <- wilcox.test(L[,2], H[,2])[["p.value"]]

  
}

rownames(Pvalue) <- names(cluster)[3:ncol(cluster)]
Pvalue$`t-test` <- as.numeric(Pvalue$`t-test`)
Pvalue$`w-test` <- as.numeric(Pvalue$`w-test`)

name <- file.path(statiDir, paste0("statis_MRI_", newDate, ".csv"))
write.csv(Pvalue, name)
