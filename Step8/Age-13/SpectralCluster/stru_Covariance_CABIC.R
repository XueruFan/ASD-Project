# Structural covariance analysis on OoS scores of each two DK regions.
# Male, ASD, <13 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list = ls())

# 加载必要的包
packages <- c("tidyverse", "stats", "dplyr", "ggplot2", "reshape2", "pheatmap",
              "magrittr", "readr", "openxlsx")
sapply(packages, require, character.only = TRUE)

statDir <- ("E:/PhDproject/CABIC/result/pred/Spect513/StruCorr")

# 自定义 cor.mtest() 函数来计算 p 值矩阵
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0  # 对角线的 p 值设为 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value  # 填入 p 值
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# 读取 Cluster 数据
cluster <- read.csv("E:/PhDproject/CABIC/result/pred/Spect513/cabic_cluster_predictions_513_240928.csv")
colnames(cluster)[45] <- "clusterID"

# 读取脑区数据
centile <- read.csv("E:/PhDproject/CABIC/result/cabic_centile_240928.csv")
centile <- subset(centile, dx == "CN" & sex == "Male")
centile <- centile[, -4:-5]
centile$clusterID <- 3

# 使用 bind_rows() 按列名匹配合并
All <- bind_rows(cluster, centile)


volumeNames <- names(centile)[which(names(centile) == "bankssts"):which(names(centile) == "insula")]
results <- list()

# 遍历每个 clusterID
for (clust in unique(All$clusterID)) {
  data_subset <- All %>% filter(clusterID == clust)
  
  # 初始化残差矩阵
  residuals_matrix <- matrix(NA, nrow = nrow(data_subset), ncol = length(volumeNames))
  colnames(residuals_matrix) <- volumeNames
  rownames(residuals_matrix) <- data_subset$participant
  
  for (v in volumeNames) {
    model <- lm(as.formula(paste(v, "~ Site + TCV")), data = data_subset, na.action = na.exclude)
    residuals_vector <- rep(NA, nrow(data_subset))
    residuals_vector[as.numeric(names(residuals(model)))] <- residuals(model)
    residuals_matrix[, v] <- residuals_vector
  }
  
  # 计算相关性矩阵并进行 Fisher z 变换
  corr_matrix <- cor(residuals_matrix, use = "complete.obs", method = "pearson")
  z_matrix <- 0.5 * log((1 + corr_matrix) / (1 - corr_matrix))  # Fisher z 变换
  
  p_matrix <- cor.mtest(residuals_matrix)
  p_matrix[upper.tri(p_matrix, diag = TRUE)] <- NA
  p_matrix[p_matrix >= 0.05] <- NA
  
  # 保存 z 变换后的相关性矩阵
  write.csv(z_matrix, file.path(statDir, paste0("cabic_Cluster_", clust, "_z_correlation.csv")))
  write.csv(p_matrix, file.path(statDir, paste0("cabic_Cluster_", clust, "_p.csv")))
  
  results[[as.character(clust)]] <- list(correlation = z_matrix, p_values = p_matrix)
}

# 两两比较
cluster_pairs <- combn(unique(All$clusterID), 2, simplify = FALSE)

for (pair in cluster_pairs) {
  clust1 <- as.character(pair[1])
  clust2 <- as.character(pair[2])
  
  cluster1_z <- results[[clust1]]$correlation
  cluster2_z <- results[[clust2]]$correlation
  
  # 计算 z 值差异
  z_diff <- (cluster1_z - cluster2_z) / 
    sqrt((1 / (nrow(All %>% filter(clusterID == pair[1])) - 3)) +
           (1 / (nrow(All %>% filter(clusterID == pair[2])) - 3)))
  
  z_diff_p <- 2 * (1 - pnorm(abs(z_diff)))  # 差异的 p 值
  z_diff_p[upper.tri(z_diff_p, diag = TRUE)] <- NA
  z_diff_p[z_diff_p >= 0.05] <- NA
  
  # 保存 Z 值和 p 值矩阵
  write.csv(z_diff, file.path(statDir, paste0("cabic_Z_diff_Cluster_", clust1, "_vs_", clust2, ".csv")))
  write.csv(z_diff_p, file.path(statDir, paste0("cabic_P_value_Cluster_", clust1, "_vs_", clust2, ".csv")))
}