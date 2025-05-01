# Structural covariance analysis on OoS scores of each two DK regions.
# Male, ASD, <13 years old, GMM Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list = ls())
packages <- c("tidyverse", "stats", "dplyr", "ggplot2", "reshape2", "pheatmap",
              "magrittr", "readr", "openxlsx")
sapply(packages, require, character.only = TRUE)

# 定义路径
abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/GMM513")
statDir <- file.path(abideDir, "Analysis/Statistic/GMM513")

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
name <- paste0("Cluster_240610.csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]

# 读取脑区数据
centile <- read.csv(file.path(phenoDir, "abide_All_centile_513_240315.csv"))

# 合并 Cluster 和 脑区数据
All <- merge(cluster, centile, by = "participant", all.x = TRUE)
cn <- subset(centile, dx == "CN" & sex == "Male")
cn$clusterID <- 3

# 使用 bind_rows() 按列名匹配合并
All <- bind_rows(All, cn)

# 修改各自站点名字的缩写
All$Site <- gsub("ABIDEII-NYU_1", "NYU", All$Site)
All$Site <- gsub("ABIDEII-NYU_2", "NYU", All$Site)
All$Site <- gsub("ABIDEII-KKI_1", "KKI", All$Site)
All$Site <- gsub("ABIDEII-SDSU_1", "SDSU", All$Site)
All$Site <- gsub("ABIDEII-UCLA_1", "UCLA", All$Site)
All$Site <- gsub("UCLA_1", "UCLA", All$Site)
All$Site <- gsub("UCLA_2", "UCLA", All$Site)
All$Site <- gsub("ABIDEII-GU_1", "GU", All$Site)
All$Site <- gsub("ABIDEII-UCD_1", "UCD", All$Site)
All$Site <- gsub("ABIDEII-EMC_1", "EMC", All$Site)
All$Site <- gsub("TRINITY", "TCD", All$Site)
All$Site <- gsub("ABIDEII-TCD_1", "TCD", All$Site)
All$Site <- gsub("ABIDEII-USM_1", "USM", All$Site)
All$Site <- gsub("ABIDEII-IU_1", "IU", All$Site)
All$Site <- gsub("ABIDEII-U_MIA_1", "UMIA", All$Site)
All$Site <- gsub("ABIDEII-ETH_1", "ETH", All$Site)
All$Site <- gsub("UM_1", "UM", All$Site)
All$Site <- gsub("UM_2", "UM", All$Site)
All$Site <- gsub("ABIDEII-OHSU_1", "OHSU", All$Site)
All$Site <- gsub("STANFORD", "SU1", All$Site)
All$Site <- gsub("ABIDEII-SU_2", "SU2", All$Site)
All$Site <- gsub("LEUVEN_2", "KUL", All$Site)
All$Site <- gsub("CALTECH", "CALT", All$Site)


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
  write.csv(z_matrix, file.path(statDir, paste0("abide_Cluster_", clust, "_z_correlation.csv")))
  write.csv(p_matrix, file.path(statDir, paste0("abide_Cluster_", clust, "_p.csv")))
  
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
  write.csv(z_diff, file.path(statDir, paste0("abide_Z_diff_Cluster_", clust1, "_vs_", clust2, ".csv")))
  write.csv(z_diff_p, file.path(statDir, paste0("abide_P_value_Cluster_", clust1, "_vs_", clust2, ".csv")))
}