# this script is used to do ABIDE I and II GMM Clustering (age 6-17.9 male only)
# Xue-Ru Fan 13 march 2024 @BNU

rm(list=ls())

packages <- c("mclust", "ggplot2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/GmmCluster")
resDate <- "240315"
newDate <- "240610"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_A_centile_dev_", resDate, ".csv")))

data_raw <- subset(abide_centile, dx == "ASD" & sex == "Male")

# 有4个ASD男性被去掉了，找到他们
# rows_with_na <- apply(data_raw, 1, function(x) any(is.na(x)))
# na_rows <- data_raw[rows_with_na, ]
# 去掉他们
data_centile <- na.omit(data_raw[, -2:-6])
data <- data_centile[, -1]
data <- data.frame(scale(data))

################################## Part 1：进行GMM聚类并简单可视化 #################################
set.seed(941205)
model <- Mclust(data)
summary(model)
clusters <- model$classification

################################# save result
data_cluster <- cbind(clusters, data_centile)
colnames(data_cluster)[1] <- "clusterID"
name <- paste0("abide_A_asd_male_dev_GMM_Cluster_", newDate, ".csv")
write.csv(data_cluster, file.path(resDir, name), row.names = F)

# 可视化
pca_result <- prcomp(data, scale. = TRUE)
data_pca <- data.frame(pca_result$x[, 1:2], cluster = as.factor(clusters))
ggplot(data_pca, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(shape = 16, size = 3) +  # shape=16为实心圆点
  scale_color_manual(values = c("#0064b5", "#ff6347")) + 
  theme_minimal() + # 使用简洁的主题
  labs(x = "Dimension 1",  y = "Dimension 2") + # 添加坐标轴标签
  theme(axis.title = element_text(face = "bold", size = 12), # 自定义X轴标题样式
        axis.text = element_text(face = "bold", size = 10), # 增大轴刻度文本大小
        panel.background = element_rect(fill = "white"), # 背景颜色
        panel.grid.major = element_line(color = "grey90"), # 主要网格线颜色
        panel.grid.minor = element_blank(), # 不显示次要网格线
        legend.position = "none") # 不显示图例

################################# save plot
name <- file.path(plotDir, paste0("abide_A_asd_male_dev_GMM_PCA_", newDate, ".png"))
ggsave(name, width = 8, height = 6, units = "in", dpi = 500)


################################## Part 2：保存模型结果 ############################################

# 提取参数
parameters <- model$parameters
means <- parameters$mean
variances <- parameters$variance
proportions <- model$parameters$pro

# 创建数据框
means_df <- as.data.frame(t(means))
names(means_df) <- names(data)
means_df$Cluster <- 1:nrow(means_df)

proportions_df <- as.data.frame(proportions)
names(proportions_df) <- c("Proportion")
proportions_df$Cluster <- 1:nrow(proportions_df)

summary_df <- merge(means_df, proportions_df, by = "Cluster")

# 保存为CSV文件
name <- paste0("abide_A_asd_male_dev_GMM_Summary_", newDate, ".csv")
write.csv(summary_df, file.path(resDir, name), row.names = F)


################################## Part 3：使用特征消除法评估每个特征的贡献 ########################
model_bic <- model$bic

feature_importance <- vector("numeric", length = ncol(data))

for (i in 1:ncol(data)) {
  data_subset <- data[, -i]  # 移除当前特征
  set.seed(941205)
  model_subset <- Mclust(data_subset)
  feature_importance[i] <- model_bic - model_subset$bic  # BIC差异
}

# 根据BIC差异评估特征重要性
names(feature_importance) <- colnames(data)
features <- names(data) # 特征名称向量
plot_data <- data.frame(Feature = features, Importance = feature_importance)

# 对数据框根据重要性评分进行排序
plot_data <- plot_data[order(plot_data$Importance), ]

# 使用ggplot2绘制条形图
ggplot(plot_data, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#66cdaa") +
  coord_flip() +  # 翻转坐标轴，使特征名称更容易阅读
  theme_minimal() +  # 使用简洁的主题
  labs(x = "Feature", y = "Importance") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜，以防重叠

################################# save plot
name <- file.path(plotDir, paste0("abide_A_asd_male_dev_GMM_RM_Rank_", newDate, ".png"))
ggsave(name, width = 7, height = 8, units = "in", dpi = 500)


################################# 保存该结果
name <- paste0("abide_A_asd_male_dev_GMM_RM_Rank_", newDate, ".csv")
write.csv(plot_data, file.path(resDir, name), row.names = F)
