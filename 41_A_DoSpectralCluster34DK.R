# this script is used to do ABIDE I and II Spectral Clustering (age 6-17.9 male only)
# 使用谱聚类算法对人群进行聚类分析，并把对分类特征的重要性评估出来
# Xue-Ru Fan 13 march 2024 @BNU

rm(list=ls())

packages <- c("kernlab", "cluster", "ggplot2", "Rtsne", "randomForest")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster")
resDate <- "240315" # 上次分析
newDate <- "240610" # 这次分析

abide_centile <- read.csv(file.path(dataDir, paste0("abide_A_centile_dev_", resDate, ".csv")))

data_raw <- subset(abide_centile, dx == "ASD" & sex == "Male")

# 有4个ASD男性被去掉了，找到他们
# rows_with_na <- apply(data_raw, 1, function(x) any(is.na(x)))
# na_rows <- data_raw[rows_with_na, ]
# 去掉他们
data_centile <- na.omit(data_raw[, -2:-13])
data <- data_centile[, -1]

################################## Part 0：标准化 ##################################################
# 为什么要标准化
# 消除量纲影响：不同特征的量纲不同，会导致在计算距离时量纲大的特征主导结果。标准化可以消除量纲的影响。
# 提高算法效率：许多算法（如梯度下降法）对特征的尺度敏感，标准化可以加快收敛速度。
# 提高准确性：在聚类等算法中，标准化可以使所有特征对相似度计算同等重要，从而提高结果的准确性。

data <- data.frame(scale(data))
################################## Part 1: 选择最佳聚类数目 ########################################
# 定义一个函数，用于计算不同聚类数目下的轮廓系数
calculate_silhouette <- function(data, max_clusters = 10) {
  silhouette_scores <- numeric(max_clusters - 1)
  
  for (k in 2:max_clusters) {
    # 执行谱聚类
    set.seed(941205)
    cluster_result <- specc(as.matrix(data), centers = k)
    cluster_membership <- cluster_result@.Data
    
    # 计算轮廓系数
    silhouette_scores[k-1] <- mean(silhouette(cluster_membership, dist(data))[, 3])
  }
  
  return(silhouette_scores)
}

# 计算轮廓系数
max_clusters <- 10 # 可以根据需要调整最大聚类数目
silhouette_scores <- calculate_silhouette(data, max_clusters)

# 找出轮廓系数最大时的聚类数
optimal_clusters <- which.max(silhouette_scores) + 1

# 可视化轮廓系数
# 创建一个数据框，用于绘图
plot_data <- data.frame(Cluster_Numbers = 2:max_clusters, Silhouette_Scores = silhouette_scores)

# 使用ggplot2进行绘图
ggplot(plot_data, aes(x = Cluster_Numbers, y = Silhouette_Scores)) +
  geom_line(color = "#00BFC4", size = 2) + # 定义线条颜色和大小
  geom_point(color = "#F8766D", size = 4) + # 定义点的颜色和大小
  theme_minimal() + # 使用简洁的主题
  labs(x = "Number of Clusters",  y = "Silhouette Scores") + # 添加坐标轴标签
  theme(axis.title = element_text(face = "bold", size = 12), # 自定义X轴标题样式
        axis.text = element_text(face = "bold", size = 10), # 增大轴刻度文本大小
        panel.background = element_rect(fill = "white"), # 背景颜色
        panel.grid.major = element_line(color = "grey90"), # 主要网格线颜色
        panel.grid.minor = element_blank(), # 不显示次要网格线
        legend.position = "none") # 不显示图例

################################# save plot
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_Silhouette_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 8, height = 6, units = "in", dpi = 500)


################################## Part 2: 根据最佳聚类数目进行谱聚类 ##############################

# 执行谱聚类
set.seed(941205)
cluster_result <- specc(as.matrix(data), centers = optimal_clusters)

# 获取聚类结果
cluster_membership <- cluster_result@.Data

################################# save result
data_cluster <- cbind(cluster_membership, data_centile)
# 结果发现标准化和不标准化特征值对聚类结果没有影响
colnames(data_cluster)[1] <- "clusterID"
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
write.csv(data_cluster, file.path(resDir, name), row.names = F)

################################## Part 3: 使用特征消除法评估每个特征的贡献 ########################
# 定义一个函数，用于计算给定数据的谱聚类轮廓系数
calculate_silhouette_for_spectral_clustering <- function(data, centers) {
  spectral_clustering_result <- specc(as.matrix(data), centers = centers)
  cluster_membership <- spectral_clustering_result@.Data
  silhouette_score <- mean(silhouette(cluster_membership, dist(data))[, 3])
  return(silhouette_score)
}

# 计算原始数据的轮廓系数（我们已经知道最佳聚类数为 optimal_clusters）
original_silhouette <- calculate_silhouette_for_spectral_clustering(data, optimal_clusters)

# 初始化一个向量，用于存储每个特征移除后的轮廓系数
silhouettes_after_removal <- numeric(length = ncol(data))

# 逐一移除每个特征，并计算轮廓系数
for (i in 1:ncol(data)) {
  data_without_feature <- data[, -i]
  set.seed(941205)
  silhouettes_after_removal[i] <- calculate_silhouette_for_spectral_clustering(data_without_feature,
                                                                               optimal_clusters)
}

feature_importance <- original_silhouette - silhouettes_after_removal # 比较移除特征前后的轮廓系数
features <- names(data) # 特征名称向量
plot_data <- data.frame(Feature = features, Importance = feature_importance)

# 对数据框根据重要性评分进行排序
plot_data <- plot_data[order(plot_data$Importance, decreasing = TRUE), ]

# 使用ggplot2绘制条形图
ggplot(plot_data, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#66cdaa") +
  coord_flip() +  # 翻转坐标轴，使特征名称更容易阅读
  theme_minimal() +  # 使用简洁的主题
  labs(x = "Feature", y = "Importance") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜，以防重叠

################################# save plot
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_RM_Rank_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 8, units = "in", dpi = 500)


################################# 保存该结果
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_RM_Rank_", newDate, ".csv")
write.csv(plot_data, file.path(resDir, name), row.names = F)