# this script is used to do ABIDE All Clustering (age 6-17.9 male only)
# 使用高斯混合模型算法对人群进行聚类分析，并把对分类特征的重要性评估出来
# 注意，这里只是用34个脑区作为分类指标
# Xue-Ru Fan 13 march 2024 @BNU
###################################################
# Part 0: 对每个特征值的centile进行标准化
# Part 1: 选择不同聚类数的轮廓系数，确定最佳聚类数目，png
# Part 2: 根据最佳聚类数目进行谱聚类，保存结果csv，保存了cluster编号和41个脑指标的原始centile
# Part 3: 使用特征消除法评估每个特征的贡献，保存排序的png和xlsx
##################################################

rm(list=ls())

packages <- c("mclust", "ggplot2", "cluster")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster")
resDate <- "240315"
newDate <- "240610"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_A_centile_dev_", resDate, ".csv")))

data_raw <- subset(abide_centile, dx == "ASD" & sex == "Male")

# 有4个ASD男性被去掉了，找到他们
# rows_with_na <- apply(data_raw, 1, function(x) any(is.na(x)))
# na_rows <- data_raw[rows_with_na, ]
# 去掉他们
data_centile <- na.omit(data_raw[, -2:-6])
data <- data_centile[, -1:-8]

################################## Part 0：标准化 ##################################################
data <- data.frame(scale(data))



################################## Part 1: 选择最佳聚类数目 ########################################
# 定义一个函数，用于计算不同聚类数目下的轮廓系数
calculate_silhouette <- function(data, max_clusters = 10) {
  silhouette_scores <- numeric(max_clusters - 1)
  
  for (k in 2:max_clusters) {
    # 执行GMM聚类
    set.seed(941205)
    gmm_model <- Mclust(as.matrix(data), G = k)
    cluster_membership <- gmm_model$classification
    
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

################################# save plot
name <- paste0("asd_male_GMM_Cluster_Silhouette_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 8, height = 6, units = "in", dpi = 500)

# 使用ggplot2进行绘图
ggplot(plot_data, aes(x = Cluster_Numbers, y = Silhouette_Scores)) +
  geom_line(color = "#00BFC4", size = 2) + # 定义线条颜色和大小
  geom_point(color = "#F8766D", size = 4) + # 定义点的颜色和大小
  theme_minimal() + # 使用简洁的主题
  labs(x = "聚类簇数",  y = "轮廓系数") + # 添加坐标轴标签
  theme(text = element_text(family = "STSong"),
    axis.title = element_text(face = "bold", size = 12), # 自定义X轴标题样式
    axis.text = element_text(face = "bold", size = 10), # 增大轴刻度文本大小
    panel.background = element_rect(fill = "white"), # 背景颜色
    panel.grid.major = element_line(color = "grey90"), # 主要网格线颜色
    panel.grid.minor = element_blank(), # 不显示次要网格线
    legend.position = "none") # 不显示图例
dev.off()


# ################################## Part 1：进行GMM聚类并简单可视化 ###############################
# set.seed(941205)
# model <- Mclust(data)
# summary(model)
# clusters <- model$classification
# 
# ################################# save result
# data_cluster <- cbind(clusters, data_centile)
# colnames(data_cluster)[1] <- "clusterID"
# name <- paste0("abide_A_asd_male_dev_GMM_Cluster_", newDate, ".csv")
# write.csv(data_cluster, file.path(resDir, name), row.names = F)
# 
# # 可视化
# pca_result <- prcomp(data, scale. = TRUE)
# data_pca <- data.frame(pca_result$x[, 1:2], cluster = as.factor(clusters))
# ggplot(data_pca, aes(x = PC1, y = PC2, color = cluster)) +
#   geom_point(shape = 16, size = 3) +  # shape=16为实心圆点
#   scale_color_manual(values = c("#0064b5", "#ff6347")) + 
#   theme_minimal() + # 使用简洁的主题
#   labs(x = "Dimension 1",  y = "Dimension 2") + # 添加坐标轴标签
#   theme(axis.title = element_text(face = "bold", size = 12), # 自定义X轴标题样式
#         axis.text = element_text(face = "bold", size = 10), # 增大轴刻度文本大小
#         panel.background = element_rect(fill = "white"), # 背景颜色
#         panel.grid.major = element_line(color = "grey90"), # 主要网格线颜色
#         panel.grid.minor = element_blank(), # 不显示次要网格线
#         legend.position = "none") # 不显示图例
# 
# ################################# save plot
# name <- file.path(plotDir, paste0("abide_A_asd_male_dev_GMM_PCA_", newDate, ".png"))
# ggsave(name, width = 8, height = 6, units = "in", dpi = 500)


################################## Part 2: 根据最佳聚类数目进行谱聚类 ##############################

# 执行谱聚类
set.seed(941205)
gmm_model <- Mclust(as.matrix(data), G = optimal_clusters)
# 获取聚类结果
cluster_membership <- gmm_model$classification

################################# save result
data_cluster <- cbind(cluster_membership, data_centile)

colnames(data_cluster)[1] <- "clusterID"
name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
write.csv(data_cluster, file.path(resDir, name), row.names = F)

# ################################## Part 2：保存模型结果 ##########################################
# 
# # 提取参数
# parameters <- model$parameters
# means <- parameters$mean
# variances <- parameters$variance
# proportions <- model$parameters$pro
# 
# # 创建数据框
# means_df <- as.data.frame(t(means))
# names(means_df) <- names(data)
# means_df$Cluster <- 1:nrow(means_df)
# 
# proportions_df <- as.data.frame(proportions)
# names(proportions_df) <- c("Proportion")
# proportions_df$Cluster <- 1:nrow(proportions_df)
# 
# summary_df <- merge(means_df, proportions_df, by = "Cluster")
# 
# # 保存为CSV文件
# name <- paste0("abide_A_asd_male_dev_GMM_Summary_", newDate, ".csv")
# write.csv(summary_df, file.path(resDir, name), row.names = F)


################################## Part 3: 使用特征消除法评估每个特征的贡献 ########################
# 计算轮廓系数的函数，使用GMM聚类
calculate_silhouette_for_gmm_clustering <- function(data, centers) {
  gmm_clustering_result <- Mclust(as.matrix(data), G = centers)
  cluster_membership <- gmm_clustering_result$classification
  silhouette_score <- mean(silhouette(cluster_membership, dist(data))[, 3])
  return(silhouette_score)
}

# 计算原始数据的轮廓系数（我们已经知道最佳聚类数为 optimal_clusters）
original_silhouette <- calculate_silhouette_for_gmm_clustering(data, optimal_clusters)

# 初始化一个向量，用于存储每个特征移除后的轮廓系数
silhouettes_after_removal <- numeric(length = ncol(data))

# 逐一移除每个特征，并计算轮廓系数
for (i in 1:ncol(data)) {
  data_without_feature <- data[, -i]
  set.seed(941205)
  silhouettes_after_removal[i] <- calculate_silhouette_for_gmm_clustering(data_without_feature,
                                                                          optimal_clusters)
}


feature_importance <- original_silhouette - silhouettes_after_removal # 比较移除特征前后的轮廓系数


features <- names(data) # 特征名称向量
# 创建对应的中文名向量
features_chinese <- c("颞上沟后侧", 
                      "前扣带尾部",
                      "额中回尾部",
                      "楔叶皮层",
                      "内嗅皮层",
                      "梭状回",
                      "顶下皮层",
                      "颞下回",
                      "扣带回峡部",
                      "外侧枕叶",
                      "外侧眶额",
                      "舌回",
                      "内侧眶额",
                      "颞中回",
                      "海马旁回",
                      "中央旁小叶",
                      "额下回盖部",
                      "额下回眶部",
                      "额下回三角部",
                      "距状沟周围皮层",
                      "中央后回",
                      "后扣带",
                      "中央前回",
                      "楔前叶",
                      "前扣带喙部",
                      "额中回喙部",
                      "额上回",
                      "顶上皮层",
                      "颞上回",
                      "缘上回",
                      "额极",
                      "颞极",
                      "颞横皮层",
                      "脑岛")


plot_data <- data.frame(Feature = features_chinese, Importance = feature_importance)

# 对数据框根据重要性评分进行排序
plot_data <- plot_data[order(plot_data$Importance, decreasing = TRUE), ]

################################# save plot
name_rank <- paste0("asd_male_GMM_Cluster_RMrank_", newDate, ".png")
CairoPNG(file.path(plotDir, name_rank), width = 7, height = 8, units = "in", dpi = 500)

# 使用ggplot2绘制条形图
ggplot(plot_data, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#66cdaa") +
  coord_flip() +  # 翻转坐标轴，使特征名称更容易阅读
  theme_minimal() +  # 使用简洁的主题
  labs(x = "DK分区", y = "人群聚类贡献度") +
  theme(text = element_text(family = "STSong"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜，以防重叠
dev.off()



################################# 保存该结果
gmm_rank <- data.frame(脑区 = features_chinese, Feature = features, 贡献度 = feature_importance)
gmm_rank <- gmm_rank[order(gmm_rank$贡献度, decreasing = TRUE), ]

name <- paste0("asd_male_GMM_Cluster_RMrank_", newDate, ".xlsx")
write.xlsx(gmm_rank, file.path(resDir, name))
