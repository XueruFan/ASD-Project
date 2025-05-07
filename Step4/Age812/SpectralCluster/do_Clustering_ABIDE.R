# Performe spectral clustering analysis using 34 regional volume OoS scores as classification
# features to subgroup male ASD individuals (5~8.9 years) from the ABIDE dataset.
# Xue-Ru Fan 13 march 2024 @BNU
###################################################
# Part 1: Choose the optimal number of clusters by computing silhouette scores for different cluster
#         numbers (ranging from 2 to 10).
# Part 2: Use the optimal cluster number determined from last part to perform spectral clustering.
# SM Part 1: Evaluate feature importance using recursive feature elimination (RFE) to identify the most
#         distinct feature for clustering.
##################################################

rm(list=ls())

packages <- c("kernlab", "cluster", "ggplot2", "Rtsne", "openxlsx")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/Spect812")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect812")
resDate <- "240315"
newDate <- "240610"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_All_centile_812_", resDate, ".csv")))

data_raw <- subset(abide_centile, dx == "ASD" & sex == "Male")
data_centile <- na.omit(data_raw[, -2:-6])
data <- data_centile[, -1:-8] # 全局脑指标不作为分类变量

################## Part 1: Choose the optimal number of clusters ###################################
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
optimal_clusters
# # 可视化轮廓系数
# # 创建一个数据框，用于绘图
# plot_data <- data.frame(Cluster_Numbers = 2:max_clusters, Silhouette_Scores = silhouette_scores)
# 
# ggplot(plot_data, aes(x = Cluster_Numbers, y = Silhouette_Scores)) +
#   geom_line(color = "#aec4dc", size = 3) + # 定义线条颜色和大小
#   geom_point(color = "#e78f89", size = 6) + # 定义点的颜色和大小
#   theme_minimal() + # 使用简洁的主题
#   labs(x = "Number of Clusters", y = "Silhouette Scores") + # 添加坐标轴标签
#   theme(
#     axis.title = element_text(face = "bold", size = 16), # 自定义轴标题样式
#     axis.text = element_text(face = "bold", size = 16), # 增大轴刻度文本大小
#     panel.background = element_rect(fill = "white", color = NA), # 移除背景外框
#     panel.grid.major = element_line(color = "grey90", size = 0.5), # 主要网格线颜色和粗细
#     panel.grid.minor = element_line(color = "grey95", size = 0.25), # 次要网格线颜色和粗细
#     panel.border = element_blank(), # 移除黑色外框
#     legend.position = "none" # 不显示图例
#     )
# 
# ################################# save plot
# name <- paste0("Silhouette_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 8, height = 6, units = "in", dpi = 500)


################## Part 2: Use the optimal cluster number to perform spectral clustering ###########

# 执行谱聚类
set.seed(941205)
cluster_result <- specc(as.matrix(data), centers = optimal_clusters)

# 获取聚类结果
cluster_membership <- cluster_result@.Data

################################# save result
data_cluster <- cbind(cluster_membership, data_centile)
# 结果发现标准化和不标准化特征值对聚类结果没有影响
colnames(data_cluster)[1] <- "clusterID"
name <- paste0("Cluster_", newDate, ".csv")
write.csv(data_cluster, file.path(resDir, name), row.names = F)
