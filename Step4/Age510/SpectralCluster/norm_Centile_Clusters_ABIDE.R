# Perform Jarque-Bera normality tests on the individual OoS centile of ABIDE brain metrics for each 
# clusters
# Xue-Ru Fan 13 march 2024 @BNU
##################################################

rm(list=ls())

packages <- c("mclust", "ggplot2", "cluster", "Cairo", "tseries", "openxlsx")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Spect510")
newDate <- "240610"

abide_centile <- read.csv(file.path(resDir, paste0("Cluster_", newDate, ".csv")))

cluster_1 <- subset(abide_centile, clusterID == 1)[, -1:-2]
cluster_2 <- subset(abide_centile, clusterID == 2)[, -1:-2]

############ cluster 1
# 检查数据框中有无缺失值
any(is.na(cluster_1))

# 对每列进行 Jarque-Bera 检验并提取 p 值
p_values <- sapply(cluster_1, function(column) jarque.bera.test(column)$p.value)

# 创建新的数据框，包含列名和对应的 p 值
datap <- data.frame(
  # Column = names(data),  # 第一列：data 数据框中的列名
  P_Value = p_values     # 第二列：对应的 p 值
)

# 将 P_Value 列转换为字符类型，保留适当的小数位数
# datap$P_Value <- sprintf("%.4f", datap$P_Value)

write.csv(datap, file.path(resDir, paste0("Cluster_1_Centile_NormTest_", newDate, ".csv")))


############ cluster 2
# 检查数据框中有无缺失值
any(is.na(cluster_2))

# 对每列进行 Jarque-Bera 检验并提取 p 值
p_values <- sapply(cluster_2, function(column) jarque.bera.test(column)$p.value)

# 创建新的数据框，包含列名和对应的 p 值
datap <- data.frame(
  # Column = names(data),  # 第一列：data 数据框中的列名
  P_Value = p_values     # 第二列：对应的 p 值
)

# 将 P_Value 列转换为字符类型，保留适当的小数位数
# datap$P_Value <- sprintf("%.4f", datap$P_Value)

write.csv(datap, file.path(resDir, paste0("Cluster_2_Centile_NormTest_", newDate, ".csv")))
