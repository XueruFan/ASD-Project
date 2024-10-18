# 对脑指标的个体偏离百分位数进行正态性检验
    # Jarque-Bera 检验：
    # 用途：基于偏度和峰度来测试数据是否来自正态分布。
    # 原理：正态分布的偏度应该为0，峰度为3
    # Jarque-Bera 检验通过计算实际偏度和峰度与正态分布的偏离程度，来检验正态性。
    # 假设：
    # 零假设（H0H_0H0）：数据服从正态分布。
    # 备择假设（H1H_1H1）：数据不服从正态分布。
    # 使用场景：Jarque-Bera 检验适合大样本数据。
# Xue-Ru Fan 13 march 2024 @BNU
###################################################
# csv
##################################################

rm(list=ls())

packages <- c("mclust", "ggplot2", "cluster", "Cairo", "tseries", "openxlsx")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster")
resDate <- "240315"
newDate <- "240610"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_All_centile_515_", resDate, ".csv")))

data_raw <- subset(abide_centile, dx == "ASD" & sex == "Male")

data_centile <- na.omit(data_raw[, -2:-6])
data <- data_centile[, -1]


# 检查数据框中有无缺失值
any(is.na(data))

# 对每列进行 Jarque-Bera 检验并提取 p 值
p_values <- sapply(data, function(column) jarque.bera.test(column)$p.value)

# 创建新的数据框，包含列名和对应的 p 值
datap <- data.frame(
  # Column = names(data),  # 第一列：data 数据框中的列名
  P_Value = p_values     # 第二列：对应的 p 值
)

# 将 P_Value 列转换为字符类型，保留适当的小数位数
datap$P_Value <- sprintf("%.4f", datap$P_Value)

write.csv(datap, file.path(resDir, paste0("Centile_NormTest_515_", newDate, ".csv")))
