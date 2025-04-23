# Perform Jarque-Bera normality tests on the individual OoS centile of CABIC brain metrics
# Xue-Ru Fan 13 march 2024 @BNU
##################################################

rm(list=ls())

packages <- c("mclust", "ggplot2", "cluster", "Cairo", "tseries", "openxlsx")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result")
newDate <- "250117"
oldDate <- "240928"

name <- paste0("cabic_centile_510_", oldDate, ".csv")
cabic_centile <- read.csv(file.path(resuDir, name))
data_raw <- subset(cabic_centile, dx == "ASD" & sex == "Male")

data_centile <- na.omit(data_raw[, -2:-5])
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
# datap$P_Value <- sprintf("%.4f", datap$P_Value)

write.csv(datap, file.path(resuDir, paste0("centile_normtest_510_", newDate, ".csv")))
