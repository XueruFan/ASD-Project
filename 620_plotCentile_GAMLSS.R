# 本代码用来可视化两个聚类人群脑指标centile的发育轨线
# GAMLSS模型
# Xue-Ru Fan 04 Jan 2024 @BNU
##################################
rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamlss",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/SpectralCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/SpectralCluster/Centile/GAMLSS")
resDate <- "240315"
newDate <- "240610"

# define variables
ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("asd_male_Spectral_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]
centile <- read.csv(file.path(phenoDir, paste0("abide_A_centile_", resDate, ".csv")))
All <- merge(cluster, centile, by = "participant", All.x = TRUE)

volumeNames <- names(centile)[c(which(names(centile) == "GMV"):which(names(centile) == "insula"))]

for (volumeName in volumeNames){
  # volumeName <- "GMV"
  studyFIT <- All[, c("clusterID", "Age", volumeName)]
  colnames(studyFIT)[2:3] <- c("x", "y")
  
  predictions <- list()
  for (cluster in unique(studyFIT$clusterID)) {
    # cluster <- "1"
    
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)
    
    # 拟合GAMLSS模型
    gamlss_model <- gamlss(y ~ pb(x), data = cluster_data, family = BE) # Beta分布
    
    # # 提取 global deviance, AIC, SBC 等值
    # global_deviance <- gamlss_model[["G.deviance"]]
    # aic_value <- gamlss_model[["aic"]]
    # sbc_value <- gamlss_model[["sbc"]] # BIC 在 GAMLSS 中也称为 SBC (Schwarz Bayesian Criterion)
    
    # 准备预测数据
    pred_data <- data.frame(x = seq(6, 18, length.out = 100))

    # 使用 predictAll 来获取所有的预测值
    preds <- predictAll(gamlss_model, newdata = pred_data)
    
    # 提取 mu（均值）的预测值，并将其添加到 pred_data
    # 对于 Beta 分布，你的 gamlss 模型会分别拟合 mu（均值）和 sigma（尺度），这里提取了 mu 的预测值
    #（pred$mu），因为在 Beta 分布中，mu 是我们感兴趣的均值预测。sigma 决定分布的尺度或形状，
    # 具体来说是数据的分散程度。sigma 越大，数据点的分散性越大；越小，数据点越集中。
    
    pred_data$y_pred <- preds$mu
    pred_data$clusterID <- cluster  # 添加分组标识
    
    # 将预测结果添加到列表中
    predictions[[cluster]] <- pred_data
    
  }
  
  # 合并所有预测结果
  all_predictions <- bind_rows(predictions)
  
  # 使用ggplot2绘制预测结果
  ggplot(all_predictions, aes(x = x, y = y_pred, color = factor(clusterID))) +
    geom_line(lwd = 2, alpha = 1) +  # 绘制预测线
    scale_color_manual(values = c("#719988", "#faa264")) +  # 自定义分组颜色
    # 添加散点图
    geom_point(data = studyFIT, aes(x = x, y = y, color = factor(clusterID)), alpha = .2, size = 2,
               shape = 16) +
    scale_x_continuous(breaks = seq(6, 18, 2),
                       label = c("6 yr", "8 yr", "10 yr", "12 yr", "14 yr", "16 yr", "18 yr")) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.position = "None")
  
  
  name <- paste0("SC_Centile_GAMLSS_", volumeName, "_", newDate,".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 10,
         height = 10, unit = "cm")
}
