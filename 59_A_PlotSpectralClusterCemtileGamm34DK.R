# 画出ASD男性两个亚型的34+7个脑指标的百分位数随龄发育GAMM曲线
# Xue-Ru Fan 04 Jan 2024 @BNU

rm(list=ls())

# load packages
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamm4",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster34DK")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster34DK/GAMMcentile")
resDate <- "240315"
newDate <- "240610"

# define variables
# id_group <- c("1", "2") # this code is for 2 clusters
ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
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
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)
    
    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(x, k = 4), data = cluster_data)
    
    # 准备预测数据
    pred_data <- data.frame(x = seq(6, 18, length.out = 100))
    
    # 获取预测值及其标准误差
    preds <- predict(gamm_model$gam, newdata = pred_data, type = "response", se = TRUE)
    
    # 向pred_data中添加预测值和置信区间
    pred_data$y_pred <- preds$fit
    pred_data$ci_lower <- preds$fit - 1.96 * preds$se.fit
    pred_data$ci_upper <- preds$fit + 1.96 * preds$se.fit
    pred_data$clusterID <- cluster  # 添加分组标识
    
    # 将预测结果添加到列表中
    predictions[[cluster]] <- pred_data
  }
  
  # 合并所有预测结果
  all_predictions <- bind_rows(predictions)
  
  
  # 使用ggplot2绘制预测结果
  ggplot(all_predictions, aes(x = x, y = y_pred, color = factor(clusterID))) +
    geom_line(lwd = 2, alpha = 1) +  # 绘制预测线
    # 绘制置信区间
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(clusterID)), alpha = 0.1,
                lwd = 0.1) +
    scale_color_manual(values = c("#719988", "#faa264")) +  # 自定义分组颜色
    scale_fill_manual(values = c("#719988", "#faa264")) +  # 自定义填充颜色
    # 添加散点图
    geom_point(data = studyFIT, aes(x = x, y = y, color = factor(clusterID)), alpha = .2, size = 2,
               shape = 16) +
    scale_x_continuous(breaks = seq(6, 18, 2),
                       label = c("6 yr", "8 yr", "10 yr", "12 yr", "14 yr", "16 yr", "18 yr")) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.position = "None")
  

  name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", volumeName, "_GAMMcentile_", newDate,".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 10,
         height = 10, unit = "cm")
}
