# Plot GAMM between OoS scores and age of each cluster for CABIC
# Xue-Ru Fan 04 Jan 2024 @BNU
##################################
rm(list=ls())

# load packages
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamm4", "openxlsx",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/510")
plotDir <- file.path(cabicDir, "result/pred/510/Plot/GAMM")
resDate <- "240928"

# define variables
ageRange <- log((seq(5, 10, 0.1)*365.245)+280)

name <- paste0("cabic_cluster_predictions_510_", resDate, ".csv")
centile <- read.csv(file.path(resuDir, name))

centile$predicted_cluster[centile$predicted_cluster == "1"] <- "H"
centile$predicted_cluster[centile$predicted_cluster == "2"] <- "1"
centile$predicted_cluster[centile$predicted_cluster == "H"] <- "2"

volumeNames <- names(centile)[c(which(names(centile) == "bankssts"):which(names(centile) == "insula"))]

for (volumeName in volumeNames) {
  # volumeName <- "bankssts"
  
  studyFIT <- centile[, c("predicted_cluster", "Age", "TCV", volumeName)]
  colnames(studyFIT)[4] <- "y"
  
  # predictions_site <- list()
  predictions_no_site <- list()
  
  # 定义颜色（黄绿色和橙色）对应每个 predicted_cluster
  cluster_colors <- c("#719988", "#faa264")
  
  # # 第一部分：包含站点固定效应，收集透明颜色线条数据
  # for (cluster in unique(studyFIT$predicted_cluster)) {
  #   # 筛选当前 predicted_cluster 的数据
  #   cluster_data <- filter(studyFIT, predicted_cluster == cluster)
  #   cluster_data$Site <- factor(cluster_data$Site)
  #   
  #   # 拟合 GAMM 模型（包含站点固定效应）
  #   gamm_model_site <- gamm4(y ~ s(x, k = 4) + Site, data = cluster_data)
  #   
  #   # 准备预测数据，增加 Site 信息
  #   pred_data_site <- expand.grid(x = seq(6, 18, length.out = 100), Site = levels(cluster_data$Site))
  #   
  #   # 获取预测值（包含站点效应）
  #   preds_site <- predict(gamm_model_site$gam, newdata = pred_data_site, type = "response")
  #   
  #   # 将预测值加入 pred_data_site
  #   pred_data_site$y_pred <- preds_site
  #   pred_data_site$predicted_cluster <- cluster  # 添加 predicted_cluster 标识
  #   
  #   # 将预测结果添加到列表中
  #   predictions_site[[cluster]] <- pred_data_site
  # }
  # 
  # # 合并所有站点的预测结果
  # all_predictions_site <- bind_rows(predictions_site)
  
  # 第二部分：忽略站点效应，收集彩色线条数据，并恢复置信区间
  for (cluster in unique(studyFIT$predicted_cluster)) {
    # cluster <- 1
    
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, predicted_cluster == cluster)
    
    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(Age, bs = "tp", k = 4) + TCV, data = cluster_data)
    
    # 准备预测数据（忽略站点）
    pred_data_no_site <- data.frame(Age = seq(5, 10, length.out = 100))
    pred_data_no_site$TCV <- mean(cluster_data$TCV)  # 使用TCV的平均值
    
    
    # 获取预测值和标准误（忽略站点效应）
    preds_no_site <- predict(gamm_model$gam, newdata = pred_data_no_site,
                             type = "response", se.fit = TRUE)
    
     # 将预测值加入 pred_data_no_site
    pred_data_no_site$y_pred <- preds_no_site$fit
    pred_data_no_site$se_fit <- preds_no_site$se.fit
    
    # 计算置信区间上下限
    pred_data_no_site$ci_lower <- pred_data_no_site$y_pred - 1.96 * pred_data_no_site$se_fit
    pred_data_no_site$ci_upper <- pred_data_no_site$y_pred + 1.96 * pred_data_no_site$se_fit
    
    # 截断置信区间使其保持在0和1之间
    pred_data_no_site$ci_lower <- pmax(pred_data_no_site$ci_lower, 0)
    pred_data_no_site$ci_upper <- pmin(pred_data_no_site$ci_upper, 1)
    
    pred_data_no_site$predicted_cluster <- cluster  # 添加 predicted_cluster 标识
    
    # 将预测结果添加到列表中
    predictions_no_site[[cluster]] <- pred_data_no_site
  }
  
  # 合并忽略站点效应的预测结果
  all_predictions_no_site <- bind_rows(predictions_no_site)
  
  # 创建 ggplot2 图形，先添加透明的彩色线条来表示站点效应
  ggplot() +
    geom_point(data = studyFIT, aes(x = Age, y = y, color = factor(predicted_cluster)),
               alpha = 0.2, size = 2, shape = 16) + # 透明原始数据点
    annotate("segment", x = 5, xend = 10, y = 0.5, yend = 0.5, color = "#e6e6e6",
             linetype = "solid", size = 3, alpha = 0.5) +
    # geom_line(data = all_predictions_site, aes(x = x, y = y_pred,
    #                                            group = interaction(Site, predicted_cluster), 
    #                                            color = factor(predicted_cluster)),
    #           lwd = 0.5, alpha = 0.2) +  # 彩色透明线条表示站点效应
    # 添加彩色的忽略站点效应的线条和置信区间
    geom_ribbon(data = all_predictions_no_site, aes(x = Age, ymin = ci_lower, ymax = ci_upper,
                                                    fill = factor(predicted_cluster)), 
                alpha = 0.2, linetype = 0) +  # 绘制彩色置信区间带
    geom_line(data = all_predictions_no_site,
              aes(x = Age, y = y_pred, color = factor(predicted_cluster)), lwd = 3) +
    scale_color_manual(values = cluster_colors) +  # 自定义颜色（黄绿色和橙色）
    scale_fill_manual(values = cluster_colors) +   # 填充颜色与线条颜色一致
    scale_x_continuous(breaks = seq(5, 10, 2),
                       labels = c("5yr", "7", "9")) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.position = "None")  # 不显示图例
  
  
  # 保存图像
  name <- paste0("gamm_Centile_", volumeName, "_", resDate, ".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 20, height = 20, unit = "cm")
}
