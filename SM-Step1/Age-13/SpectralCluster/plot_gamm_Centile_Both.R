# Plot GAMM between OoS scores and age within subgroup H for both CABIC and ABIDE
# Xue-Ru Fan 04 Jan 2024 @BNU
##################################
rm(list=ls())

# load packages
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamm4", "openxlsx",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/513")
resDate <- "240928"
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect513/GAMM_Both")
newDate <- "240610"
oldDate <- "240315"


# define variables
ageRange <- log((seq(1, 13, 0.1)*365.245)+280)

name <- paste0("cabic_cluster_predictions_513_", resDate, ".csv")
cabic <- read.csv(file.path(resuDir, name))
cabic <- subset(cabic, predicted_cluster == "2")
cabic$Data <- "CABIC"

name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]
abide <- read.csv(file.path(phenoDir, paste0("abide_A_centile_", oldDate, ".csv")))
abide <- merge(cluster, abide, by = "participant", All.x = TRUE)
abide <- subset(abide, clusterID == "2")
abide$Data <- "ABIDE"

volumeNames <- names(cabic)[c(which(names(cabic) == "bankssts"):which(names(cabic) == "insula"))]

for (volumeName in volumeNames) {
  # volumeName <- "bankssts"
  
  studyFIT_cabic <- cabic[, c("Data", "Age", "TCV", volumeName)]
  studyFIT_abide <- abide[, c("Data", "Age", "TCV", volumeName)]
  studyFIT <- rbind(studyFIT_cabic, studyFIT_abide)
  
  colnames(studyFIT)[4] <- "y"
  
  # predictions_site <- list()
  predictions_no_site <- list()
  
  # 定义颜色（黄绿色和橙色）对应每个 predicted_cluster
  cluster_colors <- c("#d9ca39", "#e29135")

  
  # 第二部分：忽略站点效应，收集彩色线条数据，并恢复置信区间
  for (site in unique(studyFIT$Data)) {

    # 筛选当前cluster的数据
    site_data <- filter(studyFIT, Data == site)
    
    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(Age, bs = "tp", k = 4) + TCV, data = site_data)
    
    # 准备预测数据（忽略站点）
    pred_data_no_site <- data.frame(Age = seq(1,13, length.out = 100))
    pred_data_no_site$TCV <- mean(site_data$TCV)  # 使用TCV的平均值
    
    
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
    
    pred_data_no_site$Data <- site  # 添加 predicted_cluster 标识
    
    # 将预测结果添加到列表中
    predictions_no_site[[site]] <- pred_data_no_site
  }
  
  
  # 合并忽略站点效应的预测结果
  all_predictions_no_site <- bind_rows(predictions_no_site)
  
  
  studyFIT$Data <- as.factor(studyFIT$Data)
  
  # 创建 ggplot2 图形，先添加透明的彩色线条来表示站点效应
  ggplot() +
    geom_point(data = studyFIT, aes(x = Age, y = y, color = Data, shape = Data, size = Data),
               alpha = 0.3) + # 透明原始数据点
    scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
    scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
    annotate("segment", x = 1, xend = 13, y = 0.5, yend = 0.5, color = "#e6e6e6",
             linetype = "solid", size = 3, alpha = 0.5) +
    # geom_line(data = all_predictions_site, aes(x = x, y = y_pred,
    #                                            group = interaction(Site, predicted_cluster), 
    #                                            color = factor(predicted_cluster)),
    #           lwd = 0.5, alpha = 0.2) +  # 彩色透明线条表示站点效应
    # 添加彩色的忽略站点效应的线条和置信区间
    geom_ribbon(data = all_predictions_no_site, aes(x = Age, ymin = ci_lower, ymax = ci_upper,
                                                    fill = Data), 
                alpha = 0.2, linetype = 0) +  # 绘制彩色置信区间带
    geom_line(data = all_predictions_no_site,
              aes(x = Age, y = y_pred, color = factor(Data)), lwd = 3) +
    scale_color_manual(values = cluster_colors) +  # 自定义颜色（黄绿色和橙色）
    scale_fill_manual(values = cluster_colors) +   # 填充颜色与线条颜色一致
    scale_x_continuous(breaks = seq(1,13, 2),
                       labels = c("1 yr", "3", "5", "7", "9", "11", "13")) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(legend.position = c(0.03, 0.06),
          legend.title = element_blank(),
          legend.key.width  = unit(1, "cm"),
          axis.text = element_text(size = 12))
  
  
  # 保存图像
  name <- paste0("gamm_Centile_", volumeName, "_", resDate, ".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 20, height = 20, unit = "cm")
}
