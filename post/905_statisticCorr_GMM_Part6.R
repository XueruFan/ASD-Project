# 本代码用来通过GAMM建模分析centile和认知的相关
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
# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster/Corr")
statDir <- file.path(abideDir, "Analysis/Statistic/GmmCluster")
resDate <- "240315"
newDate <- "240610"

# # define variables
# ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]
centile <- read.csv(file.path(phenoDir, paste0("abide_A_centile_", resDate, ".csv")))
All <- merge(cluster, centile, by = "participant", All.x = TRUE)

# 修改各自站点名字的缩写
All$Site <- gsub("ABIDEII-NYU_1", "NYU", All$Site)
All$Site <- gsub("ABIDEII-NYU_2", "NYU", All$Site)
All$Site <- gsub("ABIDEII-KKI_1", "KKI", All$Site)
All$Site <- gsub("ABIDEII-SDSU_1", "SDSU", All$Site)
All$Site <- gsub("ABIDEII-UCLA_1", "UCLA", All$Site)
All$Site <- gsub("UCLA_1", "UCLA", All$Site)
All$Site <- gsub("UCLA_2", "UCLA", All$Site)
All$Site <- gsub("ABIDEII-GU_1", "GU", All$Site)
All$Site <- gsub("ABIDEII-UCD_1", "UCD", All$Site)
All$Site <- gsub("ABIDEII-EMC_1", "EMC", All$Site)
All$Site <- gsub("TRINITY", "TCD", All$Site)
All$Site <- gsub("ABIDEII-TCD_1", "TCD", All$Site)
All$Site <- gsub("ABIDEII-USM_1", "USM", All$Site)
All$Site <- gsub("ABIDEII-IU_1", "IU", All$Site)
All$Site <- gsub("ABIDEII-U_MIA_1", "UMIA", All$Site)
All$Site <- gsub("ABIDEII-ETH_1", "ETH", All$Site)
All$Site <- gsub("UM_1", "UM", All$Site)
All$Site <- gsub("UM_2", "UM", All$Site)
All$Site <- gsub("ABIDEII-OHSU_1", "OHSU", All$Site)
All$Site <- gsub("STANFORD", "SU1", All$Site)
All$Site <- gsub("ABIDEII-SU_2", "SU2", All$Site)
All$Site <- gsub("LEUVEN_2", "KUL", All$Site)
All$Site <- gsub("CALTECH", "CALT", All$Site)


# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"

names_cog <- c("ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_RRB")
# names_cog <- c("FIQ", "VIQ", "PIQ", "ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT",
#                "ADOS_2_RRB", "SRS_AWARENESS_RAW", "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW",
#                "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW", "SRS_AWARENESS_T", "SRS_COGNITION_T",
#                "SRS_COMMUNICATION_T", "SRS_MOTIVATION_T", "SRS_MANNERISMS_T", "SRS_TOTAL_T",
#                "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
#                "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C", "VINELAND_ABC_Standard",
#                "VINELAND_COMMUNICATION_STANDARD", "VINELAND_DAILYLIVING_STANDARD",
#                "VINELAND_SOCIAL_STANDARD", "BMI")

cog <- pheno[, c("participant", names_cog)]

All <- merge(All, cog, by = "participant", All.x = TRUE)


# volumeNames <- names(centile)[c(which(names(centile) == "bankssts"):which(names(centile) == "insula"))]
volumeNames <- c("middletemporal", "insula")
# volumeNames <- c("insula")
# [1] "bankssts"                 "caudalanteriorcingulate"  "caudalmiddlefrontal"     
# [4] "cuneus"                   "entorhinal"               "fusiform"                
# [7] "inferiorparietal"         "inferiortemporal"         "isthmuscingulate"        
# [10] "lateraloccipital"         "lateralorbitofrontal"     "lingual"                 
# [13] "medialorbitofrontal"      "middletemporal"           "parahippocampal"         
# [16] "paracentral"              "parsopercularis"          "parsorbitalis"           
# [19] "parstriangularis"         "pericalcarine"            "postcentral"             
# [22] "posteriorcingulate"       "precentral"               "precuneus"               
# [25] "rostralanteriorcingulate" "rostralmiddlefrontal"     "superiorfrontal"         
# [28] "superiorparietal"         "superiortemporal"         "supramarginal"           
# [31] "frontalpole"              "temporalpole"             "transversetemporal"      
# [34] "insula"   

# 创建一个空的数据框来保存模型统计值
model_stats <- data.frame(VolumeName = character(), CogName = character(), ClusterID = integer(),
                          R_squared = numeric(), Estimate = numeric(), SE = numeric(),
                          t_value = numeric(), P_value_t = character(),
                          edf = numeric(), red_edf = numeric(),
                          F_value = numeric(), P_value_smooth = character())

for (volumeName in volumeNames){
  # volumeName <- "bankssts"
  for (cogname in names_cog) {
    # cogname <- "FIQ"
    
  studyFIT <- All[, c("clusterID", volumeName, "TCV", cogname)]
  colnames(studyFIT)[2] <- "x"
  colnames(studyFIT)[4] <- "y"
    
    for (cluster in unique(studyFIT$clusterID)) {
      # cluster <- 1
      
      # 筛选当前cluster的数据
      cluster_data <- filter(studyFIT, clusterID == cluster)
      cluster_data[cluster_data < 0] <- NA
      cluster_data <- na.omit(cluster_data) 
      # 拟合GAMM模型
      gamm_model <- gamm4(y ~ s(x, bs = "tp", k = 4) + TCV, data = cluster_data)
      
      # 提取模型统计值
      model_summary <- summary(gamm_model$gam)
      
      # 提取模型统计值
      r_squared <- model_summary$r.sq
      estimate <- model_summary$p.table[1, "Estimate"]
      se <- model_summary$p.table[1, "Std. Error"]
      t_value <- model_summary$p.table[1, "t value"]
      p_value_t <- model_summary$p.table[1, "Pr(>|t|)"]
      
      edf <- model_summary$s.table[1, "edf"]
      red_edf <- model_summary$s.table[1, "Ref.df"]
      f_value <- model_summary$s.table[1, "F"]
      p_value_smooth <- model_summary$s.table[1, "p-value"]
      
      # 保存统计值
      model_stats <- rbind(model_stats, 
                           data.frame(VolumeName = volumeName, CogName = cogname, ClusterID = cluster, 
                                      R_squared = r_squared, Estimate = estimate, SE = se,
                                      t_value = t_value, P_value_t = p_value_t,
                                      edf = edf, red_edf = red_edf,
                                      F_value = f_value, P_value_smooth = p_value_smooth))
    }
  }
}

# 对 p_value 进行多重比较校正
# model_stats$P_value_t <- p.adjust(model_stats$P_value_t, method = "BH")
# model_stats$P_value_smooth <- p.adjust(model_stats$P_value_smooth, method = "BH")

model_stats<- model_stats %>%
  arrange(P_value_smooth)


# 对除前两列（"VolumeName" 和 "ClusterID"）以外的数值型列保留4位小数
model_stats[, -c(1, 2, 3)] <- lapply(model_stats[, -c(1, 2, 3)], function(x) {
  if(is.numeric(x)) {
    sprintf("%.4f", x)  # 保留4位小数并转换为字符格式
  } else {
    x
  }
})

# 设置列名
colnames(model_stats)[4:12] <- c("R²", "估计值", "标准误", "t值", "p值", "有效自由度",
                                 "剩余自由度", "F值", "平滑项p值")

# 保存模型统计值到 Excel 文件
write.xlsx(model_stats, file.path(statDir, paste0("asd_male_dev_GC_corr_GAMM_", newDate,
                                                  ".xlsx")), rowNames = F, colNames = T)


############# 画图
# 筛选出平滑项p值小于0.05的VolumeName
filtered_volumes <- model_stats %>%
  dplyr::filter(model_stats[[12]] < 0.05) %>%
  dplyr::pull(VolumeName)  # 将结果转换为向量

volumeNames <- unique(filtered_volumes)

for (volumeName in volumeNames) {
  for (cogname in names_cog) {
  # volumeName <- "insula"
  # cogname <- "ADOS_2_SEVERITY_TOTAL"

  studyFIT <- All[, c("clusterID", volumeName, "TCV", cogname)]
  studyFIT[studyFIT < 0] <- NA
  studyFIT <- na.omit(studyFIT)
  colnames(studyFIT)[2] <- "x"
  colnames(studyFIT)[4] <- "y"

  # predictions_site <- list()
  predictions_no_site <- list()

  # 定义颜色（黄绿色和橙色）对应每个 clusterID
  cluster_colors <- c("#719988", "#faa264")


  for (cluster in 1:2) {
    # cluster <- 2

    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)

    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(x, bs = "tp", k = 4) + TCV, data = cluster_data)

    # 准备预测数据（忽略站点）
    pred_data_no_site <- data.frame(x = seq(0, 1, length.out = 100))
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
    # pred_data_no_site$ci_lower <- pmax(pred_data_no_site$ci_lower, 0)
    # pred_data_no_site$ci_upper <- pmin(pred_data_no_site$ci_upper, 1)

    pred_data_no_site$clusterID <- cluster  # 添加 clusterID 标识

    # 将预测结果添加到列表中
    predictions_no_site[[cluster]] <- pred_data_no_site
  }

  # 合并忽略站点效应的预测结果
  all_predictions_no_site <- bind_rows(predictions_no_site)

  # 创建 ggplot2 图形，先添加透明的彩色线条来表示站点效应
  ggplot() +
    geom_point(data = studyFIT, aes(x = x, y = y, color = factor(clusterID)),
               alpha = 0.2, size = 2, shape = 16) + # 透明原始数据点
    # geom_line(data = all_predictions_site, aes(x = x, y = y_pred,
    #                                            group = interaction(Site, clusterID),
    #                                            color = factor(clusterID)),
    #           lwd = 0.5, alpha = 0.2) +  # 彩色透明线条表示站点效应
    # 添加彩色的忽略站点效应的线条和置信区间
    geom_ribbon(data = all_predictions_no_site, aes(x = x, ymin = ci_lower, ymax = ci_upper,
                                                    fill = factor(clusterID)),
                alpha = 0.2, linetype = 0) +  # 绘制彩色置信区间带
    geom_line(data = all_predictions_no_site,
              aes(x = x, y = y_pred, color = factor(clusterID)), lwd = 3) +
    scale_color_manual(values = cluster_colors) +  # 自定义颜色（黄绿色和橙色）
    scale_fill_manual(values = cluster_colors) +   # 填充颜色与线条颜色一致
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.position = "None")  # 不显示图例


  # 保存图像
  name <- paste0("GC_Corr_GAMM_", volumeName, "_", cogname, "_", newDate, ".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 20, height = 20, unit = "cm")
  }
}
