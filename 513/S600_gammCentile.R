# 本代码用来建模两个聚类人群脑指标centile的发育轨线
# GAMM模型
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
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect513/GAMM")
statDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
resDate <- "240315"
newDate <- "240610"


name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]
centile <- read.csv(file.path(phenoDir, paste0("abide_All_centile_", resDate, ".csv")))
All <- merge(cluster, centile, by = "participant", All.x = TRUE)


volumeNames <- names(centile)[c(which(names(centile) == "bankssts"):which(names(centile) == "insula"))]

# 创建一个空的数据框来保存模型统计值
model_stats <- data.frame(VolumeName = character(), ClusterID = integer(),
                          R_squared = numeric(), Estimate = numeric(), SE = numeric(),
                          t_value = numeric(), P_value_t = character(),
                          edf = numeric(), red_edf = numeric(),
                          F_value = numeric(), P_value_smooth = character())

for (volumeName in volumeNames){
  # volumeName <- "bankssts"
  
  studyFIT <- All[, c("clusterID", "Age", "TCV", volumeName)]
  colnames(studyFIT)[4] <- "y"
  
  for (cluster in unique(studyFIT$clusterID)) {
    # cluster <- 1
    
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)

    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(Age, bs = "tp", k = 4) + TCV, data = cluster_data)
    
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
                         data.frame(VolumeName = volumeName, ClusterID = cluster, 
                                    R_squared = r_squared, Estimate = estimate, SE = se,
                                    t_value = t_value, P_value_t = p_value_t,
                                    edf = edf, red_edf = red_edf,
                                    F_value = f_value, P_value_smooth = p_value_smooth))
  }
}

# 对 p_value 进行校正
model_stats$P_adj <- p.adjust(model_stats$P_value_smooth, method = "fdr")

model_stats<- model_stats %>%
  arrange(P_adj)


# 对除前两列（"VolumeName" 和 "ClusterID"）以外的数值型列保留4位小数
model_stats[, -c(1, 2)] <- lapply(model_stats[, -c(1, 2)], function(x) {
  if(is.numeric(x)) {
    sprintf("%.4f", x)  # 保留4位小数并转换为字符格式
  } else {
    x
  }
})

# 设置列名
colnames(model_stats)[3:11] <- c("R²", "估计值", "标准误", "t值", "p值", "有效自由度",
                                 "剩余自由度", "F值", "平滑项p值")

# 保存模型统计值到 Excel 文件
write.xlsx(model_stats, file.path(statDir, paste0("gamm_Centile_", newDate,
                                                  ".xlsx")), rowNames = F, colNames = T)
