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
clustDir <- file.path(abideDir, "Analysis/Cluster/SpectralCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/SpectralCluster/Centile/GAMM")
statDir <- file.path(abideDir, "Analysis/Statistic")
resDate <- "240315"
newDate <- "240610"

# define variables
ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("asd_male_Spectral_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
cluster <- cluster[, c("clusterID", "participant")]
centile <- read.csv(file.path(phenoDir, paste0("abide_A_centile_", resDate, ".csv")))
All <- merge(cluster, centile, by = "participant", All.x = TRUE)

# # 认知行为
# pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
# colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# FIQ <- pheno[, c("participant", "FIQ")]
# 
# All <- merge(All, FIQ, by = "participant", All.x = TRUE)

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


volumeNames <- names(centile)[c(which(names(centile) == "GMV"):which(names(centile) == "insula"))]


################################### Part 1: 不考虑站点效应 #########################################
# 创建一个空的数据框来保存模型统计值
model_stats <- data.frame(VolumeName = character(), ClusterID = integer(),
                          R_squared = numeric(), Estimate = numeric(), SE = numeric(),
                          t_value = numeric(), P_value_t = character(),
                          edf = numeric(), red_edf = numeric(),
                          F_value = numeric(), P_value_smooth = character())

for (volumeName in volumeNames){
  # volumeName <- "GMV"
  
  studyFIT <- All[, c("clusterID", "Age", volumeName)]
  colnames(studyFIT)[2:3] <- c("x", "y")
  
  for (cluster in unique(studyFIT$clusterID)) {
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)

    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(x, k = 4), data = cluster_data)
    
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

model_stats<- model_stats %>%
  arrange(P_value_smooth)

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
write.xlsx(model_stats, file.path(statDir, paste0("asd_male_dev_SC_statis_GAMM_NegSite_", newDate,
                                                  ".xlsx")), rowNames = F, colNames = T)


################################### Part 2: 考虑站点效应 ###########################################

# 创建一个空的数据框来保存模型统计值
model_stats <- data.frame(VolumeName = character(), ClusterID = integer(),
                          R_squared = numeric(), Estimate = numeric(), SE = numeric(),
                          t_value = numeric(), P_value_t = character(),
                          edf = numeric(), red_edf = numeric(),
                          F_value = numeric(), P_value_smooth = character())

for (volumeName in volumeNames){
  # volumeName <- "GMV"
  
  studyFIT <- All[, c("clusterID", "Age", volumeName, "Site")]
  colnames(studyFIT)[2:3] <- c("x", "y")
  
  for (cluster in unique(studyFIT$clusterID)) {
    # cluster <- 1
    # 筛选当前cluster的数据
    cluster_data <- filter(studyFIT, clusterID == cluster)
    cluster_data$Site <- factor(cluster_data$Site)
    
    # 拟合GAMM模型
    gamm_model <- gamm4(y ~ s(x, k = 4) + Site, data = cluster_data)
    
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

model_stats<- model_stats %>%
  arrange(P_value_smooth)

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
write.xlsx(model_stats, file.path(statDir, paste0("asd_male_dev_SC_statis_GAMM_ConSite_", newDate,
                                                  ".xlsx")), rowNames = F, colNames = T)
