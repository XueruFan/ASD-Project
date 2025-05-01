# Correlation analysis on structural covariance and cognitive behaviors.
# Male, ASD, <13 years old, GMM Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list = ls())

# 加载必要的包
packages <- c("tidyverse", "stats", "dplyr", "ggplot2", "reshape2", "pheatmap",
              "magrittr", "readr", "openxlsx")
sapply(packages, require, character.only = TRUE)

# 定义路径
abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/GMM513")
statDir <- file.path(abideDir, "Analysis/Statistic/GMM513")
resDate <- "240315"
newDate <- "240610"

# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "GMV")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
All <- merge(cluster, pheno, by = "participant", All.x = TRUE)
All[All < 0] <- NA

# 修改各自站点名字的缩写
All$SITE_ID <- gsub("ABIDEII-NYU_1", "NYU", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-NYU_2", "NYU", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-KKI_1", "KKI", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-SDSU_1", "SDSU", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-UCLA_1", "UCLA", All$SITE_ID)
All$SITE_ID <- gsub("UCLA_1", "UCLA", All$SITE_ID)
All$SITE_ID <- gsub("UCLA_2", "UCLA", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-GU_1", "GU", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-UCD_1", "UCD", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-EMC_1", "EMC", All$SITE_ID)
All$SITE_ID <- gsub("TRINITY", "TCD", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-TCD_1", "TCD", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-USM_1", "USM", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-IU_1", "IU", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-U_MIA_1", "UMIA", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-ETH_1", "ETH", All$SITE_ID)
All$SITE_ID <- gsub("UM_1", "UM", All$SITE_ID)
All$SITE_ID <- gsub("UM_2", "UM", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-OHSU_1", "OHSU", All$SITE_ID)
All$SITE_ID <- gsub("STANFORD", "SU1", All$SITE_ID)
All$SITE_ID <- gsub("ABIDEII-SU_2", "SU2", All$SITE_ID)
All$SITE_ID <- gsub("LEUVEN_2", "KUL", All$SITE_ID)
All$SITE_ID <- gsub("CALTECH", "CALT", All$SITE_ID)


names_cog_p <- c("FIQ", "ADOS_2_SOCAFFECT" , "ADOS_2_TOTAL", "ADI_R_SOCIAL_TOTAL_A", 
                 "ADI_R_RRB_TOTAL_C", "SRS_TOTAL_RAW", "SRS_COGNITION_RAW", "SRS_AWARENESS_RAW", 
                 "SRS_COMMUNICATION_RAW", "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW")
names_cog_s <- c("ADOS_2_RRB")

cor_brain <- read.csv(file.path(statDir, "StruCova_Both.csv"))
cor_brain$ROI1 <- paste0(cor_brain$ROI1, "_centile")
cor_brain$ROI2 <- paste0(cor_brain$ROI2, "_centile")

L <- subset(All, clusterID == "1")
H <- subset(All, clusterID == "2")


################################## Part 1: Group L #################################################

# 初始化结果数据框
results_L <- data.frame(
  ROI1 = character(),
  ROI2 = character(),
  name_cog = character(),
  r = numeric(),
  p = numeric(),
  df = numeric(),
  t = numeric(),
  stringsAsFactors = FALSE
)

# 遍历显著脑区对和认知变量
for (i in 2:nrow(cor_brain)) {
  
  roi1 <- cor_brain$ROI1[i]
  roi2 <- cor_brain$ROI2[i]
  
  # 遍历认知变量 names_cog_p 进行相关分析
  for (name_cog in names_cog_p) {
  # name_cog <- names_cog_p[2]
  temp_L <- L[!is.na(L[[name_cog]]), c(roi1, roi2, "SITE_ID", "TCV_centile", name_cog)]
  
  
  # 对两个脑区分别进行回归，控制 Site 和 TCV
  x1_lm <- lm(as.formula(paste(roi1, "~ SITE_ID + TCV_centile")), data = temp_L)
  x2_lm <- lm(as.formula(paste(roi2, "~ SITE_ID + TCV_centile")), data = temp_L)
  
  # 提取残差
  temp_L$roi1_residuals <- residuals(x1_lm)
  temp_L$roi2_residuals <- residuals(x2_lm)
  
  # 计算每个被试的 partial_p
  temp_L$partial_p <- temp_L$roi1_residuals * temp_L$roi2_residuals
  
    
  # 确保非 NA 数据不少于 30 个
  if (nrow(temp_L) >= 30) {
    
    y_lm <- lm(as.formula(paste(name_cog, "~ SITE_ID + TCV_centile")), data = temp_L)
    temp_L$y_residuals <- residuals(y_lm)
    

      
      cor_test <- cor.test(temp_L$partial_p, temp_L$y_residuals, method = "pearson")
        
        # 保存结果
        results_L <- rbind(results_L, data.frame(
          ROI1 = roi1,
          ROI2 = roi2,
          name_cog = name_cog,
          r = cor_test$estimate,
          p = cor_test$p.value,
          df = df.residual(y_lm),  # 提取自由度,
          t = cor_test$statistic
        ))
    }
  }
  
  for (name_cog in names_cog_s) {
    # name_cog <- names_cog_s[1]
    temp_L <- L[!is.na(L[[roi1]]) & !is.na(L[[roi2]]) & !is.na(L[[name_cog]]),
                c(roi1, roi2, "SITE_ID", "TCV_centile", name_cog)]
    
    
    # 对两个脑区分别进行回归，控制 Site 和 TCV
    x1_lm <- lm(as.formula(paste(roi1, "~ SITE_ID + TCV_centile")), data = temp_L)
    x2_lm <- lm(as.formula(paste(roi2, "~ SITE_ID + TCV_centile")), data = temp_L)
    
    # 提取残差
    temp_L$roi1_residuals <- residuals(x1_lm)
    temp_L$roi2_residuals <- residuals(x2_lm)
    
    # 计算每个被试的 partial_p
    temp_L$partial_p <- temp_L$roi1_residuals * temp_L$roi2_residuals
    
    
    # 确保非 NA 数据不少于 30 个
    if (nrow(temp_L) >= 30) { 
      
      y_lm <- lm(as.formula(paste(name_cog, "~ SITE_ID + TCV_centile")), data = temp_L)
      temp_L$y_residuals <- residuals(y_lm)
      
        
        cor_test <- cor.test(temp_L$partial_p, temp_L$y_residuals, method = "spearman")
        
        # 保存结果
        results_L <- rbind(results_L, data.frame(
          ROI1 = roi1,
          ROI2 = roi2,
          name_cog = name_cog,
          r = cor_test$estimate,
          p = cor_test$p.value,
          df = df.residual(y_lm),  # 提取自由度
          t = cor_test$statistic
        ))
    }
  }
  
}
  



################################## Part 2: Group H #################################################

# 初始化结果数据框
results_H <- data.frame(
  ROI1 = character(),
  ROI2 = character(),
  name_cog = character(),
  r = numeric(),
  p = numeric(),
  df = numeric(),
  t = numeric(),
  stringsAsFactors = FALSE
)

# 遍历显著脑区对和认知变量
for (i in c(1)) {
  
  roi1 <- cor_brain$ROI1[i]
  roi2 <- cor_brain$ROI2[i]
  
  # 遍历认知变量 names_cog_p 进行相关分析
  for (name_cog in names_cog_p) {
    # name_cog <- names_cog_p[2]
    temp_H <- H[!is.na(H[[name_cog]]), c(roi1, roi2, "SITE_ID", "TCV_centile", name_cog)]
    
    
    # 对两个脑区分别进行回归，控制 Site 和 TCV
    x1_lm <- lm(as.formula(paste(roi1, "~ SITE_ID + TCV_centile")), data = temp_H)
    x2_lm <- lm(as.formula(paste(roi2, "~ SITE_ID + TCV_centile")), data = temp_H)
    
    # 提取残差
    temp_H$roi1_residuals <- residuals(x1_lm)
    temp_H$roi2_residuals <- residuals(x2_lm)
    
    # 计算每个被试的 partial_p
    temp_H$partial_p <- temp_H$roi1_residuals * temp_H$roi2_residuals
    
    
    # 确保非 NA 数据不少于 30 个
    if (nrow(temp_H) >= 30) {
      
      y_lm <- lm(as.formula(paste(name_cog, "~ SITE_ID + TCV_centile")), data = temp_H)
      temp_H$y_residuals <- residuals(y_lm)
      
      
      cor_test <- cor.test(temp_H$partial_p, temp_H$y_residuals, method = "pearson")
      
      # 保存结果
      results_H <- rbind(results_H, data.frame(
        ROI1 = roi1,
        ROI2 = roi2,
        name_cog = name_cog,
        r = cor_test$estimate,
        p = cor_test$p.value,
        df = df.residual(y_lm),  # 提取自由度,
        t = cor_test$statistic
      ))
    }
  }
  
  for (name_cog in names_cog_s) {
    # name_cog <- names_cog_s[1]
    temp_H <- H[!is.na(H[[name_cog]]), c(roi1, roi2, "SITE_ID", "TCV_centile", name_cog)]
    
    
    # 对两个脑区分别进行回归，控制 Site 和 TCV
    x1_lm <- lm(as.formula(paste(roi1, "~ SITE_ID + TCV_centile")), data = temp_H)
    x2_lm <- lm(as.formula(paste(roi2, "~ SITE_ID + TCV_centile")), data = temp_H)
    
    # 提取残差
    temp_H$roi1_residuals <- residuals(x1_lm)
    temp_H$roi2_residuals <- residuals(x2_lm)
    
    # 计算每个被试的 partial_p
    temp_H$partial_p <- temp_H$roi1_residuals * temp_H$roi2_residuals
    
    
    # 确保非 NA 数据不少于 30 个
    if (nrow(temp_H) >= 30) { 
      
      y_lm <- lm(as.formula(paste(name_cog, "~ SITE_ID + TCV_centile")), data = temp_H)
      temp_H$y_residuals <- residuals(y_lm)
      
      
      cor_test <- cor.test(temp_H$partial_p, temp_H$y_residuals, method = "spearman")
      
      # 保存结果
      results_H <- rbind(results_H, data.frame(
        ROI1 = roi1,
        ROI2 = roi2,
        name_cog = name_cog,
        r = cor_test$estimate,
        p = cor_test$p.value,
        df = df.residual(y_lm),  # 提取自由度
        t = cor_test$statistic
      ))
    }
  }
  
}
  

### 保存结果
name <- paste0("corr_strucova_H_", newDate, ".csv")
write.csv(results_H, file.path(statDir, name), row.names = F)
name <- paste0("corr_strucova_L_", newDate, ".csv")
write.csv(results_L, file.path(statDir, name), row.names = F)

