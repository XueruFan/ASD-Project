# Calculate mediation effect in the correlations we detected
# Male, ASD, <13 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "mediation",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
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
all_site_levels <- levels(factor(All$SITE_ID))
# 选择自变量列
names_brain <- names(cluster)[10:ncol(cluster)]

names_cog <- c("FIQ", "ADOS_2_SOCAFFECT" , "ADOS_2_TOTAL", "ADOS_2_RRB","ADI_R_SOCIAL_TOTAL_A",  
               "ADI_R_RRB_TOTAL_C", "SRS_TOTAL_RAW", "SRS_COGNITION_RAW", "SRS_AWARENESS_RAW", 
               "SRS_COMMUNICATION_RAW", "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW")
names_med <- c("FIQ", names_brain)

#### select: 控制变量、自变量、因变量
names_col <- c("clusterID", "SITE_ID", "TCV_centile", names_brain, names_cog)

H <- All[, names_col]
H[H < 0] <- NA

H <- subset(H, clusterID == "2")
H <- H[, -1]

results_H <- data.frame()


# 脑-行为
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    for (name_med in names_med) {
      if (name_brain != name_med) {
        temp <- H[, c(name_brain, "SITE_ID", "TCV_centile", name_med, name_cog)]
        temp <- na.omit(temp)
        
        # 对变量进行标准化
        temp$x <- scale(temp[[name_brain]])
        temp$med <- scale(temp[[name_med]])
        temp$y <- scale(temp[[name_cog]])
        # temp$x <- temp[[name_brain]]
        # temp$med <- temp[[name_med]]
        # temp$y <- temp[[name_cog]]
        # 
        # 步骤1: 估计路径 A (X -> M)
        model_M1 <- lm(med ~ x + SITE_ID + TCV_centile, data = temp)
        
        # 提取路径 A 的系数和标准误
        A <- coef(model_M1)["x"]
        se_A <- summary(model_M1)$coefficients["x", "Std. Error"]
        p_A <- summary(model_M1)$coefficients["x", "Pr(>|t|)"]  # 提取 A 的 p 值
        
        # 步骤2: 估计路径 B 和路径 C' (M -> Y，控制 X)
        model_Y1 <- lm(y ~ x + med + SITE_ID + TCV_centile, data = temp)
        
        # 提取路径 B 和 C' 的系数和标准误
        B <- coef(model_Y1)["med"]
        se_B <- summary(model_Y1)$coefficients["med", "Std. Error"]
        p_B <- summary(model_Y1)$coefficients["med", "Pr(>|t|)"]  # 提取 B 的 p 值
        
        # 提取直接效应 C'
        direct_effect <- coef(model_Y1)["x"]
        p_C_prime <- summary(model_Y1)$coefficients["x", "Pr(>|t|)"]  # 提取直接效应的 p 值
        
        
        # 步骤3: 估计总效应 C (X -> Y)
        model_C <- lm(y ~ x + SITE_ID + TCV_centile, data = temp)
        
        # 提取总效应 C
        total_effect <- coef(model_C)["x"]
        
        # 计算间接效应 A * B
        indirect_effect <- A * B
        
        # Sobel Test 计算间接效应的显著性
        sobel_test_stat <- (A * B) / sqrt(B^2 * se_A^2 + A^2 * se_B^2)
        p_value <- 2 * (1 - pnorm(abs(sobel_test_stat)))
        
        results_H <- rbind(results_H, data.frame(
          brain = name_brain,
          cog = name_cog,
          med = name_med,
          A_effect = A,
          A_pvalue = p_A,
          B_effect = B,
          B_pvalue = p_B,
          DirectEffect = direct_effect,
          DirectEffect_pvalue = p_C_prime,
          TotalEffect = total_effect,
          IndirectEffect = indirect_effect,
          SobelTestStat = sobel_test_stat,
          Sobel_pvalue = p_value))
      }
    }
  }
}


name <- paste0("corr_median_Part4_Select_H_brain2cog_", newDate, ".csv")
write.csv(results_H, file.path(statiDir, name), row.names = F)



# 行为-脑
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    for (name_med in names_med) {
      if (name_cog != name_med) {
        temp <- H[, c(name_brain, "SITE_ID", "TCV_centile", name_med, name_cog)]
        temp <- na.omit(temp)
        
        # 对变量进行标准化
        temp$y <- scale(temp[[name_brain]])
        temp$med <- scale(temp[[name_med]])
        temp$x <- scale(temp[[name_cog]])

        # 步骤1: 估计路径 A (X -> M)
        model_M1 <- lm(med ~ x + SITE_ID + TCV_centile, data = temp)
        
        # 提取路径 A 的系数和标准误
        A <- coef(model_M1)["x"]
        se_A <- summary(model_M1)$coefficients["x", "Std. Error"]
        p_A <- summary(model_M1)$coefficients["x", "Pr(>|t|)"]  # 提取 A 的 p 值
        
        # 步骤2: 估计路径 B 和路径 C' (M -> Y，控制 X)
        model_Y1 <- lm(y ~ x + med + SITE_ID + TCV_centile, data = temp)
        
        # 提取路径 B 和 C' 的系数和标准误
        B <- coef(model_Y1)["med"]
        se_B <- summary(model_Y1)$coefficients["med", "Std. Error"]
        p_B <- summary(model_Y1)$coefficients["med", "Pr(>|t|)"]  # 提取 B 的 p 值
        
        # 提取直接效应 C'
        direct_effect <- coef(model_Y1)["x"]
        p_C_prime <- summary(model_Y1)$coefficients["x", "Pr(>|t|)"]  # 提取直接效应的 p 值
        
        
        # 步骤3: 估计总效应 C (X -> Y)
        model_C <- lm(y ~ x + SITE_ID + TCV_centile, data = temp)
        
        # 提取总效应 C
        total_effect <- coef(model_C)["x"]
        
        # 计算间接效应 A * B
        indirect_effect <- A * B
        
        # Sobel Test 计算间接效应的显著性
        sobel_test_stat <- (A * B) / sqrt(B^2 * se_A^2 + A^2 * se_B^2)
        p_value <- 2 * (1 - pnorm(abs(sobel_test_stat)))
        
        results_H <- rbind(results_H, data.frame(
          brain = name_brain,
          cog = name_cog,
          med = name_med,
          A_effect = A,
          A_pvalue = p_A,
          B_effect = B,
          B_pvalue = p_B,
          IndirectEffect = indirect_effect,
          DirectEffect = direct_effect,
          DirectEffect_pvalue = p_C_prime,
          TotalEffect = total_effect,
          SobelTestStat = sobel_test_stat,
          Sobel_pvalue = p_value))
      }
    }
  }
}


name <- paste0("corr_median_Part4_Select_H_cog2brain_", newDate, ".csv")
write.csv(results_H, file.path(statiDir, name), row.names = F)
