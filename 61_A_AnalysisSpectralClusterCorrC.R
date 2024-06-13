# 本代码用来分析由谱聚类的方法分类的两组ASD男性的变量之间的相关系数和显著性水平
# 雪如 2024年2月28日于北师大办公室

rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster")
statiDir <- file.path(abideDir, "Analysis/Statistic")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster/Corr")
resDate <- "240315"
newDate <- "240610"

# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_All_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "GMV")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
# 根据脑指标对于人群聚类的贡献进行加权
rank <- read.csv(file.path(clustDir, paste0("abide_A_asd_male_dev_Spectral_Cluster_RM_Rank_",
                                            newDate, ".csv")))
All <- merge(cluster, pheno, by = "participant", All.x = TRUE)


# 控制变量：site、FIQ
# 自变量：全部脑指标的加权合成新指标
# 因变量：认知
# 用全部脑指标每一个分别作为自变量


# 选择自变量列
names_brain <- names(cluster)[3:ncol(cluster)]

names_cog <- c("ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT",
               "ADOS_2_RRB", "SRS_AWARENESS_RAW", "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW",
               "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW", "SRS_AWARENESS_T", "SRS_COGNITION_T",
               "SRS_COMMUNICATION_T", "SRS_MOTIVATION_T", "SRS_MANNERISMS_T", "SRS_TOTAL_T",
               "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
               "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C", "VINELAND_ABC_Standard",
               "VINELAND_COMMUNICATION_STANDARD", "VINELAND_DAILYLIVING_STANDARD",
               "VINELAND_SOCIAL_STANDARD", "BMI")



#### select: 控制变量、自变量、因变量
names_col <- c("clusterID", "SITE_ID", "FIQ", names_brain, names_cog)
All <- All[, names_col]
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

L <- subset(All, clusterID == "1")
L <- L[, -1]
H <- subset(All, clusterID == "2")
H <- H[, -1]

# 初始化空数据框来存储相关计算的结果
results_L <- data.frame()
results_H <- data.frame()

# 
# # 给rank中的脑指标按照贡献大小赋值
# rank$Feature <- paste0(rank$Feature, "_centile")
# rank$Rank <- 0 # 创建一个新的列来保存排名
# # 将正数和负数分开排名
# positive_ranks <- rank(rank$Importance[rank$Importance > 0], ties.method = "first")
# negative_ranks <- rank(-rank$Importance[rank$Importance < 0], ties.method = "first")
# # 将排名结果赋值回数据框中
# rank$Rank[rank$Importance > 0] <- positive_ranks
# rank$Rank[rank$Importance < 0] <- -negative_ranks
# rank$Importance <- rank$Rank
# rank <- rank[, -3]


# 
# importance <- rank %>% select(Feature, Importance)



################### L组

# 确保 SITE_ID 和 FIQ 是因子和数值类型
L$SITE_ID <- as.factor(L$SITE_ID)
L$FIQ <- as.numeric(L$FIQ)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(L) & name_cog %in% names(L)) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_L <- L[!is.na(y) & !is.na(x) & !is.na(L$FIQ), c(name_brain, "SITE_ID", "FIQ", name_cog)]
        temp_L$y <- y[!is.na(y) & !is.na(x) & !is.na(L$FIQ)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID 和 FIQ
            y_lm <- lm(y ~ SITE_ID + FIQ, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + FIQ")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 FIQ
            y_lm <- lm(y ~ FIQ, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ FIQ")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          }
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value
          ))
        }
      }
    }
  }
}

########### 给P值排序
L_sorted <- arrange(results_L, p_value)

### 保存下来结果
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_statis_CorrC_L_", newDate, ".csv")
write.csv(L_sorted, file.path(statiDir, name), row.names = F)

################### H组

# 确保 SITE_ID 和 FIQ 是因子和数值类型
H$SITE_ID <- as.factor(H$SITE_ID)
H$FIQ <- as.numeric(H$FIQ)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$FIQ), c(name_brain, "SITE_ID", "FIQ", name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$FIQ)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID 和 FIQ
            y_lm <- lm(y ~ SITE_ID + FIQ, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + FIQ")), data = temp_H)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 FIQ
            y_lm <- lm(y ~ FIQ, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ FIQ")), data = temp_H)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          }
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value
          ))
        }
      }
    }
  }
}

########### 给P值排序
H_sorted <- arrange(results_H, p_value)

### 保存下来结果
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_statis_CorrC_H_", newDate, ".csv")
write.csv(H_sorted, file.path(statiDir, name), row.names = F)
