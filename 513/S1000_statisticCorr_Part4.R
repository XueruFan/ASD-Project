# 本代码用来分析两组ASD发育男性的变量之间的相关系数和显著性水平
# 使用斯皮尔曼相关，脑指标（centile）分别作为自变量
# 雪如 2024年2月28日于北师大办公室
##################################
# Site是控制变量，TCV作为协变量控制掉
# 保存原始的相关系数和p值csv文件，另外，筛选p小于0.05的结果，保存csv文件
# 按照显著性水平结果，绘制相关图png

##################################

rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
# plotDir <- file.path(abideDir, "Plot/Cluster/Spect513/Corr/Part4")
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

# 选择自变量列
names_brain <- names(cluster)[10:ncol(cluster)]
# names_brain_l <- c("middletemporal_centile", "pericalcarine_centile")
# names_brain_h <- c("caudalanteriorcingulate_centile", "insula_centile")
# names_brain <- c("isthmuscingulate_centile", "entorhinal_centile", "precuneus_centile",
#                  "middletemporal_centile", "transversetemporal_centile", "postcentral_centile",
#                  "rostralmiddlefrontal_centile", "fusiform_centile", "inferiortemporal_centile",
#                  "rostralanteriorcingulate_centile")


names_cog_p <- c("FIQ", "ADOS_2_SOCAFFECT" , "ADOS_2_TOTAL", "ADI_R_SOCIAL_TOTAL_A", 
                 "ADI_R_RRB_TOTAL_C", "SRS_TOTAL_RAW", "SRS_COGNITION_RAW", "SRS_AWARENESS_RAW", 
                 "SRS_COMMUNICATION_RAW", "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW")
names_cog_s <- c("ADOS_2_RRB")

#### select: 控制变量、自变量、因变量
names_col <- c("clusterID", "SITE_ID", "TCV_centile", names_brain, names_cog_p, names_cog_s)
temp <- All[, names_col]
temp[temp < 0] <- NA

L <- subset(temp, clusterID == "1")
L <- L[, -1]
H <- subset(temp, clusterID == "2")
H <- H[, -1]

# 初始化空数据框来存储相关计算的结果
results_L <- data.frame()
results_H <- data.frame()



################### L组

# 确保 SITE_ID 和 FIQ 是因子和数值类型
L$SITE_ID <- as.factor(L$SITE_ID)

# 初始化空数据框来存储结果
results_L <- data.frame()

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_s) {
    if (name_brain %in% names(L) & name_cog %in% names(L)) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        temp_L <- L[!is.na(y) & !is.na(x) & !is.na(L$TCV_centile),
                    c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        temp_L$y <- y[!is.na(y) & !is.na(x) & !is.na(L$TCV_centile)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID 和 TCV
            y_lm <- lm(y ~ SITE_ID + TCV_centile, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + TCV_centile")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
            
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 TCV
            y_lm <- lm(y ~ TCV_centile, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV_centile")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
            
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          }
          
          # 保存 Spearman 相关系数、p 值和 R² 到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值
          ))
        }
      }
    }
  }
}

###################### 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Pearson 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_p) {
    if (name_brain %in% names(L) & name_cog %in% names(L)) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        temp_L <- L[!is.na(y) & !is.na(x) & !is.na(L$TCV_centile),
                    c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        temp_L$y <- y[!is.na(y) & !is.na(x) & !is.na(L$TCV_centile)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID 和 TCV
            y_lm <- lm(y ~ SITE_ID + TCV_centile, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + TCV_centile")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Pearson 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 TCV
            y_lm <- lm(y ~ TCV_centile, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV_centile")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Pearson 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          }
          
          # 保存 Pearson 相关系数、p 值和 R² 到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值
          ))
        }
      }
    }
  }
}

# 调整 p 值，使用 FDR 方法
results_L$P_adj <- p.adjust(results_L$p_value, method = "fdr")
########### 给P值排序
L_sorted <- arrange(results_L, p_value)

################### H组

# 确保 SITE_ID 和 FIQ 是因子和数值类型
H$SITE_ID <- as.factor(H$SITE_ID)

# 初始化空数据框来存储结果
results_H <- data.frame()

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_s) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$TCV_centile),
                    c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$TCV_centile)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE_ID)) > 1) {
            y_lm <- lm(y ~ SITE_ID + TCV_centile, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + TCV_centile")), data = temp_H)
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          } else {
            y_lm <- lm(y ~ TCV_centile, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV_centile")), data = temp_H)
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          }
          
          # 保存结果
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值
          ))
        }
      }
    }
  }
}

# 继续计算 Pearson 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_p) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$TCV_centile),
                    c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$TCV_centile)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            y_lm <- lm(y ~ SITE_ID + TCV_centile, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID + TCV_centile")), data = temp_H)
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          } else {
            y_lm <- lm(y ~ TCV_centile, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV_centile")), data = temp_H)
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            R2_value <- summary(y_lm)$r.squared
            df_value <- df.residual(y_lm)
          }
          
          # 保存结果
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值
          ))
        }
      }
    }
  }
}

# 调整 p 值，使用 FDR 方法
results_H$P_adj <- p.adjust(results_H$p_value, method = "fdr")
H_sorted <- arrange(results_H, p_value)


L_sorted <- L_sorted[L_sorted$p_value < 0.05, ] # 删除P值大于等于0.05的行
H_sorted <- H_sorted[H_sorted$p_value < 0.05, ]

### 保存结果
name <- paste0("corr_Part4_Select_H_", newDate, ".csv")
write.csv(H_sorted, file.path(statiDir, name), row.names = F)
name <- paste0("corr_Part4_Select_L_", newDate, ".csv")
write.csv(L_sorted, file.path(statiDir, name), row.names = F)

# 合并L组和H组的选择结果
L_selected <- L_sorted
H_selected <- H_sorted
L_selected$cluster <- "L"
H_selected$cluster <- "H"
selected <- rbind(L_selected, H_selected)

name <- paste0("corr_Part4_Select_LH_", newDate, ".csv")
write.csv(selected, file.path(statiDir, name), row.names = F)
