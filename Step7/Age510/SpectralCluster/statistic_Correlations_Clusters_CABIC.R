# Calculate cluster-wise correlation with their statistical significance for all variable pairs
# Male, ASD, 5~9.9 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/510/Corr")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_510_", resDate, ".csv")
cluster <- read.csv(file.path(cabicDir, "result/pred/510", name))

cabic_all <- merge(cluster, pheno, by = "participant")

names_brain <- names(cluster)[11:(ncol(cluster)-1)]

names_cog_p <- c("FIQ", "ADOS_SA", "ADOS_TOTAL", "ADIR_SOCI", "ADIR_RRB",
                 "SRS_SA", "SRS_SCOG", "SRS_SCOM", "SRS_SM", "SRS_AM", "SRS_TOTAL")
names_cog_s <- c("ADOS_RRB")


#### select: 控制变量、自变量、因变量
names_col <- c("predicted_cluster", "SITE", "TCV", "Age", names_brain, names_cog_p, names_cog_s)


temp <- cabic_all[, names_col]
temp[temp < 0] <- NA
colnames(temp)[1] <- "clusterID"

temp <- subset(temp, Age < 13)

L <- subset(temp, clusterID == "1")
L <- L[, c(-1,-4)]
H <- subset(temp, clusterID == "2")
H <- H[, c(-1,-4)]

# 初始化空数据框来存储相关计算的结果
results_L <- data.frame()
results_H <- data.frame()


################### L组

# 确保 SITE 和 FIQ 是因子和数值类型
L$SITE <- as.factor(L$SITE)


###################### 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_s) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_L <- L[!is.na(y), c(name_brain, "SITE", "TCV", name_cog)]
        temp_L$y <- y[!is.na(y)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV")), data = temp_L)
            
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
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 TCV
            y_lm <- lm(y ~ TCV, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV")), data = temp_L)
            
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
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值  # 保存 R² 值
          ))
        }
      }
    }
}
###################### 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Person 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_p) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_L <- L[!is.na(y), c(name_brain, "SITE", "TCV", name_cog)]
        temp_L$y <- y[!is.na(y)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          } else {
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 TCV
            y_lm <- lm(y ~ TCV, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          }
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值  # 保存 R² 值
          ))
        }
      }
  }
}

results_L$P_adj <- p.adjust(results_L$p_value, method = "fdr")
########### 给P值排序
L_sorted <- arrange(results_L, p_value)


################### H组

# 确保 SITE 和 FIQ 是因子和数值类型
H$SITE <- as.factor(H$SITE)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog_s) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$TCV),
                    c(name_brain, "SITE", "TCV", name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$TCV)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV")), data = temp_H)
            
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
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 FIQ
            y_lm <- lm(y ~ TCV, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV")), data = temp_H)
            
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
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值  # 保存 R² 值
          ))
        }
      }
    }
  }
}
for (name_brain in names_brain) {
  for (name_cog in names_cog_p) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$TCV),
                    c(name_brain, "SITE", "TCV", name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$TCV)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV")), data = temp_H)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          } else {
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 FIQ
            y_lm <- lm(y ~ TCV, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV")), data = temp_H)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "pearson")
            # 提取 R² 值
            R2_value <- summary(y_lm)$r.squared
            # 提取自由度
            df_value <- df.residual(y_lm)  # 提取自由度
          }
          
          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value,
            R2 = R2_value,  # 保存 R² 值
            df = df_value  # 保存自由度  # 保存 R² 值  # 保存 R² 值
          ))
        }
      }
    }
  }
}
results_H$P_adj <- p.adjust(results_H$p_value, method = "fdr")
########### 给P值排序
H_sorted <- arrange(results_H, p_value)

L_sorted <- L_sorted[L_sorted$p_value < 0.05, ] # 删除数据框中P列大于或等于0.05的行
H_sorted <- H_sorted[H_sorted$p_value < 0.05, ]

name <- paste0("corr_Part4_Select_H_", resDate, ".csv")
write.csv(H_sorted, file.path(resuDir, name), row.names = F)
name <- paste0("corr_Part4_Select_L_", resDate, ".csv")
write.csv(L_sorted, file.path(resuDir, name), row.names = F)

L_sorted$cluster <- "L"
H_sorted$cluster <- "H"

selected <- rbind(L_sorted, H_sorted)
name <- paste0("corr_Part4_Select_LH_", resDate, ".csv")
write.csv(selected, file.path(resuDir, name), row.names = F)

