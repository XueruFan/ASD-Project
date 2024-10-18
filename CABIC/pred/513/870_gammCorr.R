rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/513/Corr")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_513_", resDate, ".csv")
cluster <- read.csv(file.path(cabicDir, "result/pred/513", name))

cabic_all <- merge(cluster, pheno, by = "participant")

names_brain <- names(cluster)[11:(ncol(cluster)-1)]
# names_brain <- c("parstriangularis", "postcentral", "entorhinal", "inferiortemporal", 
# "supramarginal", "insula", "rostralanteriorcingulate", "transversetemporal")

names_cog <- c("FIQ", "ADOS_SA", "ADOS_RRB","ADOS_TOTAL", "ADIR_SOCI", "ADIR_RRB",
                 "SRS_SA", "SRS_SCOG", "SRS_SCOM", "SRS_SM", "SRS_AM", "SRS_TOTAL")


#### select: 控制变量、自变量、因变量
names_col <- c("predicted_cluster", "SITE", "TCV", "Age", names_brain, names_cog)


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

for (name_brain in names_brain) {
  for (name_cog in names_cog) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y)) >= 30) {
        # temp_L <- L[!is.na(y), c(name_brain, "SITE", "AGE_AT_SCAN", "TCV", name_cog)]
        temp_L <- L[!is.na(y), c(name_brain, "SITE", "TCV", name_cog)]
        
        colnames(temp_L)[1] <- "x"
        colnames(temp_L)[4] <- "y"
        
        if (nrow(temp_L) >= 30) {
          
          # 拟合GAMM模型
          # gamm_model <- gamm4(y ~ s(x, k = 4) + s(AGE_AT_SCAN, k = 4) + TCV + SITE, data = temp_L)
          gamm_model <- gamm4(y ~ s(x, k = 5) + TCV + SITE, data = temp_L)
          
          
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
          results_L <- rbind(results_L, 
                               data.frame(brain = name_brain, cog = name_cog,
                                          R_squared = r_squared, Estimate = estimate, SE = se,
                                          t_value = t_value, P_value_t = p_value_t,
                                          edf = edf, red_edf = red_edf,
                                          F_value = f_value, P_value_smooth = p_value_smooth))
        }
      }
  }
}

# 对 p_value 进行校正
results_L$P_adj <- p.adjust(results_L$P_value_smooth, method = "fdr")

results_L<- results_L %>%
  arrange(P_adj)


H$SITE <- as.factor(H$SITE)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) >= 30)) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        # temp_H <- H[!is.na(y) & !is.na(x), c(name_brain, "SITE", "AGE_AT_SCAN", "TCV", name_cog)]
        temp_H <- H[!is.na(y), c(name_brain, "SITE", "TCV", name_cog)]
        
        colnames(temp_H)[1] <- "x"
        colnames(temp_H)[4] <- "y"
        
        if (nrow(temp_H) >= 30) {
          
          # 拟合GAMM模型
          # gamm_model <- gamm4(y ~ s(x, k = 4) + s(AGE_AT_SCAN, k = 4) + TCV + SITE, data = temp_H)
          gamm_model <- gamm4(y ~ s(x, k = 5) + TCV + SITE, data = temp_H)          
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
          results_H <- rbind(results_H, 
                             data.frame(brain = name_brain, cog = name_cog,
                                        R_squared = r_squared, Estimate = estimate, SE = se,
                                        t_value = t_value, P_value_t = p_value_t,
                                        edf = edf, red_edf = red_edf,
                                        F_value = f_value, P_value_smooth = p_value_smooth))
        }
      }
  }
}

# 对 p_value 进行校正
results_H$P_adj <- p.adjust(results_H$P_value_smooth, method = "fdr")

results_H<- results_H %>%
  arrange(P_adj)

L_sorted <- results_L[results_L$P_value_smooth < 0.05, ] # 删除数据框中P列大于或等于0.05的行
H_sorted <- results_H[results_H$P_value_smooth < 0.05, ]


# # 对除前两列（"VolumeName" 和 "ClusterID"）以外的数值型列保留4位小数
# model_stats[, -c(1, 2)] <- lapply(model_stats[, -c(1, 2)], function(x) {
#   if(is.numeric(x)) {
#     sprintf("%.4f", x)  # 保留4位小数并转换为字符格式
#   } else {
#     x
#   }
# })
# 
# # 设置列名
# colnames(model_stats)[3:11] <- c("R²", "估计值", "标准误", "t值", "p值", "有效自由度",
#                                  "剩余自由度", "F值", "平滑项p值")
# 
name <- paste0("corr_gamm_H_", resDate, ".csv")
write.csv(H_sorted, file.path(resuDir, name), row.names = F)
name <- paste0("corr_gamm_L_", resDate, ".csv")
write.csv(L_sorted, file.path(resuDir, name), row.names = F)

L_sorted$cluster <- "L"
H_sorted$cluster <- "H"

selected <- rbind(L_sorted, H_sorted)
name <- paste0("corr_gamm_LH_", resDate, ".csv")
write.csv(selected, file.path(resuDir, name), row.names = F)