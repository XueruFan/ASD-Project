# 本代码用来建模两个聚类人群脑指标centile的发育轨线
# GAMM模型
# Xue-Ru Fan 04 Jan 2024 @BNU
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


names_cog <- c("FIQ", "ADOS_2_SOCAFFECT" , "ADOS_2_RRB", "ADOS_2_TOTAL", "ADI_R_SOCIAL_TOTAL_A", 
               "ADI_R_RRB_TOTAL_C", "SRS_TOTAL_RAW", "SRS_COGNITION_RAW", "SRS_AWARENESS_RAW", 
               "SRS_COMMUNICATION_RAW", "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW")

#### select: 控制变量、自变量、因变量
# names_col <- c("clusterID", "AGE_AT_SCAN", "SITE_ID", "TCV_centile", names_brain, names_cog)
names_col <- c("clusterID", "SITE_ID", "TCV_centile", names_brain, names_cog)

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

for (name_brain in names_brain) {
  for (name_cog in names_cog) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y)) >= 30) {
        # temp_L <- L[!is.na(y), c(name_brain, "SITE_ID", "AGE_AT_SCAN", "TCV_centile", name_cog)]
        temp_L <- L[!is.na(y), c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        
        colnames(temp_L)[1] <- "x"
        colnames(temp_L)[4] <- "y"
        
        if (nrow(temp_L) >= 30) {
          
          # 拟合GAMM模型
          # gamm_model <- gamm4(y ~ s(x, k = 4) + s(AGE_AT_SCAN, k = 4) + TCV_centile + SITE_ID, data = temp_L)
          gamm_model <- gamm4(y ~ s(x, k = 5) + TCV_centile + SITE_ID, data = temp_L)
          
          
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


H$SITE_ID <- as.factor(H$SITE_ID)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) >= 30)) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        # temp_H <- H[!is.na(y) & !is.na(x), c(name_brain, "SITE_ID", "AGE_AT_SCAN", "TCV_centile", name_cog)]
        temp_H <- H[!is.na(y), c(name_brain, "SITE_ID", "TCV_centile", name_cog)]
        
        colnames(temp_H)[1] <- "x"
        colnames(temp_H)[4] <- "y"
        
        if (nrow(temp_H) >= 30) {
          
          # 拟合GAMM模型
          # gamm_model <- gamm4(y ~ s(x, k = 4) + s(AGE_AT_SCAN, k = 4) + TCV_centile + SITE_ID, data = temp_H)
          gamm_model <- gamm4(y ~ s(x, k = 5) + TCV_centile + SITE_ID, data = temp_H)          
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
# 保存模型统计值到 Excel 文件
write.xlsx(L_sorted, file.path(statiDir, paste0("corr_gamm_L_", newDate,
                                                  ".xlsx")), rowNames = F, colNames = T)
write.xlsx(H_sorted, file.path(statiDir, paste0("corr_gamm_H_", newDate,
                                               ".xlsx")), rowNames = F, colNames = T)

L_selected <- L_sorted#[L_sorted$coef < -0.3 | L_sorted$coef > 0.3, ]
H_selected <- H_sorted#[H_sorted$coef < -0.3 | H_sorted$coef > 0.3, ]
L_selected$cluster <- "L"
H_selected$cluster <- "H"

selected <- rbind(L_selected, H_selected)
name <- paste0("corr_gamm_LH_", newDate, ".csv")
write.csv(selected, file.path(statiDir, name), row.names = F)