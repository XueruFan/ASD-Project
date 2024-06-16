# 本代码用来分析两组ASD男性的变量之间的相关系数和显著性水平
# 斯皮尔曼相关，Site作为控制变量，34个脑指标分别作为自变量,
# 保存原始的相关系数和p值csv文件，另外，筛选p小于0.01且r大于0.2的结果，保存csv文件
# 按照显著性水平结果，绘制相关图png
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
clustDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster34DK")
statiDir <- file.path(abideDir, "Analysis/Statistic")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster34DK/Corr/CorrF")
resDate <- "240315"
newDate <- "240610"

# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_All_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "bankssts")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
# 根据脑指标对于人群聚类的贡献进行加权
rank <- read.csv(file.path(clustDir, paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_RM_Rank_",
                                            newDate, ".csv")))
All <- merge(cluster, pheno, by = "participant", All.x = TRUE)


# 选择自变量列
names_brain <- names(cluster)[3:ncol(cluster)]

names_cog <- c("FIQ", "VIQ", "PIQ", "ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT",
               "ADOS_2_RRB", "SRS_AWARENESS_RAW", "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW",
               "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW", "SRS_AWARENESS_T", "SRS_COGNITION_T",
               "SRS_COMMUNICATION_T", "SRS_MOTIVATION_T", "SRS_MANNERISMS_T", "SRS_TOTAL_T",
               "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
               "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C", "VINELAND_ABC_Standard",
               "VINELAND_COMMUNICATION_STANDARD", "VINELAND_DAILYLIVING_STANDARD",
               "VINELAND_SOCIAL_STANDARD", "BMI")



#### select: 控制变量、自变量、因变量
names_col <- c("clusterID", "SITE_ID", names_brain, names_cog)
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



################### L组

# 确保 SITE_ID 是因子类型
L$SITE_ID <- as.factor(L$SITE_ID)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(L) & name_cog %in% names(L)) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列
        temp_L <- L[!is.na(y) & !is.na(x), c(name_brain, "SITE_ID", name_cog)]
        temp_L$y <-  y[!is.na(y) & !is.na(x)]
        
        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID
            y_lm <- lm(y ~ SITE_ID, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID")), data = temp_L)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID
            y_lm <- lm(y ~ 1, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ 1")), data = temp_L)
            
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_statis_CorrF_L_", newDate, ".csv")
write.csv(L_sorted, file.path(statiDir, name), row.names = F)

################### H组

H$SITE_ID <- as.factor(H$SITE_ID)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列
        temp_H <- H[!is.na(y) & !is.na(x), c(name_brain, "SITE_ID", name_cog)]
        temp_H$y <-  y[!is.na(y) & !is.na(x)]
        
        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE_ID
            y_lm <- lm(y ~ SITE_ID, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE_ID")), data = temp_H)
            
            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)
            
            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 FIQ
            y_lm <- lm(y ~ 1, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ 1")), data = temp_H)
            
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_statis_CorrF_H_", newDate, ".csv")
write.csv(H_sorted, file.path(statiDir, name), row.names = F)

#################### 筛选显著且相关大于0.2的结果
H_sorted <- H_sorted %>%
  filter(abs(coef) > 0.2, p_value < 0.01)
L_sorted <- L_sorted %>%
  filter(abs(coef) > 0.2, p_value < 0.01)
colnames(L_sorted)[3:4] <- c("Lr", "Lp")
colnames(H_sorted)[3:4] <- c("Hr", "Hp")
sorted <- full_join(L_sorted, H_sorted, by = c("name_cog", "name_brain"))

sorted <- sorted %>%
  mutate(SortValue = ifelse(is.na(Lp), abs(Hp), abs(Lp)))
# 根据 SortValue 列对数据框进行排序
sorted <- sorted %>%
  arrange(SortValue)

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_statis_CorrF_Final", newDate, ".csv")
write.csv(sorted[, -7], file.path(statiDir, name), row.names = F)


##################################### Part 2: 画图 #################################################


for (i in 1:nrow(sorted)) {
  to_plot_names <- c(sorted[i, "name_brain"], sorted[i, "name_cog"])

  plotPoint_L <- L[, to_plot_names]
  plotPoint_L <- plotPoint_L[!is.na(plotPoint_L[[2]]), ]
  plotPoint_H <- H[, to_plot_names]
  plotPoint_H <- plotPoint_H[!is.na(plotPoint_H[[2]]), ]
  
  if (nrow(plotPoint_L) < 30 & nrow(plotPoint_H) < 30) {
    next  # 如果数据点过少，跳过当前循环
  }
  
  colnames(plotPoint_L)[1:2] <- c("x","y")
  colnames(plotPoint_H)[1:2] <- c("x","y")
  
  
  if (is.na(sorted[i, "Lp"])) {
    ggplot() +
      geom_point(data = plotPoint_L, aes(x = x, y = y, color = "lightgray"),
                 alpha = .8, size = 2, shape = 16) +
      geom_smooth(data = plotPoint_L, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
                  color = "lightgray", fill = "lightgray") +
      geom_point(data = plotPoint_H, aes(x = x, y = y, color = "#faa264"),
                 alpha = .8, size = 2, shape = 16) +
      geom_smooth(data = plotPoint_H, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
                  color = "#faa264", fill = "#faa264") +
      scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
      xlab(to_plot_names[1]) +
      ylab(to_plot_names[2]) +
      theme_cowplot() +
      theme(legend.position = "none", # without legend
            axis.text.y = element_text(size = 15, face = "bold"),
            axis.text.x = element_text(size = 15, face = "bold")) +
      scale_color_identity()
  } else {
    ggplot() +
      geom_point(data = plotPoint_H, aes(x = x, y = y, color = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264")),
                 alpha = .8, size = 2, shape = 16) +
      geom_smooth(data = plotPoint_H, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
                  color = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264"),
                  fill = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264")) +
      scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
      geom_point(data = plotPoint_L, aes(x = x, y = y, color = "#719988"),
                 alpha = .8, size = 2, shape = 16) +
      geom_smooth(data = plotPoint_L, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
                  color = "#719988", fill = "#719988") +
      xlab(to_plot_names[1]) +
      ylab(to_plot_names[2]) +
      theme_cowplot() +
      theme(legend.position = "none", # without legend
            axis.text.y = element_text(size = 15, face = "bold"),
            axis.text.x = element_text(size = 15, face = "bold")) +
      scale_color_identity()
  }
  name <- paste0("CorrF_", i, "_", to_plot_names[1], "_", to_plot_names[2], "_", newDate, ".png")
  ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
}
