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


# 给rank中的脑指标按照贡献大小赋值
rank$Feature <- paste0(rank$Feature, "_centile")
rank$Rank <- 0 # 创建一个新的列来保存排名
# 将正数和负数分开排名
positive_ranks <- rank(rank$Importance[rank$Importance > 0], ties.method = "first")
negative_ranks <- rank(-rank$Importance[rank$Importance < 0], ties.method = "first")
# 将排名结果赋值回数据框中
rank$Rank[rank$Importance > 0] <- positive_ranks
rank$Rank[rank$Importance < 0] <- -negative_ranks
rank$Importance <- rank$Rank
rank <- rank[, -3]


################### Part 1 ：首先使用全部脑指标进行加权后合成的指标作为自变量 ######################

################### Part 11
# 控制变量：site、FIQ
# 自变量：全部脑指标的加权合成新指标
# 因变量：认知

################### L组

importance <- rank %>% select(Feature, Importance)

# 确保 SITE_ID 和 FIQ 是因子和数值类型
L$SITE_ID <- as.factor(L$SITE_ID)
L$FIQ <- as.numeric(L$FIQ)

# 对 L 数据框中的 names_brain 列进行加权
weighted_names_brain <- L %>%
  select(all_of(importance$Feature)) %>%
  mutate(across(everything(), ~ . * importance$Importance[match(cur_column(), importance$Feature)]))

# 将加权后的列合并成一个新的自变量
weighted_var <- rowSums(weighted_names_brain, na.rm = TRUE)
if(length(weighted_var) != nrow(L)) {
  stop("Length of weighted_var does not match number of rows in L")
}

# 将加权合成的新变量添加到 L 数据框中
L$weighted_var <- weighted_var

# 循环 names_cog 列，计算 Spearman 相关性
for (name_cog in names_cog) {
  if (name_cog %in% names(L)) {
    y <- L[[name_cog]] # 提取因变量
    # 确保 y 中非 NA 值不少于 30 个
    if (sum(!is.na(y)) >= 30) {
      # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
      temp_L <- L[!is.na(y) & !is.na(L$FIQ), c("weighted_var", "SITE_ID", "FIQ", name_cog)]
      temp_L$y <- y[!is.na(y) & !is.na(L$FIQ)]
      
      if (nrow(temp_L) >= 30) {
        if (length(unique(temp_L$SITE_ID)) > 1) { # 确保 SITE_ID 还有多个水平
          # 对 y 和 weighted_var 进行回归，控制 SITE_ID 和 FIQ
          y_lm <- lm(y ~ SITE_ID + FIQ, data = temp_L)
          weighted_var_lm <- lm(weighted_var ~ SITE_ID + FIQ, data = temp_L)
          
          # 提取残差
          y_residuals <- residuals(y_lm)
          weighted_var_residuals <- residuals(weighted_var_lm)
          
          # 使用 Spearman 相关性计算残差之间的相关性
          cor_test <- cor.test(y_residuals, weighted_var_residuals, method = "spearman")
        } else {
          # 如果 SITE_ID 只有一个水平，不控制 SITE_ID，但仍然控制 FIQ
          y_lm <- lm(y ~ FIQ, data = temp_L)
          weighted_var_lm <- lm(weighted_var ~ FIQ, data = temp_L)
          
          # 提取残差
          y_residuals <- residuals(y_lm)
          weighted_var_residuals <- residuals(weighted_var_lm)
          
          # 使用 Spearman 相关性计算残差之间的相关性
          cor_test <- cor.test(y_residuals, weighted_var_residuals, method = "spearman")
        }
        
        # 保存 Spearman 相关系数和显著性到新的数据框中
        results_L <- rbind(results_L, data.frame(
          name_cog = name_cog,
          coef = cor_test$estimate,
          p_value = cor_test$p.value
        ))
      }
    }
  }
}

########### 给P值排序
L_sorted <- arrange(results_L, p_value)

### 保存下来结果
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_statis_CorrA_L_", newDate, ".csv")
write.csv(L_sorted, file.path(statiDir, name), row.names = F)
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_statis_CorrA_H_", newDate, ".csv")
write.csv(H_sorted, file.path(statiDir, name), row.names = F)




##################################### Part 2: 画图 #################################################

##### 画L组图 

for (i in 1:nrow(L_sorted)) {
  to_plot_names <- c(L_sorted[i, 1], L_sorted[i, 2])
  
  plotPoint <- L[, to_plot_names]
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  if (nrow(plotPoint) < 40) {
    next  # 如果行数少于40，跳过此次循环的剩余部分，也就是说，不够40个的就不看了
  }
  
  # 如果行数不少于40，继续执行下面的代码
  colnames(plotPoint) <- c("x","y")
  note_p <- paste0("p = ", round(L_sorted[i, "P"], 4))
  note_r <- paste0("r = ", round(L_sorted[i, "R"], 4))
  note_n <- paste0("n = ", round(L_sorted[i, "R"], 4))
  
  ggplot(plotPoint, aes(x = x, y = y)) +
    geom_point(color = "#add8e6", alpha = .8, size = 2, shape = 16) +  # 添加散点图层
    geom_smooth(method = "lm", se = T, lwd = 2, color = "#add8e6", fill = "#add8e6") +
    theme_cowplot() +
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
    xlab(to_plot_names[1]) +
    ylab(to_plot_names[2]) +
    annotate("text", x = Inf, y = Inf, label = note_r, hjust = 1.2, vjust = 3, size = 7) +
    annotate("text", x = Inf, y = Inf, label = note_p, hjust = 1.2, vjust = 1.2, size = 7) +
    
    theme(legend.position = "none", # without legend
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"))
  
  name <- paste0("L_", i, "_", to_plot_names[1], "_", to_plot_names[2], "_", newDate, ".png")
  ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
}


####################################### 画H组图 ####################################################

for (i in 1:nrow(H_sorted)) {
  to_plot_names <- c(H_sorted[i, 1], H_sorted[i, 2])
  
  plotPoint <- H[, to_plot_names]
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  if (nrow(plotPoint) < 40) {
    next  # 如果行数少于40，跳过此次循环的剩余部分，也就是说，不够40个的就不看了
  }
  
  # 如果行数不少于40，继续执行下面的代码
  colnames(plotPoint) <- c("x","y")
  note_p <- paste0("p = ", round(H_sorted[i, "P"], 4))
  note_r <- paste0("r = ", round(H_sorted[i, "R"], 4))

  ggplot(plotPoint, aes(x = x, y = y)) +
    geom_point(color = "#ffb699", alpha = .8, size = 2, shape = 16) +  # 添加散点图层
    geom_smooth(method = "lm", se = T, lwd = 2, color = "#ffb699", fill = "#ffb699") +
    theme_cowplot() +
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
    xlab(to_plot_names[1]) +
    ylab(to_plot_names[2]) +
    annotate("text", x = Inf, y = Inf, label = note_r, hjust = 1.2, vjust = 3, size = 7) +
    annotate("text", x = Inf, y = Inf, label = note_p, hjust = 1.2, vjust = 1.2, size = 7) +
    theme(legend.position = "none", # without legend
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 15, face = "bold"))
  
  name <- paste0("H_", i, "_", to_plot_names[1], "_", to_plot_names[2], "_", newDate, ".png")
  ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
}
