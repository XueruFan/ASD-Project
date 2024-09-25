# 本代码用来可视化两个聚类人群不同的认知行为和脑指标之间的相关
# 只画上一步里相关大于0.2并且p值小于0.01的相关
# 为了使所有的线都能集中一些，画图时给y做了标准化
# Vineland量表和PIQ的得分是越高越好，因此画图时原始值*-1，这样图示中都是y越高症状越严重
# 雪如 2024年8月5日于北师大办公室
##################################
# Part 1: L组，png
# Part 2: H组，png
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
clustDir <- file.path(abideDir, "Analysis/Cluster/SpectralCluster")
statiDir <- file.path(abideDir, "Analysis/Statistic")
plotDir <- file.path(abideDir, "Plot/Cluster/SpectralCluster/Corr")
resDate <- "240315"
newDate <- "240610"

# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("asd_male_Spectral_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "GMV")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
All <- merge(cluster, pheno, by = "participant", All.x = TRUE)

L <- subset(All, clusterID == "1")
L[L < 0] <- NA

H <- subset(All, clusterID == "2")
L[L < 0] <- NA

corr <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_corr_Part3_Site_LH_", newDate,
                                            "_forPlot.csv")))
corr_L <- corr[!is.na(corr$Lr), 1:2]
corr_H <- corr[!is.na(corr$Hr), 1:2]



##################### Part 1: 画L组 ################################################################
# 创建一个空的数据框来累积所有的绘图数据
plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(corr_L)) {
  to_plot_names <- c(corr_L[i, "name_brain"], corr_L[i, "name_cog"])
  
  plotPoint_L <- L[, to_plot_names]
  plotPoint_L <- plotPoint_L[!is.na(plotPoint_L[[2]]), ]
  colnames(plotPoint_L)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint_L$type <- paste0(corr_L[i, "name_brain"], "_", corr_L[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint_L)
}

grouped_data <- group_by(plot_data, type)

# 给y值都做一下标准化，这样画出的线之间距离不远
plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

# Vineland的量表是得分越高越好
plot_data <- plot_data %>%
  mutate(y_scaled = case_when(
    type %in% c("sGMV_centile_VINELAND_ABC_Standard", "WMV_centile_VINELAND_COMMUNICATION_STANDARD")
    ~ -1 * y_scaled,
    TRUE ~ y_scaled
  ))


# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("entorhinal_centile_BMI",
                                        "parstriangularis_centile_ADI_R_NONVERBAL_TOTAL_BV",
                                        "insula_centile_ADI_R_NONVERBAL_TOTAL_BV",
                                        "lingual_centile_ADI_R_SOCIAL_TOTAL_A",
                                        "rostralmiddlefrontal_centile_SRS_AWARENESS_T",
                                        "rostralanteriorcingulate_centile_SRS_TOTAL_T",
                                        "rostralanteriorcingulate_centile_SRS_COMMUNICATION_T",
                                        "rostralanteriorcingulate_centile_SRS_MANNERISMS_T",
                                        "caudalanteriorcingulate_centile_SRS_MANNERISMS_T",
                                        "sGMV_centile_VINELAND_ABC_Standard",                  
                                        "WMV_centile_VINELAND_COMMUNICATION_STANDARD")))


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_SRS_TOTAL_T"] <- "RAC - SRS T"
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_SRS_COMMUNICATION_T"] <- "RAC - SRS C"
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_SRS_MANNERISMS_T"] <- "RAC - SRS M"
levels(plot_data$type)[levels(plot_data$type) ==
                         "caudalanteriorcingulate_centile_SRS_MANNERISMS_T"] <- "CAC - SRS M"
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralmiddlefrontal_centile_SRS_AWARENESS_T"] <- "RMF - SRS A"
levels(plot_data$type)[levels(plot_data$type) ==
                         "insula_centile_ADI_R_NONVERBAL_TOTAL_BV"] <- "insula - ADIR N"
levels(plot_data$type)[levels(plot_data$type) ==
                         "parstriangularis_centile_ADI_R_NONVERBAL_TOTAL_BV"] <- "PT - ADIR N"
levels(plot_data$type)[levels(plot_data$type) ==
                         "lingual_centile_ADI_R_SOCIAL_TOTAL_A"] <- "lingual - ADIR S"
levels(plot_data$type)[levels(plot_data$type) ==
                         "sGMV_centile_VINELAND_ABC_Standard"] <- "sGMV - VL A"
levels(plot_data$type)[levels(plot_data$type) ==
                         "WMV_centile_VINELAND_COMMUNICATION_STANDARD"] <- "WMV - VL C"
levels(plot_data$type)[levels(plot_data$type) ==
                         "entorhinal_centile_BMI"] <- "entorhinal - BMI"



ggplot(plot_data, aes(x = x, y = y_scaled, group = type, color = type)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_fixed(ratio = .7) +  # 设置坐标轴比例为 1:1
  theme_cowplot() +
  theme(legend.position = "right",
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold")) +
  scale_color_manual(values = c("RAC - SRS M" = "#FB9A99",
                                "RAC - SRS C" = "#FB9A99",
                                "RAC - SRS T" = "#FB9A99",
                                "CAC - SRS M" = "#FB9A99",
                                "RMF - SRS A" = "#FF7F00",
                                "lingual - ADIR S" = "#ad5f8a",              
                                "insula - ADIR N" = "#6A3D9A",             
                                "PT - ADIR N" = "blue",   
                                "sGMV - VL A" = "#999999",                  
                                "WMV - VL C" = "#999999",         
                                "entorhinal - BMI" = "black"))

name <- paste0("SC_corr_Site_L_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)


##################### Part 2: 画H组 ################################################################
# 创建一个空的数据框来累积所有的绘图数据
plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(corr_H)) {
  to_plot_names <- c(corr_H[i, "name_brain"], corr_H[i, "name_cog"])
  
  plotPoint_H <- H[, to_plot_names]
  plotPoint_H <- plotPoint_H[!is.na(plotPoint_H[[2]]), ]
  colnames(plotPoint_H)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint_H$type <- paste0(corr_H[i, "name_brain"], "_", corr_H[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint_H)
}

grouped_data <- group_by(plot_data, type)

# 给y值都做一下标准化，这样画出的线之间距离不远
plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

# Vineland量表、IQ是得分越高越好
plot_data <- plot_data %>%
  mutate(y_scaled = case_when(
    type %in% c("parsorbitalis_centile_VINELAND_SOCIAL_STANDARD", "medialorbitofrontal_centile_PIQ")
    ~ -1 * y_scaled,
    TRUE ~ y_scaled
  ))


# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("temporalpole_centile_BMI",
                                        "medialorbitofrontal_centile_PIQ",
                                        "rostralanteriorcingulate_centile_SRS_MANNERISMS_T",
                                        "rostralanteriorcingulate_centile_SRS_AWARENESS_T",
                                        "caudalanteriorcingulate_centile_SRS_AWARENESS_T",
                                        "parsorbitalis_centile_VINELAND_SOCIAL_STANDARD")))


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_SRS_MANNERISMS_T"] <- "RAC - SRS M"
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_SRS_AWARENESS_T"] <- "RAC - SRS A"
levels(plot_data$type)[levels(plot_data$type) ==
                         "caudalanteriorcingulate_centile_SRS_AWARENESS_T"] <- "CAC - SRS A"
levels(plot_data$type)[levels(plot_data$type) ==
                         "parsorbitalis_centile_VINELAND_SOCIAL_STANDARD"] <- "PA - VL S"
levels(plot_data$type)[levels(plot_data$type) ==
                         "medialorbitofrontal_centile_PIQ"] <- "MOF - PIQ"
levels(plot_data$type)[levels(plot_data$type) ==
                         "temporalpole_centile_BMI"] <- "TP - BMI"

ggplot(plot_data, aes(x = x, y = y_scaled, group = type, color = type)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_fixed(ratio = .7) +  # 设置坐标轴比例为 1:1
  theme_cowplot() +
  theme(legend.position = "right",
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold")) +
  scale_color_manual(values = c("RAC - SRS M" = "#FB9A99",
                                "RAC - SRS A" = "#FB9A99",
                                "CAC - SRS A" = "#FB9A99",
                                "PA - VL S" = "#E3B839",
                                "MOF - PIQ" = "#999999",         
                                "TP - BMI" = "black"))

name <- paste0("SC_corr_Site_H_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)