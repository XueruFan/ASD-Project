# 本代码用来可视化两个聚类人群不同的认知行为和脑指标之间的相关
# 只画上一步里相关大于0.4的相关
# Site Part3
# 雪如 2024年8月5日于北师大办公室
##################################
# Part 1: L组，png
# Part 2: H组，png
##################################
rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss", "Cairo")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
statiDir <- file.path(abideDir, "Analysis/Statistic/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster/Corr")
resDate <- "240315"
newDate <- "240610"

# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "GMV")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
All <- merge(cluster, pheno, by = "participant", All.x = TRUE)

L <- subset(All, clusterID == "1")
L[L < 0] <- NA

H <- subset(All, clusterID == "2")
L[L < 0] <- NA

corr <- read.csv(file.path(statiDir, paste0("asd_male_dev_GC_corr_Part3_Site_LH_", newDate,
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

# # Vineland的量表是得分越高越好
# plot_data <- plot_data %>%
#   mutate(y_scaled = case_when(
#     type %in% c("sGMV_centile_VINELAND_ABC_Standard", "WMV_centile_VINELAND_COMMUNICATION_STANDARD")
#     ~ -1 * y_scaled,
#     TRUE ~ y_scaled
#   ))
# 
# 
# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("WMV_centile_VINELAND_COMMUNICATION_STANDARD",
                                        "sGMV_centile_VINELAND_COMMUNICATION_STANDARD",
                                        "WMV_centile_VINELAND_ABC_Standard",
                                        "sGMV_centile_VINELAND_ABC_Standard",
                                        "TCV_centile_VINELAND_ABC_Standard",
                                        "totalSA2_centile_VINELAND_ABC_Standard")))


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "WMV_centile_VINELAND_ABC_Standard"] <- "WMV - 综合分"
levels(plot_data$type)[levels(plot_data$type) ==
                         "WMV_centile_VINELAND_COMMUNICATION_STANDARD"] <- "WMV - 沟通"
levels(plot_data$type)[levels(plot_data$type) ==
                         "sGMV_centile_VINELAND_ABC_Standard"] <- "sGMV - 综合分"
levels(plot_data$type)[levels(plot_data$type) ==
                         "sGMV_centile_VINELAND_COMMUNICATION_STANDARD"] <- "sGMV - 沟通"
levels(plot_data$type)[levels(plot_data$type) ==
                         "TCV_centile_VINELAND_ABC_Standard"] <- "TCV - 综合分"
levels(plot_data$type)[levels(plot_data$type) ==
                         "totalSA2_centile_VINELAND_ABC_Standard"] <- "SA - 综合分"


name <- paste0("GC_corr_Site_L_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

ggplot(plot_data, aes(x = x, y = y, group = type, color = type)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = .7) +  # 设置坐标轴比例为 1:1
  # coord_fixed(ratio = .7) +  # 设置坐标轴比例为 1:1
  theme_cowplot() +
  theme(text = element_text(family = "STSong"),
        legend.position = "right",
        legend.title = element_blank(),
        # axis.text.y = element_blank(),    # 去掉y轴文字
        # axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold")) +
  scale_color_manual(values = c("WMV - 综合分" = "#f8b500",
                                "WMV - 沟通" = "black",
                                "sGMV - 综合分" = "#007bbb",
                                "sGMV - 沟通" = "#884898",
                                "TCV - 综合分" = "#df7163",
                                "SA - 综合分" = "#00a381"))

dev.off()






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
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

# # Vineland量表、IQ是得分越高越好
# plot_data <- plot_data %>%
#   mutate(y_scaled = case_when(
#     type %in% c("parsorbitalis_centile_VINELAND_SOCIAL_STANDARD", "medialorbitofrontal_centile_PIQ")
#     ~ -1 * y_scaled,
#     TRUE ~ y_scaled
#   ))


# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("superiorfrontal_centile_BMI",
                                        "caudalmiddlefrontal_centile_BMI",
                                        "superiorparietal_centile_BMI",
                                        "isthmuscingulate_centile_BMI",
                                        "pericalcarine_centile_BMI",
                                        "temporalpole_centile_BMI")))

plot_data$y[plot_data$y < 0] <- NA
plot_data <- na.omit(plot_data)


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "superiorfrontal_centile_BMI"] <- "额上回"
levels(plot_data$type)[levels(plot_data$type) ==
                         "caudalmiddlefrontal_centile_BMI"] <- "额中回尾部"
levels(plot_data$type)[levels(plot_data$type) ==
                         "superiorparietal_centile_BMI"] <- "顶上皮层"
levels(plot_data$type)[levels(plot_data$type) ==
                         "isthmuscingulate_centile_BMI"] <- "扣带回峡部"
levels(plot_data$type)[levels(plot_data$type) ==
                         "pericalcarine_centile_BMI"] <- "距状沟周围皮层"
levels(plot_data$type)[levels(plot_data$type) ==
                         "temporalpole_centile_BMI"] <- "颞极"

name <- paste0("GC_corr_Site_H_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

ggplot(plot_data, aes(x = x, y = y, group = type, color = type)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_cowplot() +
  theme(text = element_text(family = "STSong"),
        legend.position = "right",
        legend.title = element_blank(),
        # axis.text.y = element_blank(),    # 去掉y轴文字
        # axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold")) +
  scale_color_manual(values = c("额上回" = "#f8b500",
                                "额中回尾部" = "black",
                                "顶上皮层" = "#007bbb",
                                "扣带回峡部" = "#884898",
                                "距状沟周围皮层" = "#df7163",
                                "颞极" = "#00a381"))

dev.off()
