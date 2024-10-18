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
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect513/Corr/Part4")
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
All[All < 0] <- NA
All$clusterID[All$clusterID == 1] <- "L"
All$clusterID[All$clusterID == 2] <- "H"

corr_l <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_L_", newDate, ".csv")))
corr_h <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_H_", newDate, ".csv")))

corr_l$cluster <- "L"
corr_h$cluster <- "H"

corr <- rbind(corr_l, corr_h)

corr_ados <- corr[grep("ADOS", corr$name_cog), ]
# corr_srs <- corr[grep("SRS", corr$name_cog), ]
corr_fiq <- corr[grep("FIQ", corr$name_cog), ]


# 定义自定义颜色向量
custom_colors <- c("#c85554", "#69821b", "#93ca76", "#9CB0C3", "#43676b", "#ee7800",
                   "#dccb18", "#8491c3", "#00a381", "#cc7eb1", "#e09e87")


##################### Part 1: ados-rrb ################################################################
# 创建一个空的数据框来累积所有的绘图数据

temp <- corr_ados[grep("RRB", corr_ados$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
  
  plotPoint <- subset(All, clusterID == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]
  
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  colnames(plotPoint)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, clusterID)

# 给y值都做一下标准化，这样画出的线之间距离不远
plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

# Vineland的量表是得分越高越好
# plot_data <- plot_data %>%
#   mutate(y_scaled = case_when(
#     type %in% c("sGMV_centile_VINELAND_ABC_Standard", "WMV_centile_VINELAND_COMMUNICATION_STANDARD")
#     ~ -1 * y_scaled,
#     TRUE ~ y_scaled
#   ))
unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("parstriangularis_centile_ADOS_2_RRB",
                                        "postcentral_centile_ADOS_2_RRB",
                                        "supramarginal_centile_ADOS_2_RRB",
                                        "inferiortemporal_centile_ADOS_2_RRB",
                                        "entorhinal_centile_ADOS_2_RRB")))

# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "parstriangularis_centile_ADOS_2_RRB"] <- "Pars triangularis"
levels(plot_data$type)[levels(plot_data$type) ==
                         "postcentral_centile_ADOS_2_RRB"] <- "Postcentral"
levels(plot_data$type)[levels(plot_data$type) ==
                         "entorhinal_centile_ADOS_2_RRB"] <- "Entorhinal"
levels(plot_data$type)[levels(plot_data$type) ==
                         "supramarginal_centile_ADOS_2_RRB"] <- "Supramarginal"
levels(plot_data$type)[levels(plot_data$type) ==
                         "inferiortemporal_centile_ADOS_2_RRB"] <- "Inferior temporal"



ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
                      color = type, linetype = clusterID)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
  guides(linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.86),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(6,7,10,4,3)]) +
  annotate("text", x = 0.04, y = -0.26, label = "*", size = 12, color = custom_colors[6]) +
  annotate("text", x = 0.04, y = -0.12, label = "*", size = 12, color = custom_colors[7]) +
  annotate("text", x = 0.97, y = 0.1, label = "*", size = 12, color = custom_colors[10]) +
  annotate("text", x = 0.04, y = -0.44, label = "*", size = 12, color = custom_colors[4]) +
  annotate("text", x = 0.06, y = -0.6, label = "**", size = 12, color = custom_colors[3]) +
  annotate("text", x = 0.7, y = -0.57, label = "Restricted and Repetitive Behavior", size = 6, color = "#474a4d")

name <- paste0("corr_ADOS_2_RRB_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



##################### Part 2: ados-soca ############################################################
# 创建一个空的数据框来累积所有的绘图数据

temp <- corr_ados[grep("SOCA", corr_ados$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
  
  plotPoint <- subset(All, clusterID == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]
  
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  colnames(plotPoint)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, clusterID)

# 给y值都做一下标准化，这样画出的线之间距离不远
plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("rostralanteriorcingulate_centile_ADOS_2_SOCAFFECT",
                                        "insula_centile_ADOS_2_SOCAFFECT",
                                        "transversetemporal_centile_ADOS_2_SOCAFFECT")))


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralanteriorcingulate_centile_ADOS_2_SOCAFFECT"] <- "Rostral anterior cingulate"
levels(plot_data$type)[levels(plot_data$type) ==
                         "insula_centile_ADOS_2_SOCAFFECT"] <- "Insula"
levels(plot_data$type)[levels(plot_data$type) ==
                         "transversetemporal_centile_ADOS_2_SOCAFFECT"] <- "Transverse temporal"



ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
                      color = type, linetype = clusterID)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
  guides(linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.9),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(8,5,11)]) +
  annotate("text", x = 0.08, y = -0.62, label = "*", size = 12, color = custom_colors[8]) +
  annotate("text", x = 0.07, y = -0.46, label = "**", size = 12, color = custom_colors[5]) +
  annotate("text", x = 0.2, y = -0.3, label = "*", size = 12, color = custom_colors[11]) +
  annotate("text", x = 0.9, y = -0.59, label = "Social Affect", size = 6, color = "#474a4d")

name <- paste0("corr_ADOS_2_SOCAF_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



##################### Part 3: fiq ################################################################
# 创建一个空的数据框来累积所有的绘图数据

temp <- corr_fiq[grep("FIQ", corr_fiq$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")

  plotPoint <- subset(All, clusterID == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]

  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]

  colnames(plotPoint)[1:2] <- c("x","y")

  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, clusterID)

# 给y值都做一下标准化，这样画出的线之间距离不远
plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("superiorfrontal_centile_FIQ",
                                        "rostralmiddlefrontal_centile_FIQ",
                                        "caudalanteriorcingulate_centile_FIQ")))


# modify name
levels(plot_data$type)[levels(plot_data$type) ==
                         "superiorfrontal_centile_FIQ"] <- "Superior frontal"
levels(plot_data$type)[levels(plot_data$type) ==
                         "rostralmiddlefrontal_centile_FIQ"] <- "Rostral middle frontal"
levels(plot_data$type)[levels(plot_data$type) ==
                         "caudalanteriorcingulate_centile_FIQ"] <- "Caudal anterior cingulate"



ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
                      color = type, linetype = clusterID)) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.9),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(9,1,2)]) +
  annotate("text", x = 0.98, y = -0.31, label = "**", size = 12, color = custom_colors[9]) +
  annotate("text", x = 0.04, y = -0.19, label = "**", size = 12, color = custom_colors[1]) +
  annotate("text", x = 0.96, y = -0.17, label = "*", size = 12, color = custom_colors[1]) +
  annotate("text", x = 0.06, y = -0.41, label = "*", size = 12, color = custom_colors[2]) +
  annotate("text", x = 0.8, y = -0.37, label = "FIQ", size = 6, color = "#474a4d")

name <- paste0("corr_FIQ_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)




##################### Part 3: ados-total ################################################################
# 创建一个空的数据框来累积所有的绘图数据

# temp <- corr_ados[grep("TOTAL", corr_ados$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# 
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("insula_centile_ADOS_2_TOTAL",
#                                         "transversetemporal_centile_ADOS_2_TOTAL",
#                                         "entorhinal_centile_ADOS_2_TOTAL")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "entorhinal_centile_ADOS_2_TOTAL"] <- "Entorhinal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "insula_centile_ADOS_2_TOTAL"] <- "Insula"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "transversetemporal_centile_ADOS_2_TOTAL"] <- "Transverse temporal"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.9),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(3,10,5)]) +
#   annotate("text", x = 0.06, y = -0.42, label = "*", size = 12, color = custom_colors[3]) +
#   annotate("text", x = 0.26, y = -0.42, label = "*", size = 12, color = custom_colors[10]) +
#   annotate("text", x = 0.06, y = -0.31, label = "*", size = 12, color = custom_colors[5]) +
#   annotate("text", x = 0.95, y = -0.41, label = "Total", size = 6, color = "#474a4d")
# 
# name <- paste0("corr_ADOS_2_TOTAL_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



# ##################### Part 4: srs-cog ################################################################
# # 创建一个空的数据框来累积所有的绘图数据
# 
# temp <- corr_srs[grep("COG", corr_srs$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# 
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("rostralmiddlefrontal_centile_SRS_COGNITION_RAW",
#                                         "rostralanteriorcingulate_centile_SRS_COGNITION_RAW",
#                                         "precentral_centile_SRS_COGNITION_RAW")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralmiddlefrontal_centile_SRS_COGNITION_RAW"] <- "Rostral middle frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralanteriorcingulate_centile_SRS_COGNITION_RAW"] <- "Rostral anterior cingulate"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "precentral_centile_SRS_COGNITION_RAW"] <- "Precentral"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.9),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(1,8,4)]) +
#   annotate("text", x = 0.66, y = 0.6, label = "**", size = 12, color = custom_colors[1]) +
#   annotate("text", x = 0.8, y = 0.6, label = "*", size = 12, color = custom_colors[8]) +
#   annotate("text", x = 0.8, y = 0.4, label = "*", size = 12, color = custom_colors[4]) +
#   annotate("text", x = 0.86, y = -0.37, label = "Social Cognition", size = 6, color = "#474a4d")
# 
# name <- paste0("corr_SRS_COG_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# 
# 
# 
# ##################### Part 5: srs-com ################################################################
# # 创建一个空的数据框来累积所有的绘图数据
# 
# temp <- corr_srs[grep("COM", corr_srs$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# 
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("superiorfrontal_centile_SRS_COMMUNICATION_RAW",
#                                         "rostralmiddlefrontal_centile_SRS_COMMUNICATION_RAW",
#                                         "rostralanteriorcingulate_centile_SRS_COMMUNICATION_RAW",
#                                         "postcentral_centile_SRS_COMMUNICATION_RAW",
#                                         "insula_centile_SRS_COMMUNICATION_RAW",
#                                         "fusiform_centile_SRS_COMMUNICATION_RAW")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralmiddlefrontal_centile_SRS_COMMUNICATION_RAW"] <- "Rostral middle frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralanteriorcingulate_centile_SRS_COMMUNICATION_RAW"] <- "Rostral anterior cingulate"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "postcentral_centile_SRS_COMMUNICATION_RAW"] <- "Postcentral"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "superiorfrontal_centile_SRS_COMMUNICATION_RAW"] <- "Superior frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "insula_centile_SRS_COMMUNICATION_RAW"] <- "Insula"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "fusiform_centile_SRS_COMMUNICATION_RAW"] <- "Fusiform"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.84),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         # legend.text = element_text(size = 10),  # 减小图例字体大小
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(9,1,8,2,3,6)]) +
#   annotate("text", x = 0.92, y = 0.65, label = "*", size = 12, color = custom_colors[9]) +
#   annotate("text", x = 0.7, y = 0.65, label = "**", size = 12, color = custom_colors[1]) +
#   annotate("text", x = 0.87, y = 0.56, label = "*", size = 12, color = custom_colors[8]) +
#   annotate("text", x = 0.87, y = 0.5, label = "*", size = 12, color = custom_colors[2]) +
#   annotate("text", x = 0.94, y = 0.05, label = "**", size = 12, color = custom_colors[3]) +
#   annotate("text", x = 0.94, y = 0.18, label = "**", size = 12, color = custom_colors[6]) +
#   annotate("text", x = 0.82, y = -0.42, label = "Social Communication", size = 6, color = "#474a4d")
# 
# 
# name <- paste0("corr_SRS_COM_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# 
# 
# 
# 
# ##################### Part 6: srs-man ################################################################
# # 创建一个空的数据框来累积所有的绘图数据
# 
# temp <- corr_srs[grep("MAN", corr_srs$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("parsorbitalis_centile_SRS_MANNERISMS_RAW",
#                                         "rostralmiddlefrontal_centile_SRS_MANNERISMS_RAW",
#                                         "postcentral_centile_SRS_MANNERISMS_RAW",
#                                         "isthmuscingulate_centile_SRS_MANNERISMS_RAW",
#                                         "lateraloccipital_centile_SRS_MANNERISMS_RAW")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralmiddlefrontal_centile_SRS_MANNERISMS_RAW"] <- "Rostral middle frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "isthmuscingulate_centile_SRS_MANNERISMS_RAW"] <- "Isthmus cingulate"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "postcentral_centile_SRS_MANNERISMS_RAW"] <- "Postcentral"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "lateraloccipital_centile_SRS_MANNERISMS_RAW"] <- "Lateral occipital"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "parsorbitalis_centile_SRS_MANNERISMS_RAW"] <- "Pars orbitalis"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.86),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         # legend.text = element_text(size = 10),  # 减小图例字体大小
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(15,1,2,7,13)]) +
#   annotate("text", x = 0.94, y = 0.10, label = "*", size = 12, color = custom_colors[15]) +
#   annotate("text", x = 0.81, y = 0.82, label = "**", size = 12, color = custom_colors[1]) +
#   annotate("text", x = 0.82, y = 0.62, label = "*", size = 12, color = custom_colors[2]) +
#   annotate("text", x = 0.93, y = 0.45, label = "**", size = 12, color = custom_colors[7]) +
#   annotate("text", x = 0.94, y = -0.5, label = "*", size = 12, color = custom_colors[13]) +
#   annotate("text", x = 0.84, y = -0.65, label = "Social Mannerisms", size = 6, color = "#474a4d")
# 
# name <- paste0("corr_SRS_MAN_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# 
# 
# 
# 
# ##################### Part 7: srs-mon ################################################################
# # 创建一个空的数据框来累积所有的绘图数据
# 
# temp <- corr_srs[grep("MOT", corr_srs$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# 
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("superiorfrontal_centile_SRS_MOTIVATION_RAW",
#                                         "rostralmiddlefrontal_centile_SRS_MOTIVATION_RAW",
#                                         "postcentral_centile_SRS_MOTIVATION_RAW",
#                                         "entorhinal_centile_SRS_MOTIVATION_RAW",
#                                         "fusiform_centile_SRS_MOTIVATION_RAW")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralmiddlefrontal_centile_SRS_MOTIVATION_RAW"] <- "Rostral middle frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "superiorfrontal_centile_SRS_MOTIVATION_RAW"] <- "Superior frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "postcentral_centile_SRS_MOTIVATION_RAW"] <- "Postcentral"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "entorhinal_centile_SRS_MOTIVATION_RAW"] <- "Entorhinal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "fusiform_centile_SRS_MOTIVATION_RAW"] <- "Fusiform"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.86),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         # legend.text = element_text(size = 10),  # 减小图例字体大小
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(9,1,2,5,6)]) +
#   annotate("text", x = 0.85, y = 0.63, label = "*", size = 12, color = custom_colors[9]) +
#   annotate("text", x = 0.01, y = -0.28, label = "*", size = 12, color = custom_colors[1]) +
#   annotate("text", x = 0.85, y = 0.55, label = "*", size = 12, color = custom_colors[2]) +
#   annotate("text", x = 0.94, y = 0.3, label = "*", size = 12, color = custom_colors[5]) +
#   annotate("text", x = 0.94, y = 0.15, label = "*", size = 12, color = custom_colors[6]) +
#   annotate("text", x = 0.84, y = -0.35, label = "Social Motivation", size = 6, color = "#474a4d")
# 
# name <- paste0("corr_SRS_MOT_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# 
# 
# 
# ##################### Part 8: srs-tot ################################################################
# # 创建一个空的数据框来累积所有的绘图数据
# 
# temp <- corr_srs[grep("TOTAL", corr_srs$name_cog), ]
# 
# plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
# 
# for (i in 1:nrow(temp)) {
#   to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "clusterID")
#   
#   plotPoint <- subset(All, clusterID == temp[i, "cluster"])
#   plotPoint <- plotPoint[, to_plot_names]
#   
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   colnames(plotPoint)[1:2] <- c("x","y")
#   
#   # 添加数据到累积数据框
#   plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
#   plot_data <- rbind(plot_data, plotPoint)
# }
# grouped_data <- group_by(plot_data, type, clusterID)
# 
# # 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))
# 
# # 设置 type 的顺序
# plot_data <- plot_data %>%
#   mutate(type = factor(type, levels = c("superiorfrontal_centile_SRS_TOTAL_RAW",
#                                         "rostralmiddlefrontal_centile_SRS_TOTAL_RAW",
#                                         "postcentral_centile_SRS_TOTAL_RAW",
#                                         "insula_centile_SRS_TOTAL_RAW",
#                                         "fusiform_centile_SRS_TOTAL_RAW")))
# 
# 
# # modify name
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "rostralmiddlefrontal_centile_SRS_TOTAL_RAW"] <- "Rostral middle frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "superiorfrontal_centile_SRS_TOTAL_RAW"] <- "Superior frontal"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "postcentral_centile_SRS_TOTAL_RAW"] <- "Postcentral"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "insula_centile_SRS_TOTAL_RAW"] <- "Insula"
# levels(plot_data$type)[levels(plot_data$type) ==
#                          "fusiform_centile_SRS_TOTAL_RAW"] <- "Fusiform"
# 
# 
# 
# ggplot(plot_data, aes(x = x, y = y_scaled, group = interaction(type, clusterID),
#                       color = type, linetype = clusterID)) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 2) +
#   scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#   # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
#   scale_linetype_manual(values = c("L" = "twodash", "H" = "solid")) +
#   guides(color = guide_legend(order = 1),
#          linetype = guide_legend(override.aes = list(color = "#474a4d"))) +  # 强制线形颜色为灰色
#   theme_cowplot() +
#   theme(legend.position = c(0.02, 0.86),
#         legend.title = element_blank(),
#         legend.key.width  = unit(1, "cm"), # 控制图例线条
#         # legend.text = element_text(size = 10),  # 减小图例字体大小
#         axis.text.y = element_blank(),    # 去掉y轴文字
#         axis.ticks.y = element_blank(),  # 去掉y轴刻度
#         axis.title = element_blank(),
#         axis.text.x = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors[c(9,1,2,3,6)]) +
#   annotate("text", x = 0.92, y = 0.56, label = "*", size = 12, color = custom_colors[9]) +
#   annotate("text", x = 0.72, y = 0.56, label = "**", size = 12, color = custom_colors[1]) +
#   annotate("text", x = 0.94, y = 0.65, label = "**", size = 12, color = custom_colors[2]) +
#   annotate("text", x = 0.95, y = 0.03, label = "*", size = 12, color = custom_colors[3]) +
#   annotate("text", x = 0.95, y = 0.2, label = "*", size = 12, color = custom_colors[6]) +
#   annotate("text", x = 0.95, y = -0.35, label = "Total", size = 6, color = "#474a4d")
# 
# name <- paste0("corr_SRS_TOT_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
