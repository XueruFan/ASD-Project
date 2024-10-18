
rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/513/Corr")
plotDir <- file.path(cabicDir, "result/pred/513/Plot/Corr")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_513_", resDate, ".csv")
cluster <- read.csv(file.path(cabicDir, "result/pred/513", name))

All <- merge(cluster, pheno, by = "participant")
All$predicted_cluster[All$predicted_cluster == 1] <- "L"
All$predicted_cluster[All$predicted_cluster == 2] <- "H"

corr_l <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_L_", resDate, ".csv")))
corr_h <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_H_", resDate, ".csv")))

corr_l$cluster <- "L"
corr_h$cluster <- "H"

corr <- rbind(corr_l, corr_h)
corr <- corr[c(23,25,31),]

# corr_ados <- corr[grep("ADOS", corr$name_cog), ]
# corr_srs <- corr[grep("SRS", corr$name_cog), ]


# 定义自定义颜色向量
custom_colors <- c("#c85554", "#69821b", "#93ca76", "#9CB0C3", "#43676b", "#ee7800",
                   "#dccb18", "#8491c3", "#00a381", "#cc7eb1", "#e09e87", "#83ccd2")


##################### Part 1: ados-rrb ################################################################
# 创建一个空的数据框来累积所有的绘图数据

temp <- corr[grep("RRB", corr$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "predicted_cluster")
  
  plotPoint <- subset(All, predicted_cluster == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]
  
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  colnames(plotPoint)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, predicted_cluster)

# 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("inferiortemporal_ADOS_RRB")))

levels(plot_data$type)[levels(plot_data$type) ==
                         "inferiortemporal_ADOS_RRB"] <- "Inferior temporal"



ggplot(plot_data, aes(x = x, y = y, group = interaction(type, predicted_cluster),
                      color = type, linetype = predicted_cluster)) +
  geom_point(aes(shape = predicted_cluster), size = 2.6, alpha = 0.4, color = custom_colors[c(4)]) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid"), guide = "none") +  # 去掉cluster的图例
  scale_shape_manual(values = c("H" = 16), guide = "none") +  # 去掉cluster的形状图例
  guides(color = guide_legend(order = 1)) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.98),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(4)]) +
  annotate("text", x = 0.06, y = 0.7, label = "*", size = 12, color = custom_colors[4]) +
  annotate("text", x = 0.275, y = 5.5, label = "Restricted and Repetitive Behavior", size = 6, color = "#474a4d")

name <- paste0("corr_ADOS_RRB_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



##################### Part 2: ados-soca ############################################################
temp <- corr[grep("SA", corr$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "predicted_cluster")
  
  plotPoint <- subset(All, predicted_cluster == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]
  
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  colnames(plotPoint)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, predicted_cluster)

unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("transversetemporal_ADOS_SA")))

levels(plot_data$type)[levels(plot_data$type) ==
                         "transversetemporal_ADOS_SA"] <- "Transverse temporal"



ggplot(plot_data, aes(x = x, y = y, group = interaction(type, predicted_cluster),
                      color = type, linetype = predicted_cluster)) +
  geom_point(aes(shape = predicted_cluster), size = 2.6, alpha = 0.4, color = custom_colors[c(11)]) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid"), guide = "none") +  # 去掉cluster的图例
  scale_shape_manual(values = c("H" = 16), guide = "none") +  # 去掉cluster的形状图例
  guides(color = guide_legend(order = 1)) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.98),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(11)]) +
  annotate("text", x = 0.07, y = 10.5, label = "**", size = 12, color = custom_colors[11]) +
  annotate("text", x = 0.9, y = .5, label = "Social Affect", size = 6, color = "#474a4d")

name <- paste0("corr_ADOS_SA_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



##################### Part 3: srs-am ################################################################
# 创建一个空的数据框来累积所有的绘图数据

temp <- corr[grep("SRS_AM", corr$name_cog), ]

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())

for (i in 1:nrow(temp)) {
  to_plot_names <- c(temp[i, "name_brain"], temp[i, "name_cog"], "predicted_cluster")
  
  plotPoint <- subset(All, predicted_cluster == temp[i, "cluster"])
  plotPoint <- plotPoint[, to_plot_names]
  
  plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
  
  colnames(plotPoint)[1:2] <- c("x","y")
  
  # 添加数据到累积数据框
  plotPoint$type <- paste0(temp[i, "name_brain"], "_", temp[i, "name_cog"])
  plot_data <- rbind(plot_data, plotPoint)
}

grouped_data <- group_by(plot_data, type, predicted_cluster)

# 给y值都做一下标准化，这样画出的线之间距离不远
# plot_data <- mutate(grouped_data, y_scaled = as.numeric(scale(y)))

unique(plot_data$type)

# 设置 type 的顺序
plot_data <- plot_data %>%
  mutate(type = factor(type, levels = c("isthmuscingulate_SRS_AM")))

levels(plot_data$type)[levels(plot_data$type) ==
                         "isthmuscingulate_SRS_AM"] <- "Isthmus cingulate"

# 去掉上2.5% 的极端值
plot_data <- plot_data %>%
  filter(y > quantile(y, 0.025) & y < quantile(y, 0.975))

ggplot(plot_data, aes(x = x, y = y, group = interaction(type, predicted_cluster),
                      color = type, linetype = predicted_cluster)) +
  geom_point(aes(shape = predicted_cluster), size = 2.6, alpha = 0.4, color = custom_colors[c(12)]) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # coord_fixed(ratio = 1) +  # 设置坐标轴比例为 1:1
  scale_linetype_manual(values = c("L" = "twodash", "H" = "solid"), guide = "none") +  # 去掉cluster的图例
  scale_shape_manual(values = c("H" = 16), guide = "none") +  # 去掉cluster的形状图例
  guides(color = guide_legend(order = 1)) +  # 强制线形颜色为灰色
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.98),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"), # 控制图例线条
        axis.text.y = element_blank(),    # 去掉y轴文字
        axis.ticks.y = element_blank(),  # 去掉y轴刻度
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = custom_colors[c(12)]) +
  annotate("text", x = 0.02, y = 15, label = "*", size = 12, color = custom_colors[12]) +
  annotate("text", x = 0.82, y = 4.8, label = "Autistic Mannerisms", size = 6, color = "#474a4d")

name <- paste0("corr_SRS_AM_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
