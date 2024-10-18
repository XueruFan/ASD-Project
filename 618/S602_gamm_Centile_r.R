# 本代码用来把随龄变化显著的脑区R2画在脑图上
# Xue-Ru Fan 04 Jan 2024 @BNU
##################################
rm(list=ls())

# load packages
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamm4", "openxlsx",
              "readr", "ggridges", "tidyr", "ggsegDefaultExtra")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect618")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect618/GAMM")
statDir <- file.path(abideDir, "Analysis/Statistic/Spect618")
resDate <- "240315"
newDate <- "240610"

stat <- read.xlsx(file.path(statDir, paste0("gamm_Centile_", newDate, ".xlsx")))

# stat_toplot <- stat[stat$平滑项p值 < 0.05,]

stat_L <- subset(stat, ClusterID == 1)
stat_H <- subset(stat, ClusterID == 2)

############ R2
to_plot_L <- data.frame(stat_L[, c(1,3)])
colnames(to_plot_L) <- c("label", "value")
to_plot_L$label <- paste0("lh_", to_plot_L$label)
to_plot_L$value <- as.numeric(to_plot_L$value)

ggseg(.data = to_plot_L, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#fce2c4", "#65318e"),  # 下限-中间-上限
                       limits = c(0, 0.36),  # 设置上下限
                       breaks = c(0, 0.12, 0.24, 0.36),  # 自定义区间
                       labels = c("0", "0.12", "0.24", "0.36")) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("gamm_Centile_R2_L_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)


to_plot_H<- data.frame(stat_H[, c(1,3)])
colnames(to_plot_H) <- c("label", "value")
to_plot_H$label <- paste0("lh_", to_plot_H$label)
to_plot_H$value <- as.numeric(to_plot_H$value)

ggseg(.data = to_plot_H, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#fce2c4", "#65318e"),  # 下限-中间-上限
                       limits = c(0, 0.36),  # 设置上下限
                       breaks = c(0, 0.12, 0.24, 0.36),  # 自定义区间
                       labels = c("0", "0.12", "0.24", "0.36")) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("gamm_Centile_R2_H_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)




############ R2 显著
stat_toplot <- stat[stat$平滑项p值 < 0.05,]

stat_L <- subset(stat_toplot, ClusterID == 1)
stat_H <- subset(stat_toplot, ClusterID == 2)

to_plot_L <- data.frame(stat_L[, c(1,3)])
colnames(to_plot_L) <- c("label", "value")
to_plot_L$label <- paste0("lh_", to_plot_L$label)
to_plot_L$value <- as.numeric(to_plot_L$value)

ggseg(.data = to_plot_L, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#fce2c4", "#65318e"),  # 下限-中间-上限
                       limits = c(0, 0.36),  # 设置上下限
                       breaks = c(0, 0.12, 0.24, 0.36),  # 自定义区间
                       labels = c("0", "0.12", "0.24", "0.36"),
                       na.value = "lightgray") +  # 将缺失数据的区域设置为白色
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("gamm_Centile_R2p_L_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)


to_plot_H<- data.frame(stat_H[, c(1,3)])
colnames(to_plot_H) <- c("label", "value")
to_plot_H$label <- paste0("lh_", to_plot_H$label)
to_plot_H$value <- as.numeric(to_plot_H$value)

ggseg(.data = to_plot_H, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#fce2c4", "#65318e"),  # 下限-中间-上限
                       limits = c(0, 0.36),  # 设置上下限
                       breaks = c(0, 0.12, 0.24, 0.36),  # 自定义区间
                       labels = c("0", "0.12", "0.24", "0.36"),
                       na.value = "lightgray") +  # 将缺失数据的区域设置为白色
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("gamm_Centile_R2p_H_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
