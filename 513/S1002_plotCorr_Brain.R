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

corr_l <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_L_", newDate, ".csv")))
corr_h <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_H_", newDate, ".csv")))

corr_l$name_brain <- gsub("_centile", "", corr_l$name_brain)
corr_h$name_brain <- gsub("_centile", "", corr_h$name_brain)
corr_l$label <- paste0("lh_", corr_l$name_brain)
corr_h$label <- paste0("rh_", corr_h$name_brain)

corr <- rbind(corr_l, corr_h)


corr_soca <- corr[grep("SOCA", corr$name_cog), ]
corr_rrb<- corr[grep("RRB", corr$name_cog), ]
corr_fiq <- corr[grep("FIQ", corr$name_cog), ]


############ soca
to_plot <- corr_soca[, c(6,3)]
colnames(to_plot) <- c("label", "value")
to_plot$value <- as.numeric(to_plot$value)

ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#c9171e"),  # 下限-中间-上限
                       limits = c(0.1, 0.25),  # 设置上下限
                       breaks = c(0.1, 0.25),  # 自定义区间
                       labels = c("0.18", "0.25"),
                       na.value = "#f5f5f5") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("corr_ADOS_2_SOCAF_brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


############ RRB
to_plot <- corr_rrb[, c(6,3)]
colnames(to_plot) <- c("label", "value")
to_plot$value <- as.numeric(to_plot$value)

ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradientn(colors = c("#ffffff", "#c9171e"),  # 下限-中间-上限
                       limits = c(0.1, 0.25),  # 设置上下限
                       breaks = c(0.1, 0.25),  # 自定义区间
                       labels = c("0.18", "0.25"),
                       na.value = "#f5f5f5") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("corr_ADOS_2_RRB_brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


############ fiq
# to_plot <- corr_fiq[, c(6,3)]
# colnames(to_plot) <- c("label", "value")
# to_plot$value <- as.numeric(to_plot$value)
# 
# ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
#       position = "stacked", size = 1.2) +
#   theme_void() +
#   theme(legend.title = element_blank(), legend.position = "bottom",
#         legend.key.width = unit(1, "cm")) +
#   scale_fill_gradientn(colors = c("#165e83", "#ffffff", "#c9171e"),
#                        limits = c(-0.3, 0.3),  # 设置上下限
#                        breaks = c(-0.3, -0.15, 0, 0.15, 0.3),  # 自定义区间
#                        labels = c("-0.3", "-0.15", "0", "0.15", "0.3"),
#                        na.value = "#f5f5f5") +
#   guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
# 
# name <- paste0("corr_FIQ_brain_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


# ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
#       position = "stacked", size = 1.2) +
#   theme_void() +
#   theme(legend.title = element_blank(), 
#         legend.position = "right",  # 将colorbar放在右侧
#         legend.key.width = unit(1, "cm")) +
#   scale_fill_gradientn(colors = c("#ffffff", "#c9171e"),  # 下限-中间-上限
#                        limits = c(0, 0.3),  # 设置上下限
#                        breaks = c(0, 0.1, 0.2,  0.3),  # 自定义区间
#                        labels = c("0", "0.1", "0.2", "0.3"),
#                        na.value = "#f5f5f5") +
#   guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, 
#                                 ticks = FALSE, direction = "vertical"))  # 竖着显示colorbar
# name <- paste0("corr_colorbar_", resDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
