
rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss", "ggsegDefaultExtra")
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

corr_l$label <- paste0("lh_", corr_l$name_brain)
corr_h$label <- paste0("rh_", corr_h$name_brain)

corr <- rbind(corr_l, corr_h)
corr <- corr[c(23,25,31),]


############ soca
to_plot <- corr[1, c(6,3)]
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

name <- paste0("corr_ADOS_SA_brain_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


############ RRB
to_plot <- corr[2, c(6,3)]
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

name <- paste0("corr_ADOS_RRB_brain_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


############ AM
to_plot <- corr[3, c(6,3)]
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

name <- paste0("corr_SRS_AM_brain_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
