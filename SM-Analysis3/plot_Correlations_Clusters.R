# Plot cluster-wise correlation both in ABIDE and CABIC
# Male, ASD, <13 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################
rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss", "ggsegDefaultExtra")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
cabicDir <- 'E:/PhDproject/CABIC'

phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
resuDir <- file.path(cabicDir, "result/inde/Spect-13/Corr")
plotDir <- file.path(cabicDir, "result/inde/Spect-13/Plot/Corr")
oldDate <- "240315"
newDate <- "240610"
resDate <- "240928"

# ABDIE
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", oldDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "GMV")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")
ABIDE <- merge(cluster, pheno, by = "participant", ABIDE.x = TRUE)
ABIDE[ABIDE < 0] <- NA

ABIDE$clusterID[ABIDE$clusterID == 2] <- "H"
ABIDE$clusterID[ABIDE$clusterID == 1] <- "L"
ABIDE <- ABIDE %>% rename(ADOS_RRB = ADOS_2_RRB)
ABIDE <- ABIDE %>% rename(ADOS_TOTAL = ADOS_2_TOTAL)
ABIDE <- ABIDE %>% rename(ADOS_SA = ADOS_2_SOCAFFECT)
ABIDE <- ABIDE %>% rename(SRS_AM = SRS_MANNERISMS_RAW)


# CABIC
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("Cluster_", resDate, ".csv")
cluster <- read.csv(file.path('E:/PhDproject/CABIC/result/inde/Spect-13', name))
CABIC <- merge(cluster, pheno, by = "participant")
CABIC$clusterID[CABIC$clusterID == 1] <- "H"
CABIC$clusterID[CABIC$clusterID == 2] <- "L"
names(CABIC)[which(names(CABIC) == "bankssts"):which(names(CABIC) == "insula")] <- 
  paste0(names(CABIC)[which(names(CABIC) == "bankssts"):which(names(CABIC) == "insula")], "_centile")



############################ Group H ###############################################################

corr_abide_h <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_H_", newDate, ".csv")))

corr_abide_h$name_cog[corr_abide_h$name_cog == "ADOS_2_RRB"] <- "ADOS_RRB"
corr_abide_h$name_cog[corr_abide_h$name_cog == "ADOS_2_SOCAFFECT"] <- "ADOS_SA"
corr_abide_h$name_cog[corr_abide_h$name_cog == "ADOS_2_TOTAL"] <- "ADOS_TOTAL"
corr_abide_h$name_cog[corr_abide_h$name_cog == "ADOS_2_RRB"] <- "ADOS_RRB"
corr_abide_h$name_cog[corr_abide_h$name_cog == "SRS_AWARENESS_RAW"] <- "SRS_SA"
corr_abide_h$name_cog[corr_abide_h$name_cog == "SRS_MANNERISMS_RAW"] <- "SRS_AM"
corr_abide_h$name_cog[corr_abide_h$name_cog == "SRS_COMMUNICATION_RAW"] <- "SRS_SCOM"
corr_abide_h$name_cog[corr_abide_h$name_cog == "SRS_MOTIVATION_RAW"] <- "SRS_SM"
corr_abide_h$name_cog[corr_abide_h$name_cog == "SRS_TOTAL_RAW"] <- "SRS_TOTAL"
corr_abide_h$name_cog[corr_abide_h$name_cog == "ADI_R_SOCIAL_TOTAL_A"] <- "ADIR_SOCI"

corr_cabic_h <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_H_", resDate, ".csv")))
corr_cabic_h$name_brain <- paste0(corr_cabic_h$name_brain, "_centile")

both <- inner_join(corr_cabic_h, corr_abide_h, by = c("name_brain", "name_cog"),
                   suffix = c(".cabic", ".abide"))
both <- subset(both, sign(coef.cabic) == sign(coef.abide))
name <- paste0("corr_Part4_Both_H_ABIDE_indeCABIC_", newDate, ".csv")
write.csv(both, file.path(statiDir, name), row.names = F)

both$clusterID <- "H"

# # # 定义自定义颜色向量
# custom_colors <- c("#c85554", "#69821b", "#93ca76", "#9CB0C3", "#43676b", "#ee7800",
#                    "#dccb18", "#8491c3", "#00a381", "#cc7eb1", "#e09e87")


i <- 1

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
to_plot_names <- c(both[i, "name_brain"], both[i, "name_cog"], "clusterID")

plotPoint.abide <- subset(ABIDE, clusterID == both[i, "clusterID"])
plotPoint.abide <- plotPoint.abide[, to_plot_names]
plotPoint.abide <- plotPoint.abide[!is.na(plotPoint.abide[[2]]), ]
colnames(plotPoint.abide)[1:2] <- c("x","y")
plotPoint.cabic <- subset(CABIC, clusterID == both[i, "clusterID"])
plotPoint.cabic <- plotPoint.cabic[, to_plot_names]
plotPoint.cabic <- plotPoint.cabic[!is.na(plotPoint.cabic[[2]]), ]
colnames(plotPoint.cabic)[1:2] <- c("x","y")

plotPoint.abide$Data <- "ABIDE"
plotPoint.cabic$Data <- "CABIC"

plotPoint <- rbind(plotPoint.abide[, c("x","y","Data")], plotPoint.cabic[, c("x","y","Data")])
plotPoint$Data <- as.factor(plotPoint$Data)

ggplot(plotPoint, aes(x = x, y = y, color = Data, shape = Data, size = Data)) +
  geom_point(alpha = 0.3) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24)) +# 自定义散点形状和大小
  scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
  scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
  # 颜色自定义
  scale_color_manual(values = c("ABIDE" = "#d9ca39", "CABIC" = "#e29135")) +
  # 主题设置
  theme_cowplot() +
  # 设置X轴标题
  xlab("") +
  # xlab("Transverse temporal") +
  ylab("") +
  theme(legend.position = c(0.02, 0.95),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"),
        axis.text.x = element_text(size = 12))
  # annotate("text", x = 1, y = 9, label = "*", size = 12, color = "black") +
  # annotate("text", x = 1, y = 15, label = "**", size = 12, color = "black") +
  # annotate("text", x = 0.19, y = 21, label = "ADOS  Social Affect", size = 7, color = "black")

name <- paste0("corr_Part4_ADOS_SA_transversetemporal_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



i <- 2

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
to_plot_names <- c(both[i, "name_brain"], both[i, "name_cog"], "clusterID")

plotPoint.abide <- subset(ABIDE, clusterID == both[i, "clusterID"])
plotPoint.abide <- plotPoint.abide[, to_plot_names]
plotPoint.abide <- plotPoint.abide[!is.na(plotPoint.abide[[2]]), ]
colnames(plotPoint.abide)[1:2] <- c("x","y")
plotPoint.cabic <- subset(CABIC, clusterID == both[i, "clusterID"])
plotPoint.cabic <- plotPoint.cabic[, to_plot_names]
plotPoint.cabic <- plotPoint.cabic[!is.na(plotPoint.cabic[[2]]), ]
colnames(plotPoint.cabic)[1:2] <- c("x","y")

plotPoint.abide$Data <- "ABIDE"
plotPoint.cabic$Data <- "CABIC"

plotPoint <- rbind(plotPoint.abide[, c("x","y","Data")], plotPoint.cabic[, c("x","y","Data")])
plotPoint$Data <- as.factor(plotPoint$Data)

ggplot(plotPoint, aes(x = x, y = y, color = Data, shape = Data, size = Data)) +
  geom_point(alpha = 0.3) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +# 自定义散点形状和大小
  scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
  scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
  # 颜色自定义
  scale_color_manual(values = c("ABIDE" = "#d9ca39", "CABIC" = "#e29135")) +
  # 主题设置
  theme_cowplot() +
  # 设置X轴标题
  xlab("") +
  ylab("") +
  theme(legend.position = c(0.02, 0.95),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"),
        axis.text.x = element_text(size = 12))
  # annotate("text", x = 1, y = 3.1, label = "*", size = 12, color = "black") +
  # annotate("text", x = 1, y = 1.5, label = "*", size = 12, color = "black") +
  # annotate("text", x = 0.11, y = 7, label = "ADOS  RRB", size = 7, color = "black")


name <- paste0("corr_Part4_ADOS_RRB_inferiortemporal_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



i <- 3

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
to_plot_names <- c(both[i, "name_brain"], both[i, "name_cog"], "clusterID")

plotPoint.abide <- subset(ABIDE, clusterID == both[i, "clusterID"])
plotPoint.abide <- plotPoint.abide[, to_plot_names]
plotPoint.abide <- plotPoint.abide[!is.na(plotPoint.abide[[2]]), ]
colnames(plotPoint.abide)[1:2] <- c("x","y")
plotPoint.cabic <- subset(CABIC, clusterID == both[i, "clusterID"])
plotPoint.cabic <- plotPoint.cabic[, to_plot_names]
plotPoint.cabic <- plotPoint.cabic[!is.na(plotPoint.cabic[[2]]), ]
colnames(plotPoint.cabic)[1:2] <- c("x","y")

plotPoint.abide$Data <- "ABIDE"
plotPoint.cabic$Data <- "CABIC"

plotPoint <- rbind(plotPoint.abide[, c("x","y","Data")], plotPoint.cabic[, c("x","y","Data")])
plotPoint$Data <- as.factor(plotPoint$Data)

ggplot(plotPoint, aes(x = x, y = y, color = Data, shape = Data, size = Data)) +
  geom_point(alpha = 0.3) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 44), breaks = c(0, 11, 22, 33, 44)) +# 自定义散点形状和大小
  scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
  scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
  # 颜色自定义
  scale_color_manual(values = c("ABIDE" = "#d9ca39", "CABIC" = "#e29135")) +
  # 主题设置
  theme_cowplot() +
  # 设置X轴标题
  xlab("") +
  ylab("") +
  theme(legend.position = c(0.02, 0.95),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"),
        axis.text.x = element_text(size = 12))
  # annotate("text", x = 1, y = 17, label = "*", size = 12, color = "black") +
  # annotate("text", x = 1, y = 11.5, label = "*", size = 12, color = "black") +
  # annotate("text", x = 0.11, y = 39, label = "ADOS  Total", size = 7, color = "black")


name <- paste0("corr_Part4_ADOS_Total_transversetemporal_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



i <- 4

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
to_plot_names <- c(both[i, "name_brain"], both[i, "name_cog"], "clusterID")

plotPoint.abide <- subset(ABIDE, clusterID == both[i, "clusterID"])
plotPoint.abide <- plotPoint.abide[, to_plot_names]
plotPoint.abide <- plotPoint.abide[!is.na(plotPoint.abide[[2]]), ]
colnames(plotPoint.abide)[1:2] <- c("x","y")
plotPoint.cabic <- subset(CABIC, clusterID == both[i, "clusterID"])
plotPoint.cabic <- plotPoint.cabic[, to_plot_names]
plotPoint.cabic <- plotPoint.cabic[!is.na(plotPoint.cabic[[2]]), ]
colnames(plotPoint.cabic)[1:2] <- c("x","y")

plotPoint.abide$Data <- "ABIDE"
plotPoint.cabic$Data <- "CABIC"

plotPoint <- rbind(plotPoint.abide[, c("x","y","Data")], plotPoint.cabic[, c("x","y","Data")])
plotPoint$Data <- as.factor(plotPoint$Data)

ggplot(plotPoint, aes(x = x, y = y, color = Data, shape = Data, size = Data)) +
  geom_point(alpha = 0.3) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 44), breaks = c(0, 11, 22, 33, 44)) +# 自定义散点形状和大小
  scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
  scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
  # 颜色自定义
  scale_color_manual(values = c("ABIDE" = "#d9ca39", "CABIC" = "#e29135")) +
  # 主题设置
  theme_cowplot() +
  # 设置X轴标题
  xlab("") +
  ylab("") +
  theme(legend.position = c(0.02, 0.95),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"),
        axis.text.x = element_text(size = 12))
  # annotate("text", x = 1, y = 13, label = "*", size = 12, color = "black") +
  # annotate("text", x = 0.99, y = 17.5, label = "**", size = 12, color = "black") +
  # annotate("text", x = 0.26, y = 39, label = "SRS  Autistic Mannerisms", size = 7, color = "black")


name <- paste0("corr_Part4_SRS_AM_isthmuscingulate_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)


############################ Part 2: Group L #######################################################

corr_abide_l <- read.csv(file.path(statiDir, paste0("corr_Part4_Select_L_", newDate, ".csv")))

corr_abide_l$name_cog[corr_abide_l$name_cog == "ADOS_2_RRB"] <- "ADOS_RRB"
corr_abide_l$name_cog[corr_abide_l$name_cog == "ADOS_2_SOCAFFECT"] <- "ADOS_SA"
corr_abide_l$name_cog[corr_abide_l$name_cog == "ADOS_2_TOTAL"] <- "ADOS_TOTAL"
corr_abide_l$name_cog[corr_abide_l$name_cog == "ADOS_2_RRB"] <- "ADOS_RRB"
corr_abide_l$name_cog[corr_abide_l$name_cog == "SRS_AWARENESS_RAW"] <- "SRS_SA"
corr_abide_l$name_cog[corr_abide_l$name_cog == "SRS_MANNERISMS_RAW"] <- "SRS_AM"
corr_abide_l$name_cog[corr_abide_l$name_cog == "SRS_COMMUNICATION_RAW"] <- "SRS_SCOM"
corr_abide_l$name_cog[corr_abide_l$name_cog == "SRS_MOTIVATION_RAW"] <- "SRS_SM"
corr_abide_l$name_cog[corr_abide_l$name_cog == "SRS_TOTAL_RAW"] <- "SRS_TOTAL"
corr_abide_l$name_cog[corr_abide_l$name_cog == "ADI_R_SOCIAL_TOTAL_A"] <- "ADIR_SOCI"

corr_cabic_l <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_L_", resDate, ".csv")))
corr_cabic_l$name_brain <- paste0(corr_cabic_l$name_brain, "_centile")

both <- inner_join(corr_cabic_l, corr_abide_l, by = c("name_brain", "name_cog"),
                   suffix = c(".cabic", ".abide"))
both <- subset(both, sign(coef.cabic) == sign(coef.abide))
name <- paste0("corr_Part4_Both_L_ABIDE_indeCABIC_", newDate, ".csv")
write.csv(both, file.path(statiDir, name), row.names = F)

both$clusterID <- "L"


i <- 1

plot_data <- data.frame(x = numeric(), y = numeric(), type = character())
to_plot_names <- c(both[i, "name_brain"], both[i, "name_cog"], "clusterID")

plotPoint.abide <- subset(ABIDE, clusterID == both[i, "clusterID"])
plotPoint.abide <- plotPoint.abide[, to_plot_names]
plotPoint.abide <- plotPoint.abide[!is.na(plotPoint.abide[[2]]), ]
colnames(plotPoint.abide)[1:2] <- c("x","y")
plotPoint.cabic <- subset(CABIC, clusterID == both[i, "clusterID"])
plotPoint.cabic <- plotPoint.cabic[, to_plot_names]
plotPoint.cabic <- plotPoint.cabic[!is.na(plotPoint.cabic[[2]]), ]
colnames(plotPoint.cabic)[1:2] <- c("x","y")

plotPoint.abide$Data <- "ABIDE"
plotPoint.cabic$Data <- "CABIC"

plotPoint <- rbind(plotPoint.abide[, c("x","y","Data")], plotPoint.cabic[, c("x","y","Data")])
plotPoint$Data <- as.factor(plotPoint$Data)

ggplot(plotPoint, aes(x = x, y = y, color = Data, shape = Data, size = Data)) +
  geom_point(alpha = 0.3) +  # 添加散点
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +# 自定义散点形状和大小
  scale_shape_manual(values = c("ABIDE" = 16, "CABIC" = 16)) +  # 自定义形状
  scale_size_manual(values = c("ABIDE" = 2.5, "CABIC" = 2.5)) +  # 自定义大小
  # 颜色自定义
  scale_color_manual(values = c("ABIDE" = "#ade084", "CABIC" = "#018a67")) +
  # 主题设置
  theme_cowplot() +
  # 设置X轴标题
  xlab("") +
  ylab("") +
  theme(legend.position = c(0.02, 0.95),
        legend.title = element_blank(),
        legend.key.width  = unit(1, "cm"),
        axis.text.x = element_text(size = 12))

name <- paste0("corr_Part4_ADOS_RRB_parstriangularis_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)



################画脑区示意图
regions <- c("lh_parstriangularis")
to_plot <- data.frame(label = regions, value = 1)

ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", view = "lateral", size = 1.2) +
  theme_void() +
  theme(legend.position = "none") +  # 去掉图例
  scale_fill_gradientn(colors = c("#d26b66"), na.value = "#f5f5f5")

name <- paste0("parstriangularis_on_brain_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 5, units = "in", dpi = 500)
