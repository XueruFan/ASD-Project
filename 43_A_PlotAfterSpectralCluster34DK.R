# this script is used to visualize subgroup after Spectral clustering (dev ASD Male only)
# 绘图使用谱聚类算法（这里只是用34个脑区作为分类指标）对人群进行聚类后的两类人群
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# Part 1: 画出两个分型的34个分区体积centile的概率密度分布图，png
# Part 2: 画出两个分型7个全局指标常模分（中位数+四分位距）的概率密度图，png
# Part 3: 画出两个分型的34个脑区的异常（5%和95%）流行率，png和csv
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster34DK")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster34DK")
resDate <- "240315"
newDate <- "240610"

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
cluster_result <- read.csv(file.path(resDir, name))

id_group <- c("1", "2") # cluster 1 是L， 2是H


######################### Part 1: Project 34 region centile on the brain ###########################
group1 <- (subset(cluster_result, clusterID == "1"))[, -1:-9]
group2 <- (subset(cluster_result, clusterID == "2"))[, -1:-9]
######################### cluster 1
asd_parc <- group1
median <- apply(asd_parc, 2, median, na.rm = TRUE)

label <- names(asd_parc)
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median
asd_parc_centile$label <- paste0("lh_", asd_parc_centile$label)

median <- as.numeric(asd_parc_centile$median)
label <- as.character(asd_parc_centile$label)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "#0064b5", high = "white", limits = c(min(median), max(median))) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_1_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)


######################### cluster 2
asd_parc <- group2
median <- apply(asd_parc, 2, median, na.rm = TRUE)
label <- names(asd_parc)
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median
asd_parc_centile$label <- paste0("lh_", asd_parc_centile$label)
median <- as.numeric(asd_parc_centile$median)
label <- as.character(asd_parc_centile$label)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "white", high = "#ff6347", limits = c(min(median), max(median))) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_2_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "abideDir", "dataDir", "newDate", "cluster_result", "resDir",
                     "id_group")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 2: Plot density for 7 global ########################################
# 初始化空数据框
all_data <- data.frame()
for (id in id_group){
  # id <- "1"
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", id, "))[, -1:-2]")))
  eval(parse(text = paste0("global <- group", id)))
  global <- global[which(names(global) == "GMV"):which(names(global) == "totalSA2")]
  asd_long <- gather(global, key = "Measure", value = "Centile", GMV:totalSA2, na.rm = TRUE,
                     factor_key = TRUE)
  asd_long$Group <- factor(rep(id, nrow(asd_long)))  # 添加群体标识
  all_data <- rbind(all_data, asd_long)  # 合并数据
}
  
# 修改度量标准的显示名称
all_data$Measure <- factor(all_data$Measure, levels = c("TCV", "WMV", "GMV", "sGMV", "Ventricles",
                                                        "meanCT2", "totalSA2"),
                           labels = c("Total Cerebrum Volume", "White Matter Volume",
                                      "Grey Matter Volume", "Subcortical Volume",
                                      "Ventricular Volume", "Mean Cortical Thickness",
                                      "Total Surface Area"))

ggplot(all_data, aes(x = Centile, y = Measure, fill = Measure, group = interaction(Group, Measure))) +
  geom_density_ridges(scale = 1.6, quantile_lines = TRUE, size = 1, quantiles = 2, alpha = .9) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("Grey Matter Volume" = "#4682b4",
                               "White Matter Volume" = "#dcdcdc",
                               "Subcortical Volume" = "#808080",
                               "Total Cerebrum Volume" = "#f0e68c",
                               "Mean Cortical Thickness" = "#008b8b",
                               "Total Surface Area" = "#ffa07a",
                               "Ventricular Volume" = "#ffb6c1")) +
  coord_fixed(ratio = 0.2) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_Centile_Global_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 5, units = "in", dpi = 500)


objects_to_keep <- c("plotDir", "abideDir", "dataDir", "resDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 3: Project extram percentage on the brain ###########################
# 34个脑区的异常流行率, 对于每个脑区来说有多少人的体积常模分属于极端值（位于thr以外，这里thr取5%）

thr <- 0.05

# cluster 1
# 小于5% !!!
abno <- apply(group1 <= thr, 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g1 <- data.frame(label)
abno_g1$perc <- as.numeric(abno)
# 保存一个排序
abno_g1_sorted <- abno_g1 %>% arrange(desc(perc))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_1_Abno5_Rank_", newDate, ".csv")
write.csv(abno_g1_sorted, file.path(resDir, name), row.names = F)
# 画脑图
abno_g1$label <- as.character(paste0("lh_", abno_g1$label))
ggseg(.data = abno_g1, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_1_Abno5_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# 大于95%
abno <- apply(group1 >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g1 <- data.frame(label)
abno_g1$perc <- as.numeric(abno)
abno_g1$label <- as.character(paste0("lh_", abno_g1$label))
ggseg(.data = abno_g1, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_1_Abno95_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

#  cluster 2
# 大于95% !!!
abno <- apply(group2 >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g2 <- data.frame(label)
abno_g2$perc <- as.numeric(abno)
# 保存一个排序
abno_g2_sorted <- abno_g2 %>% arrange(desc(perc))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_2_Abno95_Rank_", newDate, ".csv")
write.csv(abno_g2_sorted, file.path(resDir, name), row.names = F)
# 画脑图
abno_g2$label <- as.character(paste0("lh_", abno_g2$label))
ggseg(.data = abno_g2, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_2_Abno95_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# 小于5%
abno <- apply(group2 <= thr, 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g2 <- data.frame(label)
abno_g2$perc <- as.numeric(abno)
abno_g2$label <- as.character(paste0("lh_", abno_g2$label))
ggseg(.data = abno_g2, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_2_Abno5_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
