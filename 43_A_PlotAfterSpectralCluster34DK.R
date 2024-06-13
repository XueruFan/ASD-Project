# this script is used to visualize subgroup after Spectral clustering (dev ASD Male only)
# Xue-Ru Fan 25 April 2023 @BNU

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster")
resDate <- "240315"
newDate <- "240610"

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
cluster_result <- read.csv(file.path(resDir, name))

id_group <- c("1", "2") # cluster 1 是L， 2是H


######################### Part 1: Project 34 region centile on the brain ###########################

group1 <- (subset(cluster_result, clusterID == "1"))[, -1:-2]
group2 <- (subset(cluster_result, clusterID == "2"))[, -1:-2]

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

objects_to_keep <- c("plotDir", "abideDir", "dataDir", "newDate", "cluster_result", "group1",
                     "group2", "resDir")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 4: Project extram percentage on the brain ###########################
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
