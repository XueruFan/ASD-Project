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

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_", newDate, ".csv")
cluster_result <- read.csv(file.path(resDir, name))

id_group <- c("1", "2") # cluster 1 是L， 2是H

######################### Part 1: Plot density for 34 regions ######################################

for (id in id_group){
  # id <- "1"
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", id, "))[, -1:-2]")))
  eval(parse(text = paste0("parc <- group", id)))
  parc <- parc[which(names(parc) == "bankssts"):which(names(parc) == "insula")]
  parc_long <- gather(parc, key = "Measure", value = "Centile", bankssts:insula, na.rm = TRUE,
                      factor_key = TRUE)
  centile <- as.numeric(parc_long$Centile)
  measure <- parc_long$Measure
  measure2 <- reorder(measure, centile, FUN = median) # from small to large
  ggplot(parc_long, aes(x = centile, y = measure2, fill = measure2, alpha = .5)) +
    geom_density_ridges(scale = 1.5, quantile_lines = TRUE, size = .3, quantiles = 4) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_fill_manual(values = rep("transparent", 34)) +
    coord_fixed(ratio = 0.17) + 
    xlab("") +
    ylab("") +
    theme_ridges() + 
    theme(legend.position = "none", # without legend
          axis.text.y = element_text(size = 6, face = "bold"),
          axis.text.x = element_text(size = 6, face = "bold"))
  name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_", id, "_Centile_Regional_", newDate, ".png")
  ggsave(file.path(plotDir, name), width = 4.6, height = 6, units = "in", dpi = 500)
}

objects_to_keep <- c("plotDir", "abideDir", "dataDir", "resDir", "newDate", "cluster_result", "id_group")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 2: Plot density for 7 global ########################################

for (id in id_group){
  # id <- "1"
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", id, "))[, -1:-2]")))
  eval(parse(text = paste0("global <- group", id)))
  global <- global[which(names(global) == "GMV"):which(names(global) == "totalSA2")]
  asd_long <- gather(global, key = "Measure", value = "Centile", GMV:totalSA2, na.rm = TRUE,
                     factor_key = TRUE)
  centile <- as.numeric(asd_long$Centile)
  measure <- asd_long$Measure
  measure2 <- reorder(measure, centile, FUN = median) # from small to large
  # modify name
  levels(measure2)[levels(measure2) == "GMV"] <- "Grey Matter Volume"
  levels(measure2)[levels(measure2) == "WMV"] <- "White Matter Volume"
  levels(measure2)[levels(measure2) == "sGMV"] <- "Subcortical Volume"
  levels(measure2)[levels(measure2) == "TCV"] <- "Total Cerebrum Volume"
  levels(measure2)[levels(measure2) == "meanCT2"] <- "Mean Cortical Thickness"
  levels(measure2)[levels(measure2) == "Ventricles"] <- "Ventricular Volume"
  levels(measure2)[levels(measure2) == "totalSA2"] <- "Total Surface Area"

  ggplot(asd_long, aes(x = centile, y = measure2, fill = measure2)) +
    geom_density_ridges(scale = 1.8, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
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
  name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_", id, "_Centile_Global_", newDate, ".png")
  ggsave(file.path(plotDir, name), width = 6, height = 5, units = "in", dpi = 500)
  
}
      
objects_to_keep <- c("plotDir", "abideDir", "dataDir", "resDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 3: Project 34 region centile on the brain ###########################

group1 <- (subset(cluster_result, clusterID == "1"))[, -1:-2]
group2 <- (subset(cluster_result, clusterID == "2"))[, -1:-2]

######################### cluster 1
asd_parc <- group1
median <- apply(asd_parc, 2, median, na.rm = TRUE)

label <- names(asd_parc)[8:41]
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median[8:41]
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

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_1_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)


######################### cluster 2
asd_parc <- group2
median <- apply(asd_parc, 2, median, na.rm = TRUE)
label <- names(asd_parc)[8:41]
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median[8:41]
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

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_2_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "abideDir", "dataDir", "newDate", "cluster_result", "group1",
                     "group2", "resDir")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 4: Project extram percentage on the brain ###########################
# 34个脑区的异常流行率, 对于每个脑区来说有多少人的体积常模分属于极端值（位于thr以外，这里thr取5%）

thr <- 0.05

# cluster 1
# 小于5% !!!
abno <- apply(group1[, 8:41] <= thr, 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g1 <- data.frame(label)
abno_g1$perc <- as.numeric(abno)
# 保存一个排序
abno_g1_sorted <- abno_g1 %>% arrange(desc(perc))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_1_Abno5_Rank_", newDate, ".csv")
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_1_Abno5_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# 大于95%
abno <- apply(group1[, 8:41] >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_1_Abno95_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

#  cluster 2
# 大于95% !!!
abno <- apply(group2[, 8:41] >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g2 <- data.frame(label)
abno_g2$perc <- as.numeric(abno)
# 保存一个排序
abno_g2_sorted <- abno_g2 %>% arrange(desc(perc))
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_2_Abno95_Rank_", newDate, ".csv")
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_2_Abno95_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# 小于5%
abno <- apply(group2[, 8:41] <= thr, 2, function(x) sum(x) / length(x) * 100)
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
name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_2_Abno5_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
