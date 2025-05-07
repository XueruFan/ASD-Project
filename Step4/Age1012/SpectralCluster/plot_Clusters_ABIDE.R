# Plot OoS centile scores of each cluster for ABIDE data
# Xue-Ru Fan 13 march 2024 @BNU
###################################################
# Part 1: Project the OoS centile (median) for 34 regional volumes
# Part 2: Plot the probability density of OoS centile (median) for 7 global measures
# Part 3: Plot the probability density of OoS centile (median) for 34 regional volumes
# Part 4: Project extreme percentage (<2.5% and >97.5%) of OoS centils on the brain
##################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Spect112")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect112")
resDate <- "240315"
newDate <- "240610"

name <- paste0("Cluster_", newDate, ".csv")
cluster_result <- read.csv(file.path(resDir, name))


##### Part 1: Project the OoS centile (median) for 34 regional volumes in ASD males on brain #######
group1 <- (subset(cluster_result, clusterID == "1"))[, -1:-9]
group2 <- (subset(cluster_result, clusterID == "2"))[, -1:-9]

######################### cluster 1
asd_parc <- group1
median <- apply(asd_parc, 2, median, na.rm = TRUE)

label <- names(asd_parc)
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median
asd_parc_centile$label <- paste0("lh_", asd_parc_centile$label)
name <- paste0("Cluster_1_Centile_", newDate, ".xlsx")
write.xlsx(asd_parc_centile, file.path(resDir, name))

median <- as.numeric(asd_parc_centile$median)
label <- as.character(asd_parc_centile$label)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  # scale_fill_gradient(low = "#7c9d97", high = "white", limits = c(min(median), max(median))) +
  scale_fill_gradient(low = "#5378ac", high = "white", limits = c(0.24, 0.44),
                      breaks = seq(0.24, 0.44, by = 0.05)) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("Cluster_1_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)


######################### cluster 2
asd_parc <- group2
median <- apply(asd_parc, 2, median, na.rm = TRUE)
label <- names(asd_parc)
asd_parc_centile <- data.frame(label)
asd_parc_centile$median <- median
asd_parc_centile$label <- paste0("lh_", asd_parc_centile$label)
name <- paste0("Cluster_2_Centile_", newDate, ".xlsx")
write.xlsx(asd_parc_centile, file.path(resDir, name))
median <- as.numeric(asd_parc_centile$median)
label <- as.character(asd_parc_centile$label)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "white", high = "#d26b66",
                       limits = c(0.54, 0.79),
                       breaks = seq(0.54, 0.79, by = 0.05)) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("Cluster_2_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "abideDir", "newDate", "cluster_result", "resDir", "id_group")
rm(list = (setdiff(ls(), objects_to_keep)))


##### Part 2: Plot the probability density of OoS centile (median) for 7 global measures ###########
# 初始化空数据框
all_data <- data.frame()
id_group <- c("1", "2")

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
all_data$Measure <- factor(all_data$Measure, levels = rev(c("TCV", "WMV", "GMV", "sGMV", "Ventricles",
                                                        "meanCT2", "totalSA2")),
                           # labels = rev(c("脑总容积", "白质总体积",
                           #            "皮层灰质总体积", "皮层下灰质总体积",
                           #            "脑脊液总体积", "平均皮层厚度",
                           #            "皮层总表面积")))
                           labels = rev(c("TCV", "WMV",
                                      "GMV", "sGMV",
                                      "CSF", "mCT",
                                      "tSA")))

name_global <- paste0("Cluster_Centile_Global_", newDate, ".png")
# CairoPNG(file.path(plotDir, name_global), width = 6, height = 5, units = "in", dpi = 500)

ggplot(all_data, aes(x = Centile, y = Measure, fill = Group, group = interaction(Group, Measure))) +
  geom_density_ridges(scale = 1.6, quantile_lines = TRUE, size = 1, quantiles = 2, alpha = .9) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("1" = "#b4d4c7", "2" = "#facaae")) +  # 此处为两组人设置不同颜色
  coord_fixed(ratio = 0.2) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  # theme(text = element_text(family = "STSong"),
  theme(
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name_global), width = 6, height = 5, units = "in", dpi = 500)


objects_to_keep <- c("plotDir", "abideDir", "resDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


##### Part 3: Plot the probability density of OoS centile (median) for 34 regional volumes #########
# 初始化空数据框
all_data <- data.frame()
id_group <- c("1", "2") # cluster 1 是L， 2是H

for (id in id_group){
  # id <- "1"
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", id, "))[, -1:-2]")))
  eval(parse(text = paste0("regional <- group", id)))
  regional <- regional[which(names(regional) == "bankssts"):which(names(regional) == "insula")]
  asd_long <- gather(regional, key = "Measure", value = "Centile", bankssts:insula, na.rm = TRUE,
                     factor_key = TRUE)
  asd_long$Group <- factor(rep(id, nrow(asd_long)))  # 添加群体标识
  all_data <- rbind(all_data, asd_long)  # 合并数据
}

# 修改度量标准的显示名称
# cn_labels <- rev(c("额上回","额中回喙部","额中回尾部",
#                         "额下回盖部","额下回三角部","额下回眶部","额极",
#                         "外侧眶额","内侧眶额",
#                         "前扣带喙部","前扣带尾部",
#                         "中央前回","中央旁小叶","中央后回",
#                         "缘上回","后扣带","扣带回峡部",
#                         "楔前叶","顶上皮层","顶下皮层",
#                         "颞横皮层","颞上沟后侧",
#                         "颞上回","颞中回","颞下回",
#                         "梭状回","海马旁回","内嗅皮层","颞极",
#                         "外侧枕叶","舌回","距状沟周围皮层","楔叶皮层",
#                         "脑岛"))
en_labels <- rev(c("superiorfrontal", "rostralmiddlefrontal", "caudalmiddlefrontal",
                   "parsopercularis", "parstriangularis", "parsorbitalis", "frontalpole",
                   "lateralorbitofrontal", "medialorbitofrontal",
                   "rostralanteriorcingulate", "caudalanteriorcingulate",
                   "precentral", "paracentral", "postcentral",
                   "supramarginal", "posteriorcingulate", "isthmuscingulate",
                   "precuneus", "superiorparietal", "inferiorparietal",
                   "transversetemporal", "bankssts",
                   "superiortemporal", "middletemporal", "inferiortemporal",
                   "fusiform", "parahippocampal", "entorhinal", "temporalpole",
                   "lateraloccipital", "lingual", "pericalcarine", "cuneus",
                   "insula"))
cn_labels <- rev(c("Superior frontal", "Rostral middle frontal", "Caudal middle frontal",
                   "Pars opercularis", "Pars triangularis", "Pars orbitalis", "Frontal pole",
                   "Lateral orbital frontal", "Medial orbital frontal",
                   "Rostral anterior cingulate", "Caudal anterior cingulate",
                   "Precentral", "Paracentral", "Postcentral",
                   "Supramarginal", "Posterior cingulate", "Isthmus cingulate",
                   "Precuneus", "Superior parietal", "Inferior parietal",
                   "Transverse temporal", "Banks superior temporal",
                   "Superior temporal", "Middle temporal", "Inferior temporal",
                   "Fusiform", "Parahippocampal", "Entorhinal", "Temporal pole",
                   "Lateral occipital", "Lingual", "Pericalcarine", "Cuneus",
                   "Insula"))


all_data$Measure <- factor(all_data$Measure, levels = en_labels, labels = cn_labels)

name_regional <- paste0("Cluster_Centile_Regional_", newDate, ".png")
# CairoPNG(file.path(plotDir, name_regional), width = 4, height = 8, units = "in", dpi = 500)

# color <- c("#95483f","#a22041","#b7282e",
#            "#d3381c","#ea5506","#ec6800","#f08300",
#            "#e17b34","#bf783a",
#            "#c39143","#b48a76",
#            "#c37854","#bfa46f","#a8c97f",
#            "#38b48b","#7c9d97","#316745",
#            "#00552e","#2f5d50","#006e54",
#            "#8aaee6","#0091ff",
#            "#2ca9e1","#0075c2","#165e83",
#            "#192f60","#19448e","#4169e1","#778899",
#            "#4a488e","#674598","#522f60","#460e44",
#            "white")

ggplot(all_data, aes(x = Centile, y = Measure, fill = Group, group = interaction(Group, Measure))) +
  # geom_density_ridges(scale = 1.7, quantile_lines = TRUE, size = 0.5, quantiles = 2, alpha = .94) +
  geom_density_ridges(scale = 2, quantile_lines = TRUE, size = 0.75, quantiles = 2, alpha = .9) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("1" = "#b4d4c7", "2" = "#facaae")) +
  # coord_fixed(ratio = 0.15) + 
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6))

# dev.off()
# ggsave(file.path(plotDir, name_regional), width = 4, height = 8, units = "in", dpi = 500)
ggsave(file.path(plotDir, name_regional), width = 5, height = 8, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "abideDir", "resDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


##### Part 4: Part 4: Project extreme percentage (<2.5% and >97.5%) of OoS centils on the brain ####
# 34个脑区的异常流行率, 对于每个脑区来说有多少人的体积常模分属于极端值
group1 <- (subset(cluster_result, clusterID == "1"))[, -1:-9]
group2 <- (subset(cluster_result, clusterID == "2"))[, -1:-9]

label_map <- c("额上回" = "Superior frontal",
               "额中回喙部" = "Rostral middle frontal",
               "额中回尾部" = "Caudal middle frontal",
               "额下回盖部" = "Pars opercularis",
               "额下回三角部" = "Pars triangularis",
               "额下回眶部" = "Pars orbitalis",
               "额极" = "Frontal pole",
               "外侧眶额" = "Lateral orbital frontal",
               "内侧眶额" = "Medial orbital frontal",
               "前扣带喙部" = "Rostral anterior cingulate",
               "前扣带尾部" = "Caudal anterior cingulate",
               "中央前回" = "Precentral",
               "中央旁小叶" = "Paracentral",
               "中央后回" = "Postcentral",
               "缘上回" = "Supramarginal",
               "后扣带" = "Posterior cingulate",
               "扣带回峡部" = "Isthmus cingulate",
               "楔前叶" = "Precuneus",
               "顶上皮层" = "Superior ",
               "顶下皮层" = "Inferior parietal",
               "颞横皮层" = "Transverse temporal",
               "颞上沟后侧" = "Banks superior temporal",
               "颞上回" = "Superior temporal",
               "颞中回" = "Middle temporal",
               "颞下回" = "Inferior temporal",
               "梭状回" = "Fusiform",
               "海马旁回" = "Parahippocampal",
               "内嗅皮层" = "Entorhinal",
               "颞极" = "Temporal pole",
               "外侧枕叶" = "Lateral occipital",
               "舌回" = "Lingual",
               "距状沟周围皮层" = "Pericalcarine",
               "楔叶皮层" = "Cuneus",
               "脑岛" = "Insula")
label_map_df <- data.frame(
  EnglishLabel = names(label_map),
  fulllabel = as.character(label_map),  # 确保中文标签是字符类型
  stringsAsFactors = FALSE
)

label <- c("superiorfrontal", "rostralmiddlefrontal", "caudalmiddlefrontal",
                   "parsopercularis", "parstriangularis", "parsorbitalis", "frontalpole",
                   "lateralorbitofrontal", "medialorbitofrontal",
                   "rostralanteriorcingulate", "caudalanteriorcingulate",
                   "precentral", "paracentral", "postcentral",
                   "supramarginal", "posteriorcingulate", "isthmuscingulate",
                   "precuneus", "superiorparietal", "inferiorparietal",
                   "transversetemporal", "bankssts",
                   "superiortemporal", "middletemporal", "inferiortemporal",
                   "fusiform", "parahippocampal", "entorhinal", "temporalpole",
                   "lateraloccipital", "lingual", "pericalcarine", "cuneus",
                   "insula")
label_map_df <- cbind(label_map_df, label)

###########
thr <- 0.025
# thr <- 0.01

# cluster 1
abno <- apply(group1 < thr, 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g1 <- data.frame(label)
abno_g1$perc <- as.numeric(abno)
abno_g1 <- merge(abno_g1, label_map_df, by = "label")
# 保存一个排序
abno_g1_sorted <- abno_g1 %>% arrange(desc(perc))
# 为 abno_g1_sorted 数据框添加中文标签列


name <- paste0("Cluster_1_Ab025rank_", newDate, ".xlsx")
# name <- paste0("Cluster_1_Ab01rank_", newDate, ".xlsx")
write.xlsx(abno_g1_sorted, file.path(resDir, name))

# 画脑图
abno_g1$label <- as.character(paste0("lh_", abno_g1$label))
ggseg(.data = abno_g1, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "white", high = "#d26b66", limits = c(0, 16),
                      breaks = seq(0, 16, by = 4),  # 自定义区间
                      labels = c("0", "4%", "8%", "12%", "16%")) +
  # scale_fill_gradient(low = "white", high = "#d26b66", limits = c(0, 8),
  #                     breaks = seq(0, 8, by=2),  # 自定义区间
  #                     labels = c("0", "2%", "4%", "6%", "8%")) +
  # scale_fill_viridis_c(option = "inferno", limits = c(0, 18),  # 设置上下限
  #                      breaks = seq(0, 18, by=6),  # 自定义区间
  #                      labels = c("0", "6", "12", "18")) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("Cluster_1_Ab025brain_V2_", newDate, ".png")
# name <- paste0("Cluster_1_Ab01brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# # 大于95%
# abno <- apply(group1 >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
# label <- names(abno)
# abno_g1 <- data.frame(label)
# abno_g1$perc <- as.numeric(abno)
# abno_g1$label <- as.character(paste0("lh_", abno_g1$label))
# ggseg(.data = abno_g1, mapping = aes(fill = perc), color = "black", atlas = dk,
#       position = "stacked", hemisphere = "left", size = 1.2) +
#   theme_void() +
#   theme(legend.title = element_blank(), legend.position = "bottom",
#         legend.key.width = unit(1, "cm")) +
#   scale_fill_viridis_c(option = "inferno") +
#   guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
# name <- paste0("Cluster_1_Ab975brain_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

############  cluster 2
# 大于95% !!!
abno <- apply(group2 >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g2 <- data.frame(label)
abno_g2$perc <- as.numeric(abno)
abno_g2 <- merge(abno_g2, label_map_df, by = "label")
abno_g2_sorted <- abno_g2 %>% arrange(desc(perc))
name <- paste0("Cluster_2_Ab975rank_", newDate, ".xlsx")
# name <- paste0("Cluster_2_Ab99rank_", newDate, ".xlsx")
write.xlsx(abno_g2_sorted, file.path(resDir, name))

# 画脑图
abno_g2$label <- as.character(paste0("lh_", abno_g2$label))
ggseg(.data = abno_g2, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "white", high = "#d26b66", limits = c(0, 16),
                      breaks = seq(0, 16, by = 4),  # 自定义区间
                      labels = c("0", "4%", "8%", "12%", "16%")) +
  # scale_fill_gradient(low = "white", high = "#d26b66", limits = c(0, 5),
  #                     # breaks = seq(0, 18, by=6),  # 自定义区间
  #                     breaks = c(0, 2.5, 5),  # 自定义区间
  #                     labels = c("0", "2.5%", "5%")) 
  # scale_fill_viridis_c(option = "inferno", limits = c(0, 18),  # 设置上下限
  #                      breaks = seq(0,18, by=6),  # 自定义区间
  #                      labels = c("0", "6", "12", "18")) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("Cluster_2_Ab975brain_V2_", newDate, ".png")
# name <- paste0("Cluster_2_Ab99_brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
# # 小于5%
# abno <- apply(group2 <= thr, 2, function(x) sum(x) / length(x) * 100)
# label <- names(abno)
# abno_g2 <- data.frame(label)
# abno_g2$perc <- as.numeric(abno)
# abno_g2$label <- as.character(paste0("lh_", abno_g2$label))
# ggseg(.data = abno_g2, mapping = aes(fill = perc), color = "black", atlas = dk,
#       position = "stacked", hemisphere = "left", size = 1.2) +
#   theme_void() +
#   theme(legend.title = element_blank(), legend.position = "bottom",
#         legend.key.width = unit(1, "cm")) +
#   scale_fill_viridis_c(option = "inferno") +
#   guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
# name <- paste0("Spect_Cluster_2_Ab025brain_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)



##### SM Part 1: Arrange extreme percentage ########################################################
# merge
AbP <- full_join(abno_g1_sorted, abno_g2_sorted, by = "label", suffix = c(".L", ".H"))
AbP <- AbP[, c(3,4,1,2,5)]
colnames(AbP) <- c("脑区", "FullName", "Region","L", "H")
# 保留小数点后2位
AbP$L2dig <- sprintf("%.2f", AbP$L)
AbP$H2dig <- sprintf("%.2f", AbP$H)

# 保存文件
name <- paste0("Cluster_AbPerc_", newDate, ".xlsx")
write.xlsx(AbP, file.path(resDir, name))
