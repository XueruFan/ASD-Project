# this script is used to visualize subgroup after GMM clustering (dev ASD Male only)
# 绘图使用高斯混合聚类算法（这里只是用34个脑区作为分类指标）对人群进行聚类后的两类人群
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# 6-17.9岁
# Part 1: 画出两个分型34个分区体积centile中位数投射出的脑图，png
# Part 2: 画出两个分型7个全局指标常模分（中位数+四分位距）的概率密度图，png
# Part 3: 画出两个分型34个分区体积常模分（中位数+四分位距）的概率密度分布图，png
# Part 4: 画出两个分型34个脑区的异常（5%和95%）流行率，png和xlsx
# Part 5: 把两组人按照L组的异常率从大到小排列脑区，xlsx
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/inde/515")
plotDir <- file.path(resuDir, "Plot/Cluster")
newDate <- "240928"

name <- paste0("Cluster_", newDate, ".csv")
cluster_result <- read.csv(file.path(resuDir, name))


######################### Part 1: Project 34 region centile on the brain ###########################
group1 <- (subset(cluster_result, clusterID == "2"))[, which(names(cluster_result) == "bankssts"):
                                                               which(names(cluster_result) == "insula")]
group2 <- (subset(cluster_result, clusterID == "1"))[, which(names(cluster_result) == "bankssts"):
                                                               which(names(cluster_result) == "insula")]
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
  # scale_fill_gradient(low = "#47885e", high = "white", limits = c(min(median), max(median))) +
  scale_fill_gradient2(low = "#47885e", mid = "white", high = "#db8449", midpoint = 0.5,
                       limits = c(min(asd_parc_centile$median), max(asd_parc_centile$median))) +
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
median <- as.numeric(asd_parc_centile$median)
label <- as.character(asd_parc_centile$label)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient(low = "white", high = "#db8449", limits = c(min(median), max(median))) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- paste0("Cluster_2_Centile_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "cluster_result", "resuDir")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 2: Plot density for 7 global ########################################
# 初始化空数据框
all_data <- data.frame()
id_group <- c("1", "2") # cluster 1 是H， 2是L

for (id in id_group){
  # id <- "1"
  ID = 3-as.numeric(id)
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", ID, "))[, -1:-2]")))
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
                           labels = c("脑总容积", "白质总体积",
                                      "皮层灰质总体积", "皮层下灰质总体积",
                                      "脑脊液总体积", "平均皮层厚度",
                                      "皮层总表面积"))

name_global <- paste0("Cluster_Centile_Global_", newDate, ".png")
CairoPNG(file.path(plotDir, name_global), width = 6, height = 5, units = "in", dpi = 500)

ggplot(all_data, aes(x = Centile, y = Measure, fill = Group, group = interaction(Group, Measure))) +
  geom_density_ridges(scale = 1.6, quantile_lines = TRUE, size = 0.7, quantiles = 2, alpha = .75) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("1" = "#a6c8b2", "2" = "#f0cfa0")) +  # 此处为两组人设置不同颜色
  coord_fixed(ratio = 0.2) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"))
dev.off()


objects_to_keep <- c("plotDir", "resuDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 3: Plot density for 34 regional ########################################
# 初始化空数据框
all_data <- data.frame()
id_group <- c("1", "2") # cluster 1 是L， 2是H

for (id in id_group){
  # id <- "1"
  ID = 3-as.numeric(id)
  eval(parse(text = paste0("group", id, " <- (subset(cluster_result, clusterID == ", ID, "))[, -1:-2]")))
  eval(parse(text = paste0("regional <- group", id)))
  regional <- regional[which(names(regional) == "bankssts"):which(names(regional) == "insula")]
  asd_long <- gather(regional, key = "Measure", value = "Centile", bankssts:insula, na.rm = TRUE,
                     factor_key = TRUE)
  asd_long$Group <- factor(rep(id, nrow(asd_long)))  # 添加群体标识
  all_data <- rbind(all_data, asd_long)  # 合并数据
}

# 修改度量标准的显示名称
cn_labels <- rev(c("额上回","额中回喙部","额中回尾部",
                        "额下回盖部","额下回三角部","额下回眶部","额极",
                        "外侧眶额","内侧眶额",
                        "前扣带喙部","前扣带尾部",
                        "中央前回","中央旁小叶","中央后回",
                        "缘上回","后扣带","扣带回峡部",
                        "楔前叶","顶上皮层","顶下皮层",
                        "颞横皮层","颞上沟后侧",
                        "颞上回","颞中回","颞下回",
                        "梭状回","海马旁回","内嗅皮层","颞极",
                        "外侧枕叶","舌回","距状沟周围皮层","楔叶皮层",
                        "脑岛"))
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

all_data$Measure <- factor(all_data$Measure, levels = en_labels, labels = cn_labels)

name_regional <- paste0("Cluster_Centile_Regional_", newDate, ".png")
CairoPNG(file.path(plotDir, name_regional), width = 4, height = 8, units = "in", dpi = 500)

# color <- c("#95483f","#a22041","#b7282e",
#            "#d3381c","#ea5506","#ec6800","#f08300",
#            "#e17b34","#bf783a",
#            "#c39143","#b48a76",
#            "#c37854","#bfa46f","#a8c97f",
#            "#38b48b","#47885e","#316745",
#            "#00552e","#2f5d50","#006e54",
#            "#8aaee6","#0091ff",
#            "#2ca9e1","#0075c2","#165e83",
#            "#192f60","#19448e","#4169e1","#778899",
#            "#4a488e","#674598","#522f60","#460e44",
#            "white")

ggplot(all_data, aes(x = Centile, y = Measure, fill = Group, group = interaction(Group, Measure))) +
  geom_density_ridges(scale = 1.7, quantile_lines = TRUE, size = 0.5, quantiles = 2, alpha = .75) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("1" = "#a6c8b2", "2" = "#f0cfa0")) +
  coord_fixed(ratio = 0.15) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))

dev.off()

objects_to_keep <- c("plotDir", "resuDir", "newDate", "cluster_result")
rm(list = (setdiff(ls(), objects_to_keep)))


######################### Part 4: Project extram percentage on the brain ###########################
# 34个脑区的异常流行率, 对于每个脑区来说有多少人的体积常模分属于极端值（位于thr以外，这里thr取5%）
group1 <- (subset(cluster_result, clusterID == "2"))[, which(names(cluster_result) == "bankssts"):
                                                               which(names(cluster_result) == "insula")]
group2 <- (subset(cluster_result, clusterID == "1"))[, which(names(cluster_result) == "bankssts"):
                                                               which(names(cluster_result) == "insula")]

label_map <- c("额上回" = "superiorfrontal",
               "额中回喙部" = "rostralmiddlefrontal",
               "额中回尾部" = "caudalmiddlefrontal",
               "额下回盖部" = "parsopercularis",
               "额下回三角部" = "parstriangularis",
               "额下回眶部" = "parsorbitalis",
               "额极" = "frontalpole",
               "外侧眶额" = "lateralorbitofrontal",
               "内侧眶额" = "medialorbitofrontal",
               "前扣带喙部" = "rostralanteriorcingulate",
               "前扣带尾部" = "caudalanteriorcingulate",
               "中央前回" = "precentral",
               "中央旁小叶" = "paracentral",
               "中央后回" = "postcentral",
               "缘上回" = "supramarginal",
               "后扣带" = "posteriorcingulate",
               "扣带回峡部" = "isthmuscingulate",
               "楔前叶" = "precuneus",
               "顶上皮层" = "superiorparietal",
               "顶下皮层" = "inferiorparietal",
               "颞横皮层" = "transversetemporal",
               "颞上沟后侧" = "bankssts",
               "颞上回" = "superiortemporal",
               "颞中回" = "middletemporal",
               "颞下回" = "inferiortemporal",
               "梭状回" = "fusiform",
               "海马旁回" = "parahippocampal",
               "内嗅皮层" = "entorhinal",
               "颞极" = "temporalpole",
               "外侧枕叶" = "lateraloccipital",
               "舌回" = "lingual",
               "距状沟周围皮层" = "pericalcarine",
               "楔叶皮层" = "cuneus",
               "脑岛" = "insula")
label_map_df <- data.frame(
  EnglishLabel = names(label_map),
  label = as.character(label_map),  # 确保中文标签是字符类型
  stringsAsFactors = FALSE
)

###########
thr <- 0.025

# cluster 1
# 小于5% !!!
abno <- apply(group1 <= thr, 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g1 <- data.frame(label)
abno_g1$perc <- as.numeric(abno)
abno_g1 <- merge(abno_g1, label_map_df, by = "label")
# 保存一个排序
abno_g1_sorted <- abno_g1 %>% arrange(desc(perc))
# 为 abno_g1_sorted 数据框添加中文标签列


name <- paste0("Cluster_1_Ab025rank_", newDate, ".xlsx")
write.xlsx(abno_g1_sorted, file.path(resuDir, name))

# 画脑图
abno_g1$label <- as.character(paste0("lh_", abno_g1$label))
ggseg(.data = abno_g1, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("Cluster_1_Ab025brain_", newDate, ".png")
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
name <- paste0("Cluster_1_Ab975brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)

############  cluster 2
# 大于95% !!!
abno <- apply(group2 >= (1-thr), 2, function(x) sum(x) / length(x) * 100)
label <- names(abno)
abno_g2 <- data.frame(label)
abno_g2$perc <- as.numeric(abno)
abno_g2 <- merge(abno_g2, label_map_df, by = "label")
abno_g2_sorted <- abno_g2 %>% arrange(desc(perc))
name <- paste0("Cluster_2_Ab975rank_", newDate, ".xlsx")
write.xlsx(abno_g2_sorted, file.path(resuDir, name))

# 画脑图
abno_g2$label <- as.character(paste0("lh_", abno_g2$label))
ggseg(.data = abno_g2, mapping = aes(fill = perc), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_viridis_c(option = "inferno") +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))
name <- paste0("Cluster_2_Ab975brain_", newDate, ".png")
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
name <- paste0("GMM_Cluster_2_Ab025brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)



######################### Part 5: Agrrange extram percentage #######################################
# merge
AbP <- full_join(abno_g1_sorted, abno_g2_sorted, by = "label", suffix = c(".L", ".H"))
AbP <- AbP[, c(3,1,2,4)]
colnames(AbP) <- c("脑区","Region","L", "H")
# 保留小数点后四位
AbP$L <- round(AbP$L, digits = 4)
AbP$H <- round(AbP$H, digits = 4)

# 保存文件
name <- paste0("Cluster_AbPerc_", newDate, ".xlsx")
write.xlsx(AbP, file.path(resuDir, name))
