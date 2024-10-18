# this script is used to visualize ABIDE centile (6-17.9 only)
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# Part 1: ASD男性34个分区体积centile（中位数）的概率密度图，png
# Part 2: ASD男性34个分区体积centile的中位数投射出的脑图，png
# Part 3: ASD男性7个全局指标常模分（中位数）的概率密度图，png
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "Cairo")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
plotDir <- file.path(abideDir, "Plot/Density")
resdate <- "240315"
newDate <- "240610"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_All_centile_513_", resdate, ".csv")))
asd_centile <- subset(abide_centile, dx == "ASD" & sex == "Male")


######################### Part 1: Plot density for 34 regions ######################################
# ASD 34个分区体积Centile分（中位数）的概率密度图

asd_parc <- asd_centile[, which(names(asd_centile) == "bankssts"):ncol(asd_centile)]


# transfer to long data
all_data <- gather(asd_parc, key = "Measure", value = "Centile", bankssts:insula, na.rm = TRUE,
                   factor_key = TRUE)

# 修改度量标准的显示名称
# cn_labels <- rev(c("额上回","额中回喙部","额中回尾部",
#                    "额下回盖部","额下回三角部","额下回眶部","额极",
#                    "外侧眶额","内侧眶额",
#                    "前扣带喙部","前扣带尾部",
#                    "中央前回","中央旁小叶","中央后回",
#                    "缘上回","后扣带","扣带回峡部",
#                    "楔前叶","顶上皮层","顶下皮层",
#                    "颞横皮层","颞上沟后侧",
#                    "颞上回","颞中回","颞下回",
#                    "梭状回","海马旁回","内嗅皮层","颞极",
#                    "外侧枕叶","舌回","距状沟周围皮层","楔叶皮层",
#                    "脑岛"))
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

name_regional <- paste0("centile_regional_density_513_", newDate, ".png")
# CairoPNG(file.path(plotDir, name_regional), width = 4, height = 8, units = "in", dpi = 500)

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

ggplot(all_data, aes(x = Centile, y = Measure, fill = Measure)) +
  # geom_density_ridges(scale = 1.7, quantile_lines = TRUE, size = 0.5, quantiles = 2, alpha = .95) +
  geom_density_ridges(scale = 2, quantile_lines = TRUE, size = 0.75, quantiles = 2, alpha = .94) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = rep("white", length(levels(all_data$Measure)))) +  # 使用提供的颜色列表
  coord_fixed(ratio = 0.15) + 
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6))
ggsave(file.path(plotDir, name_regional), width = 5, height = 8, units = "in", dpi = 500)
# dev.off()

objects_to_keep <- c("plotDir", "abideDir", "resDir", "newDate", "cluster_result", "asd_centile",
                     "asd_parc")
rm(list = (setdiff(ls(), objects_to_keep)))



######################### Part 2: Project 34 region centile on the brain ###########################
# ASD患者34个分区体积常模分的中位数的脑图

# calculate median
asd_parc_centile <- data.frame(paste0("lh_", names(asd_parc)))
colnames(asd_parc_centile) = "label"
asd_parc_centile$median <- apply(asd_parc, 2, median, na.rm = TRUE)
name <- paste0("Cluster_Centile_", newDate, ".xlsx")
write.xlsx(asd_parc_centile, file.path(resDir, name))

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient2(low = "#126cb5", mid = "white", high = "#a52a2a", midpoint = 0.5) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- file.path(plotDir, paste0("centile_regional_median_513_", newDate, ".png"))
ggsave(name, width = 7.8, height = 3, units = "in", dpi = 500)


######################### Part3: Plot density for 7 measures #######################################
# ASD组7个全局指标常模分（中位数+四分位距）的概率密度图

asd_global <- asd_centile[which(names(asd_centile) == "GMV"):which(names(asd_centile) == "totalSA2")]

# transfer to long data
all_data <- gather(asd_global, key = "Measure", value = "Centile", GMV:totalSA2, na.rm = TRUE,
                   factor_key = TRUE)

# 修改度量标准的显示名称
all_data$Measure <- factor(all_data$Measure, levels = c("TCV", "WMV", "GMV", "sGMV", "Ventricles",
                                                        "meanCT2", "totalSA2"),
                           # labels = c("脑总容积", "白质总体积",
                           #            "皮层灰质总体积", "皮层下灰质总体积",
                           #            "脑脊液总体积", "平均皮层厚度",
                           #            "皮层总表面积"))
                           labels = c("TCV", "WMV",
                                      "GMV", "sGMV",
                                      "CSF", "mCT",
                                      "tSA"))

name_global <- paste0("centile_global_density_513_", newDate, ".png")
# CairoPNG(file.path(plotDir, name_global), width = 6, height = 5, units = "in", dpi = 500)

ggplot(all_data, aes(x = Centile, y = Measure, fill = Measure)) +
  # geom_density_ridges(scale = 1.6, quantile_lines = TRUE, size = 0.7, quantiles = 2, alpha = .88) +
  geom_density_ridges(scale = 1.6, quantile_lines = TRUE, size = 1, quantiles = 2, alpha = .94) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = rep("white", length(levels(all_data$Measure)))) +
  coord_fixed(ratio = 0.2) + 
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name_global), width = 6, height = 5, units = "in", dpi = 500)
