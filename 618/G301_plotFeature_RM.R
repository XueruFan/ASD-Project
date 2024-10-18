# this script is used to plot GMM cluster top areas on brain
# Xue-Ru Fan 9 August 2024 @BNU
###################################################
# png
###################################################
rm(list=ls())
packages <- c("ggseg", "ggplot2", "ggsegDefaultExtra", "showtext", "extrafont")
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE' # winds
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster")
newDate <- "240610"


####################### plot specific area #########################################################

# 定义脑区名称、颜色和对应的中文名称
region <- c("fusiform","supramarginal","precuneus","lateral orbitofrontal","transverse temporal",
            "postcentral","paracentral","isthmus cingulate","cuneus","frontal pole")

chinese = c("梭状回","缘上回","楔前叶",
            "外侧眶额","颞横皮层","中央后回","中央旁小叶",
            "扣带回峡部","楔叶皮层","额极")

color <- rev(c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786",
           "#D8576B", "#ED7953", "#FB9F3A", "#FDC328", "#F0F921"))

brain <- data.frame(region, chinese, color)

brain$region <- factor(brain$region, levels = region)

colors_named <- setNames(brain$color, brain$region)
labels_named <- setNames(brain$chinese, brain$region)

ggseg(.data = brain, atlas = dkextra, mapping = aes(fill = region), color = "black",
      hemisphere = "left", size = 0.8) +
  scale_fill_manual(values = colors_named, labels = labels_named, breaks = region,
                    na.value = "lightgray") +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        # legend.position = "right",            # 图例位置，可以是 "right", "left", "top", "bottom"
        legend.key.size = unit(0.7, "cm"),      # 图例大小
        legend.text = element_text(size = 8, family = "STSong", face = "bold"))   # 图例文字大小

name <- paste0("asd_male_GMM_Cluster_RMrank_Brain_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 3.5, units = "in", dpi = 500)
