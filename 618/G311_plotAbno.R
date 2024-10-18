# 绘图使用高斯混合聚类算法聚类后的脑区异常流行率
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# png
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "openxlsx", "dplyr", "Cairo", "ggradar")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Gmm618")
plotDir <- file.path(abideDir, "Plot/Cluster/Gmm618")
resDate <- "240315"
newDate <- "240610"

name <- paste0("Cluster_AbPerc_", newDate, ".xlsx")
abno <- read.xlsx(file.path(resDir, name))

rownames(abno) <- abno[,1]
abno <- abno[, 3:4]

abno <- abno/100
abno <- data.frame(t(abno))
abno[3,] <- abs(abno[1,]-abno[2,])

# 提取第三行的数值为向量
third_row_values <- as.numeric(abno[3, ])
# 使用第三行的数值对列排序
abno <- abno[, order(third_row_values, na.last = TRUE)]
abno <- abno[1:2,]

abno$人群 <- rownames(abno)
abno <- abno[, c(35, 34:1)]

name <- paste0("Cluster_AbnoRadar_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 10, height = 10, units = "in", dpi = 500)
ggradar(plot.data = abno,
        font.radar = "STSong",
        values.radar = c(NA, "5%", "10%"),
        grid.min = 0, grid.mid = 0.05, grid.max = 0.1,
        background.circle.colour="lightgray",#雷达图背景色的颜色
        gridline.mid.colour="lightgray",#轴值中位圈的颜色
        gridline.max.colour="lightgray",#最大轴值圈的颜色
        gridline.min.linetype = 1,#轴值最小圈的线条类型
        gridline.mid.linetype = 1,#轴值中位圈的的线条类型
        gridline.max.linetype = 1,#最大轴值圈的线条类型
        group.line.width = 2,
        axis.label.size = 3.5,
        axis.label.offset = 1.05,       # 可选，轴标签与中心点的距离比例
        grid.label.size = 4,            # 可选，网格标签的字体大小
        group.point.size = 4,           # 可选，组点的大小
        group.colours = c("#db8449", "#47885e"),
        fill = T, fill.alpha = 0.2,
        legend.text.size = 14,
        legend.position = "bottom")
dev.off()
