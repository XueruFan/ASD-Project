# 把对分类特征的权重（SVM-RFECVSVM-RFECV）画出来
# ASD 男性 6-17.9岁
# 注意，这里只是用34个脑区作为分类指标
# Xue-Ru Fan 13 march 2024 @BNU
###################################################
# png
##################################################

rm(list=ls())

packages <- c("mclust", "ggplot2", "cluster", "Cairo", "openxlsx", "ggseg", "ggsegDefaultExtra")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
dataDir <- file.path(abideDir, "Preprocessed")
resDir <- file.path(abideDir, "Analysis/Cluster/Gmm515")
plotDir <- file.path(abideDir, "Plot/Cluster/Gmm515")
resDate <- "240315"
newDate <- "240610"

feature_weight <- read.xlsx(file.path(resDir, "select_region_sorted.xlsx"))

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

# 创建英文和中文标签的对应表
label_map <- data.frame(Region = en_labels, Chinese_Label = cn_labels)

# 使用 merge 函数，将 feature_weight 和 label_map 合并
feature_weight_with_cn <- merge(feature_weight, label_map, by = "Region", all.x = TRUE)

# 查看结果
print(feature_weight_with_cn)

# 如果需要按照权重排序
feature_weight_with_cn <- feature_weight_with_cn[order(-feature_weight_with_cn$Weight), ]


plot_data <- data.frame(Feature = feature_weight_with_cn$Chinese_Label, 
                        Importance = feature_weight_with_cn$Weight)

# 对数据框根据重要性评分进行排序
plot_data <- plot_data[order(plot_data$Importance, decreasing = TRUE), ]

################################# save plot
name_rank <- paste0("Cluster_SVMRFECV_Feature_Rank_", newDate, ".png")
CairoPNG(file.path(plotDir, name_rank), width = 7, height = 8, units = "in", dpi = 500)

# 使用ggplot2绘制条形图
ggplot(plot_data, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#66cdaa") +
  coord_flip() +  # 翻转坐标轴，使特征名称更容易阅读
  theme_minimal() +  # 使用简洁的主题
  labs(x = "", y = "特征权重") +
  theme(text = element_text(family = "STSong"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜，以防重叠
dev.off()



################################# 保存该结果
rank <- feature_weight_with_cn[order(feature_weight_with_cn$Weight, decreasing = TRUE), ]

name <- paste0("Cluster_SVMRFECV_Feature_Rank_", newDate, ".xlsx")
write.xlsx(rank, file.path(resDir, name))




################### 画在脑子上
to_plot<- data.frame(rank[, c(1,2)])
# to_plot <- to_plot[abs(to_plot$value) > 1, ] # 筛选权重>1的脑区
colnames(to_plot) <- c("label", "value")
to_plot$label <- paste0("lh_", to_plot$label)
to_plot$value <- as.numeric(to_plot$value)

# ggseg(.data = to_plot, atlas = dkextra, mapping = aes(fill = value), color = "black",
#       hemisphere = "left", size = 0.8) +
ggseg(.data = to_plot, mapping = aes(fill = value), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  scale_fill_gradient2(low = "#006a6c", mid = "white", high = "#bd6856", midpoint = 0) +  # 设置颜色范围，0为白色
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  # scale_fill_gradientn(colors = c("#ffffff", "#fce2c4", "#65318e"),  # 下限-中间-上限
  #                      limits = c(0, 0.6),  # 设置上下限
  #                      breaks = c(0, 0.2, 0.4, 0.6),  # 自定义区间
  #                      labels = c("0", "0.2", "0.4", "0.6")) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))


name <- paste0("Cluster_SVMRFECV_Feature_Weight_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7.8, height = 3, units = "in", dpi = 500)
