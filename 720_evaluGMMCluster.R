# 从人数上分析脑区异常对于人群的预测准确率
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# Part 1: 根据聚类贡献排名前10的脑区来预测，png
# Part 2: 根据异常流行率排名靠前的脑区来预测，png
# Part 3: 把1和2的图画在一起，png
# Part 4：梭状回、颞中回、中央后回、额极一起预测
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx", "reshape2")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster")
resDate <- "240315"
newDate <- "240610"

name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
cluster_result <- read.csv(file.path(resDir, name))



################################# Part 1 ###########################################################

cluster_less_1 <- cluster_result[with(cluster_result, fusiform < 0.5), ]
cluster_less_2 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5), ]
cluster_less_3 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5), ]
cluster_less_4 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5), ]
cluster_less_5 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5), ]
cluster_less_6 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5 & postcentral < 0.5), ]
cluster_less_7 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5 & postcentral < 0.5 & 
                                        paracentral < 0.5), ]
cluster_less_8 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5 & postcentral < 0.5 & 
                                        paracentral < 0.5 & isthmuscingulate < 0.5), ]
cluster_less_9 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5 & postcentral < 0.5 & 
                                        paracentral < 0.5 & isthmuscingulate < 0.5 & 
                                        cuneus < 0.5), ]
cluster_less_10 <- cluster_result[with(cluster_result, fusiform < 0.5 & supramarginal < 0.5 & 
                                        precuneus < 0.5 & lateralorbitofrontal < 0.5 & 
                                        transversetemporal < 0.5 & postcentral < 0.5 & 
                                        paracentral < 0.5 & isthmuscingulate < 0.5 & 
                                        cuneus < 0.5 & temporalpole < 0.5), ]

cluster_more_1 <- cluster_result[with(cluster_result, fusiform > 0.5), ]
cluster_more_2 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5), ]
cluster_more_3 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5), ]
cluster_more_4 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5), ]
cluster_more_5 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                        transversetemporal > 0.5), ]
cluster_more_6 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                        transversetemporal > 0.5 & postcentral > 0.5), ]
cluster_more_7 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                        transversetemporal > 0.5 & postcentral > 0.5 & 
                                        paracentral > 0.5), ]
cluster_more_8 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                        transversetemporal > 0.5 & postcentral > 0.5 & 
                                        paracentral > 0.5 & isthmuscingulate > 0.5), ]
cluster_more_9 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                        precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                        transversetemporal > 0.5 & postcentral > 0.5 & 
                                        paracentral > 0.5 & isthmuscingulate > 0.5 & 
                                        cuneus > 0.5), ]
cluster_more_10 <- cluster_result[with(cluster_result, fusiform > 0.5 & supramarginal > 0.5 & 
                                         precuneus > 0.5 & lateralorbitofrontal > 0.5 & 
                                         transversetemporal > 0.5 & postcentral > 0.5 & 
                                         paracentral > 0.5 & isthmuscingulate > 0.5 & 
                                         cuneus > 0.5 & frontalpole > 0.5), ]


#######################
evalu <- data.frame(A = rep(NA, 10),
                    H = rep(NA, 10),
                    L = rep(NA, 10))


for (i in 1:10) {
  # ＞0.5的人中H
  hn <- eval(parse(text = paste0("cluster_more_", i, "[with(cluster_more_", i, ", clusterID == 2), ]")))
  # ＜0.5的人中L
  ln <- eval(parse(text = paste0("cluster_less_", i, "[with(cluster_less_", i, ", clusterID == 1), ]")))
  
  eval(parse(text = paste0("evalu[", i, ", 'H'] <- nrow(hn) / nrow(cluster_more_", i, ")")))
  eval(parse(text = paste0("evalu[", i, ", 'L'] <- nrow(ln) / nrow(cluster_less_", i, ")")))
  eval(parse(text = paste0("evalu[", i, ", 'A'] <- (nrow(hn) + nrow(ln)) / (nrow(cluster_more_", 
                           i, ") + nrow(cluster_less_", i, "))")))
}

write.xlsx(evalu, file.path(resDir, paste0("asd_male_GMM_Cluster_Evalue_", newDate,
                                                  ".xlsx")), rowNames = T, colNames = T)


evalu$Top <- as.numeric(rownames(evalu))


name <- paste0("asd_male_GMM_Cluster_Evalue_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

# 使用 ggplot 绘制折线图
ggplot(data = evalu, aes(x = Top)) + 
  geom_line(aes(y = H, colour = "H组"), size = 2) +  # 为H列绘制折线
  geom_line(aes(y = L, colour = "L组"), size = 2) +  # 为L列绘制折线
  geom_line(aes(y = A, colour = "两组"), size = 2, alpha = .8) +  # 为A列绘制折线
  labs(x = "", y = "") +
  theme_minimal() +  # 使用简洁的主题
  scale_colour_manual(values = c("两组" = "black", "H组" = "#faa264", "L组" = "#719988")) +  # 定义线条颜色
  scale_x_continuous(breaks = seq(0, 10, by = 2),
                     labels = c(seq(0, 10, by = 2))) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1),
                     labels = c("50%", "60%", "70%", "80%", "90%", "100%")) +
  theme(text = element_text(family = "STSong"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "#e5e4e6", size = 0.5),
        panel.grid.minor = element_line(colour = "#f3f3f3", size = 0.01),
        legend.position = c(0.95, 0.1),  # 将图例放置在右下角
        legend.justification = c(1, 0),   # 图例的对齐方式，右下
        legend.text = element_text(family = "STSong", size = 15, face = "bold"),
        legend.key.width = unit(1, "cm"),  # 调整图例键的宽度
        legend.spacing = unit(0.5, "cm"))
dev.off()




################################# Part 2 ###########################################################

cluster_less_1 <- cluster_result[with(cluster_result, middletemporal < 0.5), ]

cluster_less_2 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5), ]

cluster_less_3 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5), ]

cluster_less_4 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5), ]

cluster_less_5 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5), ]

cluster_less_6 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5 & lingual < 0.5), ]

cluster_less_7 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5 & lingual < 0.5 & 
                                        medialorbitofrontal < 0.5), ]

cluster_less_8 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5 & lingual < 0.5 & 
                                        medialorbitofrontal < 0.5 & fusiform < 0.5), ]

cluster_less_9 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5 & lingual < 0.5 & 
                                        medialorbitofrontal < 0.5 & fusiform < 0.5 &
                                        pericalcarine < 0.5), ]

cluster_less_10 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                        inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                        inferiortemporal < 0.5 & lingual < 0.5 & 
                                        medialorbitofrontal < 0.5 & fusiform < 0.5 &
                                        pericalcarine < 0.5 & parsopercularis < 0.5), ]

cluster_less_11 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                         inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                         inferiortemporal < 0.5 & lingual < 0.5 & 
                                         medialorbitofrontal < 0.5 & fusiform < 0.5 &
                                         pericalcarine < 0.5 & parsopercularis < 0.5 &
                                         postcentral < 0.5), ]

cluster_less_12 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                         inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                         inferiortemporal < 0.5 & lingual < 0.5 & 
                                         medialorbitofrontal < 0.5 & fusiform < 0.5 &
                                         pericalcarine < 0.5 & parsopercularis < 0.5 &
                                         postcentral < 0.5 & insula < 0.5), ]

cluster_less_13 <- cluster_result[with(cluster_result, middletemporal < 0.5 & temporalpole < 0.5 & 
                                         inferiorparietal < 0.5 & frontalpole < 0.5 & 
                                         inferiortemporal < 0.5 & lingual < 0.5 & 
                                         medialorbitofrontal < 0.5 & fusiform < 0.5 &
                                         pericalcarine < 0.5 & parsopercularis < 0.5 &
                                         postcentral < 0.5 & insula < 0.5 & precentral < 0.5), ]



cluster_more_1 <- cluster_result[with(cluster_result, middletemporal > 0.5), ]

cluster_more_2 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5), ]

cluster_more_3 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5), ]

cluster_more_4 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5), ]

cluster_more_5 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                        inferiortemporal > 0.5), ]

cluster_more_6 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                        inferiortemporal > 0.5 & lingual > 0.5), ]

cluster_more_7 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                        inferiortemporal > 0.5 & lingual > 0.5 & 
                                        medialorbitofrontal > 0.5), ]

cluster_more_8 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                        inferiortemporal > 0.5 & lingual > 0.5 & 
                                        medialorbitofrontal > 0.5 & fusiform > 0.5), ]

cluster_more_9 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                        inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                        inferiortemporal > 0.5 & lingual > 0.5 & 
                                        medialorbitofrontal > 0.5 & fusiform > 0.5 &
                                        pericalcarine > 0.5), ]

cluster_more_10 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                         inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                         inferiortemporal > 0.5 & lingual > 0.5 & 
                                         medialorbitofrontal > 0.5 & fusiform > 0.5 &
                                         pericalcarine > 0.5 & parsopercularis > 0.5), ]

cluster_more_11 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                         inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                         inferiortemporal > 0.5 & lingual > 0.5 & 
                                         medialorbitofrontal > 0.5 & fusiform > 0.5 &
                                         pericalcarine > 0.5 & parsopercularis > 0.5 &
                                         postcentral > 0.5), ]

cluster_more_12 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                         inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                         inferiortemporal > 0.5 & lingual > 0.5 & 
                                         medialorbitofrontal > 0.5 & fusiform > 0.5 &
                                         pericalcarine > 0.5 & parsopercularis > 0.5 &
                                         postcentral > 0.5 & insula > 0.5), ]

cluster_more_13 <- cluster_result[with(cluster_result, middletemporal > 0.5 & temporalpole > 0.5 & 
                                         inferiorparietal > 0.5 & frontalpole > 0.5 & 
                                         inferiortemporal > 0.5 & lingual > 0.5 & 
                                         medialorbitofrontal > 0.5 & fusiform > 0.5 &
                                         pericalcarine > 0.5 & parsopercularis > 0.5 &
                                         postcentral > 0.5 & insula > 0.5 & precentral > 0.5), ]

#######################
evalu_abno <- data.frame(A = rep(NA, 13),
                    H = rep(NA, 13),
                    L = rep(NA, 13))


for (i in 1:13) {
  # ＞0.5的人中H
  hn <- eval(parse(text = paste0("cluster_more_", i, "[with(cluster_more_", i, ", clusterID == 2), ]")))
  # ＜0.5的人中L
  ln <- eval(parse(text = paste0("cluster_less_", i, "[with(cluster_less_", i, ", clusterID == 1), ]")))
  
  eval(parse(text = paste0("evalu_abno[", i, ", 'H'] <- nrow(hn) / nrow(cluster_more_", i, ")")))
  eval(parse(text = paste0("evalu_abno[", i, ", 'L'] <- nrow(ln) / nrow(cluster_less_", i, ")")))
  eval(parse(text = paste0("evalu_abno[", i, ", 'A'] <- (nrow(hn) + nrow(ln)) / (nrow(cluster_more_", 
                           i, ") + nrow(cluster_less_", i, "))")))
}

write.xlsx(evalu_abno, file.path(resDir, paste0("asd_male_GMM_Cluster_Evalue_Abno_", newDate,
                                           ".xlsx")), rowNames = T, colNames = T)


evalu_abno$Top <- as.numeric(rownames(evalu_abno))


name <- paste0("asd_male_GMM_Cluster_Evalue_Abno_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

# 使用 ggplot 绘制折线图
ggplot(data = evalu_abno, aes(x = Top)) + 
  geom_line(aes(y = H, colour = "H组"), size = 2) +  # 为H列绘制折线
  geom_line(aes(y = L, colour = "L组"), size = 2) +  # 为L列绘制折线
  geom_line(aes(y = A, colour = "两组"), size = 2, alpha = .8) +  # 为A列绘制折线
  labs(x = "", y = "") +
  theme_minimal() +  # 使用简洁的主题
  scale_colour_manual(values = c("两组" = "black", "H组" = "#faa264", "L组" = "#719988")) +  # 定义线条颜色
  scale_x_continuous(breaks = seq(0, 13, by = 2),
                     labels = c(seq(0, 13, by = 2))) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1),
                     labels = c("50%", "60%", "70%", "80%", "90%", "100%")) +
  theme(text = element_text(family = "STSong"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "#e5e4e6", size = 0.5),
        panel.grid.minor = element_line(colour = "#f3f3f3", size = 0.01),
        legend.position = c(0.95, 0.1),  # 将图例放置在右下角
        legend.justification = c(1, 0),   # 图例的对齐方式，右下
        legend.text = element_text(family = "STSong", size = 15, face = "bold"),
        legend.key.width = unit(1, "cm"),  # 调整图例键的宽度
        legend.spacing = unit(0.5, "cm"))
dev.off()


######################################### Part 3 ###################################################

colnames(evalu_abno)[1:3] <- c("Aa", "Ha", "La")
Ev <- merge(evalu, evalu_abno, by = "Top", all = TRUE)

name <- paste0("asd_male_GMM_Cluster_Evalue_All_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

# 使用 ggplot 绘制折线图
ggplot(data = Ev, aes(x = Top)) + 
  geom_line(aes(y = H, colour = "H组", linetype = "贡献度预测"), size = 2) +  # H组实线
  geom_line(aes(y = L, colour = "L组", linetype = "贡献度预测"), size = 2) +  # L组实线
  geom_line(aes(y = A, colour = "两组", linetype = "贡献度预测"), size = 2) +  # 两组实线
  geom_line(aes(y = Ha, colour = "H组", linetype = "异常率预测"), size = 2) +  # H组虚线
  geom_line(aes(y = La, colour = "L组", linetype = "异常率预测"), size = 2) +  # L组虚线
  geom_line(aes(y = Aa, colour = "两组", linetype = "异常率预测"), size = 2) +  # 两组虚线
  labs(x = "", y = "") +
  theme_minimal() +  # 使用简洁的主题
  scale_colour_manual(values = c("两组" = "gray", "H组" = "#faa264", "L组" = "#719988")) +  # 定义线条颜色
  scale_linetype_manual(values = c("贡献度预测" = 1, "异常率预测" = 4)) +  # 定义线型
  scale_x_continuous(breaks = seq(0, 13, by = 2),
                     labels = c(seq(0, 13, by = 2))) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1),
                     labels = c("50%", "60%", "70%", "80%", "90%", "100%")) +
  theme(text = element_text(family = "STSong"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "#e5e4e6", size = 0.5),
        panel.grid.minor = element_line(colour = "#f3f3f3", size = 0.01),
        legend.position = c(0.95, 0.1),  # 将图例放置在右下角
        legend.justification = c(1, 0),   # 图例的对齐方式，右下
        legend.text = element_text(family = "STSong", size = 15, face = "bold"),
        legend.key.width = unit(1, "cm"),  # 调整图例键的宽度
        legend.spacing = unit(0.5, "cm"))
dev.off()


############################ Part 4 ################################################################

# 定义你想使用的脑区
regions <- c("postcentral", "fusiform", "frontalpole", "middletemporal")

# 创建一个列表，用于存储不同组合的预测结果
evalu <- data.frame(A = rep(NA, 15),  # 共15种组合 (C(4,1) + C(4,2) + C(4,3) + C(4,4) = 15)
                    H = rep(NA, 15),
                    L = rep(NA, 15))

# 生成1到4个脑区的所有组合
combinations <- list()
counter <- 1
for (i in 1:4) {
  combn_result <- combn(regions, i, simplify = FALSE)
  for (comb in combn_result) {
    combinations[[counter]] <- comb
    counter <- counter + 1
  }
}

# 遍历每个组合，并计算预测率
for (i in 1:length(combinations)) {
  region_combination <- combinations[[i]]
  
  # 创建筛选条件：小于0.5
  condition_less <- paste(region_combination, "< 0.5", collapse = " & ")
  # 创建筛选条件：大于0.5
  condition_more <- paste(region_combination, "> 0.5", collapse = " & ")
  
  # 小于0.5的人
  eval(parse(text = paste0("cluster_less_", i, " <- cluster_result[with(cluster_result, ", condition_less, "), ]")))
  # 大于0.5的人
  eval(parse(text = paste0("cluster_more_", i, " <- cluster_result[with(cluster_result, ", condition_more, "), ]")))
  
  # 计算H组和L组的预测率
  hn <- eval(parse(text = paste0("cluster_more_", i, "[with(cluster_more_", i, ", clusterID == 2), ]")))
  ln <- eval(parse(text = paste0("cluster_less_", i, "[with(cluster_less_", i, ", clusterID == 1), ]")))
  
  eval(parse(text = paste0("evalu[", i, ", 'H'] <- nrow(hn) / nrow(cluster_more_", i, ")")))
  eval(parse(text = paste0("evalu[", i, ", 'L'] <- nrow(ln) / nrow(cluster_less_", i, ")")))
  eval(parse(text = paste0("evalu[", i, ", 'A'] <- (nrow(hn) + nrow(ln)) / (nrow(cluster_more_", 
                           i, ") + nrow(cluster_less_", i, "))")))
}


# # 添加组合编号（Top）
evalu$Combination <- 1:nrow(evalu)

# 为每个组合生成标签（脑区名称组合）
evaluDk <- evalu
evaluDk$Dk <- sapply(combinations, function(x) paste(x, collapse = " + "))
# 保存结果到Excel文件
write.xlsx(evaluDk, file.path(resDir, paste0("asd_male_GMM_Cluster_Evalue_4dk_", newDate,
                                           ".xlsx")), rowNames = T, colNames = T)

# 添加每个组合的脑区个数
evalu$RegionCount <- sapply(combinations, length)
evalu[,1:3] <- evalu[,1:3] - 0.5

evalu <- evalu[, -1]

# 将数据框转换为长格式
evalu_long <- melt(evalu, id.vars = c("Combination", "RegionCount"), 
                   variable.name = "Group", value.name = "PredictionRate")

# 保存并绘制图形
name <- paste0("asd_male_GMM_Cluster_Evalue_4dk_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)

# 使用 ggplot 绘制直方图，x 轴显示脑区组合名称，并按不同脑区个数排序
ggplot(data = evalu_long, aes(x = reorder(Combination, RegionCount), y = PredictionRate, fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 1, alpha = 0.85) +  # 通过width参数设置每个条形图的宽度
  labs(x = "", y = "") + # 修改x和y轴标签
  theme_minimal() +  # 使用简洁的主题
  scale_fill_manual(values = c("H" = "#faa264", "L" = "#719988")) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
  labels = c("50%", "60%", "70%", "80%", "90%", "100%")) +
  theme(text = element_text(family = "STSong"),
        axis.text = element_text(size = 12, face = "bold", hjust = 1),  # 将脑区组合名称旋转，以避免重叠
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "lightgray", size = 0.1),  # 调整主网格线颜色和粗细
        panel.grid.minor = element_line(colour = "lightgray", size = 0.1), 
        legend.position = "",
        # legend.justification = c(1, 0),   # 图例的对齐方式，右下
        legend.text = element_text(family = "STSong", size = 15, face = "bold"),
        legend.key.width = unit(1, "cm"),  # 调整图例键的宽度
        legend.spacing = unit(0.5, "cm")) +
  facet_wrap(~ RegionCount, scales = "free_x", nrow = 1)  # 按 RegionCount 分组

dev.off()

