# 本代码用来分析由两组ASD男性的变量之间的相关系数和显著性水平
# 雪如 2024年2月28日于北师大办公室
##################################
rm(list=ls())

# load packages
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "gamm4", "openxlsx",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
phenoDir <- file.path(abideDir, "Preprocessed")
clustDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster/Corr")
statDir <- file.path(abideDir, "Analysis/Statistic/GmmCluster")
resDate <- "240315"
newDate <- "240610"


# 认知行为
pheno <- read.csv(file.path(phenoDir, paste0("abide_A_All_", resDate, ".csv")))
colnames(pheno)[which(names(pheno) == "Participant")] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
start <- which(names(cluster) == "bankssts")
colnames(cluster)[start:ncol(cluster)] <- paste0(colnames(cluster)[start:ncol(cluster)], "_centile")

All <- merge(cluster, pheno, by = "participant", All.x = TRUE)


# 选择自变量列
names_brain <- c("middletemporal_centile", "insula_centile")
# names_brain <- names(cluster)[10:ncol(cluster)]

# # 选择因变量列
# names_cog <- c("FIQ", "VIQ", "PIQ", "ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT",
#                "ADOS_2_RRB", "SRS_AWARENESS_RAW", "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW",
#                "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW", "SRS_AWARENESS_T", "SRS_COGNITION_T",
#                "SRS_COMMUNICATION_T", "SRS_MOTIVATION_T", "SRS_MANNERISMS_T", "SRS_TOTAL_T",
#                "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
#                "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C", "VINELAND_ABC_Standard",
#                "VINELAND_COMMUNICATION_STANDARD", "VINELAND_DAILYLIVING_STANDARD",
#                "VINELAND_SOCIAL_STANDARD", "BMI")
names_cog <- c("ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT","ADOS_2_RRB")

names_col <- c("clusterID", names_brain, names_cog)
All <- All[, names_col]
All[All < 0] <- NA

L <- subset(All, clusterID == "1")
L <- L[, -1]
H <- subset(All, clusterID == "2")
H <- H[, -1]


##################################### Part 1: 计算变量之间的相关 ###################################

# L组
results_L <- data.frame(x = character(),
                        y = character(),
                        R = numeric(),
                        P = numeric(),
                        stringsAsFactors = FALSE)

for (i in 1:length(names_brain)) { # 自变量
  for (j in (length(names_brain) +1):ncol(L)) { # 因变量
    cor_test_result <- cor.test(L[, i], L[, j], method = "spearman")
    results_L <- rbind(results_L, data.frame(x = names(L)[i],
                                             y = names(L)[j],
                                             R = cor_test_result$estimate,
                                             P = cor_test_result$p.value))
  }
}
# results_L$P <- p.adjust(results_L$P, method = "BH")
# H组
results_H <- data.frame(x = character(),
                        y = character(),
                        R = numeric(),
                        P = numeric(),
                        stringsAsFactors = FALSE)

for (i in 1:length(names_brain)) { # 自变量
  for (j in (length(names_brain) +1):ncol(H)) { # 因变量
    cor_test_result <- cor.test(H[, i], H[, j], method = "spearman")
    results_H <- rbind(results_H, data.frame(x = names(H)[i],
                                             y = names(H)[j],
                                             R = cor_test_result$estimate,
                                             P = cor_test_result$p.value))
  }
}
# results_H$P <- p.adjust(results_H$P, method = "BH")

# ########### 给P值排序
# L_sorted <- arrange(results_L, P)
# L_sorted <- L_sorted[abs(L_sorted$R) >= 0.2, ] # 删除R列绝对值小于0.2的行
# L_sorted <- L_sorted[L_sorted$P < 0.05, ] # 删除数据框中P列大于或等于0.05的行
# 
# H_sorted <- arrange(results_H, P)
# H_sorted <- H_sorted[abs(H_sorted$R) >= 0.2, ]
# H_sorted <- H_sorted[H_sorted$P < 0.05, ]

# ### 保存下来结果
# name <- paste0("abide_A_asd_male_dev_GMM_Cluster_statis_Corr_L_", resDate, ".csv")
# write.csv(L_sorted, file.path(statiDir, name), row.names = F)
# name <- paste0("abide_A_asd_male_dev_GMM_Cluster_statis_Corr_H_", resDate, ".csv")
# write.csv(H_sorted, file.path(statiDir, name), row.names = F)



# ##################################### Part 2: 画图 #################################################
# 
# ##### 画L组图 
# 
# for (i in 1:nrow(L_sorted)) {
#   to_plot_names <- c(L_sorted[i, 1], L_sorted[i, 2])
#   
#   plotPoint <- L[, to_plot_names]
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   if (nrow(plotPoint) < 40) {
#     next  # 如果行数少于40，跳过此次循环的剩余部分，也就是说，不够40个的就不看了
#   }
#   
#   # 如果行数不少于40，继续执行下面的代码
#   colnames(plotPoint) <- c("x","y")
#   note_p <- paste0("p = ", round(L_sorted[i, "P"], 4))
#   note_r <- paste0("r = ", round(L_sorted[i, "R"], 4))
#   
#   ggplot(plotPoint, aes(x = x, y = y)) +
#     geom_point(color = "#add8e6", alpha = .8, size = 2, shape = 16) +  # 添加散点图层
#     geom_smooth(method = "lm", se = T, lwd = 2, color = "#add8e6", fill = "#add8e6") +
#     theme_cowplot() +
#     scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#     xlab(to_plot_names[1]) +
#     ylab(to_plot_names[2]) +
#     annotate("text", x = Inf, y = Inf, label = note_r, hjust = 1.2, vjust = 3, size = 7) +
#     annotate("text", x = Inf, y = Inf, label = note_p, hjust = 1.2, vjust = 1.2, size = 7) +
#     
#     theme(legend.position = "none", # without legend
#           axis.text.y = element_text(size = 15, face = "bold"),
#           axis.text.x = element_text(size = 15, face = "bold"))
#   
#   name <- paste0("L_", i, "_", to_plot_names[1], "_", to_plot_names[2], "_", resDate, ".png")
#   ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# }
# 
# 
# ####################################### 画H组图 ####################################################
# 
# for (i in 1:nrow(H_sorted)) {
#   to_plot_names <- c(H_sorted[i, 1], H_sorted[i, 2])
#   
#   plotPoint <- H[, to_plot_names]
#   plotPoint <- plotPoint[!is.na(plotPoint[[2]]), ]
#   
#   if (nrow(plotPoint) < 40) {
#     next  # 如果行数少于40，跳过此次循环的剩余部分，也就是说，不够40个的就不看了
#   }
#   
#   # 如果行数不少于40，继续执行下面的代码
#   colnames(plotPoint) <- c("x","y")
#   note_p <- paste0("p = ", round(H_sorted[i, "P"], 4))
#   note_r <- paste0("r = ", round(H_sorted[i, "R"], 4))
# 
#   ggplot(plotPoint, aes(x = x, y = y)) +
#     geom_point(color = "#ffb699", alpha = .8, size = 2, shape = 16) +  # 添加散点图层
#     geom_smooth(method = "lm", se = T, lwd = 2, color = "#ffb699", fill = "#ffb699") +
#     theme_cowplot() +
#     scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#     xlab(to_plot_names[1]) +
#     ylab(to_plot_names[2]) +
#     annotate("text", x = Inf, y = Inf, label = note_r, hjust = 1.2, vjust = 3, size = 7) +
#     annotate("text", x = Inf, y = Inf, label = note_p, hjust = 1.2, vjust = 1.2, size = 7) +
#     theme(legend.position = "none", # without legend
#           axis.text.y = element_text(size = 15, face = "bold"),
#           axis.text.x = element_text(size = 15, face = "bold"))
#   
#   name <- paste0("H_", i, "_", to_plot_names[1], "_", to_plot_names[2], "_", resDate, ".png")
#   ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# }
