# 本代码用来分析两组ASD发育男性的变量之间的相关系数和显著性水平
# 使用斯皮尔曼相关，脑指标（centile）分别作为自变量
# 雪如 2024年2月28日于北师大办公室
##################################
# Site是控制变量，TCV,AGE作为协变量控制掉
# 保存原始的相关系数和p值csv文件，另外，筛选p小于0.05的结果，保存csv文件
# 按照显著性水平结果，绘制相关图png

##################################


rm(list=ls())
packAGEs <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packAGEs,instAll.packAGEs,character.only=TRUE)
sapply(packAGEs, require, character.only = TRUE)

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/inde/515/Corr")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("Cluster_", resDate, ".csv")
cluster <- read.csv(file.path(cabicDir, "result/inde/515", name))


cabic_all <- merge(cluster, pheno, by = "participant")


names_brain <- names(cluster)[4:(ncol(cluster)-1)]

names_cog <- names(pheno)[9:ncol(pheno)]
names_cog <- names_cog[-2]

#### select: 控制变量、自变量、因变量
names_col <- c("clusterID", "SITE", "TCV", "AGE", names_brain[-4], names_cog)
temp <- cabic_all[, names_col]
temp[temp < 0] <- NA
colnames(temp)[1] <- "clusterID"

temp <- subset(temp, AGE < 15 & AGE >= 5)

L <- subset(temp, clusterID == "1")
L <- L[, -1]
H <- subset(temp, clusterID == "2")
H <- H[, -1]

# 初始化空数据框来存储相关计算的结果
results_L <- data.frame()
results_H <- data.frame()


################### L组

# 确保 SITE 和 FIQ 是因子和数值类型
L$SITE <- as.factor(L$SITE)


# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(L) & name_cog %in% names(L)) {
      y <- L[[name_cog]] # 提取因变量
      x <- L[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_L <- L[!is.na(y) & !is.na(x) & !is.na(L$TCV),
                    c(name_brain, "SITE", "TCV", "AGE",name_cog)]
        temp_L$y <- y[!is.na(y) & !is.na(x) & !is.na(L$TCV)]

        if (nrow(temp_L) >= 30) {
          if (length(unique(temp_L$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV + AGE, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV + AGE")), 
                       data = temp_L)

            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)

            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 TCV
            y_lm <- lm(y ~ TCV + AGE, data = temp_L)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV + AGE")), data = temp_L)

            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)

            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          }

          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_L <- rbind(results_L, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value
          ))
        }
      }
    }
  }
}
results_L$P_adj <- p.adjust(results_L$p_value, method = "bonferroni")
########### 给P值排序
L_sorted <- arrange(results_L, p_value)


################### H组

# 确保 SITE 和 FIQ 是因子和数值类型
H$SITE <- as.factor(H$SITE)

# 循环 names_brain 列，计算每个自变量与 names_cog 因变量之间的 Spearman 相关性
for (name_brain in names_brain) {
  for (name_cog in names_cog) {
    if (name_brain %in% names(H) & name_cog %in% names(H)) {
      y <- H[[name_cog]] # 提取因变量
      x <- H[[name_brain]] # 提取自变量
      # 确保 y 和 x 中非 NA 值不少于 30 个
      if (sum(!is.na(y) & !is.na(x)) >= 30) {
        # 创建临时数据框，包含做相关分析需要的列，并过滤掉 FIQ 为 NA 的行
        temp_H <- H[!is.na(y) & !is.na(x) & !is.na(H$TCV),
                    c(name_brain, "SITE", "TCV", "AGE",name_cog)]
        temp_H$y <- y[!is.na(y) & !is.na(x) & !is.na(H$TCV)]

        if (nrow(temp_H) >= 30) {
          if (length(unique(temp_H$SITE)) > 1) { # 确保 SITE 还有多个水平
            # 对 y 和 x 进行回归，控制 SITE 和 FIQ
            y_lm <- lm(y ~ SITE + TCV + AGE, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ SITE + TCV + AGE")),
                       data = temp_H)

            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)

            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          } else {
            # 如果 SITE 只有一个水平，不控制 SITE，但仍然控制 FIQ
            y_lm <- lm(y ~ TCV + AGE, data = temp_H)
            x_lm <- lm(as.formula(paste(name_brain, "~ TCV + AGE")), data = temp_H)

            # 提取残差
            y_residuals <- residuals(y_lm)
            x_residuals <- residuals(x_lm)

            # 使用 Spearman 相关性计算残差之间的相关性
            cor_test <- cor.test(y_residuals, x_residuals, method = "spearman")
          }

          # 保存 Spearman 相关系数和显著性到新的数据框中
          results_H <- rbind(results_H, data.frame(
            name_cog = name_cog,
            name_brain = name_brain,
            coef = cor_test$estimate,
            p_value = cor_test$p.value
          ))
        }
      }
    }
  }
}
results_H$P_adj <- p.adjust(results_H$p_value, method = "bonferroni")
########### 给P值排序
H_sorted <- arrange(results_H, p_value)

L_sorted <- L_sorted[L_sorted$p_value < 0.05, ] # 删除数据框中P列大于或等于0.05的行
H_sorted <- H_sorted[H_sorted$p_value < 0.05, ]


name <- paste0("cabic_Part6_H_", resDate, ".csv")
write.csv(H_sorted, file.path(resuDir, name), row.names = F)
### 保存下来结果
name <- paste0("cabic_Part6_L_", resDate, ".csv")
write.csv(L_sorted, file.path(resuDir, name), row.names = F)

# L_selected <- L_sorted[L_sorted$coef < -0.3 | L_sorted$coef > 0.3, ]
# H_selected <- H_sorted[H_sorted$coef < -0.3 | H_sorted$coef > 0.3, ]
# L_selected$cluster <- "L"
# H_selected$cluster <- "H"
# 
# selected <- rbind(L_selected, H_selected)
# name <- paste0("corr_Part6_MAIN_", newDate, ".csv")
# write.csv(selected, file.path(statiDir, name), row.names = F)

# 
# #################### 筛选显著且相关大于0.2的结果
# H_sorted <- H_sorted %>%
#   # filter(abs(coef) > 0.2, p_value < 0.05)
#   filter(p_value < 0.05)
# L_sorted <- L_sorted %>%
#   # filter(abs(coef) > 0.2, p_value < 0.05)
#   filter(p_value < 0.05)
# colnames(L_sorted)[3:4] <- c("Lr", "Lp")
# colnames(H_sorted)[3:4] <- c("Hr", "Hp")
# sorted <- full_join(L_sorted, H_sorted, by = c("name_cog", "name_brain"))
# 
# sorted <- sorted %>%
#   mutate(SortValue = ifelse(is.na(Lp), abs(Hp), abs(Lp)))
# # 根据 SortValue 列对数据框进行排序
# sorted <- sorted %>%
#   arrange(SortValue)
# 
# sorted$Ln <- NA
# sorted$Hn <- NA
# 
# for (i in 1:nrow(sorted)) {
#   temp <- L[, c(sorted[i, "name_brain"], sorted[i, "name_cog"])]
#   temp <- temp[!is.na(temp[[2]]),]
#   sorted[i, "Ln"] <- nrow(temp)
#   temp <- H[, c(sorted[i, "name_brain"], sorted[i, "name_cog"])]
#   temp <- temp[!is.na(temp[[2]]),]
#   sorted[i, "Hn"] <- nrow(temp)
# }
# 
# name <- paste0("asd_male_dev_GC_corr_Part3_SiteFIQ_LH_", newDate, ".csv")
# write.csv(sorted[, -7], file.path(statiDir, name), row.names = F)
# 
# 
# #################### 画图
# 
# 
# for (i in 1:nrow(sorted)) {
#   to_plot_names <- c(sorted[i, "name_brain"], sorted[i, "name_cog"])
# 
#   plotPoint_L <- L[, to_plot_names]
#   plotPoint_L <- plotPoint_L[!is.na(plotPoint_L[[2]]), ]
#   plotPoint_H <- H[, to_plot_names]
#   plotPoint_H <- plotPoint_H[!is.na(plotPoint_H[[2]]), ]
# 
#   if (nrow(plotPoint_L) < 30 & nrow(plotPoint_H) < 30) {
#     next  # 如果数据点过少，跳过当前循环
#   }
# 
#   colnames(plotPoint_L)[1:2] <- c("x","y")
#   colnames(plotPoint_H)[1:2] <- c("x","y")
# 
# 
#   if (is.na(sorted[i, "Lp"])) {
#     ggplot() +
#       geom_point(data = plotPoint_L, aes(x = x, y = y, color = "lightgray"),
#                  alpha = .8, size = 2, shape = 16) +
#       geom_smooth(data = plotPoint_L, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
#                   color = "lightgray", fill = "lightgray") +
#       geom_point(data = plotPoint_H, aes(x = x, y = y, color = "#faa264"),
#                  alpha = .8, size = 2, shape = 16) +
#       geom_smooth(data = plotPoint_H, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
#                   color = "#faa264", fill = "#faa264") +
#       scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#       xlab(to_plot_names[1]) +
#       ylab(to_plot_names[2]) +
#       theme_cowplot() +
#       theme(legend.position = "none", # without legend
#             axis.text.y = element_text(size = 15, face = "bold"),
#             axis.text.x = element_text(size = 15, face = "bold")) +
#       scale_color_identity()
#   } else {
#     ggplot() +
#       geom_point(data = plotPoint_H, aes(x = x, y = y, color = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264")),
#                  alpha = .8, size = 2, shape = 16) +
#       geom_smooth(data = plotPoint_H, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
#                   color = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264"),
#                   fill = ifelse(is.na(sorted[i, "Hp"]), "lightgray", "#faa264")) +
#       scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#       geom_point(data = plotPoint_L, aes(x = x, y = y, color = "#719988"),
#                  alpha = .8, size = 2, shape = 16) +
#       geom_smooth(data = plotPoint_L, aes(x = x, y = y), method = "lm", se = TRUE, lwd = 2,
#                   color = "#719988", fill = "#719988") +
#       xlab(to_plot_names[1]) +
#       ylab(to_plot_names[2]) +
#       theme_cowplot() +
#       theme(legend.position = "none", # without legend
#             axis.text.y = element_text(size = 15, face = "bold"),
#             axis.text.x = element_text(size = 15, face = "bold")) +
#       scale_color_identity()
#   }
#   name <- paste0("GC_corr_Part3_SiteFIQ_", to_plot_names[1], "_", to_plot_names[2], "_", newDate, ".png")
#   ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# }

