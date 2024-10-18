# 本代码用来分析高斯混合聚类的两组ASD男性的人口学和认知行为之间的差异（Part A）
# 雪如 2024年2月27日于北师大办公室
################################
# Part 01: 站点
# Part 02: 机型和厂家
# Part 03: 类型
# Part 1: 年龄
# Part 2: IQ
# Part 3: ADOS_G
# Part 4: ADOS_2
# Part 5: SRS
# Part 6: ADI_R
# Part 7: abc
# Part 8: BMI
################
# 以上每部分都会保存统计的csv文件和绘图的png文件
################
# Part Z：保存P值文件csv，之后需要手动excel，筛选出P值显著（＜0.05）的位置，保存一个xlsx文件
################################


rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2", "Cairo")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/Diff")
plotDir <- file.path(cabicDir, "result/Plot/Diff")
resDate <- "240928"

# 认知行为
pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "participant"
# 聚类信息、脑形态测量百分位数
name <- paste0("cabic_cluster_predictions_", resDate, ".csv")
cluster <- read.csv(file.path(cabicDir, "result", name))

All <- merge(cluster, pheno, by = "participant")
All <- subset(All, SEX == "M" & GROUP == "ASD")
All <- subset(All, AGE >= 6 & AGE < 18)

All[which(All$predicted_cluster == "1"), 'clusterID'] = "L"
All[which(All$predicted_cluster == "2"), 'clusterID'] = "H"
All$clusterID <- factor(All$clusterID)

evalu <- c("Median", "Mean", "SD") # 计算哪些统计值

# 新建一个空数据框用来存p值
Pvalue <- data.frame(matrix(ncol = 1, nrow = 3))
rownames(Pvalue) <- c("t-test", "w-test", "f-test")
colnames(Pvalue) <- "variable"


################################# Part 01: 站点 ####################################################

var <- All[, c("clusterID", "SITE")]

var$SITE <- as.factor(var$SITE)
var$clusterID <- as.factor(var$clusterID)

freq_var <- table(var$SITE, var$clusterID)
# 过滤掉那些频数和小于10的行
filtered_freq <- freq_var
# fisher检验
site_pvalue <- fisher.test(filtered_freq, simulate.p.value = TRUE, B = 1e5)

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

L_site_counts <- as.data.frame(table(L$SITE)) # 计算SITE的每个水平的频数，并转换为数据框
H_site_counts <- as.data.frame(table(H$SITE))
site_count <- merge(L_site_counts, H_site_counts, by = "Var1")
site_counts_sorted <- arrange(site_count, desc(Freq.x)) # 将结果按Count列从大到小排序
colnames(site_counts_sorted) <- c("Site", "L", "H")

# 将p值绑定到数据框的右边
site_counts_sorted_to_save <- site_counts_sorted
site_counts_sorted_to_save[1,4] <- paste0("p = ", site_pvalue[["p.value"]])

name <- paste0("statis_site_", resDate, ".csv")
write.csv(site_counts_sorted_to_save, file.path(resuDir, name), row.names = F)

data_long <- melt(site_counts_sorted, id.vars = "Site")

sorted_Site <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Site)

data_long$Site <- factor(data_long$Site, levels = unique(sorted_Site))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# save plot
name <- paste0("differ_site_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)
ggplot(data_long, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("#f0cfa0", "#a6c8b2")) +
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1))
dev.off()

objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))



################################# Part 1: 年龄 #####################################################
var <- All[, c("clusterID", "AGE")]
colnames(var)[2] <- "variable"

# save plot
name <- paste0("differ_age_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 5, height = 5, units = "in", dpi = 500)
ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = seq(6, 12, by = 2), labels = c("6岁", "8", "10", "12")) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  # coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 10))
dev.off()

# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- c("H", "L")
sta_ana[1,1] <- round(median(H$variable, na.rm = T), 2)
sta_ana[1,2] <- round(mean(H$variable, na.rm = T), 2)
sta_ana[1,3] <- round(sd(H$variable, na.rm = T), 2)
sta_ana[1,4] <- round(length(H$clusterID))
sta_ana[2,1] <- round(median(L$variable, na.rm = T), 2)
sta_ana[2,2] <- round(mean(L$variable, na.rm = T), 2)
sta_ana[2,3] <- round(sd(L$variable, na.rm = T), 2)
sta_ana[2,4] <- round(length(L$clusterID))
colnames(sta_ana)[4] <- "Count"

name <- paste0("statis_age_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

Pvalue["t-test","age"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["w-test","age"] <- wilcox.test(L$variable, H$variable)[["p.value"]]


objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))




################################# Part 5: SRS ###################################################
# 几个分量表的原始分
var <- All[, c("clusterID", "SRS_SA", "SRS_SCOG", "SRS_SCOM", "SRS_SM", "SRS_AM")]
var <- var[!apply(var[, 2:6], 1, function(row) all(is.na(row))), ]

colnames(var)[2:6] <- c("SRS_Aware", "SRS_Cogni", "SRS_Commu", "SRS_Motiv", "SRS_Manne")

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
                   na.rm = TRUE, factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}
name <- paste0("differ_SRS_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L SRS_Aware" = "#a6c8b2", "H SRS_Aware" = "#f0cfa0",
                               "L SRS_Cogni" = "#a6c8b2", "H SRS_Cogni" = "#f0cfa0",
                               "L SRS_Commu" = "#a6c8b2", "H SRS_Commu" = "#f0cfa0",
                               "L SRS_Motiv" = "#a6c8b2", "H SRS_Motiv" = "#f0cfa0",
                               "L SRS_Manne" = "#a6c8b2", "H SRS_Manne" = "#f0cfa0")) +
  scale_y_discrete(labels = c("L SRS_Aware" = "社交意识", "H SRS_Aware" = "",
                              "L SRS_Cogni" = "社交认知", "H SRS_Cogni" = "",
                              "L SRS_Commu" = "社交表达", "H SRS_Commu" = "",
                              "L SRS_Motiv" = "社交动机", "H SRS_Motiv" = "",
                              "L SRS_Manne" = "社交举止", "H SRS_Manne" = "")) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15)) +
  # xlim(50, 150) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
dev.off()



# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_SRS_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}


# 总分
var <- All[, c("clusterID", "SRS_TOTAL")]
colnames(var)[2] <- c("SRS_Total")
var[, 2][var[, 2] < 0] <- NA
var <- na.omit(var)


group <- c("L", "H")
varia <- names(var)[2]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

ggplot(var, aes(x = SRS_Total, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  scale_x_continuous(breaks = seq(30, 150, by = 30)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 15))

name <- paste0("differ_SRS_total_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_SRS_total_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))



################################# Part 7: abc ###################################################
abc_columns <- grep("ABC", names(All), value = TRUE)

# 上面是为了输出这些列名，下面手动贴上需要的
abc_standard_columns <- c("ABC_SB", "ABC_SR", "ABC_BO", "ABC_LC", "ABC_SA")
var <- All[, c("clusterID", abc_standard_columns)]
var <- var[!apply(var[, 2:6], 1, function(row) all(is.na(row))), ]
var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
                   na.rm = TRUE, factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

name <- paste0("differ_ABC_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ABC_SB" = "#a6c8b2", "H ABC_SB" = "#f0cfa0",
                               "L ABC_SR" = "#a6c8b2", "H ABC_SR" = "#f0cfa0",
                               "L ABC_BO" = "#a6c8b2", "H ABC_BO" = "#f0cfa0",
                               "L ABC_LC" = "#a6c8b2", "H ABC_LC" = "#f0cfa0",
                               "L ABC_SA" = "#a6c8b2", "H ABC_SA" = "#f0cfa0")) +
  # scale_y_discrete(labels = c("L VIN_Commu_t" = "沟通", "H VIN_Commu_t" = "",
  #                             "L VIN_Daily_t" = "日常生活", "H VIN_Daily_t" = "",
  #                             "L VIN_Socia_t" = "社交", "H VIN_Socia_t" = "",
  #                             "L VIN_ABC_t" = "适应行为综合得分", "H VIN_ABC_t" = "")) +
  # scale_x_continuous(breaks = seq(40, 110, by = 20)) +
  xlim(NA,50) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
dev.off()


# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_ABC_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}



# 总分
var <- All[, c("clusterID", "ABC_TOTAL")]
var[, 2][var[, 2] < 0] <- NA
var <- na.omit(var)


group <- c("L", "H")
varia <- names(var)[2]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

ggplot(var, aes(x = ABC_TOTAL, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  # scale_x_continuous(breaks = seq(30, 150, by = 30)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 15))

name <- paste0("differ_ABC_total_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_ABC_total_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 8: BMI ######################################################
var <- All[, c("clusterID", "BMI")]
colnames(var)[2] <- "variable"
var$variable[var$variable < 0] <- NA
var <- na.omit(var)

# save plot
name <- paste0("differ_bmi_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)
ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30), limits = c(NA, 30)) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  # coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 10))
dev.off()



# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- c("H", "L")
sta_ana[1,1] <- round(median(H$variable, na.rm = T), 2)
sta_ana[1,2] <- round(mean(H$variable, na.rm = T), 2)
sta_ana[1,3] <- round(sd(H$variable, na.rm = T), 2)
sta_ana[1,4] <- round(length(H$clusterID))
sta_ana[2,1] <- round(median(L$variable, na.rm = T), 2)
sta_ana[2,2] <- round(mean(L$variable, na.rm = T), 2)
sta_ana[2,3] <- round(sd(L$variable, na.rm = T), 2)
sta_ana[2,4] <- round(length(L$clusterID))
colnames(sta_ana)[4] <- "Count"

name <- paste0("statis_bmi_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

Pvalue["t-test","bmi"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["w-test","bmi"] <- wilcox.test(L$variable, H$variable)[["p.value"]]

objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 9: RBSR #####################################################
pick_columns <- grep("RBS", names(All), value = TRUE)

var <- All[, c("clusterID", pick_columns[-length(pick_columns)])]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
var <- var[!apply(var[, 2:ncol(var)], 1, function(row) all(is.na(row))), ]

var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
                   na.rm = TRUE, factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]


group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

name <- paste0("differ_RBS_", resDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L RBS_SI" = "#a6c8b2", "H RBS_SI" = "#f0cfa0",
                               "L RBS_STB" = "#a6c8b2", "H RBS_STB" = "#f0cfa0",
                               "L RBS_CB" = "#a6c8b2", "H RBS_CB" = "#f0cfa0",
                               "L RBS_RIB" = "#a6c8b2", "H RBS_RIB" = "#f0cfa0",
                               "L RBS_SAB" = "#a6c8b2", "H RBS_SAB" = "#f0cfa0",
                               "L RBS_REB" = "#a6c8b2", "H RBS_REB" = "#f0cfa0")) +
  # scale_y_discrete(labels = c("L SRS_Aware" = "社交意识", "H SRS_Aware" = "",
  #                             "L SRS_Cogni" = "社交认知", "H SRS_Cogni" = "",
  #                             "L SRS_Commu" = "社交表达", "H SRS_Commu" = "",
  #                             "L SRS_Motiv" = "社交动机", "H SRS_Motiv" = "",
  #                             "L SRS_Manne" = "社交举止", "H SRS_Manne" = "")) +
  scale_x_continuous(limits = c(NA, 15), breaks = seq(0, 15, by = 5)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
dev.off()


# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_RBS_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}


# 总分
var <- All[, c("clusterID", "RBS_TOTAL")]
colnames(var)[2] <- c("RBS_Total")
var[, 2][var[, 2] < 0] <- NA
var <- na.omit(var)


group <- c("L", "H")
varia <- names(var)[2]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

ggplot(var, aes(x = RBS_Total, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  scale_x_continuous(limits = c(NA, 75), breaks = seq(0, 75, by = 25)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 15))

name <- paste0("differ_RBS_total_", resDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)


# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
                             ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_RBS_total_", resDate, ".csv")
write.csv(sta_ana, file.path(resuDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}


objects_to_keep <- c("plotDir", "resDate", "resuDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part Z：保存P值文件 ##############################################
Pvalue <- Pvalue[, -1]

name <- paste0("statis_Pvalue_", resDate, ".csv")
write.csv(Pvalue, file.path(resuDir, name))
