# 本代码用来分析分析高斯混合聚类的两组ASD男性的人口学和认知行为之间的差异
#（Part B，这部分是只有ABIDE2才有）
# 雪如 2024年2月27日于北师大办公室
################################
# Part 9: RBSR
# Part 10: MASC
# Part 11: BRIEF
# Part 12: CBCL，保存统计的csv文件
# Part 13: WIAT
# Part 14: CPRS
################
# Part Z：保存P值文件csv，之后需要手动excel，筛选出P值显著（＜0.05）的位置，保存一个xlsx文件
################################

rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect515")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect515")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect515")
resDate <- "240315"
newDate <- "240610"

pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[1] <- "participant"

name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
colnames(cluster)[3:ncol(cluster)] <- paste0(colnames(cluster)[3:ncol(cluster)], "_centile")

All <- merge(cluster, pheno, by = "participant", all.x = TRUE)
All[which(All$clusterID == "1"), 'clusterID'] = "L"
All[which(All$clusterID == "2"), 'clusterID'] = "H"
All$clusterID <- factor(All$clusterID)

evalu <- c("Median", "Mean", "SD") # 计算哪些统计值

# 新建一个空数据框用来存p值
Pvalue <- data.frame(matrix(ncol = 1, nrow = 2))
rownames(Pvalue) <- c("t-test", "w-test")
colnames(Pvalue) <- "variable"


################################# Part 9: RBSR #####################################################
pick_columns <- grep("RBSR", names(All), value = TRUE)

var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
nona <- colSums(!is.na(var[, -1])) # 从第二列到最后一列计算非NA值的数量

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}


################################# Part 10: MASC ###################################################
pick_columns <- grep("MASC", names(All), value = TRUE)
pick_columns
# 上面是为了输出这些列名，下面手动贴上需要的
pick_columns <- pick_columns[-1]
var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
nona <- colSums(!is.na(var[, -1])) # 从第二列到最后一列计算非NA值的数量
nona

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}


################################# Part 11: BRIEF ###################################################
pick_columns <- grep("BRIEF", names(All), value = TRUE)
pick_columns
# 上面是为了输出这些列名，下面手动贴上需要的
pick_columns <- pick_columns[-1:-2]
pick_columns
var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
var_long <- var[, c(1,9)]
colnames(var_long)[2] <- "variable"

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

name <- paste0("differ_BRIEF_ORG_", newDate, ".png")
CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = seq(30, 80, by = 25)) +
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

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}

##### mean and sd
L <- data.frame(L[,9])
colnames(L) <- "variable"
H <- data.frame(H[,9])
colnames(H) <- "variable"
sta_ana <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- c("H", "L")
sta_ana[1,1] <- round(median(H$variable, na.rm = T), 2)
sta_ana[1,2] <- round(mean(H$variable, na.rm = T), 2)
sta_ana[1,3] <- round(sd(H$variable, na.rm = T), 2)
sta_ana[1,4] <- round(sum(!is.na(H[[1]])))
sta_ana[2,1] <- round(median(L$variable, na.rm = T), 2)
sta_ana[2,2] <- round(mean(L$variable, na.rm = T), 2)
sta_ana[2,3] <- round(sd(L$variable, na.rm = T), 2)
sta_ana[2,4] <- round(sum(!is.na(L[[1]])))
colnames(sta_ana)[4] <- "Count"

name <- paste0("statis_BRIEF_ORG_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

################################# Part 12: CBCL ###################################################
pick_columns <- grep("CBCL_6.18", names(All), value = TRUE)
pick_columns

var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
nona <- colSums(!is.na(var[, -1])) # 从第二列到最后一列计算非NA值的数量
nona

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}


##### mean and sd
var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)], na.rm = TRUE,
                   factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]

sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu
rownames(sta_ana) <- rev(rnames)
sta_ana$Count <- NA
sum <- summary(as.factor(var_long$measure))
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v, ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v, ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v, ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("statis_CBCL_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))


################################# Part 13: WIAT ###################################################
pick_columns <- grep("WIAT", names(All), value = TRUE)
pick_columns

var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
nona <- colSums(!is.na(var[, -1])) # 从第二列到最后一列计算非NA值的数量
nona

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}


################################# Part 14: CPRS ###################################################
pick_columns <- grep("CPRS", names(All), value = TRUE)
pick_columns

var <- All[, c("clusterID", pick_columns)]

var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
nona <- colSums(!is.na(var[, -1])) # 从第二列到最后一列计算非NA值的数量
nona

group <- c("L", "H")
varia <- names(var)[2:ncol(var)]

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v, ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v, ")[['p.value']]")))
}



################################# 保存P值文件 ######################################################
Pvalue <- Pvalue[, -1]

name <- paste0("statis_Pvalue_Part2_", newDate, ".csv")
write.csv(Pvalue, file.path(statiDir, name))
