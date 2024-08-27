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
# Part 7: VINELAND
# Part 8: BMI
################
# 以上每部分都会保存统计的csv文件和绘图的png文件
################
# Part Z：保存P值文件csv，之后需要手动excel，筛选出P值显著（＜0.05）的位置，保存一个xlsx文件
################################


rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
statiDir <- file.path(abideDir, "Analysis/Statistic/GmmCluster")
clustDir <- file.path(abideDir, "Analysis/Cluster/GmmCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/GmmCluster")
resDate <- "240315"
newDate <- "240610"

pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[1] <- "participant"

name <- paste0("asd_male_GMM_Cluster_", newDate, ".csv")
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


################################# Part 01: 站点 ####################################################

var <- All[, c("clusterID", "SITE_ID")]
# 修改各自站点名字的缩写
var$SITE_ID <- gsub("ABIDEII-NYU_1", "NYU", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-NYU_2", "NYU", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-KKI_1", "KKI", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-SDSU_1", "SDSU", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-UCLA_1", "UCLA", var$SITE_ID)
var$SITE_ID <- gsub("UCLA_1", "UCLA", var$SITE_ID)
var$SITE_ID <- gsub("UCLA_2", "UCLA", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-GU_1", "GU", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-UCD_1", "UCD", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-EMC_1", "EMC", var$SITE_ID)
var$SITE_ID <- gsub("TRINITY", "TCD", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-TCD_1", "TCD", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-USM_1", "USM", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-IU_1", "IU", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-U_MIA_1", "UMIA", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-ETH_1", "ETH", var$SITE_ID)
var$SITE_ID <- gsub("UM_1", "UM", var$SITE_ID)
var$SITE_ID <- gsub("UM_2", "UM", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-OHSU_1", "OHSU", var$SITE_ID)
var$SITE_ID <- gsub("STANFORD", "SU1", var$SITE_ID)
var$SITE_ID <- gsub("ABIDEII-SU_2", "SU2", var$SITE_ID)
var$SITE_ID <- gsub("LEUVEN_2", "KUL", var$SITE_ID)
var$SITE_ID <- gsub("CALTECH", "CALT", var$SITE_ID)

var$SITE_ID <- as.factor(var$SITE_ID)
var$clusterID <- as.factor(var$clusterID)

freq_var <- table(var$SITE_ID, var$clusterID)
# 过滤掉那些频数和小于10的行
filtered_freq <- freq_var[rowSums(freq_var) >= 10, ]
# fisher检验
site_pvalue <- fisher.test(filtered_freq, simulate.p.value = TRUE, B = 1e5)

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

L_site_counts <- as.data.frame(table(L$SITE_ID)) # 计算SITE_ID的每个水平的频数，并转换为数据框
H_site_counts <- as.data.frame(table(H$SITE_ID))
site_count <- merge(L_site_counts, H_site_counts, by = "Var1")
site_counts_sorted <- arrange(site_count, desc(Freq.x)) # 将结果按Count列从大到小排序
colnames(site_counts_sorted) <- c("Site", "L", "H")
# 只看样本量大于等于10的站点
site_counts_sorted$All <- site_counts_sorted$L + site_counts_sorted$H
site_counts_sorted <- subset(site_counts_sorted, site_counts_sorted$All >= 10)
site_counts_sorted <- site_counts_sorted[, -4]

# 将p值绑定到数据框的右边
site_counts_sorted_to_save <- site_counts_sorted
site_counts_sorted_to_save[1,4] <- paste0("p = ", site_pvalue[["p.value"]])

name <- paste0("asd_male_dev_GC_statis_site_", newDate, ".csv")
write.csv(site_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(site_counts_sorted, id.vars = "Site")

sorted_Site <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Site)

data_long$Site <- factor(data_long$Site, levels = unique(sorted_Site))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# plot
ggplot(data_long, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#f0cfa0", "#a6c8b2")) +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1))
# save plot
name <- paste0("asd_male_dev_GC_differ_site_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 02: 机型 ####################################################

var <- All[, c("clusterID", "SITE_ID")]
# 添加上机型
var <- var %>%
  mutate(
    Scan = case_when(
      SITE_ID == 'ABIDEII-NYU_1' ~ 'S Allegra',
      SITE_ID == 'ABIDEII-NYU_2' ~ 'S Allegra',
      SITE_ID == 'NYU' ~ 'S Allegra',
      SITE_ID == 'ABIDEII-KKI_1' ~ 'P Achieva',
      SITE_ID == 'ABIDEII-SDSU_1' ~ 'G MR 750',
      SITE_ID == 'ABIDEII-UCLA_1' ~ 'S TrioTim',
      SITE_ID == 'UCLA_1' ~ 'S TrioTim',
      SITE_ID == 'UCLA_2' ~ 'S TrioTim',
      SITE_ID == 'ABIDEII-GU_1' ~ 'S TrioTim',
      SITE_ID == 'ABIDEII-UCD_1' ~ 'S TrioTim',
      SITE_ID == 'ABIDEII-EMC_1' ~ 'G MR 750',
      SITE_ID == 'TRINITY' ~ 'P Achieva',
      SITE_ID == 'ABIDEII-USM_1' ~ 'S TrioTim',
      SITE_ID == 'ABIDEII-IU_1' ~ 'S TrioTim',
      SITE_ID == 'ABIDEII-U_MIA_1' ~ 'G Healthcare',
      SITE_ID == 'ABIDEII-ETH_1' ~ 'P Achieva',
      SITE_ID == 'ABIDEII-TCD_1' ~ 'P Achieva',
      SITE_ID == 'UM_1' ~ 'G Signa',
      SITE_ID == 'UM_2' ~ 'G Signa',
      SITE_ID == 'ABIDEII-OHSU_1' ~ 'S TrioTim',
      SITE_ID == 'STANFORD' ~ 'G Signa',
      SITE_ID == 'ABIDEII-SU_2' ~ 'G Signa',
      SITE_ID == 'LEUVEN_2' ~ 'P Intera',
      SITE_ID == 'CALTECH' ~ 'S TrioTim',
      SITE_ID == 'PITT' ~ 'S Allegra',
      SITE_ID == 'OLIN' ~ 'S Allegra',
      SITE_ID == 'OHSU' ~ 'S TrioTim',
      SITE_ID == 'SDSU' ~ 'G MR 750',
      SITE_ID == 'USM' ~ 'S TrioTim',
      SITE_ID == 'YALE' ~ 'S TrioTim',
      SITE_ID == 'KKI' ~ 'P Achieva',
      TRUE ~ '其他' # 默认值
    ),
    Manu = case_when(
      SITE_ID == 'ABIDEII-NYU_1' ~ 'S',
      SITE_ID == 'ABIDEII-NYU_2' ~ 'S',
      SITE_ID == 'NYU' ~ 'S',
      SITE_ID == 'ABIDEII-KKI_1' ~ 'P',
      SITE_ID == 'ABIDEII-SDSU_1' ~ 'G',
      SITE_ID == 'ABIDEII-UCLA_1' ~ 'S',
      SITE_ID == 'UCLA_1' ~ 'S',
      SITE_ID == 'UCLA_2' ~ 'S',
      SITE_ID == 'ABIDEII-GU_1' ~ 'S',
      SITE_ID == 'ABIDEII-UCD_1' ~ 'S',
      SITE_ID == 'ABIDEII-EMC_1' ~ 'G',
      SITE_ID == 'TRINITY' ~ 'P',
      SITE_ID == 'ABIDEII-USM_1' ~ 'S',
      SITE_ID == 'ABIDEII-IU_1' ~ 'S',
      SITE_ID == 'ABIDEII-U_MIA_1' ~ 'G',
      SITE_ID == 'ABIDEII-ETH_1' ~ 'P',
      SITE_ID == 'ABIDEII-TCD_1' ~ 'P',
      SITE_ID == 'UM_1' ~ 'G',
      SITE_ID == 'UM_2' ~ 'G',
      SITE_ID == 'ABIDEII-OHSU_1' ~ 'S',
      SITE_ID == 'STANFORD' ~ 'G',
      SITE_ID == 'ABIDEII-SU_2' ~ 'G',
      SITE_ID == 'LEUVEN_2' ~ 'P',
      SITE_ID == 'CALTECH' ~ 'S',
      SITE_ID == 'PITT' ~ 'S',
      SITE_ID == 'OLIN' ~ 'S',
      SITE_ID == 'OHSU' ~ 'S',
      SITE_ID == 'SDSU' ~ 'G',
      SITE_ID == 'USM' ~ 'S',
      SITE_ID == 'YALE' ~ 'S',
      SITE_ID == 'KKI' ~ 'P',
      TRUE ~ '其他' # 默认值
    )
  )

var$SITE_ID <- as.factor(var$SITE_ID)
var$clusterID <- as.factor(var$clusterID)
var$Scan <- as.factor(var$Scan)
var$Manu <- as.factor(var$Manu)

################## 先看机型有没有差异

freq_var <- table(var$Scan, var$clusterID)
# 过滤掉那些频数和小于30（其实也就是10）的机型
filtered_freq <- freq_var[rowSums(freq_var) >= 30, ]
# fisher检验
scan_pvalue <- fisher.test(filtered_freq, simulate.p.value = TRUE, B = 1e5)

L <- subset(var[, c(1,3)], clusterID == "L")
H <- subset(var[, c(1,3)], clusterID == "H")

L_scan_counts <- as.data.frame(table(L$Scan))
H_scan_counts <- as.data.frame(table(H$Scan))
scan_count <- merge(L_scan_counts, H_scan_counts, by = "Var1")
scan_counts_sorted <- arrange(scan_count, desc(Freq.x))
colnames(scan_counts_sorted) <- c("Scan", "L", "H")
# 只看样本量大于等于30的机型
scan_counts_sorted$All <- scan_counts_sorted$L + scan_counts_sorted$H
scan_counts_sorted <- subset(scan_counts_sorted, scan_counts_sorted$All >= 30)
scan_counts_sorted <- scan_counts_sorted[, -4]

# 将p值绑定到数据框的右边
scan_counts_sorted_to_save <- scan_counts_sorted
scan_counts_sorted_to_save[1,4] <- paste0("p = ", scan_pvalue[["p.value"]])

name <- paste0("asd_male_dev_GC_statis_scaner_", newDate, ".csv")
write.csv(scan_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(scan_counts_sorted, id.vars = "Scan")

sorted_Scan <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Scan)

data_long$Scan <- factor(data_long$Scan, levels = unique(sorted_Scan))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# plot
ggplot(data_long, aes(x = Scan, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#f0cfa0", "#a6c8b2")) +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1))
# save plot
name <- paste0("asd_male_dev_GC_differ_scanner_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

################## 再看厂家有没有差异

freq_var <- table(var$Manu, var$clusterID)
# fisher检验
manu_pvalue <- fisher.test(filtered_freq, simulate.p.value = TRUE, B = 1e5)

L <- subset(var[, c(1,4)], clusterID == "L")
H <- subset(var[, c(1,4)], clusterID == "H")

L_manu_counts <- as.data.frame(table(L$Manu))
H_manu_counts <- as.data.frame(table(H$Manu))
manu_count <- merge(L_manu_counts, H_manu_counts, by = "Var1")
manu_counts_sorted <- arrange(manu_count, desc(Freq.x))
colnames(manu_counts_sorted) <- c("Manu", "L", "H")

# 将p值绑定到数据框的右边
manu_counts_sorted_to_save <- manu_counts_sorted
manu_counts_sorted_to_save[1,4] <- paste0("p = ", manu_pvalue[["p.value"]])

name <- paste0("asd_male_dev_GC_statis_manu_", newDate, ".csv")
write.csv(manu_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(manu_counts_sorted, id.vars = "Manu")

sorted_Manu <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Manu)

data_long$Manu <- factor(data_long$Manu, levels = unique(sorted_Manu))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# plot
ggplot(data_long, aes(x = Manu, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#f0cfa0", "#a6c8b2")) +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))

# save plot
name <- paste0("asd_male_dev_GC_differ_manu_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 03: 类型 #####################################################
var <- All[, c("clusterID", "PDD_DSM_IV_TR")]
colnames(var)[2] <- "variable"
var$variable[var$variable < 0] <- NA
var <- na.omit(var)
var <- subset(var, variable != "0")

var[which(var$variable == "1"), "variable"] <- "Autism"
var[which(var$variable == "2"), "variable"] <- "Aspergers"
var[which(var$variable == "3"), "variable"] <- "PDD-NOS"

var$variable <- as.factor(var$variable)

freq_var <- table(var$variable, var$clusterID)

# fisher检验
type_pvalue <- fisher.test(freq_var, simulate.p.value = TRUE, B = 1e5)

L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

L_type_counts <- as.data.frame(table(L$variable))
H_type_counts <- as.data.frame(table(H$variable))
type_count <- merge(L_type_counts, H_type_counts, by = "Var1")
type_counts_sorted <- arrange(type_count, desc(Freq.x)) # 将结果按Count列从大到小排序
colnames(type_counts_sorted) <- c("type", "L", "H")

# 将p值绑定到数据框的右边
type_counts_sorted_to_save <- type_counts_sorted
type_counts_sorted_to_save[1,4] <- paste0("p = ", type_pvalue[["p.value"]])

name <- paste0("asd_male_dev_GC_statis_type_", newDate, ".csv")
write.csv(type_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(type_counts_sorted, id.vars = "type")

sorted_type <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(type)

data_long$type <- factor(data_long$type, levels = unique(sorted_type))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# plot
ggplot(data_long, aes(x = type, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#f0cfa0", "#a6c8b2")) +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
# save plot
name <- paste0("asd_male_dev_GC_differ_type_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))



################################# Part 1: 年龄 #####################################################
var <- All[, c("clusterID", "AGE_AT_SCAN")]
colnames(var)[2] <- "variable"

ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = seq(6, 18, by = 3), labels = c("6 yrs", "9", "12", "15", "18")) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 10))
# save plot
name <- paste0("asd_male_dev_GC_differ_age_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)

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

name <- paste0("asd_male_dev_GC_statis_age_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

Pvalue["t-test","age"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["w-test","age"] <- wilcox.test(L$variable, H$variable)[["p.value"]]

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 2: IQ #######################################################
var <- All[, c("clusterID", "FIQ", "VIQ", "PIQ")]
var[, 2:4][var[, 2:4] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", FIQ:PIQ, na.rm = TRUE,
                   factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]

ggplot(var_long, aes(x = variable, y = factor(measure, levels = c("L FIQ", "H FIQ", "L VIQ",
                                                                  "H VIQ", "L PIQ", "H PIQ")),
                     fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L FIQ" = "#a6c8b2", "H FIQ" = "#f0cfa0",
                               "L VIQ" = "#a6c8b2", "H VIQ" = "#f0cfa0",
                               "L PIQ" = "#a6c8b2", "H PIQ" = "#f0cfa0")) +
  xlim(50, 160) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 5))

# save plot
name <- paste0("asd_male_dev_GC_differ_iq_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 6, units = "in", dpi = 500)

# 统计检验
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")

group <- c("L", "H")
varia <- c("FIQ", "VIQ", "PIQ")

##### mean and sd
sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
colnames(sta_ana) <- evalu

rnames <- NULL
for (v in varia) {
  for (g in group) {
    rname <- paste0(g, " ", v)
    rnames <- c(rnames, rname)
  }
}
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
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v, ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}

name <- paste0("asd_male_dev_GC_statis_iq_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 3: ADOS_G ###################################################
var <- All[, c("clusterID", "ADOS_G_TOTAL", "ADOS_G_COMM", "ADOS_G_SOCIAL", "ADOS_G_STEREO_BEHAV")]
colnames(var)[2:5] <- c("ADOS_G_Total", "ADOS_G_Commu", "ADOS_G_Socia", "ADOS_G_Stere")
var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", names(var[2]):names(var[ncol(var)]),
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ADOS_G_Total" = "#a6c8b2", "H ADOS_G_Total" = "#f0cfa0",
                               "L ADOS_G_Commu" = "#a6c8b2", "H ADOS_G_Commu" = "#f0cfa0",
                               "L ADOS_G_Socia" = "#a6c8b2", "H ADOS_G_Socia" = "#f0cfa0",
                               "L ADOS_G_Stere" = "#a6c8b2", "H ADOS_G_Stere" = "#f0cfa0")) +
  coord_cartesian(xlim = c(NA, max(var_long$variable))) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_ADOS_G_", newDate, ".png")
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
    eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v, ", na.rm = T), 2)")))
    eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
  }
}
sta_ana$Median <- round(sta_ana$Median)

name <- paste0("asd_male_dev_GC_statis_ADOS_G_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 4: ADOS_2 ###################################################
var <- All[, c("clusterID", "ADOS_2_SEVERITY_TOTAL", "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT",
               "ADOS_2_RRB")]
colnames(var)[2:5] <- c("ADOS_2_Sever", "ADOS_2_Total", "ADOS_2_Socia", "ADOS_2_Restr")
var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
                   na.rm = TRUE,  factor_key = TRUE)
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ADOS_2_Sever" = "#a6c8b2", "H ADOS_2_Sever" = "#f0cfa0",
                               "L ADOS_2_Total" = "#a6c8b2", "H ADOS_2_Total" = "#f0cfa0",
                               "L ADOS_2_Socia" = "#a6c8b2", "H ADOS_2_Socia" = "#f0cfa0",
                               "L ADOS_2_Restr" = "#a6c8b2", "H ADOS_2_Restr" = "#f0cfa0")) +
  coord_cartesian(xlim = c(NA, 25)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_ADOS_2_", newDate, ".png")
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
name <- paste0("asd_male_dev_GC_statis_ADOS_2_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 5: SRS ###################################################
# 几个分量表的原始分
var <- All[, c("clusterID", "SRS_AWARENESS_RAW", "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW",
               "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW")]
colnames(var)[2:6] <- c("SRS_Aware_r", "SRS_Cogni_r", "SRS_Commu_r", "SRS_Motiv_r", "SRS_Manne_r")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L SRS_Aware_r" = "#a6c8b2", "H SRS_Aware_r" = "#f0cfa0",
                               "L SRS_Cogni_r" = "#a6c8b2", "H SRS_Cogni_r" = "#f0cfa0",
                               "L SRS_Commu_r" = "#a6c8b2", "H SRS_Commu_r" = "#f0cfa0",
                               "L SRS_Motiv_r" = "#a6c8b2", "H SRS_Motiv_r" = "#f0cfa0",
                               "L SRS_Manne_r" = "#a6c8b2", "H SRS_Manne_r" = "#f0cfa0")) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_SRS_scale_raw_", newDate, ".png")
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

name <- paste0("asd_male_dev_GC_statis_SRS_scale_raw_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}


# 几个分量表的标准分
var <- All[, c("clusterID", "SRS_AWARENESS_T", "SRS_COGNITION_T", "SRS_COMMUNICATION_T",
               "SRS_MOTIVATION_T", "SRS_MANNERISMS_T")]
colnames(var)[2:6] <- c("SRS_Aware_t", "SRS_Cogni_t", "SRS_Commu_t", "SRS_Motiv_t", "SRS_Manne_t")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L SRS_Aware_t" = "#a6c8b2", "H SRS_Aware_t" = "#f0cfa0",
                               "L SRS_Cogni_t" = "#a6c8b2", "H SRS_Cogni_t" = "#f0cfa0",
                               "L SRS_Commu_t" = "#a6c8b2", "H SRS_Commu_t" = "#f0cfa0",
                               "L SRS_Motiv_t" = "#a6c8b2", "H SRS_Motiv_t" = "#f0cfa0",
                               "L SRS_Manne_t" = "#a6c8b2", "H SRS_Manne_t" = "#f0cfa0")) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_SRS_scale_trans_", newDate, ".png")
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

name <- paste0("asd_male_dev_GC_statis_SRS_scale_trans_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}


# 总分的原始分和标准分
var <- All[, c("clusterID", "SRS_TOTAL_T", "SRS_TOTAL_RAW")]
colnames(var)[2:3] <- c("SRS_Total_t", "SRS_Total_r")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L SRS_Total_t" = "#a6c8b2", "H SRS_Total_t" = "#f0cfa0",
                               "L SRS_Total_r" = "#a6c8b2", "H SRS_Total_r" = "#f0cfa0")) +
  coord_cartesian(xlim = c(NA, 200)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 15))

name <- paste0("asd_male_dev_GC_differ_SRS_total_", newDate, ".png")
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

name <- paste0("asd_male_dev_GC_statis_SRS_total_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 6: ADI_R ####################################################
var <- All[, c("clusterID", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
                    "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C")]
colnames(var)[2:5] <- c("ADI_R_Socia", "ADI_R_Verba", "ADI_R_Nonve", "ADI_R_Restr")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ADI_R_Socia" = "#a6c8b2", "H ADI_R_Socia" = "#f0cfa0",
                               "L ADI_R_Verba" = "#a6c8b2", "H ADI_R_Verba" = "#f0cfa0",
                               "L ADI_R_Nonve" = "#a6c8b2", "H ADI_R_Nonve" = "#f0cfa0",
                               "L ADI_R_Restr" = "#a6c8b2", "H ADI_R_Restr" = "#f0cfa0")) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_ADIR_", newDate, ".png")
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
name <- paste0("asd_male_dev_GC_statis_ADIR_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 7: VINELAND ###################################################
vineland_columns <- grep("VINELAND", names(All), value = TRUE)
vineland_columns <- vineland_columns[-length(vineland_columns)]

# 先看五个标准分， 注意motor数据太少了只有9个人，这里去掉不看
vineland_standard_columns <- grep("standard", vineland_columns, ignore.case = TRUE, value = TRUE)
# 上面是为了输出这些列名，下面手动贴上需要的
vineland_standard_columns <- c("VINELAND_ABC_Standard", "VINELAND_COMMUNICATION_STANDARD",
                               "VINELAND_DAILYLIVING_STANDARD", "VINELAND_SOCIAL_STANDARD")
var <- All[, c("clusterID", vineland_standard_columns)]
colnames(var)[2:5] <- c("VIN_ABC_t", "VIN_Commu_t", "VIN_Daily_t", "VIN_Socia_t")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L VIN_Commu_t" = "#a6c8b2", "H VIN_Commu_t" = "#f0cfa0",
                               "L VIN_Daily_t" = "#a6c8b2", "H VIN_Daily_t" = "#f0cfa0",
                               "L VIN_Socia_t" = "#a6c8b2", "H VIN_Socia_t" = "#f0cfa0",
                               "L VIN_ABC_t" = "#a6c8b2", "H VIN_ABC_t" = "#f0cfa0")) +
  xlim(40, 120) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_VIN_trans_", newDate, ".png")
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

name <- paste0("asd_male_dev_GC_statis_VIN_trans_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}



# 再看量表分， 注意motor数据太少了只有9个人，这里去掉不看
vineland_scaled_columns <- grep("scaled", vineland_columns, ignore.case = TRUE, value = TRUE)
# 上面是为了输出这些列名，下面手动贴上需要的
vineland_scaled_columns <- c("VINELAND_RECEPTIVE_V_SCALED", "VINELAND_EXPRESSIVE_V_SCALED",
                             "VINELAND_WRITTEN_V_SCALED",  "VINELAND_PERSONAL_V_SCALED",
                             "VINELAND_DOMESTIC_V_SCALED", "VINELAND_COMMUNITY_V_SCALED",
                             "VINELAND_INTERPERSONAL_V_SCALED", "VINELAND_PLAY_V_SCALED",
                             "VINELAND_COPING_V_SCALED")
var <- All[, c("clusterID", vineland_scaled_columns)]
colnames(var)[2:ncol(var)] <- c("VIN_Recep", "VIN_Expre", "VIN_Write", "VIN_Perso", "VIN_Domes",
                                "VIN_Commu", "VIN_Inter", "VIN_Play", "VIN_Copin")
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

ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L VIN_Recep" = "#a6c8b2", "H VIN_Recep" = "#f0cfa0",
                               "L VIN_Expre" = "#a6c8b2", "H VIN_Expre" = "#f0cfa0",
                               "L VIN_Write" = "#a6c8b2", "H VIN_Write" = "#f0cfa0",
                               "L VIN_Perso" = "#a6c8b2", "H VIN_Perso" = "#f0cfa0",
                               "L VIN_Domes" = "#a6c8b2", "H VIN_Domes" = "#f0cfa0",
                               "L VIN_Commu" = "#a6c8b2", "H VIN_Commu" = "#f0cfa0",
                               "L VIN_Inter" = "#a6c8b2", "H VIN_Inter" = "#f0cfa0",
                               "L VIN_Play" = "#a6c8b2", "H VIN_Play" = "#f0cfa0",
                               "L VIN_Copin" = "#a6c8b2", "H VIN_Copin" = "#f0cfa0")) +
  xlim(2, 22) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))

name <- paste0("asd_male_dev_GC_differ_VIN_scale_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)


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

name <- paste0("asd_male_dev_GC_statis_VIN_scale_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
}



################################# Part 8: BMI ######################################################
var <- All[, c("clusterID", "BMI")]
colnames(var)[2] <- "variable"
var$variable[var$variable < 0] <- NA
var <- na.omit(var)

ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = c(18.5, 24.9, 29.9), limits = c(NA, 30)) +
  scale_fill_manual(values = c("L" = "#a6c8b2", "H" = "#f0cfa0")) +
  # coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 10))
# save plot
name <- paste0("asd_male_dev_GC_differ_bmi_", newDate, ".png")
ggsave(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)

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

name <- paste0("asd_male_dev_GC_statis_bmi_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

Pvalue["t-test","bmi"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["w-test","bmi"] <- wilcox.test(L$variable, H$variable)[["p.value"]]

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part Z：保存P值文件 ##############################################
Pvalue <- Pvalue[, -1]

name <- paste0("asd_male_dev_GC_statis_Pvalue_Part1_", newDate, ".csv")
write.csv(Pvalue, file.path(statiDir, name))
