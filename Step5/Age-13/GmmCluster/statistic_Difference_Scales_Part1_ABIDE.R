# Analyze demographic and cognitive behavioral differences between the clusters (Part A)
# Male, ASD, <13 years old, GMM Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################
# Part 01: Site
# Part 02: Manufacturer and scan
# Part 03: DSM tpye
# Part 1: Age
# Part 2: IQ
# SM Part 1: ADOS_G
# Part 3: ADOS_2
# Part 4: SRS
# Part 5: ADI_R
# SM Part 2: VINELAND
# Part 6: BMI
# Part Z：Save all the p values
################################


rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2", "Cairo")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
statiDir <- file.path(abideDir, "Analysis/Statistic/GMM513")
clustDir <- file.path(abideDir, "Analysis/Cluster/GMM513")
plotDir <- file.path(abideDir, "Plot/Cluster/GMM513")
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
Pvalue <- data.frame(matrix(ncol = 1, nrow = 7))
rownames(Pvalue) <- c("t-test", "t-df", "w-test", "cohend", "F-value", "F-value_p","f-test")
colnames(Pvalue) <- "variable"


################################# Part 01: Site ####################################################

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
# fisher精确检验
# Fisher 精确检验的一个特点是它并不依赖于卡方检验中的自由度概念，
# 因此 Fisher 精确检验没有明确的自由度。
# 这是因为 Fisher 精确检验主要用于精确计算小样本情况下的确切概率，而不依赖卡方检验中的大样本近似。
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
# 计算 Cramér's V 效应量
cramer_v_result <- effectsize::cramers_v(filtered_freq)
site_counts_sorted_to_save[1,5] <- paste0("Cramers_v_adjusted = ", as.numeric(cramer_v_result$Cramers_v_adjusted))

name <- paste0("statis_site_", newDate, ".csv")
write.csv(site_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(site_counts_sorted, id.vars = "Site")

sorted_Site <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Site)

data_long$Site <- factor(data_long$Site, levels = unique(sorted_Site))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# save plot
name <- paste0("differ_site_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)
ggplot(data_long, aes(x = Site, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("#facaae", "#b4d4c7")) +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 02: Manufacturer and scan ###################################

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
# 卡方检验
scan_pvalue <- chisq.test(filtered_freq)

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
# 计算 Cramér's V 效应量
cramer_v_result <- effectsize::cramers_v(filtered_freq)
scan_counts_sorted_to_save[1,5] <- paste0("Cramers_v_adjusted = ", as.numeric(cramer_v_result$Cramers_v_adjusted))
scan_counts_sorted_to_save[1,6] <- paste0("df = ", scan_pvalue$parameter)  # Fisher 检验的自由度

name <- paste0("statis_scaner_", newDate, ".csv")
write.csv(scan_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(scan_counts_sorted, id.vars = "Scan")

sorted_Scan <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Scan)

data_long$Scan <- factor(data_long$Scan, levels = unique(sorted_Scan))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))


name <- paste0("differ_scanner_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)
ggplot(data_long, aes(x = Scan, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("#facaae", "#b4d4c7")) +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

################## 再看厂家有没有差异

freq_var <- table(var$Manu, var$clusterID)
# 卡方检验
manu_pvalue <- chisq.test(freq_var)

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
# 计算 Cramér's V 效应量
cramer_v_result <- effectsize::cramers_v(freq_var)
manu_counts_sorted_to_save[1,5] <- paste0("Cramers_v_adjusted = ", as.numeric(cramer_v_result$Cramers_v_adjusted))
manu_counts_sorted_to_save[1,6] <- paste0("df = ", manu_pvalue$parameter)  # Fisher 检验的自由度

name <- paste0("statis_manu_", newDate, ".csv")
write.csv(manu_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(manu_counts_sorted, id.vars = "Manu")

sorted_Manu <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(Manu)

data_long$Manu <- factor(data_long$Manu, levels = unique(sorted_Manu))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

name <- paste0("differ_manu_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)
ggplot(data_long, aes(x = Manu, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("#facaae", "#b4d4c7")) +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 03: DSM type ################################################
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

# 卡方检验
type_pvalue <- chisq.test(freq_var)

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
# 计算 Cramér's V 效应量
cramer_v_result <- effectsize::cramers_v(freq_var)
type_counts_sorted_to_save[1,5] <- paste0("Cramers_v_adjusted = ", as.numeric(cramer_v_result$Cramers_v_adjusted))
type_counts_sorted_to_save[1,6] <- paste0("df = ", type_pvalue$parameter)  # Fisher 检验的自由度

name <- paste0("statis_type_", newDate, ".csv")
write.csv(type_counts_sorted_to_save, file.path(statiDir, name), row.names = F)

data_long <- melt(type_counts_sorted, id.vars = "type")

sorted_type <- data_long %>%
  filter(variable == "L") %>%
  arrange(desc(value)) %>%
  pull(type)

data_long$type <- factor(data_long$type, levels = unique(sorted_type))
data_long$variable <- factor(data_long$variable, levels = c("H", "L"))

# save plot
name <- paste0("differ_type_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)
ggplot(data_long, aes(x = type, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("#facaae", "#b4d4c7")) +
  # theme(text = element_text(family = "STSong"),
  theme(
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))
# dev.off()
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))



################################# Part 1: Age ######################################################
var <- All[, c("clusterID", "AGE_AT_SCAN")]
colnames(var)[2] <- "variable"

# save plot
name <- paste0("differ_age_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 5, height = 5, units = "in", dpi = 500)
ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = seq(5, 13, by = 2), labels = c("5", "7", "9", "11", "13")) +
  scale_fill_manual(values = c("L" = "#b4d4c7", "H" = "#facaae")) +
  # coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 12))
# dev.off()
ggsave(file.path(plotDir, name), width = 5, height = 4, units = "in", dpi = 500)

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

name <- paste0("statis_age_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

Pvalue["t-test","age"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["t-df", "age"] <- t.test(L$variable, H$variable)$parameter
Pvalue["w-test","age"] <- wilcox.test(L$variable, H$variable)[["p.value"]]
cohen_d_value <- effectsize::cohens_d(L$variable, H$variable)
Pvalue["cohend","age"] <- cohen_d_value$Cohens_d
Pvalue["F-value","age"] <- var.test(L$variable, H$variable)$statistic
Pvalue["F-value_p","age"] <- var.test(L$variable, H$variable)$p.value

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 2: IQ #######################################################
var <- All[, c("clusterID", "FIQ", "VIQ", "PIQ")]
var[, 2:4][var[, 2:4] < 0] <- NA

var_long <- gather(var, key = "measure", value = "variable", FIQ:PIQ, na.rm = TRUE,
                   factor_key = TRUE)
var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
var_long <- var_long[, -1]

# save plot
name <- paste0("differ_iq_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 6, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = c("L FIQ", "H FIQ", "L VIQ",
                                                                  "H VIQ", "L PIQ", "H PIQ")),
                     fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L FIQ" = "#b4d4c7", "H FIQ" = "#facaae",
                               "L VIQ" = "#b4d4c7", "H VIQ" = "#facaae",
                               "L PIQ" = "#b4d4c7", "H PIQ" = "#facaae")) +
  scale_x_continuous(breaks = seq(50, 150, by = 25)) +
  xlim(50, 150) +
  xlab("") +
  ylab("") +
  scale_y_discrete(labels = c("L FIQ" = "FIQ", "H FIQ" = "", 
                              "L VIQ" = "VIQ", "H VIQ" = "", 
                              "L PIQ" = "PIQ", "H PIQ" = "")) +  # 自定义y轴标签
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 5))
# dev.off()
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

name <- paste0("statis_iq_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['t-df', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")$parameter")))
  eval(parse(text = paste0("cohen_d_value <- effectsize::cohens_d(L$", v, ", H$", v, ")")))
  eval(parse(text = paste0("Pvalue['cohend', '", v, "'] <- cohen_d_value$Cohens_d")))
  eval(parse(text = paste0("Pvalue['F-value', '", v, "'] <- var.test(L$", v, ", H$", v, ")$statistic")))
  eval(parse(text = paste0("Pvalue['F-value_p', '", v, "'] <- var.test(L$", v, ", H$", v, ")$p.value")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))



################################# SM Part 1: ADOS_G ################################################
# var <- All[, c("clusterID", "ADOS_G_TOTAL", "ADOS_G_COMM", "ADOS_G_SOCIAL", "ADOS_G_STEREO_BEHAV")]
# colnames(var)[2:5] <- c("ADOS_G_Total", "ADOS_G_Commu", "ADOS_G_Socia", "ADOS_G_Stere")
# var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
# 
# var_long <- gather(var, key = "measure", value = "variable", names(var[2]):names(var[ncol(var)]),
#                    na.rm = TRUE, factor_key = TRUE)
# var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
# var_long <- var_long[, -1]
# 
# group <- c("L", "H")
# varia <- names(var)[2:ncol(var)]
# 
# rnames <- NULL
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     rnames <- c(rnames, rname)
#   }
# }
# 
# name <- paste0("differ_ADOS_G_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
# ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
#   geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
#   scale_fill_manual(values = c("L ADOS_G_Total" = "#b4d4c7", "H ADOS_G_Total" = "#facaae",
#                                "L ADOS_G_Commu" = "#b4d4c7", "H ADOS_G_Commu" = "#facaae",
#                                "L ADOS_G_Socia" = "#b4d4c7", "H ADOS_G_Socia" = "#facaae",
#                                "L ADOS_G_Stere" = "#b4d4c7", "H ADOS_G_Stere" = "#facaae")) +
#   coord_cartesian(xlim = c(NA, max(var_long$variable))) +
#   xlab("") +
#   ylab("") +
#   theme_ridges() +
#   theme(text = element_text(family = "STSong"),
#         legend.position = "none", # without legend
#         axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
#         axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# 
# 
# 
# # 统计检验
# L <- subset(var, clusterID == "L")
# H <- subset(var, clusterID == "H")
# 
# ##### mean and sd
# sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
# colnames(sta_ana) <- evalu
# rownames(sta_ana) <- rev(rnames)
# sta_ana$Count <- NA
# sum <- summary(as.factor(var_long$measure))
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v, ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
#   }
# }
# sta_ana$Median <- round(sta_ana$Median)
# 
# name <- paste0("statis_ADOS_G_", newDate, ".csv")
# write.csv(sta_ana, file.path(statiDir, name))
# 
# for (v in varia) {
#   # v = "ADOS_G_Total"
#   eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
#   eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
#   # 创建一个新的数据框，将 L 和 H 的数据合并，并加上 group 列
#   combined_data <- data.frame(
#     group = c(rep("L", length(eval(parse(text = paste0("L$", v))))), 
#               rep("H", length(eval(parse(text = paste0("H$", v)))))),
#     variable = c(eval(parse(text = paste0("L$", v))), 
#                  eval(parse(text = paste0("H$", v))))
#   )
#   
#   # 构建列联表
#   table_data <- table(combined_data$group, combined_data$variable)
#   
#   # 进行 Fisher 检验
#   eval(parse(text = paste0("Pvalue['f-test', '", v, "'] <- fisher.test(table_data, 
#                            simulate.p.value = TRUE, B = 1e5)[['p.value']]")))
# }
# 
# objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
# rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 3: ADOS_2 ###################################################
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

name <- paste0("differ_ADOS_2_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ADOS_2_Sever" = "#b4d4c7", "H ADOS_2_Sever" = "#facaae",
                               "L ADOS_2_Total" = "#b4d4c7", "H ADOS_2_Total" = "#facaae",
                               "L ADOS_2_Socia" = "#b4d4c7", "H ADOS_2_Socia" = "#facaae",
                               "L ADOS_2_Restr" = "#b4d4c7", "H ADOS_2_Restr" = "#facaae")) +
  coord_cartesian(xlim = c(NA, 25)) +
  # scale_y_discrete(labels = c("L ADOS_2_Sever" = "严重程度", "H ADOS_2_Sever" = "",
  #                             "L ADOS_2_Total" = "总分", "H ADOS_2_Total" = "",
  #                             "L ADOS_2_Socia" = "社交互动", "H ADOS_2_Socia" = "",
  #                             "L ADOS_2_Restr" = "刻板行为", "H ADOS_2_Restr" = "")) +
  scale_y_discrete(labels = c("L ADOS_2_Sever" = "Severity", "H ADOS_2_Sever" = "",
                              "L ADOS_2_Total" = "Total", "H ADOS_2_Total" = "",
                              "L ADOS_2_Socia" = "Social Affect", "H ADOS_2_Socia" = "",
                              "L ADOS_2_Restr" = "RRB", "H ADOS_2_Restr" = "")) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# dev.off()
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
name <- paste0("statis_ADOS_2_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['t-df', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")$parameter")))
  eval(parse(text = paste0("cohen_d_value <- effectsize::cohens_d(L$", v, ", H$", v, ")")))
  eval(parse(text = paste0("Pvalue['cohend', '", v, "'] <- cohen_d_value$Cohens_d")))
  eval(parse(text = paste0("Pvalue['F-value', '", v, "'] <- var.test(L$", v, ", H$", v, ")$statistic")))
  eval(parse(text = paste0("Pvalue['F-value_p', '", v, "'] <- var.test(L$", v, ", H$", v, ")$p.value")))
  # 创建一个新的数据框，将 L 和 H 的数据合并，并加上 group 列
  combined_data <- data.frame(
    group = c(rep("L", length(eval(parse(text = paste0("L$", v))))), 
              rep("H", length(eval(parse(text = paste0("H$", v)))))),
    variable = c(eval(parse(text = paste0("L$", v))), 
                 eval(parse(text = paste0("H$", v))))
  )
  
  # 构建列联表
  table_data <- table(combined_data$group, combined_data$variable)
  
  # 进行 Fisher 检验
  eval(parse(text = paste0("Pvalue['f-test', '", v, "'] <- fisher.test(table_data, 
                           simulate.p.value = TRUE, B = 1e5)[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 4: SRS ######################################################
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
name <- paste0("differ_SRS_scale_raw_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L SRS_Aware_r" = "#b4d4c7", "H SRS_Aware_r" = "#facaae",
                               "L SRS_Cogni_r" = "#b4d4c7", "H SRS_Cogni_r" = "#facaae",
                               "L SRS_Commu_r" = "#b4d4c7", "H SRS_Commu_r" = "#facaae",
                               "L SRS_Motiv_r" = "#b4d4c7", "H SRS_Motiv_r" = "#facaae",
                               "L SRS_Manne_r" = "#b4d4c7", "H SRS_Manne_r" = "#facaae")) +
  # scale_y_discrete(labels = c("L SRS_Aware_r" = "社交意识", "H SRS_Aware_r" = "",
  #                             "L SRS_Cogni_r" = "社交认知", "H SRS_Cogni_r" = "",
  #                             "L SRS_Commu_r" = "社交表达", "H SRS_Commu_r" = "",
  #                             "L SRS_Motiv_r" = "社交动机", "H SRS_Motiv_r" = "",
  #                             "L SRS_Manne_r" = "社交举止", "H SRS_Manne_r" = "")) +
  scale_y_discrete(labels = c("L SRS_Aware_r" = "Awareness", "H SRS_Aware_r" = "",
                              "L SRS_Cogni_r" = "Cognition", "H SRS_Cogni_r" = "",
                              "L SRS_Commu_r" = "Communication", "H SRS_Commu_r" = "",
                              "L SRS_Motiv_r" = "Motivation", "H SRS_Motiv_r" = "",
                              "L SRS_Manne_r" = "Mannerisms", "H SRS_Manne_r" = "")) +
  scale_x_continuous(breaks = seq(0, 60, by = 15)) +
  # xlim(50, 150) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# dev.off()
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

name <- paste0("statis_SRS_scale_raw_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['t-df', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")$parameter")))
  eval(parse(text = paste0("cohen_d_value <- effectsize::cohens_d(L$", v, ", H$", v, ")")))
  eval(parse(text = paste0("Pvalue['cohend', '", v, "'] <- cohen_d_value$Cohens_d")))
  eval(parse(text = paste0("Pvalue['F-value', '", v, "'] <- var.test(L$", v, ", H$", v, ")$statistic")))
  eval(parse(text = paste0("Pvalue['F-value_p', '", v, "'] <- var.test(L$", v, ", H$", v, ")$p.value")))
}


# # 几个分量表的标准分
# var <- All[, c("clusterID", "SRS_AWARENESS_T", "SRS_COGNITION_T", "SRS_COMMUNICATION_T",
#                "SRS_MOTIVATION_T", "SRS_MANNERISMS_T")]
# colnames(var)[2:6] <- c("SRS_Aware_t", "SRS_Cogni_t", "SRS_Commu_t", "SRS_Motiv_t", "SRS_Manne_t")
# var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
# 
# var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
#                    na.rm = TRUE, factor_key = TRUE)
# var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
# var_long <- var_long[, -1]
# 
# group <- c("L", "H")
# varia <- names(var)[2:ncol(var)]
# 
# rnames <- NULL
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     rnames <- c(rnames, rname)
#   }
# }
# 
# name <- paste0("differ_SRS_scale_trans_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
# ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
#   geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
#   scale_fill_manual(values = c("L SRS_Aware_t" = "#b4d4c7", "H SRS_Aware_t" = "#facaae",
#                                "L SRS_Cogni_t" = "#b4d4c7", "H SRS_Cogni_t" = "#facaae",
#                                "L SRS_Commu_t" = "#b4d4c7", "H SRS_Commu_t" = "#facaae",
#                                "L SRS_Motiv_t" = "#b4d4c7", "H SRS_Motiv_t" = "#facaae",
#                                "L SRS_Manne_t" = "#b4d4c7", "H SRS_Manne_t" = "#facaae")) +
#   scale_y_discrete(labels = c("L SRS_Aware_t" = "社交意识", "H SRS_Aware_t" = "",
#                               "L SRS_Cogni_t" = "社交认知", "H SRS_Cogni_t" = "",
#                               "L SRS_Commu_t" = "社交表达", "H SRS_Commu_t" = "",
#                               "L SRS_Motiv_t" = "社交动机", "H SRS_Motiv_t" = "",
#                               "L SRS_Manne_t" = "社交举止", "H SRS_Manne_t" = "")) +
#   scale_x_continuous(breaks = seq(40, 110, by = 20)) +
#   xlim(40,110) +
#   xlab("") +
#   ylab("") +
#   theme_ridges() +
#   theme(text = element_text(family = "STSong"),
#         legend.position = "none", # without legend
#         axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
#         axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# dev.off()
# 
# 
# # 统计检验
# L <- subset(var, clusterID == "L")
# H <- subset(var, clusterID == "H")
# 
# ##### mean and sd
# sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
# colnames(sta_ana) <- evalu
# rownames(sta_ana) <- rev(rnames)
# sta_ana$Count <- NA
# sum <- summary(as.factor(var_long$measure))
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
#   }
# }
# 
# name <- paste0("statis_SRS_scale_trans_", newDate, ".csv")
# write.csv(sta_ana, file.path(statiDir, name))
# 
# for (v in varia) {
#   eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
#   eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
# }
# 

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
  scale_fill_manual(values = c("L SRS_Total_t" = "#b4d4c7", "H SRS_Total_t" = "#facaae",
                               "L SRS_Total_r" = "#b4d4c7", "H SRS_Total_r" = "#facaae")) +
  coord_cartesian(xlim = c(NA, 200)) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 15))


name <- paste0("differ_SRS_total_", newDate, ".png")
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

name <- paste0("statis_SRS_total_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['t-df', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")$parameter")))
  eval(parse(text = paste0("cohen_d_value <- effectsize::cohens_d(L$", v, ", H$", v, ")")))
  eval(parse(text = paste0("Pvalue['cohend', '", v, "'] <- cohen_d_value$Cohens_d")))
  eval(parse(text = paste0("Pvalue['F-value', '", v, "'] <- var.test(L$", v, ", H$", v, ")$statistic")))
  eval(parse(text = paste0("Pvalue['F-value_p', '", v, "'] <- var.test(L$", v, ", H$", v, ")$p.value")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 5: ADI_R ####################################################
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

name <- paste0("differ_ADIR_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
  geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_fill_manual(values = c("L ADI_R_Socia" = "#b4d4c7", "H ADI_R_Socia" = "#facaae",
                               "L ADI_R_Verba" = "#b4d4c7", "H ADI_R_Verba" = "#facaae",
                               "L ADI_R_Nonve" = "#b4d4c7", "H ADI_R_Nonve" = "#facaae",
                               "L ADI_R_Restr" = "#b4d4c7", "H ADI_R_Restr" = "#facaae")) +
  # scale_y_discrete(labels = c("L ADI_R_Socia" = "社交互动", "H ADI_R_Socia" = "",
  #                             "L ADI_R_Verba" = "语言沟通", "H ADI_R_Verba" = "",
  #                             "L ADI_R_Nonve" = "非语言沟通", "H ADI_R_Nonve" = "",
  #                             "L ADI_R_Restr" = "刻板行为", "H ADI_R_Restr" = "")) +
  scale_y_discrete(labels = c("L ADI_R_Socia" = "Social", "H ADI_R_Socia" = "",
                              "L ADI_R_Verba" = "Verbal", "H ADI_R_Verba" = "",
                              "L ADI_R_Nonve" = "Nonverbal", "H ADI_R_Nonve" = "",
                              "L ADI_R_Restr" = "RRB", "H ADI_R_Restr" = "")) +
  # scale_x_continuous(breaks = seq(40, 110, by = 20)) +
  # xlim(40,110) +
  xlab("") +
  ylab("") +
  theme_ridges() +
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# dev.off()
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
name <- paste0("statis_ADIR_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

for (v in varia) {
  eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
                           ")[['p.value']]")))
  eval(parse(text = paste0("Pvalue['t-df', '", v, "'] <- t.test(L$", v, ", H$", v,
                           ")$parameter")))
  eval(parse(text = paste0("cohen_d_value <- effectsize::cohens_d(L$", v, ", H$", v, ")")))
  eval(parse(text = paste0("Pvalue['cohend', '", v, "'] <- cohen_d_value$Cohens_d")))
  eval(parse(text = paste0("Pvalue['F-value', '", v, "'] <- var.test(L$", v, ", H$", v, ")$statistic")))
  eval(parse(text = paste0("Pvalue['F-value_p', '", v, "'] <- var.test(L$", v, ", H$", v, ")$p.value")))
  # 创建一个新的数据框，将 L 和 H 的数据合并，并加上 group 列
  combined_data <- data.frame(
    group = c(rep("L", length(eval(parse(text = paste0("L$", v))))), 
              rep("H", length(eval(parse(text = paste0("H$", v)))))),
    variable = c(eval(parse(text = paste0("L$", v))), 
                 eval(parse(text = paste0("H$", v))))
  )
  
  # 构建列联表
  table_data <- table(combined_data$group, combined_data$variable)
  
  # 进行 Fisher 检验
  eval(parse(text = paste0("Pvalue['f-test', '", v, "'] <- fisher.test(table_data, 
                           simulate.p.value = TRUE, B = 1e5)[['p.value']]")))
}

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# SM Part 2: VINELAND ##############################################
# vineland_columns <- grep("VINELAND", names(All), value = TRUE)
# vineland_columns <- vineland_columns[-length(vineland_columns)]
# 
# # 先看五个标准分， 注意motor数据太少了只有9个人，这里去掉不看
# vineland_standard_columns <- grep("standard", vineland_columns, ignore.case = TRUE, value = TRUE)
# # 上面是为了输出这些列名，下面手动贴上需要的
# vineland_standard_columns <- c("VINELAND_ABC_Standard", "VINELAND_COMMUNICATION_STANDARD",
#                                "VINELAND_DAILYLIVING_STANDARD", "VINELAND_SOCIAL_STANDARD")
# var <- All[, c("clusterID", vineland_standard_columns)]
# colnames(var)[2:5] <- c("VIN_ABC_t", "VIN_Commu_t", "VIN_Daily_t", "VIN_Socia_t")
# var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
# 
# var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
#                    na.rm = TRUE, factor_key = TRUE)
# var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
# var_long <- var_long[, -1]
# 
# group <- c("L", "H")
# varia <- names(var)[2:ncol(var)]
# 
# rnames <- NULL
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     rnames <- c(rnames, rname)
#   }
# }
# 
# name <- paste0("differ_VIN_trans_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 7, height = 6, units = "in", dpi = 500)
# ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
#   geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
#   scale_fill_manual(values = c("L VIN_Commu_t" = "#b4d4c7", "H VIN_Commu_t" = "#facaae",
#                                "L VIN_Daily_t" = "#b4d4c7", "H VIN_Daily_t" = "#facaae",
#                                "L VIN_Socia_t" = "#b4d4c7", "H VIN_Socia_t" = "#facaae",
#                                "L VIN_ABC_t" = "#b4d4c7", "H VIN_ABC_t" = "#facaae")) +
#   scale_y_discrete(labels = c("L VIN_Commu_t" = "沟通", "H VIN_Commu_t" = "",
#                               "L VIN_Daily_t" = "日常生活", "H VIN_Daily_t" = "",
#                               "L VIN_Socia_t" = "社交", "H VIN_Socia_t" = "",
#                               "L VIN_ABC_t" = "适应行为综合得分", "H VIN_ABC_t" = "")) +
#   # scale_x_continuous(breaks = seq(40, 110, by = 20)) +
#   xlim(40,120) +
#   xlab("") +
#   ylab("") +
#   theme_ridges() +
#   theme(text = element_text(family = "STSong"),
#         legend.position = "none", # without legend
#         axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
#         axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# dev.off()
# 
# 
# # 统计检验
# L <- subset(var, clusterID == "L")
# H <- subset(var, clusterID == "H")
# 
# ##### mean and sd
# sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
# colnames(sta_ana) <- evalu
# rownames(sta_ana) <- rev(rnames)
# sta_ana$Count <- NA
# sum <- summary(as.factor(var_long$measure))
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
#   }
# }
# 
# name <- paste0("statis_VIN_trans_", newDate, ".csv")
# write.csv(sta_ana, file.path(statiDir, name))
# 
# for (v in varia) {
#   eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
#   eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
# }
# 
# 
# 
# # 再看量表分， 注意motor数据太少了只有9个人，这里去掉不看
# vineland_scaled_columns <- grep("scaled", vineland_columns, ignore.case = TRUE, value = TRUE)
# # 上面是为了输出这些列名，下面手动贴上需要的
# vineland_scaled_columns <- c("VINELAND_RECEPTIVE_V_SCALED", "VINELAND_EXPRESSIVE_V_SCALED",
#                              "VINELAND_WRITTEN_V_SCALED",  "VINELAND_PERSONAL_V_SCALED",
#                              "VINELAND_DOMESTIC_V_SCALED", "VINELAND_COMMUNITY_V_SCALED",
#                              "VINELAND_INTERPERSONAL_V_SCALED", "VINELAND_PLAY_V_SCALED",
#                              "VINELAND_COPING_V_SCALED")
# var <- All[, c("clusterID", vineland_scaled_columns)]
# colnames(var)[2:ncol(var)] <- c("VIN_Recep", "VIN_Expre", "VIN_Write", "VIN_Perso", "VIN_Domes",
#                                 "VIN_Commu", "VIN_Inter", "VIN_Play", "VIN_Copin")
# var[, 2:ncol(var)][var[, 2:ncol(var)] < 0] <- NA
# 
# var_long <- gather(var, key = "measure", value = "variable", names(var)[2]:names(var)[ncol(var)],
#                    na.rm = TRUE, factor_key = TRUE)
# var_long$measure <- paste0(var_long$clusterID, " ", var_long$measure)
# var_long <- var_long[, -1]
# 
# group <- c("L", "H")
# varia <- names(var)[2:ncol(var)]
# 
# rnames <- NULL
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     rnames <- c(rnames, rname)
#   }
# }
# 
# ggplot(var_long, aes(x = variable, y = factor(measure, levels = rnames), fill = measure)) +
#   geom_density_ridges(scale = 1.3, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
#   scale_fill_manual(values = c("L VIN_Recep" = "#b4d4c7", "H VIN_Recep" = "#facaae",
#                                "L VIN_Expre" = "#b4d4c7", "H VIN_Expre" = "#facaae",
#                                "L VIN_Write" = "#b4d4c7", "H VIN_Write" = "#facaae",
#                                "L VIN_Perso" = "#b4d4c7", "H VIN_Perso" = "#facaae",
#                                "L VIN_Domes" = "#b4d4c7", "H VIN_Domes" = "#facaae",
#                                "L VIN_Commu" = "#b4d4c7", "H VIN_Commu" = "#facaae",
#                                "L VIN_Inter" = "#b4d4c7", "H VIN_Inter" = "#facaae",
#                                "L VIN_Play" = "#b4d4c7", "H VIN_Play" = "#facaae",
#                                "L VIN_Copin" = "#b4d4c7", "H VIN_Copin" = "#facaae")) +
#   xlim(2, 22) +
#   xlab("") +
#   ylab("") +
#   theme_ridges() +
#   theme(legend.position = "none", # without legend
#         axis.text.y = element_text(size = 14, face = "bold", hjust = 0),
#         axis.text.x = element_text(size = 14, face = "bold", vjust = 5))
# 
# name <- paste0("differ_VIN_scale_", newDate, ".png")
# ggsave(file.path(plotDir, name), width = 7, height = 7, units = "in", dpi = 500)
# 
# 
# # 统计检验
# L <- subset(var, clusterID == "L")
# H <- subset(var, clusterID == "H")
# 
# ##### mean and sd
# sta_ana <- data.frame(matrix(ncol = length(evalu), nrow = length(group)*length(varia)))
# colnames(sta_ana) <- evalu
# rownames(sta_ana) <- rev(rnames)
# sta_ana$Count <- NA
# sum <- summary(as.factor(var_long$measure))
# for (v in varia) {
#   for (g in group) {
#     rname <- paste0(g, " ", v)
#     eval(parse(text = paste0("sta_ana[rname, 'Median'] <- round(median(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Mean'] <- round(mean(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'SD'] <- round(sd(", g, "$", v,
#                              ", na.rm = T), 2)")))
#     eval(parse(text = paste0("sta_ana[rname, 'Count'] <- sum[rname]")))
#   }
# }
# 
# name <- paste0("statis_VIN_scale_", newDate, ".csv")
# write.csv(sta_ana, file.path(statiDir, name))
# 
# for (v in varia) {
#   eval(parse(text = paste0("Pvalue['t-test', '", v, "'] <- t.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
#   eval(parse(text = paste0("Pvalue['w-test', '", v, "'] <- wilcox.test(L$", v, ", H$", v,
#                            ")[['p.value']]")))
# }



################################# Part 6: BMI ######################################################
var <- All[, c("clusterID", "BMI")]
colnames(var)[2] <- "variable"
var$variable[var$variable < 0] <- NA
var <- na.omit(var)

# save plot
name <- paste0("differ_bmi_", newDate, ".png")
# CairoPNG(file.path(plotDir, name), width = 6, height = 4, units = "in", dpi = 500)
ggplot(var, aes(x = variable, y = factor(clusterID, levels = c("L", "H")), fill = clusterID)) +
  geom_density_ridges(scale = 1.2, quantile_lines = TRUE, size = 1, quantiles = 4) +
  scale_x_continuous(breaks = c(13, 18, 23, 28, 33), limits = c(NA, 36)) +
  scale_fill_manual(values = c("L" = "#b4d4c7", "H" = "#facaae")) +
  # coord_fixed(ratio = 6) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(
  # theme(text = element_text(family = "STSong"),
        legend.position = "none", # without legend
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", vjust = 10))
# dev.off()
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

name <- paste0("statis_bmi_", newDate, ".csv")
write.csv(sta_ana, file.path(statiDir, name))

Pvalue["t-test","bmi"] <- t.test(L$variable, H$variable)[["p.value"]]
Pvalue["t-df", "bmi"] <- t.test(L$variable, H$variable)$parameter
Pvalue["w-test","bmi"] <- wilcox.test(L$variable, H$variable)[["p.value"]]
cohen_d_value <- effectsize::cohens_d(L$variable, H$variable)
Pvalue["cohend","bmi"] <- cohen_d_value$Cohens_d
Pvalue["F-value","bmi"] <- var.test(L$variable, H$variable)$statistic
Pvalue["F-value_p","bmi"] <- var.test(L$variable, H$variable)$p.value

objects_to_keep <- c("plotDir", "newDate", "statiDir", "All", "evalu", "Pvalue")
rm(list = (setdiff(ls(), objects_to_keep)))


################################# Part 7：Save all the p values ####################################
Pvalue <- Pvalue[, -1]

name <- paste0("statis_Pvalue_Part1_", newDate, ".csv")
write.csv(Pvalue, file.path(statiDir, name))
