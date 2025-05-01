# Perform Shapiro-Wilk test on the scales bwtween clusters
# Male, ASD, 5~9.9 years old, Spectral Clustering
# Xue-Ru Fan 13 march 2024 @BNU
##################################################

rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2", "Cairo")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
phenoDir <- file.path(abideDir, "Preprocessed")
statiDir <- file.path(abideDir, "Analysis/Statistic/Spect510")
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect510")
resDate <- "240315"
newDate <- "240610"

pheno <- read.csv(file.path(phenoDir, paste0("abide_A_all_", resDate, ".csv")))
colnames(pheno)[1] <- "participant"

name <- paste0("Cluster_", newDate, ".csv")
cluster <- read.csv(file.path(clustDir, name))
colnames(cluster)[3:ncol(cluster)] <- paste0(colnames(cluster)[3:ncol(cluster)], "_centile")

All <- merge(cluster, pheno, by = "participant", all.x = TRUE)
All[which(All$clusterID == "1"), 'clusterID'] = "H"
All[which(All$clusterID == "2"), 'clusterID'] = "L"
All$clusterID <- factor(All$clusterID)

var <- All[, c("clusterID", "AGE_AT_SCAN", "FIQ", "VIQ", "PIQ", "ADOS_2_SEVERITY_TOTAL",
               "ADOS_2_TOTAL", "ADOS_2_SOCAFFECT", "ADOS_2_RRB", "SRS_AWARENESS_RAW",
               "SRS_COGNITION_RAW", "SRS_COMMUNICATION_RAW","SRS_MOTIVATION_RAW",
               "SRS_MANNERISMS_RAW", "SRS_TOTAL_RAW", "ADI_R_SOCIAL_TOTAL_A",
               "ADI_R_VERBAL_TOTAL_BV", "ADI_R_NONVERBAL_TOTAL_BV", "ADI_R_RRB_TOTAL_C")]
var[, -1][var[, -1] < 0] <- NA
L <- subset(var, clusterID == "L")
H <- subset(var, clusterID == "H")


result_L <- data.frame(matrix(ncol = 4, nrow = length(var)-1))
colnames(result_L) <- c("Cog", "clusterID", "W", "p-value")

for (n in 1:(ncol(var)-1)) {
  # n <- "AGE_AT_SCAN"
  temp <- L[, n+1]
  Result <- shapiro.test(temp)
  result_L[n, "clusterID"] <- "L"
  result_L[n, "Cog"] <- names(var[n+1])
  result_L[n,"W"] <- Result[["statistic"]][["W"]]
  result_L[n,"p-value"] <- Result[["p.value"]]
}


result_H <- data.frame(matrix(ncol = 4, nrow = length(var)-1))
colnames(result_H) <- c("Cog", "clusterID", "W", "p-value")

for (n in 1:(ncol(var)-1)) {
  # n <- "AGE_AT_SCAN"
  temp <- H[, n+1]
  Result <- shapiro.test(temp)
  result_H[n, "clusterID"] <- "H"
  result_H[n, "Cog"] <- names(var[n+1])
  result_H[n,"W"] <- Result[["statistic"]][["W"]]
  result_H[n,"p-value"] <- Result[["p.value"]]
}

result <- rbind(result_L, result_H)

result[, 3:4] <- round(result[, 3:4], digits = 4)

name <- paste0("statis_Shapiro_Wilk_Test_", newDate, ".csv")
write.csv(result, file.path(statiDir, name), row.names = F)
