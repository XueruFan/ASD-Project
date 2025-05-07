# Compare the cluster ID of ABIDE classifier
# Male, ASD, <13 years old
# Xue-Ru Fan, 12 MAY 2025 @BNU
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
abideDir <- 'E:/PhDproject/ABIDE' # Windows
resDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
newDate <- "250117"
oldDate <- "240610"

nyu <- read.csv(file.path(resDir, paste0("Cluster_NYU_only_", oldDate, ".csv")))
all <- read.csv(file.path(resDir, paste0("Cluster_", oldDate, ".csv")))

nyu <- nyu[, c(2,1)]
all <- all[, c(2,1)]

merged_data <- merge(nyu, all, by = "participant")
colnames(merged_data)[2:3] <- c("nyu", "all")

compare <- merged_data %>%
  summarise(
    Match = sum(nyu == all),
    Mismatch = sum(nyu != all),
    Match_Percent = mean(nyu == all) * 100,
    Mismatch_Percent = mean(nyu != all) * 100
  )

write.csv(compare, file.path(resDir, "compare_ClusterID_NYU_ALL.csv"), row.names = FALSE)
write.csv(merged_data, file.path(resDir, "ClusterID_NYU_ALL.csv"), row.names = FALSE)
