# Compare the cluster indix of the same participant from both spectral or gmm clustering
# (6~9.9 years)
# Xue-Ru Fan, 12 MAY 2025 @BNU
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

abide_spect <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]
abide_gmm <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/GMM513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_spect, abide_gmm, by = "participant")

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_513_Spect_Gmm.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "SpectClusterID", "ClusterID_GMM")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_513_Spect_Gmm.csv", 
          row.names = FALSE)
