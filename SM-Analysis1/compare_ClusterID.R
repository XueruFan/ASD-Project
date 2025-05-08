# Compare the cluster indix of the same participant from both <13 analysis or a narrow age analysis
# (6~9.9 years)
# Xue-Ru Fan, 12 MAY 2025 @BNU
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

abide_narow <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect509/Cluster_240610.csv")[, c(1,2)]
abide_broad <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_narow, abide_broad, by = "participant", all.abide_narow = T)

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_AgeRange_509_513.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "NarrowClusterID", "BroadClusterID")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_AgeRange_509_513.csv", 
          row.names = FALSE)



abide_narow <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect610/Cluster_240610.csv")[, c(1,2)]
abide_broad <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_narow, abide_broad, by = "participant", all.abide_narow = T)

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_AgeRange_610_513.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "NarrowClusterID", "BroadClusterID")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_AgeRange_610_513.csv", 
          row.names = FALSE)



abide_narow <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect711/Cluster_240610.csv")[, c(1,2)]
abide_broad <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_narow, abide_broad, by = "participant", all.abide_narow = T)

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_AgeRange_711_513.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "NarrowClusterID", "BroadClusterID")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_AgeRange_711_513.csv", 
          row.names = FALSE)


abide_narow <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect812/Cluster_240610.csv")[, c(1,2)]
abide_broad <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_narow, abide_broad, by = "participant", all.abide_narow = T)

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_AgeRange_812_513.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "NarrowClusterID", "BroadClusterID")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_AgeRange_812_513.csv", 
          row.names = FALSE)



abide_narow <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect913/Cluster_240610.csv")[, c(1,2)]
abide_broad <- read.csv("E:/PhDproject/ABIDE/Analysis/Cluster/Spect513/Cluster_240610.csv")[, c(1,2)]

merged_data <- merge(abide_narow, abide_broad, by = "participant", all.abide_narow = T)
merged_data$clusterID.x[merged_data$clusterID.x == "1"] <- "H"
merged_data$clusterID.x[merged_data$clusterID.x == "2"] <- "1"
merged_data$clusterID.x[merged_data$clusterID.x == "H"] <- "2"

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID.x == clusterID.y),
    Mismatch = sum(clusterID.x != clusterID.y),
    Match_Percent = mean(clusterID.x == clusterID.y) * 100,
    Mismatch_Percent = mean(clusterID.x != clusterID.y) * 100
  )

write.csv(compare, "E:/PhDproject/ABIDE/Analysis/Cluster/compare_ClusterID_AgeRange_913_513.csv", 
          row.names = FALSE)
compare
colnames(merged_data) <- c("Participant", "NarrowClusterID", "BroadClusterID")
write.csv(merged_data, "E:/PhDproject/ABIDE/Analysis/Cluster/ABIDE_ClusterID_AgeRange_913_513.csv", 
          row.names = FALSE)
