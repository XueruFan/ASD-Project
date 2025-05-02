# Compare the cluster ID of CABIC dataset between predicting ID with ABIDE classifier and do the
# clustering on CABIC independently 
# Male, ASD, <13 years old
# Xue-Ru Fan, 12 MAY 2025 @BNU
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr", "do", "dplyr", "Cairo", "openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # MAC
cabicDir <- 'E:/PhDproject/CABIC'
predDir <- file.path(cabicDir, "result/pred/Spect513")
indeDir <- file.path(cabicDir, "result/inde/Spect-13")
newDate <- "250117"
oldDate <- "240928"

pred <- read.csv(file.path(predDir, paste0("cabic_cluster_predictions_513_", oldDate, ".csv")))
inde <- read.csv(file.path(indeDir, paste0("Cluster_", oldDate, ".csv")))

# inde的结果里，cluster1 是H，2是L，pred里cluster2 是H，1是L

inde$clusterID[inde$clusterID == "1"] <- "H"
inde$clusterID[inde$clusterID == "2"] <- "1"
inde$clusterID[inde$clusterID == "H"] <- "2"

pred <- pred[, c(1,45)]
inde <- inde[, c(2,1)]

merged_data <- merge(inde, pred, by = "participant")

compare <- merged_data %>%
  summarise(
    Match = sum(clusterID == predicted_cluster),
    Mismatch = sum(clusterID != predicted_cluster),
    Match_Percent = mean(clusterID == predicted_cluster) * 100,
    Mismatch_Percent = mean(clusterID != predicted_cluster) * 100
  )

write.csv(compare, file.path(cabicDir, "compare_ClusterID.csv"), row.names = FALSE)
