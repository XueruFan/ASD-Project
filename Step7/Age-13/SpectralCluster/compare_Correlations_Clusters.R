# Compare the correlations results between ABIDE and CABIC
# Male, ASD, <13 years old, Spectral Clustering
# Xue-Ru Fan 24 Oct 2023 @BNU
################################

rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot", "openxlsx",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
cabicDir <- 'E:/PhDproject/CABIC/result/pred/513/Corr'
resuDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
resDate <- "240928"
newDate <- "240610"


cabic_h <- read.csv(file.path(cabicDir, paste0("corr_Part4_Select_H_", resDate, ".csv")))
cabic_l <- read.csv(file.path(cabicDir, paste0("corr_Part4_Select_L_", resDate, ".csv")))

abide_h <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_H_", newDate, ".csv")))
abide_l <- read.csv(file.path(resuDir, paste0("corr_Part4_Select_L_", newDate, ".csv")))

cabic_h$name_brain <- paste0(cabic_h$name_brain, "_centile")
h_all <- merge(cabic_h, abide_h, by = "name_brain")
h_all <- subset(h_all, sign(coef.x) == sign(coef.y))

cabic_l$name_brain <- paste0(cabic_l$name_brain, "_centile")
l_all <- merge(cabic_l, abide_l, by = "name_brain")
l_all <- subset(l_all, sign(coef.x) == sign(coef.y))

write.xlsx(l_all, file.path(resuDir, "l_corr_part4_bothABIDEandCABIC.xlsx"), rowNames = F)
write.xlsx(h_all, file.path(resuDir, "h_corr_part4_bothABIDEandCABIC.xlsx"), rowNames = F)



