rm(list=ls())
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr", "stats", "gamlss")
# sapply(packages,instAll.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
cabicDir <- 'E:/PhDproject/CABIC/result/pred/513/Corr'
resuDir <- file.path(abideDir, "Analysis/Statistic/Spect513")
resDate <- "240928"
newDate <- "240610"




cabic_h <- read.csv(file.path(cabicDir, paste0("corr_Part6_Select_H_", resDate, ".csv")))
cabic_l <- read.csv(file.path(cabicDir, paste0("corr_Part6_Select_L_", resDate, ".csv")))

abide_h <- read.csv(file.path(resuDir, paste0("corr_Part6_Select_H_", newDate, ".csv")))
abide_l <- read.csv(file.path(resuDir, paste0("corr_Part6_Select_L_", newDate, ".csv")))

cabic_h$name_brain <- paste0(cabic_h$name_brain, "_centile")
h_all <- merge(cabic_h, abide_h, by = "name_brain")
h_all <- subset(h_all, sign(coef.x) == sign(coef.y))

cabic_l$name_brain <- paste0(cabic_l$name_brain, "_centile")
l_all <- merge(cabic_l, abide_l, by = "name_brain")
l_all <- subset(l_all, sign(coef.x) == sign(coef.y))





cabic_h <- read.csv(file.path(cabicDir, paste0("corr_gamm_H_", resDate, ".csv")))
cabic_l <- read.csv(file.path(cabicDir, paste0("corr_gamm_L_", resDate, ".csv")))

abide_h <- read.xlsx(file.path(resuDir, paste0("corr_gamm_H_", newDate, ".xlsx")))
abide_l <- read.xlsx(file.path(resuDir, paste0("corr_gamm_L_", newDate, ".xlsx")))

cabic_h$brain <- paste0(cabic_h$brain, "_centile")
h_all <- merge(cabic_h, abide_h, by = "brain")

cabic_l$brain <- paste0(cabic_l$brain, "_centile")
l_all <- merge(cabic_l, abide_l, by = "brain")
l_all <- l_all[, c(1,2,13, 3,4,11, 14,15,22)]
