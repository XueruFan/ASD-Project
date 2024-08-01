# 画出谱聚类分出的ASD男性两个亚型的34个脑指标的发育轨线+散点图
# 不画常模，画ABIDE经过矫正后的组平均线，和两个cluster的原始值散点及拟合的线
# Xue-Ru Fan 04 Jan 2024 @BNU

rm(list=ls())

# load packages
# for windows
source("E:/PhDproject/LBCC/lbcc/920.calc-novel-wo-subset-function.r")
source("E:/PhDproject/LBCC/lbcc/R_rainclouds.R")
# # for mac
# source("/Volumes/Xueru/PhDproject/LBCC/lbcc/920.calc-novel-wo-subset-function.r")
# source("/Volumes/Xueru/PhDproject/LBCC/lbcc/R_rainclouds.R")
packages <- c("tidyverse", "mgcv", "stringr", "reshape2", "magrittr", "ggplot2", "dplyr", "readxl",
              "stringr", "ggseg", "patchwork", "effectsize", "pwr", "cowplot",
              "readr", "ggridges", "tidyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
# abideDir <- '/Volumes/Xueru/PhDproject/ABIDE' # mac
abideDir <- 'E:/PhDproject/ABIDE' # winds
clustDir <- file.path(abideDir, "Analysis/Cluster/Cluster_A/SpectralCluster34DK")
plotDir <- file.path(abideDir, "Plot/Cluster/Cluster_A/SpectralCluster34DK/GAMM")
resDate <- "240315"
newDate <- "240610"

# define variables
id_group <- c("1", "2") # this code is for 2 clusters
ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", newDate, ".csv")
exampleFIT <- read.csv(file.path(clustDir, name))

volumeNames <- names(exampleFIT)[c(which(names(exampleFIT) == "bankssts"):which(names(exampleFIT) ==
                                                                                  "insula"))]

# load lbcc
lifeSpanDir <- "E:/PhDproject/LBCC/lbcc"
# lifeSpanDir <- "/Volumes/Xueru/PhDproject/LBCC/lbcc"
setwd(lifeSpanDir)
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")

for (volumeName in volumeNames){
  # volumeName <- "GMV"
  FIT <- readRDS(paste0("Model/FIT_",volumeName,".rds"))
  # wre代表脑图表值randem effect，transformed是原始值除以10000，Transformed.normalised是校正后的值
  ynames <- paste0(volumeName, c("Transformed.m500.wre", "Transformed", "Transformed.normalised"))
  
  
  ########## load abide data
  studyFIT <- read.csv(file.path(abideDir, "Centile_A", paste0(volumeName,".csv")))
  studyFIT <- merge(exampleFIT, studyFIT, by = "participant", all.x = T)
  studyFIT1 <- subset(studyFIT, clusterID == "1")
  studyFIT2 <- subset(studyFIT, clusterID == "2")
  
  
  ############################# 常模校正后的ABIDE 黑色线
  studyFIT$y <- studyFIT[, ynames[1]]
  newFit <- expand.grid(AgeTransformed = ageRange)
  fm <- gam(y ~ s(AgeTransformed, k = 4), data = studyFIT)
  newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  plotData_ABIDE <- newFit
  
  
  ############################# asd个体的常模预测值
  
  # cluster1
  studyFIT1$y <- studyFIT1[, ynames[2]]
  newFit <- expand.grid(AgeTransformed = ageRange)
  fm <- gam(y ~ s(AgeTransformed, k = 4), data = studyFIT1)
  newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  plotData_1 <- newFit
  
  # cluster2
  studyFIT2$y <- studyFIT2[, ynames[2]]
  newFit <- expand.grid(AgeTransformed = ageRange)
  fm <- gam(y ~ s(AgeTransformed, k = 4), data = studyFIT2)
  newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  plotData_2 <- newFit
  
  C1_PointData <- studyFIT1[, c("AgeTransformed", ynames[2])]
  colnames(C1_PointData)[2] <- "y"
  C1_PointData$cluster <- "1"
  C2_PointData <- studyFIT2[, c("AgeTransformed", ynames[2])]
  colnames(C2_PointData)[2] <- "y"
  C2_PointData$cluster <- "2"
  
  
  # 设置年龄区间
  ageTicks <- log((c(6,8,10,12,14,16,18)*365.245)+280)
  ageLimits <- log((c(6,18)*365.245)+280)
  
  # 不画常模，画ABIDE经过矫正后的组平均线，和两个cluster的原始值散点及拟合的线
  ggplot(plotData_ABIDE, aes(x = AgeTransformed, y = y)) +
    geom_line(lwd = 2, alpha = .5, color = "gray20") +
    # 添加散点图
    geom_point(data = C1_PointData, color = "#0064b5", aes(x = AgeTransformed, y = y),
               alpha = .2, size = 2, shape = 16) +
    geom_point(data = C2_PointData, color = "#ff6347", aes(x = AgeTransformed, y = y),
               alpha = .2, size = 2, shape = 16) +
    # 添加实际拟合的线
    geom_line(data = plotData_1, lwd = 2, alpha = 1, color = "#0064b5",
              aes(x = AgeTransformed, y = y)) +
    geom_line(data = plotData_2, lwd = 2, alpha = 1, color = "#ff6347",
              aes(x = AgeTransformed, y = y)) +

    theme_cowplot() +
    scale_x_continuous(limits = ageLimits, breaks = ageTicks,
                       label = c("6 yr", "8 yr", "10 yr", "12 yr", "14 yr", "16 yr", "18 yr")) +
    xlab("") + ylab("") +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_text(size = 8)) +
    theme(legend.position = "None")

  name <- paste0("abide_A_asd_male_dev_Spectral_Cluster_34DK_", volumeName, "_GAMLSS_",
                 newDate,".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 10, height = 10, unit = "cm")
  
  
}
