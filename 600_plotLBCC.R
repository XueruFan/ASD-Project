# 本代码用来可视化两个聚类人群脑指标的发育轨线
# 画LBCC和两个cluster的校正值散点及拟合的线，不画ABIDE经过矫正后的组平均线
# 画线的时候使用GAM拟合的散点
# Xue-Ru Fan 04 Jan 2024 @BNU
##################################
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
clustDir <- file.path(abideDir, "Analysis/Cluster/SpectralCluster")
plotDir <- file.path(abideDir, "Plot/Cluster/SpectralCluster/LBCC")
newDate <- "240610"

# define variables
id_group <- c("1", "2") # this code is for 2 clusters
ageRange <- log((seq(6, 18, 0.1)*365.245)+280)

name <- paste0("asd_male_Spectral_Cluster_", newDate, ".csv")
exampleFIT <- read.csv(file.path(clustDir, name))

volumeNames <- names(exampleFIT)[c(which(names(exampleFIT) == "GMV"):which(names(exampleFIT) == "insula"))]

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
  ynames <- paste0(volumeName, c("Transformed.m500.wre", "Transformed.normalised"))
  

  ########## load lbcc curve data 灰色的线
  POP.CURVE.LIST <- list(AgeTransformed = seq(log(120), log((365.25*86)), length.out = 2^8),
                         sex = c("Female","Male"))
  POP.CURVE.RAW <- do.call(what = expand.grid, args = POP.CURVE.LIST)
  CURVE <- Apply.Param(NEWData = POP.CURVE.RAW, FITParam = FIT$param)
  # 这里会报错，fs_version的数据缺失，校正的轨线不准确
  CURVE$age <- (exp(CURVE$AgeTransformed) - 280) / 365.25
  curveData <- CURVE %>% dplyr::select(AgeTransformed, age, sex, PRED.l025.pop, PRED.l250.pop,
                                       PRED.m500.pop, PRED.u750.pop, PRED.u975.pop)
  curveDataLong <- melt(curveData, id = c("AgeTransformed", "age", "sex"),
                        measure.vars = c("PRED.l025.pop", "PRED.l250.pop", "PRED.m500.pop",
                                         "PRED.u750.pop", "PRED.u975.pop"))
  curveDataLongMale <- curveDataLong %>% filter(sex == 'Male')
  
  
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
  fm <- gam(y ~ s(AgeTransformed, k = 2), data = studyFIT1)
  newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  plotData_1 <- newFit
  
  # cluster2
  studyFIT2$y <- studyFIT2[, ynames[2]]
  newFit <- expand.grid(AgeTransformed = ageRange)
  fm <- gam(y ~ s(AgeTransformed, k = 2), data = studyFIT2)
  newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  plotData_2 <- newFit
  
  C1_PointData <- studyFIT1[, c("AgeTransformed", ynames[2])]
  colnames(C1_PointData)[2] <- "y"
  C1_PointData$cluster <- "1"
  C2_PointData <- studyFIT2[, c("AgeTransformed", ynames[2])]
  colnames(C2_PointData)[2] <- "y"
  C2_PointData$cluster <- "2"
  # plotPoint <- rbind(C1_PointData, C2_PointData)
  
  
  # 设置年龄区间
  ageTicks <- log((c(6,8,10,12,14,16,18)*365.245)+280)
  ageLimits <- log((c(6,18)*365.245)+280)
  
  # plot
  lbcc_male <- subset(curveDataLongMale, variable == "PRED.m500.pop")
  
  
  # 画常模和两个cluster的校正值散点及拟合的线，不画ABIDE经过矫正后的组平均线
  ggplot(plotData_ABIDE, aes(x = AgeTransformed, y = y)) +
    # geom_line(lwd = 2, alpha = .5, color = "gray20") +
    # 添加散点图
    geom_point(data = C1_PointData, color = "#719988", aes(x = AgeTransformed, y = y),
               alpha = .2, size = 2, shape = 16) +
    geom_point(data = C2_PointData, color = "#faa264", aes(x = AgeTransformed, y = y),
               alpha = .2, size = 2, shape = 16) +
    # 添加实际拟合的线
    geom_line(data = plotData_1, lwd = 2, alpha = 1, color = "#719988",
              aes(x = AgeTransformed, y = y)) +
    geom_line(data = plotData_2, lwd = 2, alpha = 1, color = "#faa264",
              aes(x = AgeTransformed, y = y)) +
    
    # 绘制LBCC常模
    geom_line(data = lbcc_male, color = "gray",
              lwd = 2, alpha = .5, aes(x = AgeTransformed, y = value, group = sex, linetype = sex)) +

    theme_cowplot() +
    
    scale_x_continuous(limits = ageLimits, breaks = ageTicks,
                       label = c("6 yr", "8 yr", "10 yr", "12 yr", "14 yr", "16 yr", "18 yr")) +
    xlab("") + ylab("") +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_text(size = 8)) +
    theme(legend.position = "None")
  
  name <- paste0("SC_LBCC_", volumeName, "_", newDate,".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 10, height = 10, unit = "cm")
  
}
