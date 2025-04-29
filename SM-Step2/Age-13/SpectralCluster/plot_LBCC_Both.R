# Plot LBCC-corrected 34 regional volumes with GAM smoothing
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

cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result/pred/513")
resDate <- "240928"
abideDir <- 'E:/PhDproject/ABIDE' # winds
clustDir <- file.path(abideDir, "Analysis/Cluster/Spect513")
plotDir <- file.path(abideDir, "Plot/Cluster/Spect513/LBCC_Both")
newDate <- "240610"

# define variables
id_group <- c("2") # this code is for 2 clusters
ageRange <- log((seq(1, 13, 0.1) * 365.245) + 280)

name <- paste0("cabic_cluster_predictions_513_", resDate, ".csv")
cabic <- read.csv(file.path(resuDir, name))
cabic <- subset(cabic, predicted_cluster == "2")
cabic$Data <- "CABIC"

name <- paste0("Cluster_", newDate, ".csv")
abide <- read.csv(file.path(clustDir, name))
abide <- subset(abide, clusterID == "2")
abide$Data <- "ABIDE"

volumeNames <- names(cabic)[c(which(names(cabic) == "bankssts"):
                                     which(names(cabic) == "middletemporal"))]

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
  # wre代表脑图表值random effect，transformed是原始值除以10000，Transformed.normalised是校正后的值
  ynames <- paste0(volumeName, c("Transformed.m500.wre", "Transformed.normalised"))
  
  
  ########## load lbcc curve data 灰色的线
  POP.CURVE.LIST <- list(AgeTransformed = seq(log(120), log((365.25*86)), length.out = 2^8),
                         sex = c("Female","Male"))
  POP.CURVE.RAW <- do.call(what = expand.grid, args = POP.CURVE.LIST)
  CURVE <- Apply.Param(NEWData = POP.CURVE.RAW, FITParam = FIT$param)
  CURVE$age <- (exp(CURVE$AgeTransformed) - 280) / 365.25
  curveData <- CURVE %>% dplyr::select(AgeTransformed, age, sex, PRED.l025.pop, PRED.l250.pop,
                                       PRED.m500.pop, PRED.u750.pop, PRED.u975.pop)
  curveDataLong <- melt(curveData, id = c("AgeTransformed", "age", "sex"),
                        measure.vars = c("PRED.l025.pop", "PRED.l250.pop", "PRED.m500.pop",
                                         "PRED.u750.pop", "PRED.u975.pop"))
  curveDataLongMale <- curveDataLong %>% filter(sex == 'Male')
  
  studyFIT_abide <- read.csv(file.path(abideDir, "Centile", paste0(volumeName,".csv")))
  studyFIT_abide <- merge(abide, studyFIT_abide, by = "participant", all.x = T)
  
  studyFIT_cabic <- read.csv(file.path(cabicDir, "result/Centile", paste0(volumeName,".csv")))
  studyFIT_cabic <- merge(cabic, studyFIT_cabic, by = "participant", all.x = T)
  
  studyFIT1 <- studyFIT_abide
  studyFIT2 <- studyFIT_cabic
  
  
  # ############################# 常模校正后的ABIDE 黑色线
  # studyFIT$y <- studyFIT[, ynames[1]]
  # newFit <- expand.grid(AgeTransformed = ageRange)
  # fm <- gam(y ~ s(AgeTransformed, k = 4), data = studyFIT)
  # newFit$y <- predict(fm, newdata = newFit, se.fit = F)
  # plotData_ABIDE <- newFit
  
  
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
  C1_PointData$Data <- "ABIDE"
  C2_PointData <- studyFIT2[, c("AgeTransformed", ynames[2])]
  colnames(C2_PointData)[2] <- "y"
  C2_PointData$cluster <- "CABIC"
  
  # 设置年龄区间
  ageTicks <- seq(1, 13, by = 2)  # 原始年龄间隔为2年
  ageLimits <- c(1, 13)  # 限制绘图的年龄范围为5到13岁
  
  # plot
  lbcc_male <- subset(curveDataLongMale, variable == "PRED.m500.pop")
  
  
  # 画常模和两个cluster的校正值散点及拟合的线，不画ABIDE经过矫正后的组平均线
  ggplot() +
    # 添加散点图
    geom_point(data = C1_PointData, color = "#d9ca39", aes(x = (exp(AgeTransformed) - 280) / 365.25, y = y),
               alpha = .2, size = 2, shape = 16) +
    geom_point(data = C2_PointData, color = "#e29135", aes(x = (exp(AgeTransformed) - 280) / 365.25, y = y),
               alpha = .2, size = 2, shape = 16) +
    # 添加实际拟合的线
    geom_line(data = plotData_1, lwd = 2, alpha = 1, color = "#d9ca39",
              aes(x = (exp(AgeTransformed) - 280) / 365.25, y = y)) +
    geom_line(data = plotData_2, lwd = 2, alpha = 1, color = "#e29135",
              aes(x = (exp(AgeTransformed) - 280) / 365.25, y = y)) +
    
    # 绘制LBCC常模
    geom_line(data = lbcc_male, color = "gray",
              lwd = 2, alpha = .5, aes(x = (exp(AgeTransformed) - 280) / 365.25, y = value, group = sex, linetype = sex)) +
    
    theme_cowplot() +
    
    # 修改x轴为原始年龄，间隔为2年
    scale_x_continuous(limits = c(1, 13), breaks = seq(1, 13, by = 2),
                       labels = c("1 yr", "3", "5", "7", "9", "11", "13")) +
    xlab("") + ylab("") +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_text(size = 8)) +
    theme(legend.position = "None")
  
  
  
  
  name <- paste0("lbcc_", volumeName, "_", resDate,".png")
  ggsave(file.path(plotDir, name), dpi = 300, width = 10, height = 10, unit = "cm")
  
}
