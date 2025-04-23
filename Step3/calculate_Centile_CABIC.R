# Calculate OoS centile scores for all the brain measurements for CABIC data
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# Part 1: Run the code from LBCC to calculate centile scores
# Part 2: Extract the centile scores from results
###################################################

rm(list=ls())
source("E:/PhDproject/LBCC/lbcc/920.calc-novel-wo-subset-function.r")
packages <- c("tidyverse","mgcv","stringr","reshape2","magrittr","ggplot2","dplyr","readxl",
              "stringr","ggseg","patchwork","effectsize","pwr","coin", "gamlss", "parallel")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
# define filefolder
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result")
resultsDir <- file.path(resuDir, "Centile")
lifeSpanDir <- "E:/PhDproject/LBCC/lbcc"
resDate <- "240928"

# get metric names
cabic_forCentile <- read.csv(file.path(resuDir, paste0("cabic_forCentile_", resDate, ".csv")))
cabic_analysis <- cabic_forCentile[, -2] # 去掉Site列
volumeNames <- names(cabic_analysis)[c(which(names(cabic_analysis) == "GMV"):ncol(cabic_analysis))]
volumeNames <- volumeNames[-(which(volumeNames == "TCV_real"))] # 去掉真实的TCV用LBCC方法算出的TCV

######################### Part 1: Run the code from LBCC to calculate centile scores ###############
# 这段是在计算centile分，注意如果有以前的旧结果文件，需要先删掉，不然会跳过不计算哦

for (aparc in volumeNames){
  cat(aparc, "\n")
  fitRds <- file.path(lifeSpanDir, paste0("Model/FIT_", aparc, ".rds"))
  FitObj <- readRDS(file = fitRds)
  bootRds <-  file.path(lifeSpanDir, paste0("Model/BOOT_", aparc, ".rds"))
  BootObj <- readRDS(file = bootRds)
  resultCSV<- file.path(resultsDir, paste0(aparc, ".csv"))
  resultRDS <- file.path(resultsDir, paste0(aparc, ".rds"))
  # clean data
  cabic_info <- cabic_analysis[, 1:11]
  com_str <- paste0("cabic_info$", aparc, " <- cabic_analysis$", aparc)
  eval(parse(text = com_str))
  # deleta NA row
  index <- which(is.na(cabic_info[, 12]) == 1)
  if (length(index) != 0) {
    cabic_info <- cabic_info[-c(index), ]
  }
  # run!
  if (!file.exists(resultRDS)){
    RESULT <- Calc.Novel(cabic_info, Fit = FitObj, BootFit = BootObj, Apply = TRUE, NumberCores = 1)
    data <- RESULT$data
    write.csv(data[, c(1, 16:ncol(data))], resultCSV, row.names = F)
    saveRDS(RESULT, resultRDS)
    
  }
}

######################### Part 2: Extract the centile scores from results ##########################
# 提取出计算的Centile分，分年龄段保存
objects_to_keep <- c("cabicDir", "resuDir", "resDate", "cabic_forCentile", "volumeNames") # 列出要保留的对象名称
rm(list = (setdiff(ls(), objects_to_keep))) # 删除不在保留列表中的所有对象

# make a empty data frame to save data
FeatureCentile <- data.frame(cabic_forCentile$participant)
colnames(FeatureCentile) <- "participant"

for (feature in volumeNames) {
  centile <- paste0(feature, "Transformed.q.wre")
  data <- read.csv(file.path(resuDir, "Centile", paste0(feature,".csv")))
  data <- data[, c("participant", centile)]
  FeatureCentile <- merge(FeatureCentile, data, by = "participant", all.x = TRUE)
}

colnames(FeatureCentile)[2:length(FeatureCentile)] <- volumeNames

# add pheno
cabic_centile <- merge(cabic_forCentile[, c(1:3,6,12)], FeatureCentile, by = "participant", all = T)

# save result CABIC没有
name <- paste0("cabic_centile_", resDate, ".csv")
write.csv(cabic_centile, file.path(resuDir, name), row.names = F)


# ############## save the results of participants aging (1)~9.9
# cabic_centile_110 <- subset(cabic_centile, Age < 10)
# name <- paste0("cabic_centile_110_", resDate, ".csv")
# write.csv(cabic_centile_110, file.path(resuDir, name), row.names = F)

############## save the results of participants aging 5~9.9
cabic_centile_510 <- subset(cabic_centile, Age < 10 & Age >= 5)
name <- paste0("cabic_centile_510_", resDate, ".csv")
write.csv(cabic_centile_510, file.path(resuDir, name), row.names = F)
