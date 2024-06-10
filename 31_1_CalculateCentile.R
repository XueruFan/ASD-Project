# this script is used to calculate ABIDE I centile
# Xue-Ru Fan 25 April 2023 @BNU

rm(list=ls())
source("E:/PhDproject/LBCC/lbcc/920.calc-novel-wo-subset-function.r")
packages <- c("tidyverse","mgcv","stringr","reshape2","magrittr","ggplot2","dplyr","readxl",
              "stringr","ggseg","patchwork","effectsize","pwr","coin", "gamlss", "parallel")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
resultsDir <- file.path(abideDir, "Centile_1")
lifeSpanDir <- "E:/PhDproject/LBCC/lbcc"
resDate <- "240315"

abide_forCentile <- read.csv(file.path(dataDir, paste0("abide_1_forCentile_", resDate, ".csv")))
abide_analysis <- abide_forCentile[, -3:-2] # 去掉Site和Release列
volumeNames <- names(abide_analysis)[c(which(names(abide_analysis) == "GMV"):ncol(abide_analysis))]
volumeNames <- volumeNames[-(which(volumeNames == "TCV_real"))] 

################################# 报错分析 #########################################################
# volumeNames <- c("entorhinal", "middletemporal", "parsorbitalis") #没有计算出？？？？？？？
# 因为ABIDE1里的50818被试，在只分析abide1时，这三个分区的体积做完combat后变成了负数，
# 他是对照组，决定把ta删掉
abide_analysis <- abide_analysis[-which(abide_analysis[, "participant"] == "50818"), ]
########################################### Part1: Run Lbcc ########################################
# 这段是在计算centile分，注意如果有以前的旧结果文件，需要先删掉，不然会跳过不计算哦

for (aparc in volumeNames){
  # aparc <- "entorhinal","middletemporal","parsorbitalis" 
  cat(aparc, "\n")
  fitRds <- file.path(lifeSpanDir, paste0("Model/FIT_", aparc, ".rds"))
  FitObj <- readRDS(file = fitRds)
  bootRds <-  file.path(lifeSpanDir, paste0("Model/BOOT_", aparc, ".rds"))
  BootObj <- readRDS(file = bootRds)
  resultCSV<- file.path(resultsDir, paste0(aparc, ".csv"))
  resultRDS <- file.path(resultsDir, paste0(aparc, ".rds"))
  # clean data
  abide_info <- abide_analysis[, 1:11]
  com_str <- paste0("abide_info$", aparc, " <- abide_analysis$", aparc)
  eval(parse(text = com_str))
  # deleta NA row
  index <- which(is.na(abide_info[, 12]) == 1)
  if (length(index) != 0) {
    abide_info <- abide_info[-c(index), ]
  }
  # run!
  if (!file.exists(resultRDS)){
    RESULT <- Calc.Novel(abide_info, Fit = FitObj, BootFit = BootObj, Apply = TRUE, NumberCores = 1)
    data <- RESULT$data
    write.csv(data[, c(1, 16:ncol(data))], resultCSV, row.names = F)
    saveRDS(RESULT, resultRDS)
  }
}

########################################### Part2: Extract centile #################################
# 列出要保留的对象名称
objects_to_keep <- c("abideDir", "dataDir", "resDate", "abide_forCentile", "volumeNames")
rm(list = (setdiff(ls(), objects_to_keep))) # 删除不在保留列表中的所有对象

# make a empty data frame to save data
FeatureCentile <- data.frame(abide_forCentile$participant)
colnames(FeatureCentile) <- "participant"

for (feature in volumeNames) {
  centile <- paste0(feature, "Transformed.q.wre")
  data <- read.csv(file.path(abideDir, "Centile_1", paste0(feature,".csv")))
  data <- data[, c("participant", centile)]
  FeatureCentile <- merge(FeatureCentile, data, by = "participant", all.x = TRUE)
}

colnames(FeatureCentile)[2:length(FeatureCentile)] <- volumeNames

# add pheno
abide_centile <- merge(abide_forCentile[, c(1:4,7,13)], FeatureCentile, by = "participant", all = T)

# save result
write.csv(abide_centile, file.path(dataDir, paste0("abide_1_centile_", resDate, ".csv")), row.names = F)

############## save the results of participants aging 0~18
abide_centile_dev <- subset(abide_centile, Age < 18 & Age >= 6)
write.csv(abide_centile_dev, file.path(dataDir, paste0("abide_1_centile_dev_", resDate, ".csv")),
          row.names = F)
