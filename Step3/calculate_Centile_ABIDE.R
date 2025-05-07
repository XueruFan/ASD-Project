# Calculate OoS centile scores for all the brain measurements for ABIDE data
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
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
resultsDir <- file.path(abideDir, "Centile")
lifeSpanDir <- "E:/PhDproject/LBCC/lbcc"
resDate <- "240315"

# get metric names
abide_forCentile <- read.csv(file.path(dataDir, paste0("abide_All_forCentile_", resDate, ".csv")))
abide_analysis <- abide_forCentile[, -3:-2] # 去掉Site和Release列
volumeNames <- names(abide_analysis)[c(which(names(abide_analysis) == "GMV"):ncol(abide_analysis))]
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

######################### Part 2: Extract the centile scores from results ##########################
# 提取出计算的Centile分，分年龄段保存
objects_to_keep <- c("abideDir", "dataDir", "resDate", "abide_forCentile", "volumeNames") # 列出要保留的对象名称
rm(list = (setdiff(ls(), objects_to_keep))) # 删除不在保留列表中的所有对象

# make a empty data frame to save data
FeatureCentile <- data.frame(abide_forCentile$participant)
colnames(FeatureCentile) <- "participant"

for (feature in volumeNames) {
  centile <- paste0(feature, "Transformed.q.wre")
  data <- read.csv(file.path(abideDir, "Centile", paste0(feature,".csv")))
  data <- data[, c("participant", centile)]
  FeatureCentile <- merge(FeatureCentile, data, by = "participant", all.x = TRUE)
}

colnames(FeatureCentile)[2:length(FeatureCentile)] <- volumeNames

# add pheno
abide_centile <- merge(abide_forCentile[, c(1:4,7,13)], FeatureCentile, by = "participant", all = T)

# save result
name <- paste0("abide_All_centile_", resDate, ".csv")
write.csv(abide_centile, file.path(dataDir, name), row.names = F)

############## save the results of participants aging (5)~12.9
abide_centile_513 <- subset(abide_centile, Age < 13)
name <- paste0("abide_All_centile_513_", resDate, ".csv")
write.csv(abide_centile_513, file.path(dataDir, name), row.names = F)

############## save the results of participants aging (5)~9.9
abide_centile_510 <- subset(abide_centile, Age < 10)
name <- paste0("abide_All_centile_510_", resDate, ".csv")
write.csv(abide_centile_510, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 6~9.9
abide_centile_610 <- subset(abide_centile, Age < 10 & Age >= 6 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_610_", resDate, ".csv")
write.csv(abide_centile_610, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 10~12.9
abide_centile_112 <- subset(abide_centile, Age >= 10 & Age < 12 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_112_", resDate, ".csv")
write.csv(abide_centile_112, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 7~10.9
abide_centile_711 <- subset(abide_centile, Age < 11 & Age >= 7 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_711_", resDate, ".csv")
write.csv(abide_centile_711, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 8~11.9
abide_centile_812 <- subset(abide_centile, Age < 12 & Age >= 8 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_812_", resDate, ".csv")
write.csv(abide_centile_812, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 9~12.9
abide_centile_913 <- subset(abide_centile, Age < 13 & Age >= 9 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_913_", resDate, ".csv")
write.csv(abide_centile_913, file.path(dataDir, name), row.names = F)

############## save the results of participants aging 5~8.9
abide_centile_509 <- subset(abide_centile, Age < 9 & Age >= 5 & dx == "ASD" & sex == "Male")
name <- paste0("abide_All_centile_509_", resDate, ".csv")
write.csv(abide_centile_509, file.path(dataDir, name), row.names = F)
