# this script is used to do ComBat analysis for ABIDE data
# Xue-Ru Fan 16 May 2023 @BNU
###################################################
# Part 1: 把ABIDE 的1和2合并一起进行后续的分析，csv
# Part 2: 只分析ABIDE I的数据，csv
# Part 3：只分析ABIDE II的数据，csv
###################################################
rm(list=ls())
packages <- c("neuroCombat")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)
# 如果没有neuroCombat包，需要上github安装
# reference:https://doi.org/10.1016/j.neuroimage.2017.11.024
# library(devtools)
# install_github("jfortin1/neuroCombatData")
# install_github("jfortin1/neuroCombat_Rpackage")

abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
resDate <- "240315"

abide_all <- read.csv(file.path(dataDir,  paste0("abide_A_forComBat_", resDate, ".csv")))


################################## Part 1 ##################################
# 把ABIDE 的1和2合并一起进行后续的分析

abide <- abide_all
dat <- as.matrix(t(abide[, which(names(abide) == "GMV"):ncol(abide)]))
batch <- c(t(abide$Site)) # Numeric or character vector
age <- c(t(abide$age_days)) # Continuous variable
disease <- as.factor(c(t(abide$dx))) # Categorical variable
sex <- as.factor(c(t(abide$sex)))
mod <- model.matrix(~ age * sex + disease) # create model matrix for these two biological covariates

data.harmonized <- neuroCombat(dat = dat, batch = batch, mod = mod, parametric = FALSE,
                              eb = FALSE, mean.only = TRUE) # do combat

# get harmonized data
data_combated <- t(data.harmonized$dat.combat)
data_combated <- as.data.frame(data_combated)
# arrange for centile calculation
abide_forCentile <- cbind(abide[, 1:which(names(abide) == "dx")], data_combated)

# save result
name <-  paste0("abide_A_forCentile_", resDate, ".csv")
write.csv(abide_forCentile, file.path(dataDir, name), row.names = F)


################################## Part 2  ######################################
# 只分析ABIDE I的数据
abide <- subset(abide_all, Release == "1")
dat <- as.matrix(t(abide[, which(names(abide) == "GMV"):ncol(abide)]))
batch <- c(t(abide$Site)) # Numeric or character vector
age <- c(t(abide$age_days)) # Continuous variable
disease <- as.factor(c(t(abide$dx))) # Categorical variable
sex <- as.factor(c(t(abide$sex)))
mod <- model.matrix(~ age * sex + disease)
data.harmonized <- neuroCombat(dat = dat, batch = batch, mod = mod, parametric = FALSE,
                               eb = FALSE, mean.only = TRUE)
data_combated <- t(data.harmonized$dat.combat)
data_combated <- as.data.frame(data_combated)
abide_forCentile <- cbind(abide[, 1:which(names(abide) == "dx")], data_combated)
name <- paste0("abide_1_forCentile_", resDate, ".csv")
write.csv(abide_forCentile, file.path(dataDir, name), row.names = F)


################################## Part 3  #####################################
# 只分析ABIDE II的数据
abide <- subset(abide_all, Release == "2")
dat <- as.matrix(t(abide[, which(names(abide) == "GMV"):ncol(abide)]))
batch <- c(t(abide$Site)) # Numeric or character vector
age <- c(t(abide$age_days)) # Continuous variable
disease <- as.factor(c(t(abide$dx))) # Categorical variable
sex <- as.factor(c(t(abide$sex)))
mod <- model.matrix(~ age * sex + disease)
data.harmonized <- neuroCombat(dat = dat, batch = batch, mod = mod, parametric = FALSE,
                               eb = FALSE, mean.only = TRUE)
data_combated <- t(data.harmonized$dat.combat)
data_combated <- as.data.frame(data_combated)
abide_forCentile <- cbind(abide[, 1:which(names(abide) == "dx")], data_combated)
name <- paste0("abide_2_forCentile_", resDate, ".csv")
write.csv(abide_forCentile, file.path(dataDir, name), row.names = F)