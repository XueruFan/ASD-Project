# Do ComBat analysis for CABIC data
# Xue-Ru Fan 16 May 2023 @BNU
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

# define filefolder
cabicDir <- 'E:/PhDproject/CABIC'
resuDir <- file.path(cabicDir, "result")
resDate <- "240928"

cabic <- read.csv(file.path(resuDir, paste0("cabic_forComBat_", resDate, ".csv")))

dat <- as.matrix(t(cabic[, which(names(cabic) == "GMV"):ncol(cabic)]))
batch <- c(t(cabic$Site)) # Numeric or character vector
age <- c(t(cabic$age_days)) # Continuous variable
disease <- as.factor(c(t(cabic$dx))) # Categorical variable
sex <- as.factor(c(t(cabic$sex)))
mod <- model.matrix(~ age * sex + disease) # create model matrix for these two biological covariates

data.harmonized <- neuroCombat(dat = dat, batch = batch, mod = mod, parametric = FALSE,
                              eb = FALSE, mean.only = TRUE) # do combat

# get harmonized data
data_combated <- t(data.harmonized$dat.combat)
data_combated <- as.data.frame(data_combated)
# arrange for centile calculation
cabic_forCentile <- cbind(cabic[, 1:which(names(cabic) == "dx")], data_combated)

# save result
name <-  paste0("cabic_forCentile_", resDate, ".csv")
write.csv(cabic_forCentile, file.path(resuDir, name), row.names = F)
