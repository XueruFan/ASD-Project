# 本代码继上一步410和411后，汇总整理两组人各认知行为统计结果
# 雪如 2024年7月31日于北师大办公室
################################

################################
rm(list=ls())
packages <- c("ggplot2", "ggridges", "tidyr", "bestNormalize", "dplyr", "reshape2")
# sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

abideDir <- 'E:/PhDproject/ABIDE'
statiDir <- file.path(abideDir, "Analysis/Statistic")
newDate <- "240610"

# load all the differ files
VIN_trans <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_VIN_trans_", newDate, ".csv")))
VIN_scale <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_VIN_scale_", newDate, ".csv")))
type <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_type_", newDate, ".csv")))
SRS_total <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_SRS_total_", newDate, ".csv")))
SRS_scale_trans <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_SRS_scale_trans_", newDate, ".csv")))
SRS_scale_raw <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_SRS_scale_raw_", newDate, ".csv")))
site <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_site_", newDate, ".csv")))
scaner <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_scaner_", newDate, ".csv")))
manu <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_manu_", newDate, ".csv")))
iq <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_iq_", newDate, ".csv")))
bmi <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_bmi_", newDate, ".csv")))
age <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_age_", newDate, ".csv")))
ADOS_G <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_ADOS_G_", newDate, ".csv")))
ADOS_2 <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_ADOS_2_", newDate, ".csv")))
ADIR <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_ADIR_", newDate, ".csv")))
CBCL <- read.csv(file.path(statiDir, paste0("asd_male_dev_SC_statis_CBCL_", newDate, ".csv")))

# merge
data_frames <- list(age, bmi, iq, ADOS_2, ADOS_G, ADIR, SRS_scale_raw, SRS_scale_trans,
                    SRS_total, VIN_scale, VIN_trans, CBCL)

data <- do.call(rbind, data_frames)

# arrange
data[, 6] <- paste0(data$Mean, "±", data$SD)
colnames(data)[6] <- "Mean±SD"
data <- data[, c(1,5,2,6)]
colnames(data)[1:2] <- c("Variable", "N")
data[1,1] <- "H Age"
data[2,1] <- "L Age"
data[3,1] <- "H BMI"
data[4,1] <- "L BMI"

name <- paste0("asd_male_dev_SC_statis_Summary_", newDate, ".xlsx")
write.xlsx(data, file.path(statiDir, name))
