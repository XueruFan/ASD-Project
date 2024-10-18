# 整合人口信息和MRI指标测量数据，并且按照LBCC要求的格式进行整理，便于之后进行Combat分析，以及紧接着
# 的个体偏离百分位数的计算
# Xue-Ru Fan 11 May 2023 @BNU
###################################################
# Part 1: 整合abide1和2的人口信息、行为数据和mri指标测量，csv
# Part 2: 整理成做combat需要的格式，csv
###################################################
rm(list=ls())
packages <- c("plyr", "dplyr", "readxl")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
cabicDir <- 'E:/PhDproject/CABIC'
dataDir <- file.path(cabicDir, "Freesurfer", "ALL")
resuDir <- file.path(cabicDir, "result")
resDate <- "240928"

############################ Part 1: Sum All #######################################################
# 整合人口信息、行为数据和mri指标测量

global <- read.csv(file.path(resuDir, paste0("cabic_global_", resDate, ".csv")))
regional <- read.csv(file.path(resuDir, paste0("cabic_regional_", resDate, ".csv")))
cabic <- na.omit(merge(global, regional, by = "Participant"))

pheno <- read_excel(file.path(cabicDir, "CABIC_Subjects_info.xls"))
colnames(pheno)[3] <- "Participant"
cabic_all <- merge(pheno, cabic, by = "Participant")

# 有一个mri但没有pheno的被试
par_info <- pheno$Participant
par_mri <- cabic$Participant
par_mis <- setdiff(par_mri, par_info)
par_mis # 是28793号，这个被试属于GU站点，查找Phenotypic File（ABIDEII-GU_1.csv）发现，这个被试的行为
# 数据录错位了，估计是在汇总的时候，项目组给删掉了，不过这个被试是控制组的，删了就删了吧，影响不大

write.csv(cabic_all, file.path(resuDir, paste0("cabic_", resDate, ".csv")), row.names = F)


############################ Part 2: Arrange data for ComBat #######################################
# 整理成做combat需要的格式
################## arrange part A: Basic
cabic_analysis <- cabic_all[1:2]
colnames(cabic_analysis) <- c("participant", "Site")
cabic_analysis$Age <- cabic_all$AGE
cabic_analysis$age_days <- as.numeric(cabic_analysis$Age * 365.245)
cabic_analysis$age_months <- as.numeric(round(cabic_analysis$age_days / 30.4375))
cabic_analysis$age_days <- cabic_analysis$age_days + 280
cabic_analysis$sex <- cabic_all$SEX
cabic_analysis[which(cabic_analysis$sex == "F"), 'sex'] = "Female"
cabic_analysis[which(cabic_analysis$sex == "M"), 'sex'] = "Male"
cabic_analysis$study <- "XRF"
cabic_analysis$fs_version <- "FS6_T1"
cabic_analysis$country <- "China"
cabic_analysis$run <- 1
cabic_analysis$session <- 1
cabic_analysis$dx <- cabic_all$GROUP
cabic_analysis[which(cabic_analysis$dx == "TDC"), 'dx'] = "CN"

################## arrange part B: global and regional
cabic_analysis$GMV <- cabic_all$Total.cortical.gray.matter.volume
cabic_analysis$WMV <- cabic_all$Total.cerebral.white.matter.volume
cabic_analysis$sGMV <- cabic_all$Subcortical.gray.matter.volume
cabic_analysis$TCV <- cabic_analysis$GMV + cabic_analysis$WMV # 因为常模里是这么算的
cabic_analysis$TCV_real <- cabic_all$Estimated.Total.Intracranial.Volume
cabic_analysis$Ventricles <- cabic_all$Volume.of.ventricles.and.choroid.plexus
cabic_analysis$meanCT2 <- (as.numeric(cabic_all$lhMeanThickness) * as.numeric(cabic_all$lhNumVert) +
                             as.numeric(cabic_all$rhMeanThickness) * as.numeric(cabic_all$rhNumVert)
) / (as.numeric(cabic_all$lhNumVert) + as.numeric(cabic_all$rhNumVert))
cabic_analysis$totalSA2 <- as.numeric(cabic_all$lhWhiteSurfArea) + as.numeric(cabic_all$
                                                                               rhWhiteSurfArea)

################## arrange part C: 34 parcels
# 随便拷来的一个结果文件，用来提取分区名字
name_dk <- data.frame(read.table(file.path(cabicDir, "toGetName.stats")))
for (p in 1:nrow(name_dk)) {
  com_str <- paste0("cabic_analysis$", name_dk[p, 1], " <- (as.numeric(cabic_all$lh", name_dk[p, 1],
                    ") + as.numeric(cabic_all$rh", name_dk[p, 1], "))/2")
  eval(parse(text = com_str))
}

# save result
name <- paste0("cabic_forComBat_", resDate, ".csv")
write.csv(cabic_analysis, file.path(resuDir, name), row.names = F)
