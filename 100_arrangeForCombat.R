# 整合人口信息和MRI指标测量数据，并且按照LBCC要求的格式进行整理，便于之后进行Combat分析，以及紧接着
# 的个体偏离百分位数的计算
# Xue-Ru Fan 11 May 2023 @BNU
###################################################
# Part 1: 整合abide1和2的人口信息、行为数据和mri指标测量，csv
# Part 2: 整理成做combat需要的格式，csv
###################################################
rm(list=ls())
packages <- c("plyr", "dplyr")
#sapply(packages,install.packages,character.only=TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
resDate <- "240315"

############################ Part 1: Sum All #######################################################
# 整合abide1和2的人口信息、行为数据和mri指标测量

################## abide 1
abide1_global <- read.csv(file.path(dataDir, "abide_1_global_afterqc.csv"))
abide1_regional <- read.csv(file.path(dataDir, "abide_1_regional_afterqc.csv"))
abide1 <- merge(abide1_global, abide1_regional, by = "Participant")
# load pheno
abide1_pheno <- read.csv(file.path(abideDir, "ABIDE_info", "Phenotypic_V1_0b.csv"))
colnames(abide1_pheno)[2] <- "Participant"
# combine
abide1_all <- merge(abide1_pheno, abide1, by = "Participant")
# save a copy
write.csv(abide1_all, file.path(dataDir, paste0("abide_1_all_", resDate, ".csv")), row.names = F)

################## abide 2
abide2_global <- read.csv(file.path(dataDir, paste0("abide_2_global_afterqc_", resDate, ".csv")))
abide2_regional <- read.csv(file.path(dataDir, paste0("abide_2_regional_afterqc_", resDate, ".csv")))
abide2 <- merge(abide2_global, abide2_regional, by = "Participant")
abide2_pheno <- read.csv(file.path(abideDir, "ABIDE_info", "ABIDEII_Composite_Phenotypic.csv"))
colnames(abide2_pheno)[2] <- "Participant"
abide2_all <- merge(abide2_pheno, abide2, by = "Participant")

# 这里发现有一个被试有mri数据（abide2）但没有pheno数据（abide2_pheno）,找到这个被试
par_info <- abide2_pheno$Participant
par_mri <- abide2$Participant
par_mis <- setdiff(par_mri, par_info)
par_mis # 是28793号，这个被试属于GU站点，查找Phenotypic File（ABIDEII-GU_1.csv）发现，这个被试的行为
# 数据录错位了，估计是在汇总的时候，项目组给删掉了，不过这个被试是控制组的，删了就删了吧，影响不大

write.csv(abide2_all, file.path(dataDir, paste0("abide_2_all_", resDate, ".csv")), row.names = F)


################## combine abide 1 and 2
abide1_all$Release <- "1"
abide2_all$Release <- "2"

### abide1中有些列名和2不一致，这里统一按照abide2的命名修改
abide1_all <- abide1_all %>% rename(PDD_DSM_IV_TR = DSM_IV_TR)
abide1_all <- abide1_all %>% rename(ADOS_G_TOTAL = ADOS_TOTAL)
abide1_all <- abide1_all %>% rename(ADOS_G_COMM = ADOS_COMM)
abide1_all <- abide1_all %>% rename(ADOS_G_SOCIAL = ADOS_SOCIAL)
abide1_all <- abide1_all %>% rename(ADOS_G_STEREO_BEHAV = ADOS_STEREO_BEHAV)
abide1_all <- abide1_all %>% rename(ADOS_2_SOCAFFECT = ADOS_GOTHAM_SOCAFFECT)
abide1_all <- abide1_all %>% rename(ADOS_2_RRB = ADOS_GOTHAM_RRB)
abide1_all <- abide1_all %>% rename(ADOS_2_TOTAL = ADOS_GOTHAM_TOTAL)
abide1_all <- abide1_all %>% rename(ADOS_2_SEVERITY_TOTAL = ADOS_GOTHAM_SEVERITY)
abide1_all <- abide1_all %>% rename(SRS_TOTAL_RAW = SRS_RAW_TOTAL)
abide1_all <- abide1_all %>% rename(SRS_AWARENESS_RAW = SRS_AWARENESS)
abide1_all <- abide1_all %>% rename(SRS_COGNITION_RAW = SRS_COGNITION)
abide1_all <- abide1_all %>% rename(SRS_COMMUNICATION_RAW = SRS_COMMUNICATION)
abide1_all <- abide1_all %>% rename(SRS_MOTIVATION_RAW = SRS_MOTIVATION)
abide1_all <- abide1_all %>% rename(SRS_MANNERISMS_RAW = SRS_MANNERISMS)
abide1_all <- abide1_all %>% rename(NONASD_PSYDX_LABEL = COMORBIDITY)
abide1_all <- abide1_all %>% rename(CURRENT_MEDICATION_NAME = MEDICATION_NAME)
abide1_all <- abide1_all %>% rename(ADI_R_RRB_TOTAL_C = ADI_RRB_TOTAL_C)
abide1_all <- abide1_all %>% rename(VINELAND_DAILYLIVING_STANDARD = VINELAND_DAILYLVNG_STANDARD)
abide1_all <- abide1_all %>% rename(VINELAND_ABC_Standard = VINELAND_ABC_STANDARD)

abide_all <- rbind.fill(abide2_all, abide1_all)

write.csv(abide_all, file.path(dataDir, paste0("abide_A_all_", resDate, ".csv")), row.names = F)

objects_to_keep <- c("abide_all", "abideDir", "dataDir", "resDate") # 列出要保留的对象名称
rm(list = (setdiff(ls(), objects_to_keep))) # 删除不在保留列表中的所有对象


############################ Part 2: Arrange data for ComBat #######################################
# 整理成做combat需要的格式

################## arrange part A: Basic
abide_analysis <- abide_all[1:2]
colnames(abide_analysis) <- c("participant", "Site")
abide_analysis$Release <- abide_all$Release
abide_analysis$Age <- abide_all$AGE_AT_SCAN
abide_analysis$age_days <- as.numeric(abide_analysis$Age * 365.245)
abide_analysis$age_months <- as.numeric(round(abide_analysis$age_days / 30.4375))
abide_analysis$age_days <- abide_analysis$age_days + 280
abide_analysis$sex <- abide_all$SEX
abide_analysis[which(abide_analysis$sex == 2), 'sex'] = "Female"
abide_analysis[which(abide_analysis$sex == 1), 'sex'] = "Male"
abide_analysis$study <- "XRF"
abide_analysis$fs_version <- "FS6_T1"
abide_analysis$country <- NA
abide_analysis$run <- 1
abide_analysis$session <- 1
abide_analysis$dx <- abide_all$DX_GROUP
abide_analysis[which(abide_analysis$dx == 2), 'dx'] = "CN"
abide_analysis[which(abide_analysis$dx == 1), 'dx'] = "ASD"

################## arrange part B: global and regional
abide_analysis$GMV <- abide_all$Total.cortical.gray.matter.volume
abide_analysis$WMV <- abide_all$Total.cerebral.white.matter.volume
abide_analysis$sGMV <- abide_all$Subcortical.gray.matter.volume
abide_analysis$TCV <- abide_analysis$GMV + abide_analysis$WMV # 因为常模里是这么算的
abide_analysis$TCV_real <- abide_all$Estimated.Total.Intracranial.Volume
abide_analysis$Ventricles <- abide_all$Volume.of.ventricles.and.choroid.plexus
abide_analysis$meanCT2 <- (as.numeric(abide_all$lhMeanThickness) * as.numeric(abide_all$lhNumVert) +
                             as.numeric(abide_all$rhMeanThickness) * as.numeric(abide_all$rhNumVert)
) / (as.numeric(abide_all$lhNumVert) + as.numeric(abide_all$rhNumVert))
abide_analysis$totalSA2 <- as.numeric(abide_all$lhWhiteSurfArea) + as.numeric(abide_all$
                                                                               rhWhiteSurfArea)

################## arrange part C: 34 parcels
# 随便拷来的一个结果文件，用来提取分区名字
name_dk <- data.frame(read.table(file.path(dataDir, "toGetName.stats")))
for (p in 1:nrow(name_dk)) {
  com_str <- paste0("abide_analysis$", name_dk[p, 1], " <- as.numeric(abide_all$lh", name_dk[p, 1],
                    ") + as.numeric(abide_all$rh", name_dk[p, 1], ")")
  eval(parse(text = com_str))
}

# save result
name <- paste0("abide_A_forComBat_", resDate, ".csv")
write.csv(abide_analysis, file.path(dataDir, name), row.names = F)