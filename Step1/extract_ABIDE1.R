# Extracte 7 global measures and 34 regional gray matter volumes
# from the preprocessing ABIDE I data. 
# please note, this works for freesurfer6 results only
# Xue-Ru Fan 2 Jan 2024 @BNU
###################################################
# Part 0: Untar Files
# Part 1: Get sub list after qc
# Part 2: Extract 7 global measurements
# Part 3：Extract 34 regional volume, and calculate CT, SA, Vertex
###################################################

rm(list=ls())
# install.packages("R.utils")
# library(R.utils)

## define filefolder
abideDir <- 'E:/PhDproject/ABIDE'


############################ Part 0: Untar Files ###################################################

dataDir <- file.path(abideDir, "Preprocessed/ABIDE_I")
tarDir <- file.path(dataDir, "SiteWise_tar")
untarDir <- file.path(dataDir, "SubWise")

setwd(tarDir)

## get site_id
SITE_list <- list.files(tarDir)

for (site in SITE_list) {
  tarDir_site <- file.path(tarDir, site)
  setwd(tarDir_site)

  ## get sub_id
  SUB_file <- list.files(pattern = "*.tar.gz")

  ## untar files
  for (s in SUB_file) {
    untar(s, exdir = untarDir)
  }
}


############################ Part 1: Get sub list after qc #################################################

library(openxlsx)

############## load in list of subjects after manual qc
abide_qc <- read.xlsx(file.path(abideDir, "ABIDE_T1qc/abide_qc_result_afterVisualQc.xlsx"))

# 筛选出需要再看分割情况的数据
sub_qc <- subset(abide_qc, Rate == "2" | Rate == "1")

sub_qc_1 <- subset(sub_qc, ABIDE == "1") # 提取abide1的结果

sub_qc_1 <- unique(sub_qc_1$SUB_ID) # 注意：有些被试有重复个run，这里只保留一次被试编码


## get sub_id
dataDir <- file.path(abideDir, "Preprocessed/ABIDE_I/SubWise")

SUB_file_all <- list.files(dataDir) # 有预处理结果的所有被试
SUB_file <- intersect(sub_qc_1, SUB_file_all) # 取通过了qc的被试编号和有预处理结果的被试编号的交集

# 保存一份用于后续研究分析的被试清单
write.table(SUB_file, file.path(abideDir, "Preprocessed/abide_1_sub4analysis.csv"), row.names = F,
            col.names = F)


############################ Part 2: Extract 7 global measurements #################################
# 从freesurfer结果文件中提取全局指标的结果

## make an empty dataframe to save global data
abide1_global <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(abide1_global) <- c("Participant",
                             "Left hemisphere cortical gray matter volume",
                             "Right hemisphere cortical gray matter volume",
                             "Total cortical gray matter volume",
                             "Left hemisphere cerebral white matter volume",
                             "Right hemisphere cerebral white matter volume",
                             "Total cerebral white matter volume",
                             "Subcortical gray matter volume",
                             "Volume of ventricles and choroid plexus",
                             "Estimated Total Intracranial Volume")


## extra data one by one
for (s in 1:length(SUB_file)) { 
  setwd(file.path(dataDir, SUB_file[s], "stats"))
  ## sub_id
  abide1_global[s, 1] <- SUB_file[s] 
  ## global
  note <- data.frame(readLines("aseg.stats"))
  global <- data.frame(do.call('rbind', strsplit(as.character(note[14:35,]), ", ")))
  note <- data.frame(read.table("aseg.stats"))
  
  names <- names(abide1_global)[2:10]
  for (n in names) {
    where <- which(apply(global, 1, function(x) any(grep(paste0("\\b", n, "\\b"), x, value = F))))
    com_str <- paste0("abide1_global[s, '", n, "'] <- global[where, 'X4']")
    eval(parse(text = com_str))
  }
}

# save result
write.csv(abide1_global, file.path(abideDir, "Preprocessed/abide_1_global_afterqc.csv"),
          row.names = F)


############################ Part 3: Extract 34 regional volume, and calculate CT，SA，Vertex ######
# 从freesurfer结果文件中提取分区体积的结果

abide1_regional <- data.frame()

## extra data one by one
for (s in 1:length(SUB_file)) {
  setwd(file.path(dataDir, SUB_file[s], "stats"))
  
  # sub_id
  abide1_regional[s, 1] <- SUB_file[s] 
  
  ###################################  lh
  lh <- data.frame(readLines("lh.aparc.stats"))
  lh_dk <- data.frame(read.table("lh.aparc.stats"))
  # Vertex
  whereV <- which(apply(lh, 1, function(x) any(grep("NumVert, Number of Vertices", x, value = F))))
  # SA
  whereW <- which(apply(lh, 1, function(x) any(grep("WhiteSurfArea, White Surface Total Area", x,
                                                    value = F))))
  # CT
  whereM <- which(apply(lh, 1, function(x) any(grep("MeanThickness, Mean Thickness", x,
                                                    value = F))))
  
  whereL <- c(whereV, whereW, whereM)
  regional_lh <- data.frame(do.call('rbind', strsplit(as.character(lh[whereL, ]), ", ")))
  
  ###################################  rh
  rh <- data.frame(readLines("rh.aparc.stats"))
  rh_dk <- data.frame(read.table("rh.aparc.stats"))
  whereV <- which(apply(rh, 1, function(x) any(grep("NumVert, Number of Vertices", x, value = F))))
  whereW <- which(apply(rh, 1, function(x) any(grep("WhiteSurfArea, White Surface Total Area", x,
                                                    value = F))))
  whereM <- which(apply(rh, 1, function(x) any(grep("MeanThickness, Mean Thickness", x,
                                                    value = F))))
  whereR <- c(whereV, whereW, whereM)
  regional_rh <- data.frame(do.call('rbind', strsplit(as.character(rh[whereR, ]), ", ")))
  
  
  names <-  c("NumVert", "WhiteSurfArea", "MeanThickness")
  
  for (n in names) {
    ######## lh
    where <- which(apply(regional_lh, 1, function(x) any(grep(n, x, value = F))))
    com_str <- paste0("abide1_regional[s, 'lh", n, "'] <- regional_lh[where, 'X4']")
    eval(parse(text = com_str))
    ######## rh
    where <- which(apply(regional_rh, 1, function(x) any(grep(n, x, value = F))))
    com_str <- paste0("abide1_regional[s, 'rh", n, "'] <- regional_rh[where, 'X4']")
    eval(parse(text = com_str))
  }
  
  for (d in lh_dk$V1) {
    ######## lh
    where <- which(apply(lh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
    if (length(where) != 0) {
      com_str <- paste0("abide1_regional[s, 'lh", d, "'] <- lh_dk[where, 'V4']")
      eval(parse(text = com_str))
    } 
    ######## rh
    where <- which(apply(rh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
    if (length(where) != 0) {
      com_str <- paste0("abide1_regional[s, 'rh", d, "'] <- rh_dk[where, 'V4']")
      eval(parse(text = com_str))
    }
  }
}

colnames(abide1_regional)[1] <- "Participant"

# save result
write.csv(abide1_regional, file.path(abideDir, "Preprocessed/abide_1_regional_afterqc.csv"),
          row.names = F)
