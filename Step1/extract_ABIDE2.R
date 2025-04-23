# extract 7 global measures and 34 regional gray matter volumes
# from the preprocessing ABIDE II data. 
# please note, this works for freesurfer6 results only
# Xue-Ru Fan 15 May 2023 @BNU
###################################################
# Part 0: Untar Files
# Part 1: Get sub list after qc
# Part 2: Extract 7 global measurements
# Part 3：Extract 34 regional volume, and calculate CT, SA, Vertex
###################################################

rm(list=ls())
packages <- c("openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

## define filefolder
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed/ABIDE_II")
resDate <- "240315"

############################ Part 0: Untar Files ###################################################
# 解压freesurfer预处理完的tar.gz文件
tarDir <- file.path(dataDir, "SiteWise")
untarDir <- file.path(dataDir, "SubWise")

## get sub_id
SUB_file <- list.files(tarDir)

## untar files
setwd(tarDir)
for (s in SUB_file) {
  untar(s, exdir = untarDir)
}

############################ Part 1: Get sub list after qc #########################################
# 导入做完人工质控后的被试清单，与有预处理结果的被试清单取交集，整合出最终用于后续分析的被试清单
abide_qc <- read.xlsx(file.path(abideDir, "ABIDE_T1qc/abide_qc_result_afterVisualQc.xlsx"))

# 筛选出非0数据
sub_qc <- subset(abide_qc, Rate == "2" | Rate == "1")

sub_qc_2 <- subset(sub_qc, ABIDE == "2") # 提取abide2的结果

sub_qc_2 <- unique(sub_qc_2$SUB_ID) # 注意：有些被试有重复个run，这里只保留一次被试编码


## get sub_id
procDir <- file.path(dataDir, "SubWise")
SUB_file_all <- list.files(procDir) # 有预处理结果的所有被试
SUB_file <- intersect(sub_qc_2, SUB_file_all) # 取通过了qc的被试编号和有预处理结果的被试编号的交集

# 保存一份用于后续研究分析的被试清单
name <- paste0("abide_2_sub4analysis_", resDate, ".csv")
write.csv(SUB_file, file.path(abideDir, "Preprocessed", name), row.names = F, col.names = F)


############################ Part 2: Extract 7 global measurements #################################
# 从freesurfer结果文件中提取7个全局指标的结果

## make an empty dataframe to save global data
abide2_global <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(abide2_global) <- c("Participant",
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
  setwd(file.path(dataDir, "SubWise", SUB_file[s], "stats"))
  ## sub_id
  abide2_global[s, 1] <- SUB_file[s] 
  ## global
  note <- data.frame(readLines("aseg.stats"))
  global <- data.frame(do.call('rbind', strsplit(as.character(note[14:35,]), ", ")))
  note <- data.frame(read.table("aseg.stats"))
  
  names <- names(abide2_global)[2:10]
  for (n in names) {
    where <- which(apply(global, 1, function(x) any(grep(paste0("\\b", n, "\\b"), x, value = F))))
    com_str <- paste0("abide2_global[s, '", n, "'] <- global[where, 'X4']")
    eval(parse(text = com_str))
  }
}

# save result
name <- paste0("abide_2_global_afterqc_", resDate, ".csv")
write.csv(abide2_global, file.path(abideDir, "Preprocessed", name), row.names = F)


############################ Part 3: Extract 34 regional volume, and calculate CT，SA，Vertex ######
# 从freesurfer结果文件中提取34个DK分区体积的结果，并计算出CT，SA，Vertex 

abide2_regional <- data.frame()

## extra data one by one
for (s in 1:length(SUB_file)) {
  setwd(file.path(dataDir, "SubWise", SUB_file[s], "stats"))
  
  # sub_id
  abide2_regional[s, 1] <- SUB_file[s] 
  
  ###################################  lh
  lh <- data.frame(readLines("lh.aparc.stats"))
  lh_dk <- data.frame(read.table("lh.aparc.stats"))
  # Vertex
  whereV <- which(apply(lh, 1, function(x) any(grep("NumVert, Number of Vertices", x, value = F))))
  # SA
  whereW <- which(apply(lh, 1, function(x) any(grep("WhiteSurfArea, White Surface Total Area", x,
                                                    value = F))))
  # CT
  whereM <- which(apply(lh, 1, function(x) any(grep("MeanThickness, Mean Thickness", x, value = F))))
  
  whereL <- c(whereV, whereW, whereM)
  regional_lh <- data.frame(do.call('rbind', strsplit(as.character(lh[whereL, ]), ", ")))
  
  ###################################  rh
  rh <- data.frame(readLines("rh.aparc.stats"))
  rh_dk <- data.frame(read.table("rh.aparc.stats"))
  whereV <- which(apply(rh, 1, function(x) any(grep("NumVert, Number of Vertices", x, value = F))))
  whereW <- which(apply(rh, 1, function(x) any(grep("WhiteSurfArea, White Surface Total Area", x,
                                                    value = F))))
  whereM <- which(apply(rh, 1, function(x) any(grep("MeanThickness, Mean Thickness", x, value = F))))
  whereR <- c(whereV, whereW, whereM)
  regional_rh <- data.frame(do.call('rbind', strsplit(as.character(rh[whereR, ]), ", ")))
  
  
  names <-  c("NumVert", "WhiteSurfArea", "MeanThickness")
  
  for (n in names) {
    ######## lh
    where <- which(apply(regional_lh, 1, function(x) any(grep(n, x, value = F))))
    com_str <- paste0("abide2_regional[s, 'lh", n, "'] <- regional_lh[where, 'X4']")
    eval(parse(text = com_str))
    ######## rh
    where <- which(apply(regional_rh, 1, function(x) any(grep(n, x, value = F))))
    com_str <- paste0("abide2_regional[s, 'rh", n, "'] <- regional_rh[where, 'X4']")
    eval(parse(text = com_str))
  }
  
  for (d in lh_dk$V1) {
    ######## lh
    where <- which(apply(lh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
    if (length(where) != 0) {
      com_str <- paste0("abide2_regional[s, 'lh", d, "'] <- lh_dk[where, 'V4']")
      eval(parse(text = com_str))
    } 
    ######## rh
    where <- which(apply(rh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
    if (length(where) != 0) {
      com_str <- paste0("abide2_regional[s, 'rh", d, "'] <- rh_dk[where, 'V4']")
      eval(parse(text = com_str))
    }
  }
}

colnames(abide2_regional)[1] <- "Participant"

# save result
name <- paste0("abide_2_regional_afterqc_", resDate, ".csv")
write.csv(abide2_regional, file.path(abideDir, "Preprocessed", name), row.names = F)