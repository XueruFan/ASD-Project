# this script is used to extracted 7 global measures and 34 regional gray matter volumes
# from the preprocessing CABIC data. 
# Xue-Ru Fan 15 May 2023 @BNU
###################################################
# Part 1: Extract 7 global measurements
# Part 2：Extract 34 regional volume, and calaulate CT，SA，Vertex
###################################################

rm(list=ls())
packages <- c("openxlsx")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

## define filefolder
cabicDir <- 'E:/PhDproject/CABIC'
dataDir <- file.path(cabicDir, "Freesurfer", "ALL")
resuDir <- file.path(cabicDir, "result")
resDate <- "240928"

SUB_ID <- list.files(path = dataDir)


############################ Part 1: Extract 7 global measurements ###################
# 从freesurfer结果文件中提取7个全局指标的结果

## make an empty dataframe to save global data
cabic_global <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(cabic_global) <- c("Participant",
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
for (s in 1:length(SUB_ID)) { 
  cabic_global[s, 1] <- SUB_ID[s] 
  folder_path <- file.path(dataDir, SUB_ID[s], "stats")
  if (dir.exists(folder_path)) {
    setwd(folder_path)
    file_path <- "./aseg.stats"
    if (file.exists(file_path)) {
      ## global
      note <- data.frame(readLines("aseg.stats"))
      global <- data.frame(do.call('rbind', strsplit(as.character(note[1:40,]), ", ")))
      note <- data.frame(read.table("aseg.stats"))
      
      names <- names(cabic_global)[2:10]
      for (n in names) {
        where <- which(apply(global, 1, function(x) any(grep(paste0("\\b", n, "\\b"), x, value = F))))
        com_str <- paste0("cabic_global[s, '", n, "'] <- global[where, 'X4']")
        eval(parse(text = com_str))
      }
    } else {
      next
    }
    
  } else {
    next
  }
}

# save result
name <- paste0("cabic_global_", resDate, ".csv")
write.csv(cabic_global, file.path(resuDir, name), row.names = F)


############################ Part 2: Extract 34 regional volume, and calaulate CT，SA，Vertex ######
# 从freesurfer结果文件中提取34个DK分区体积的结果，并计算出CT，SA，Vertex 

cabic_regional <- data.frame()

## extra data one by one
for (s in 1:length(SUB_ID)) { 
  cabic_regional[s, 1] <- SUB_ID[s] 
  folder_path <- file.path(dataDir, SUB_ID[s], "stats")
  if (dir.exists(folder_path)) {
    setwd(folder_path)
    file_path <- "./lh.aparc.stats"
    if (file.exists(file_path)) {
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
        com_str <- paste0("cabic_regional[s, 'lh", n, "'] <- regional_lh[where, 'X4']")
        eval(parse(text = com_str))
        ######## rh
        where <- which(apply(regional_rh, 1, function(x) any(grep(n, x, value = F))))
        com_str <- paste0("cabic_regional[s, 'rh", n, "'] <- regional_rh[where, 'X4']")
        eval(parse(text = com_str))
      }
      
      for (d in lh_dk$V1) {
        ######## lh
        where <- which(apply(lh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
        if (length(where) != 0) {
          com_str <- paste0("cabic_regional[s, 'lh", d, "'] <- lh_dk[where, 'V4']")
          eval(parse(text = com_str))
        } 
        ######## rh
        where <- which(apply(rh_dk, 1, function(x) any(grep(paste0("\\b", d, "\\b"), x, value = F))))
        if (length(where) != 0) {
          com_str <- paste0("cabic_regional[s, 'rh", d, "'] <- rh_dk[where, 'V4']")
          eval(parse(text = com_str))
        }
      }
      
    } else {
      next
    }
  } else {
    next
  }
}
  

colnames(cabic_regional)[1] <- "Participant"

# save result
name <- paste0("cabic_regional_", resDate, ".csv")
write.csv(cabic_regional, file.path(resuDir, name), row.names = F)
