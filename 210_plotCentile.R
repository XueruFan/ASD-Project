# this script is used to visualize ABIDE centile (6-18 only)
# Xue-Ru Fan 25 April 2023 @BNU
###################################################
# Part 1: ASD男性34个分区体积centile（中位数+四分位距）的概率密度图，png
# Part 2: ASD男性34个分区体积centile的中位数投射出的脑图，png
# Part 3: ASD男性7个全局指标常模分（中位数+四分位距）的概率密度图，png
###################################################

rm(list=ls())
packages <- c("ggplot2", "ggseg", "ggridges", "tidyr")
# sapply(packages, install.packages, character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# define filefolder
abideDir <- 'E:/PhDproject/ABIDE'
dataDir <- file.path(abideDir, "Preprocessed")
plotDir <- file.path(abideDir, "Plot/Density")
date <- "240225"

abide_centile <- read.csv(file.path(dataDir, paste0("abide_centile_dev_", date, ".csv")))
asd_centile <- subset(abide_centile, dx == "ASD" & sex == "Male")

######################### Part 1: Plot density for 34 regions ######################################
# ASD 34个分区体积Centile分（中位数+四分位距）的概率密度图

asd_parc <- asd_centile[, which(names(asd_centile) == "bankssts"):ncol(asd_centile)]

# transfer to long data
asd_long <- gather(asd_parc, key = "Measure", value = "Centile", bankssts:insula, na.rm = TRUE,
                   factor_key = TRUE)

# modify data class
centile <- as.numeric(asd_long$Centile)
measure <- asd_long$Measure
measure2 <- reorder(measure, centile, FUN = median) # from small to large

asd_plot <- ggplot(asd_long, aes(x = centile, y = measure2, fill = measure2, alpha = .5)) +
  geom_density_ridges(scale = 1.5, quantile_lines = TRUE, size = .3, quantiles = 4) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = rep("transparent", 34)) +
  coord_fixed(ratio = 0.17) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(size = 6, face = "bold"))

asd_plot

# save plot
name <- file.path(plotDir, paste0("asd_male_centile_regional_density_", date, ".png"))
ggsave(name, width = 4.6, height = 6, units = "in", dpi = 500)


######################### Part 2: Project 34 region centile on the brain ###########################
# ASD患者34个分区体积常模分的中位数的脑图

# calculate median
asd_parc_centile <- data.frame(paste0("lh_", names(asd_parc)))
colnames(asd_parc_centile) = "label"
asd_parc_centile$median <- apply(asd_parc, 2, median, na.rm = TRUE)

ggseg(.data = asd_parc_centile, mapping = aes(fill = median), color = "black", atlas = dk,
      position = "stacked", hemisphere = "left", size = 1.2) +
  theme_void() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1, "cm")) +
  scale_fill_gradient2(low = "#126cb5", mid = "white", high = "#a52a2a", midpoint = 0.5) +
  guides(fill = guide_colourbar(frame.colour = "black", frame.linewidth = 1, ticks = FALSE))

name <- file.path(plotDir, paste0("asd_male_centile_regional_median_", date, ".png"))
ggsave(name, width = 7.8, height = 3, units = "in", dpi = 500)


######################### Part3: Plot density for 7 measures #######################################
# ASD组7个全局指标常模分（中位数+四分位距）的概率密度图

asd_global <- asd_centile[which(names(asd_centile) == "GMV"):which(names(asd_centile) == "totalSA2")]

# transfer to long data
asd_long <- gather(asd_global, key = "Measure", value = "Centile", GMV:totalSA2, na.rm = TRUE,
                   factor_key = TRUE)

# modify data class
centile <- as.numeric(asd_long$Centile)
measure <- asd_long$Measure
measure2 <- reorder(measure, centile, FUN = median) # from small to large

# modify name
levels(measure2)[levels(measure2) == "GMV"] <- "Grey Matter Volume"
levels(measure2)[levels(measure2) == "WMV"] <- "White Matter Volume"
levels(measure2)[levels(measure2) == "sGMV"] <- "Subcortical Volume"
levels(measure2)[levels(measure2) == "TCV"] <- "Total Cerebrum Volume"
levels(measure2)[levels(measure2) == "meanCT2"] <- "Mean Cortical Thickness"
levels(measure2)[levels(measure2) == "Ventricles"] <- "Ventricular Volume"
levels(measure2)[levels(measure2) == "totalSA2"] <- "Total Surface Area"

# plot 
asd_plot <- ggplot(asd_long, aes(x = centile, y = measure2, fill = measure2)) +
  geom_density_ridges(scale = 1.8, quantile_lines = TRUE, size = 0.9, quantiles = 4) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = c("Grey Matter Volume" = "#4682b4",
                               "White Matter Volume" = "#dcdcdc",
                               "Subcortical Volume" = "#808080",
                               "Total Cerebrum Volume" = "#f0e68c",
                               "Mean Cortical Thickness" = "#008b8b",
                               "Total Surface Area" = "#ffa07a",
                               "Ventricular Volume" = "#ffb6c1")) +
  coord_fixed(ratio = 0.2) + 
  xlab("") +
  ylab("") +
  theme_ridges() + 
  theme(legend.position = "none", # without legend
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"))

asd_plot

# save plot
name <- file.path(plotDir, paste0("asd_male_centile_global_density_", date, ".png"))
ggsave(name, width = 6, height = 5, units = "in", dpi = 500)
