library(tidyverse)
library(EMCluster)
# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#==============================================================#
# Part 1: Process Sorter data
#==============================================================#
# Define a vector of your experiment directories. Can load more
dirs <- c("data/raw/20220420_dauerProtocol2")

# Read in the data using easysorter: https://github.com/AndersenLab/easysorter
raw <- easysorter::read_data(dirs)

# look at score data for each color channel
score_proc <- as.data.frame(raw[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE", !is.na(strain)) %>%
  dplyr::mutate(well = paste(plate,row,col, sep = "")) #%>%
  #tidyr::pivot_longer(cols = c(green, yellow, red))

# # plot Fl by TOF for all channels
# a1_channels <- ggplot(score_proc %>% dplyr::filter(!is.na(condition) & !is.na(strain))) +
#   aes(x = TOF, y = value, fill = condition) +
#   geom_point(size=1, alpha = 0.3, shape = 21) +
#   #facet_wrap(~condition) +
#   labs(y = "Fluorescence (AU)") +
#   facet_grid(strain~condition~name) +
#   theme_bw() +
#   theme(legend.position = "none")
# a1_channels

#cowplot::ggsave2(a1_channels, filename = "plots/dauerProtocol1_sorter_FLchannels.png")
#========================================================#
# USE Daehan's Data to get dauer classes
#========================================================#
raw_dauersGWA <- easysorter::read_data("data/raw/20170917_dauerGWA3A") # FOR NOW USE THIS ASSAY FOR SETTINGS
  
proc_dauersGWA <- as.data.frame(raw_dauersGWA[1]) %>%
  dplyr::mutate(well = paste(plate,row,col, sep = ""))

# EM clustering in DI condition
df_raw <- proc_dauersGWA %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(strain %in% c("JU346", "ECA1089", "ECA1119", "ECA1121", "NIC166", "ECA1091", "ECA1120", "ECA1122")) %>%
  dplyr::filter(condition == "ascr5.800nM") %>%
  dplyr::filter(TOF>100)

y <- df_raw[c("TOF", "green", "EXT")]
ret.em_K4 <- EMCluster::init.EM(y, nclass = 4, method = "em.EM")
y_merge_K4 <- data.frame(df_raw, ret.em_K4[6])

# plot to define classes
ggplot(y_merge_K4) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "GWA3a 800nM ascr#5") +
  theme_bw() +
  ylim(0, 300) +
  xlim(0, 600)

# Caution!!! factor label changes everytime!!!
# labels=c("crap1", "crap2", "non-dauer", "dauer")
# labels=c("crap1", "dauer",  "crap2", "non-dauer")
# labels=c("non-dauer", "crap1", "crap2", "dauer")
y_merge_K4$class <- factor(y_merge_K4$class, levels=c(3,1,4,2), labels=c("dauer", "non-dauer", "other1", "other2"))

# replot with new factor names
gwa3_plot <- y_merge_K4 %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  ggplot(.) +
  aes(x = TOF, 
      y = green, 
      scale = T, color = class) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(~condition) +
  labs(title = "GWA3a 800nM ascr#5", color = "") +
  theme_bw() +
  ylim(0, 300) +
  xlim(0, 600)
gwa3_plot
cowplot::ggsave2(gwa3_plot, filename = "plots/GWA3a_classes.png", width =5, height = 5)
#========================================================#
# Try to apply those classes to our data
#========================================================#

# applying model from DI condition to all our data
ours <- score_proc[c("TOF","green", "EXT")]
ret.em_K4_ours <- EMCluster::assign.class(ours, ret.em_K4)
y_merge_K4_ours <- data.frame(score_proc, ret.em_K4_ours[6])
# CHANGES EVERY TIME!
y_merge_K4_ours$class <- factor(y_merge_K4_ours$class, levels=c(3,1,4,2), labels=c("dauer", "non-dauer", "other1", "other2"))
# plot to define classes
dauer_assay2 <- ggplot(y_merge_K4_ours) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  #scale_color_manual(values = ) # can use to set colors consistenetly
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "dauer assay 2", color = "") +
  facet_wrap(~condition) +
  theme_bw() +
  ylim(0, 300) +
  xlim(0, 600)
dauer_assay2

cowplot::ggsave2(dauer_assay2, filename = "plots/dauer_assay2_clusters_by_condition.png", width = 5, height = 5)

#====================================================================#
# dauer assay Protocol1
#====================================================================#
# Define a vector of your experiment directories. Can load more
dir1 <- c("data/raw/20220407_dauerProtocol1")

# Read in the data using easysorter: https://github.com/AndersenLab/easysorter
raw1 <- easysorter::read_data(dir1)
# look at score data for each color channel
score_proc1 <- as.data.frame(raw1[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE", !is.na(strain)) %>%
  dplyr::mutate(well = paste(plate,row,col, sep = ""))

# applying model from DI condition to all our data
ours1 <- score_proc1[c("TOF","green", "EXT")]
ret.em_K4_ours1 <- EMCluster::assign.class(ours1, ret.em_K4)
y_merge_K4_ours1 <- data.frame(score_proc1, ret.em_K4_ours1[6])
# CHANGES EVERY TIME!
y_merge_K4_ours1$class <- factor(y_merge_K4_ours1$class, levels=c(3,1,4,2), labels=c("dauer", "non-dauer", "other1", "other2"))

# plot to define classes
dauer_assay1 <- ggplot(y_merge_K4_ours1) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  #scale_color_manual(values = ) # can use to set colors consistenetly
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "dauer assay 1", color = "") +
  facet_grid(~condition) +
  theme_bw() +
  ylim(0, 300) +
  xlim(0, 600)
dauer_assay1

cowplot::ggsave2(dauer_assay1, filename = "plots/dauer_assay1_clusters_by_condition.png", width = 5, height = 5)
#cowplot::ggsave2(dauer_assay1, filename = "plots/dauer_assay1_clusters.png", width = 5, height = 5)
#===============================================================================#
# make box plots
#===============================================================================#
da1_sum <- y_merge_K4_ours1 %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(well_tot_n = n(),
                well_dauer_n = sum(class == "dauer"),
                well_nondauer_n = sum(class == "non-dauer"),
                well_other1_n = sum(class == "other1"),
                dauerFrac1 = well_dauer_n/(well_nondauer_n + well_dauer_n),
                dauerFrac2 = well_dauer_n/(well_other1_n + well_dauer_n))

da2_sum <- y_merge_K4_ours %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(well_tot_n = n(),
                well_dauer_n = sum(class == "dauer"),
                well_nondauer_n = sum(class == "non-dauer"),
                well_other1_n = sum(class == "other1"),
                dauerFrac1 = well_dauer_n/(well_nondauer_n + well_dauer_n),
                dauerFrac2 = well_dauer_n/(well_other1_n + well_dauer_n))

# box plots
bp1 <- ggplot(da1_sum %>% dplyr::distinct(well, .keep_all = T)) +
  aes(x = strain, 
      y = dauerFrac2, fill=strain,
      labs=strain) +
  ylim(0, 1)+
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  theme_bw() + 
  facet_wrap(~condition, ncol = 1) +
  labs(y = "dauer fraction (d/(d+other1))", title = "Assay1") +
  theme(legend.position = "none")
bp1

bp1b <- ggplot(da1_sum %>% dplyr::distinct(well, .keep_all = T)) +
  aes(x = strain, 
      y = dauerFrac1, fill=strain,
      labs=strain) +
  ylim(0, 1)+
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  theme_bw() + 
  facet_wrap(~condition, ncol = 1) +
  labs(y = "dauer fraction (d/(d+nond))", title = "Assay1") +
  theme(legend.position = "none")
bp1b

bp2a <- ggplot(da2_sum %>% dplyr::distinct(well, .keep_all = T)) +
  aes(x = strain, 
      y = dauerFrac1, fill=strain,
      labs=strain) +
  ylim(0, 1)+
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  theme_bw() + 
  facet_wrap(~condition, ncol = 2) +
  labs(y = "dauer fraction (d/(d+nond))", title = "Assay2") +
  theme(legend.position = "none")
bp2a

bp2b <- ggplot(da2_sum %>% dplyr::distinct(well, .keep_all = T)) +
  aes(x = strain, 
      y = dauerFrac2, fill=strain,
      labs=strain) +
  ylim(0, 1)+
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  theme_bw() + 
  facet_wrap(~condition, ncol = 2) +
  labs(y = "dauer fraction (d/(d+other1))", title = "Assay2") +
  theme(legend.position = "none")
bp2b

cowplot::ggsave2(bp1, filename = "plots/assay1_dauerF_boxplot1.png", width = 5, height = 5)
cowplot::ggsave2(bp1b, filename = "plots/assay1_dauerF_boxplot2.png", width = 5, height = 5)
cowplot::ggsave2(bp2a, filename = "plots/assay2_dauerF_boxplot1.png", width = 5, height = 5)
cowplot::ggsave2(bp2b, filename = "plots/assay2_dauerF_boxplot2.png", width = 5, height = 5)
