library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#==============================================================#
# Part 1: Process CellProfiler data
#==============================================================#
# use easyXpress to process raw images - run once
#easyXpress::tidyProject(project_dir = "/Volumes/ECA_Image/Tim/20220407_dauerProtocol1")

# read in the CellProfiler output from the process_dauer.cpproj run
now2 <- data.table::fread("~/repos/dauer_assay_protocol/CellProfiler/20220407_dauerProtocol1/output/20220414_run6/NonOverlappingWorms.csv")

# lets just plot the integrated RFP intensity / area for NonOverlappingWorms
now_proc2 <- now2 %>%
  dplyr::select(assay = Metadata_Assay, date = Metadata_Date, plate = Metadata_Plate,
                well = Metadata_Well, Parent_WormObjects, ObjectNumber, now_area = AreaShape_Area, 
                AreaShape_Center_X, AreaShape_Center_Y, now_length = Worm_Length,
                ii_RFP = Intensity_IntegratedIntensity_MaskedRFP,
                Intensity_MaxIntensity_MaskedRFP,
                Intensity_MeanIntensity_MaskedRFP,
                Intensity_MedianIntensity_MaskedRFP,
                Intensity_StdIntensity_MaskedRFP,
                Intensity_MaxIntensity_MaskedRFP,
                sd_int2 = Worm_StdIntensity_StraightenedImage_T1of1_L2of3,
                mean_int2 = Worm_MeanIntensity_StraightenedImage_T1of1_L2of3) %>%
  dplyr::mutate(rescale_ii_RFP = ii_RFP * 65536,
                cv_int = Intensity_StdIntensity_MaskedRFP/Intensity_MeanIntensity_MaskedRFP,
                log_cv_int = log(cv_int),
                cv_int2 = sd_int2/mean_int2,
                log_cv_int2 = log(cv_int2))

# plot all the intensity traits - we could do PCA to find feeding worms?
# corrlation of eigen values with env parameters
trait_cor <- round(cor(dplyr::select_if(now_proc2, is.double), use = "complete.obs", method = "pearson"), 4)
trait_heatmap <- ggplotify::as.ggplot(pheatmap::pheatmap(trait_cor,
                                            cutree_rows = NA,
                                            cutree_cols = NA,
                                            display_numbers = trait_cor))
trait_heatmap

cowplot::ggsave2(trait_heatmap, filename = "plots/all_measures_cor_heatmap.png", width = 12.5, height = 12.5)

# make a distribution for each trait
hist_df <- now_proc2 %>%
  tidyr::pivot_longer(cols = names(dplyr::select_if(now_proc2, is.double)))

# plot the histograms
all_trait_hist <- ggplot(hist_df) +
  geom_histogram(aes(x = value), bins = 100) +
  facet_wrap(~name, scales = "free", ncol = 4) +
  theme_bw()
all_trait_hist

log_cv_int_hist <- ggplot(hist_df %>% dplyr::filter(name == "log_cv_int" | name == "log_cv_int2")) +
  geom_histogram(aes(x = value), bins = 100) +
  theme_bw() +
  geom_vline(xintercept = -2.9, linetype = "dashed", color = "red") +
  facet_wrap(~name, scales = "free", ncol = 4)
log_cv_int_hist

cowplot::ggsave2(all_trait_hist, filename = "plots/all_measures_histogram.png", width = 7.5, height = 7.5)
cowplot::ggsave2(log_cv_int_hist, filename = "plots/log_cv_int_measures_histogram.png", width = 7.5, height = 7.5)

#===================================================================================#
# Try polarizing to feeding vs not feeding with log_cv_int
#===================================================================================#
now_proc3 <- now_proc2 %>%
  dplyr::mutate(fed = ifelse(log_cv_int2 > -2.9, T, F)) # could impirically find this line with azide treated then bead incubated worms?

# Read the rescaled RFP PNG
img <- png::readPNG(glue::glue("{getwd()}/CellProfiler/20220407_dauerProtocol1/output/20220414_run6/20220407-dauerProtocol1-p002-m2X_A01_w2_NonOverlappingWorms_RFP_mask.png"))
img2 <- png::readPNG(glue::glue("{getwd()}/CellProfiler/20220407_dauerProtocol1/output/20220414_run6/20220407-dauerProtocol1-p002-m2X_A01_w1_overlay.png"))
img3 <- png::readPNG(glue::glue("{getwd()}/CellProfiler/20220407_dauerProtocol1/output/20220414_run6/20220407-dauerProtocol1-p002-m2X_G07_w2_NonOverlappingWorms_RFP_mask.png"))
img4 <- png::readPNG(glue::glue("{getwd()}/CellProfiler/20220407_dauerProtocol1/output/20220414_run6/20220407-dauerProtocol1-p002-m2X_G07_w1_overlay.png"))

h<-dim(img)[1] # image height
w<-dim(img)[2] # image width

#plot rescaled RFP PNG with intensity traits from CP
well_img <- now_proc3 %>% dplyr::filter(well == "A01") %>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = log_cv_int2, shape = fed) +
  ggplot2::annotation_custom(grid::rasterGrob(img, width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::scale_shape_manual(values = c(21,24)) +
  ggplot2::geom_point(alpha = 0.75, size=4) +
  ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
  ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()
well_img

well_img2 <- now_proc3 %>% dplyr::filter(well == "A01") %>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = log_cv_int2, shape = fed) +
  ggplot2::annotation_custom(grid::rasterGrob(img2, width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::scale_shape_manual(values = c(21,24)) +
  ggplot2::geom_point(alpha = 0.75, size=4) +
  ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
  ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()
well_img2

well_img3 <- now_proc3 %>% dplyr::filter(well == "G07") %>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = log_cv_int2, shape = fed) +
  ggplot2::annotation_custom(grid::rasterGrob(img3, width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::scale_shape_manual(values = c(21,24)) +
  ggplot2::geom_point(alpha = 0.75, size=4) +
  ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
  ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()

well_img4 <- now_proc3 %>% dplyr::filter(well == "G07") %>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, fill = log_cv_int2, shape = fed) +
  ggplot2::annotation_custom(grid::rasterGrob(img4, width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::scale_shape_manual(values = c(21,24)) +
  ggplot2::geom_point(alpha = 0.75, size=4) +
  ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
  ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()

# save well
cowplot::ggsave2(well_img, filename = "plots/A01_RFP_fed.png", width = 7.5, height = 7.5)
cowplot::ggsave2(well_img2, filename = "plots/A01_BF_fed.png", width = 7.5, height = 7.5)
cowplot::ggsave2(well_img3, filename = "plots/G07_BF_fed.png", width = 7.5, height = 7.5)
cowplot::ggsave2(well_img4, filename = "plots/G07_BF_fed.png", width = 7.5, height = 7.5)

# well histogram
a01_hist <- ggplot(now_proc2 %>% dplyr::filter(well == "A01")) +
  geom_histogram(aes(x = log_cv_int2)) +
  geom_vline(xintercept = -2.9, linetype = "dashed", color = "red") +
  theme_bw()

g07_hist <- ggplot(now_proc2 %>% dplyr::filter(well == "G07")) +
  geom_histogram(aes(x = log_cv_int2)) +
  geom_vline(xintercept = -2.9, linetype = "dashed", color = "red") +
  theme_bw()

# save the histogram
cowplot::ggsave2(a01_hist, filename = "plots/A01_well_hist_RFP.png", width = 7.5, height = 7.5)
cowplot::ggsave2(g07_hist, filename = "plots/G07_well_hist_RFP.png", width = 7.5, height = 7.5)

#==============================================================#
# Part 2: Process Sorter data
#==============================================================#
# Define a vector of your experiment directories. Can load more
dirs <- c("data/raw/20220407_dauerProtocol1")

# Read in the data using easysorter: https://github.com/AndersenLab/easysorter
raw <- easysorter::read_data(dirs)

# look at score data for each color channel
score_proc <- as.data.frame(raw[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE") %>%
  tidyr::pivot_longer(cols = c(green, yellow, red))

# plot GFP by TOF
a1_channels <- ggplot(score_proc %>% dplyr::filter(!is.na(condition))) +
  aes(x = TOF, y = value, fill = condition) +
  geom_point(size=1, alpha = 0.3, shape = 21) +
  facet_grid(~condition) +
  labs(y = "Fluorescence (AU)") +
  facet_grid(~condition + name) +
  theme_bw() +
  theme(legend.position = "none")
a1_channels

cowplot::ggsave2(a1_channels, filename = "plots/dauerProtocol1_sorter_FLchannels.png")

# look at the RFP for well A01 to see if it matches our expectations
a1_channels_a01 <- ggplot(score_proc %>% dplyr::filter(col == "1" & row == "A")) +
  aes(x = TOF, y = value, fill = condition) +
  geom_point(shape = 21) +
  facet_grid(~condition + name) +
  theme_bw() +
  theme(legend.position = "none")
a1_channels_a01

cowplot::ggsave2(a1_channels_a01, filename = "plots/dauerProtocol1_sorter_FLchannels_A01.png")

# EM clustering in DI condition
score_proc_di <- as.data.frame(raw[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE") %>%
  dplyr::filter(!is.na(strain)) %>%
  dplyr::filter(condition == "ascr5.800nM") %>%
  dplyr::filter(TOF>100)
y = score_proc_di[c("TOF", "green", "EXT")]
ret.em_K4 <- EMCluster::init.EM(y, nclass = 4, method = "em.EM")
y_merge_K4 <- data.frame(score_proc_di, ret.em_K4[6])

# plot to define classes
ggplot(y_merge_K4) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0, 500) +
  xlim(0, 1000)
#--------------------------------------------------------#
# Daehan's EM clustering EXAMPLE WITH EMcluster
#--------------------------------------------------------#
raw_dauersrgCRISPR1A <- easysorter::read_data("~/Dropbox/HTA/Results/20181122_dauersrgCRISPR1A/") %>%
  dplyr::mutate(well = paste(plate,row,col, sep = "")) %>%
  dplyr::filter(well != "4G8") ### contamination

# EM clustering in DI condition
df_raw <- raw_dauersrgCRISPR1A %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(strain %in% c("JU346", "ECA1089", "ECA1119", "ECA1121", "NIC166", "ECA1091", "ECA1120", "ECA1122")) %>%
  dplyr::filter(condition == "ascr5.800nM") %>%
  dplyr::filter(TOF>100)

y <- df_raw[c("TOF", "green", "EXT")]
ret.em_K4 <- EMCluster::init.EM(y, nclass = 4, method = "em.EM")
y_merge_K4 <- data.frame(df_raw, ret.em_K4[6])

#------------------------#
# plot to define classes
#------------------------#
ggplot(y_merge_K4) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0, 500) +
  xlim(0, 1000)
# Caution!!! factor label changes everytime!!!
# labels=c("crap1", "crap2", "non-dauer", "dauer")
# labels=c("crap1", "dauer",  "crap2", "non-dauer")
# labels=c("non-dauer", "crap1", "crap2", "dauer")

y_merge_K4$class <- factor(y_merge_K4$class, levels=c(1,3,2,4), labels=c("non-dauer", "crap1", "crap2", "dauer"))

y_merge_K4 %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  ggplot(.) +
  aes(x = TOF, 
      y = green, 
      scale = T, color = class) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0, 500) +
  xlim(0, 1000)

# dauer clustering in EtOH control - applying model from DI condition
df_raw_EtOH <- raw_dauersrgCRISPR1A %>%
  dplyr::filter(!is.na(strain)) %>%
  dplyr::filter(condition == "EtOH")

y_et <- df_raw_EtOH[c("TOF","green", "EXT")]
ret.em_K4_et <- EMCluster::assign.class(y_et, ret.em_K4)
y_merge_K4_et <- data.frame(df_raw_EtOH, ret.em_K4_et[6])

#------------------------#
# plot to define classes
#------------------------#
ggplot(y_merge_K4_et) +
  aes(x = TOF, y = green, scale = T, color = as.factor(class)) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0, 500) +
  xlim(0, 1000)
# Caution!!! factor label changes everytime!!!
# labels=c("crap1", "crap2", "non-dauer", "dauer")
# labels=c("crap1", "dauer",  "crap2", "non-dauer")
# labels=c("non-dauer", "crap1", "crap2", "dauer")

y_merge_K4_et$class <- factor(y_merge_K4_et$class, levels=c(1,3,2,4), labels=c("non-dauer", "crap1", "crap2", "dauer"))

y_merge_K4_et %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  ggplot(.) +
  aes(x = TOF, 
      y = green, 
      scale = T, color = class) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0,500) +
  xlim(0,1000)

## merge and dauer fraction quanitify
df_dauersrgCRISPR1A <- rbind(y_merge_K4, y_merge_K4_et)  %>%
  dplyr::filter(green>6 & TOF < 1000 & green < 500) 

df_dauersrgCRISPR1A %>%
  dplyr::filter(!is.na(strain)) %>%
  dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  ggplot(.) +
  aes(x = TOF, 
      y = green, 
      scale = T, color = class) +
  geom_point(size=1, alpha = 0.3)+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1)) +
  theme_bw() +
  ylim(0,500) +
  xlim(0,1000) +
  facet_wrap(assay~condition)

df_dauersrgCRISPR1A %>%
  dplyr::filter(!is.na(strain)) %>%
  #dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  ggplot(.) +
  aes(x = class) +
  geom_bar() +
  facet_wrap(plate~condition)

df_dauersrgCRISPR1A_dauer <- df_dauersrgCRISPR1A %>%
  dplyr::filter(class %in% c("dauer", "non-dauer")) %>%
  dplyr::mutate(dauer = ifelse((class == "dauer"), 1, 0), worm = 1) %>%
  dplyr::select(date, experiment, round, assay, plate, condition, control, strain, row, col, dauer, worm) %>%
  dplyr::group_by(date, experiment, round, assay, plate, condition, control, strain, row, col) %>%
  dplyr::summarise(dauer_fraction = mean(dauer), n = sum(worm)) %>%
  dplyr::filter(condition != "NA") %>%
  dplyr::ungroup()

df_dauersrgCRISPR1A_dauer$condition <- factor(df_dauersrgCRISPR1A_dauer$condition, levels = c("EtOH", "ascr5.800nM"), labels = c("Control", "ascr#5_800nM"))
df_dauersrgCRISPR1A_dauer$strain <- factor(df_dauersrgCRISPR1A_dauer$strain, levels = c("JU346", "ECA1119", "ECA1089", "ECA1121", "NIC166", "ECA1120", "ECA1091", "ECA1122", "PB303", "ECA1118", "ECA1126", "ED3017", "ECA1115", "JU2106", "ECA1116", "ECA1124"), labels = c("JU346(++)", "JU346(0+)", "JU346(+0)", "JU346(00)", "NIC166(++)", "NIC166(0+)", "NIC166(+0)", "NIC166(00)", "PB303", "ECA1118", "ECA1126", "ED3017", "ECA1115", "JU2106", "ECA1116", "ECA1124"))

df_dauersrgCRISPR1A_dauer %>%
  dplyr::filter(n>=20 & n<=80) %>%
  dplyr::filter(strain %in% c("JU346(++)", "JU346(0+)", "JU346(+0)", "JU346(00)", "NIC166(++)", "NIC166(0+)", "NIC166(+0)", "NIC166(00)")) %>%
  ggplot(.) +
  aes(x = strain, 
      y = dauer_fraction, fill=condition,
      labs=strain) +
  ylim(0, 1)+
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7)+ 
  scale_fill_brewer(palette = 'Set1') +
  theme_bw() +
  theme(text = element_text(size=11, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size=13, color = "black"), axis.text = element_text(size=11, color = "black"), axis.text.x = element_text(hjust=1, vjust=0.5, angle=90)) +
  labs(x="Strain", y="Dauer fraction")

# facet by condition

## for JU346 and NIC166

df_dauersrgCRISPR1A_dauer %>%
  dplyr::filter(n>=20 & n<=80) %>%
  dplyr::filter(strain %in% c("JU346(++)", "JU346(0+)", "JU346(+0)", "JU346(00)", "NIC166(++)", "NIC166(0+)", "NIC166(+0)", "NIC166(00)")) %>%
  ggplot(.) +
  aes(x = strain, 
      y = dauer_fraction, fill=condition,
      labs=strain) +
  geom_point(aes(fill=condition), position=position_jitterdodge(), alpha = .7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7)+ 
  scale_fill_brewer(palette = 'Set1') +
  theme_bw() +
  theme(text = element_text(size=11, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size=13, color = "black"), axis.text = element_text(size=11, color = "black"), axis.text.x = element_text(hjust=1, vjust=0.5, angle=90), legend.position = 'none') +
  labs(x="Strain", y="Dauer fraction") +
  facet_wrap(~condition, ncol=1)
