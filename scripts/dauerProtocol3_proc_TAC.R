library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#==============================================================#
# Process CellProfiler data
#==============================================================#
# read in the CellProfiler output from the process_dauer.cpproj run
load("CellProfiler/20220427_dauerProtocol3/cp_data/20220427_dauerProtocol3_Analysis-20220505.RData")

# read in design and format
design <- data.table::fread("CellProfiler/20220427_dauerProtocol3/dauer_protocol3_design.csv") %>%
  dplyr::mutate(well = case_when(col %in% c(1:9) ~ paste0(row, "0", as.character(col)),
                                 !(col %in% c(1:9)) ~ paste0(row, as.character(col))),
                plate = case_when(plate %in% c(1:9) ~  paste0("p00", as.character(plate)),
                                  plate %in% c(10:99) ~  paste0("p0", as.character(plate)),
                                  !(plate %in% c(1:99)) ~  paste0("p", as.character(plate))),
                id = paste0(date, plate, well)) %>%
  tidyr::separate(condition, into = c("bead_time_min", "ascr5_um"), sep = "_") %>%
  dplyr::mutate(bead_time_min = case_when(bead_time_min == "4hr" ~ 240,
                                          bead_time_min != "4hr" ~ 20),
                bead_time_min = as.numeric(bead_time_min),
                ascr5_um = stringr::str_remove(ascr5_um, pattern = "uMascr"),
                ascr5_um = ifelse(ascr5_um == "control", 0, ascr5_um),
                ascr5_um = as.numeric(ascr5_um)) %>%
  dplyr::select(id, strain, bead_time_min, ascr5_um)

# join the data from the models
df_join1 <- dplyr::bind_rows(dauerMod_NonOverlappingWorms.model.outputs, nondauerMod_NonOverlappingWorms.model.outputs)

# try using easyXpress model selection function
df_msel <- easyXpress::modelSelection(df_join1)

# lets join the design
df_full <- df_msel %>%
  dplyr::mutate(id = paste0(Metadata_Date, Metadata_Plate, Metadata_Well)) %>%
  dplyr::left_join(design)

# edge_flag
edge_flagged <- easyXpress::edgeFlag(df_full)

# remove cluster flags
flag_rm <- edge_flagged %>%
  dplyr::filter(cluster_flag != T) %>%
  dplyr::mutate(model_select = ifelse(model_select == "dauerMod_NonOverlappingWorms.model.outputs", "dauer", "non-dauer")) # fix long model name

# caluclate dauer classification traits
dauer_df <- flag_rm %>%
  dplyr::mutate(sd_int_gut = ifelse(model_select == "dauer", Worm_StdIntensity_dauerMod_StraightenedImage_T1of1_L2of3, Worm_StdIntensity_nondauerMod_StraightenedImage_T1of1_L2of3),
                mean_int_gut = ifelse(model_select == "dauer", Worm_MeanIntensity_dauerMod_StraightenedImage_T1of1_L2of3, Worm_MeanIntensity_nondauerMod_StraightenedImage_T1of1_L2of3),
                cv_int_gut = sd_int_gut/mean_int_gut,
                log_cv_int_geut = log(cv_int_gut),
                length_um = Worm_Length * 3.2937,
                length = Worm_Length,
                width_um = 2 * AreaShape_MaximumRadius * 3.2937,
                width = AreaShape_MaximumRadius * 2,
                area = AreaShape_Area,
                area_um = AreaShape_Area * 3.2937,
                x = AreaShape_Center_X,
                y = AreaShape_Center_Y) %>%
  dplyr::select(model, Metadata_Plate:FileName_RawRFP, strain:y)

#===========================================+#
# WORKING HERE
#=============================================#
# summarize
proc_df <- flag_rm %>%
  dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
  dplyr::mutate(median_legnth = median(Worm_Length)) %>%
  dplyr::distinct(Metadata_Plate, Metadata_Well, .keep_all = T)


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