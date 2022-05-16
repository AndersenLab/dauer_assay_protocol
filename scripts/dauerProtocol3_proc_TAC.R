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
  dplyr::filter(cluster_flag != T & well_edge_flag != T) %>%
  dplyr::mutate(model_select = ifelse(model_select == "dauerMod_NonOverlappingWorms.model.outputs", "dauer", "non-dauer")) # fix long model name

# calculate dauer classification traits
dauer_df <- flag_rm %>%
  dplyr::mutate(sd_int_gut = ifelse(model_select == "dauer", Worm_StdIntensity_dauerMod_StraightenedImage_T1of1_L2of3, Worm_StdIntensity_nondauerMod_StraightenedImage_T1of1_L2of3),
                mean_int_gut = ifelse(model_select == "dauer", Worm_MeanIntensity_dauerMod_StraightenedImage_T1of1_L2of3, Worm_MeanIntensity_nondauerMod_StraightenedImage_T1of1_L2of3),
                cv_int_gut = sd_int_gut/mean_int_gut,
                log_cv_int_gut = log(cv_int_gut),
                intInt_worm = Intensity_IntegratedIntensity_MaskedRFP,
                intInt_upperq_worm = Intensity_UpperQuartileIntensity_MaskedRFP,
                length_um = Worm_Length * 3.2937,
                length = Worm_Length,
                width_um = 2 * AreaShape_MaximumRadius * 3.2937,
                width = AreaShape_MaximumRadius * 2,
                area = AreaShape_Area,
                area_um = AreaShape_Area * 3.2937,
                x = AreaShape_Center_X,
                y = AreaShape_Center_Y) %>%
  dplyr::select(model, Metadata_Plate:FileName_RawRFP, strain:y) %>%
  dplyr::mutate(plate = Metadata_Plate,
                well = Metadata_Well,
                row = stringr::str_extract(Metadata_Well, pattern = "[A-Z]"),
                col = as.numeric(stringr::str_extract(Metadata_Well, pattern = "[0-9][0-9]"))) %>%
  dplyr::group_by(plate, well) %>%
  dplyr::arrange(y, x) %>%
  dplyr::mutate(w_id = row_number())

# look at the effect of dose on worm characteristics - multivariate multiple regression: https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/
m1 <- stats::lm(cbind(length_um, width_um, area_um) ~ ascr5_um + bead_time_min, data = dauer_df)
summary(m1)
car::Anova(m1)
lh.1.out <- car::linearHypothesis(m1, hypothesis.matrix = c("bead_time_min = 0"))
lh.1.out
# results suggest bead_time_min do not improve the model fit for size metrics

# test with intensity metrics
m2 <- stats::lm(cbind(log_cv_int_gut, cv_int_gut, intInt_upperq_worm, intInt_worm) ~ ascr5_um + bead_time_min, data = dauer_df)
summary(m2)
car::Anova(m2)
lh.2.out <- car::linearHypothesis(m2, hypothesis.matrix = c("ascr5_um = 0"))
lh.2.out
# both the bead time and the ascr concentration influence intensity metrics. Good, beads only influence intensity not size!

#==================================================#
# call dauers with a classifier
#==================================================#
# try using caret package short for Classification And REgression Training - http://topepo.github.io/caret/index.html 
# need to make a truth set - Use imageJ 

#==================================================#
# Make overlays for imageJ
#==================================================#
# set wells to make overlays for
wells <- dauer_df %>%
  dplyr::filter(col > 6, row %in% c("E", "F", "G", "H")) %>%
  dplyr::distinct(well) %>%
  dplyr::pull(well)

# make a function to plot worm ids and the centroid so I can use them in imageJ
d_overlay <- function(data, wells, well_radius = 825){
  data <- data 
  for(i in unique(wells)) {
    img <- png::readPNG(glue::glue("CellProfiler/20220427_dauerProtocol3/output/Analysis-20220505/20220427-dauerProtocol3-p001-m2X_{i}_w1_overlay.png"))
    #img <- png::readPNG(glue::glue("CellProfiler/20220427_dauerProtocol3/output/Analysis-20220505/20220427-dauerProtocol3-p001-m2X_E07_w1_overlay.png"))
    
    h <- dim(img)[1] # image height
    w <- dim(img)[2] # image width
    well_radius = well_radius
    
    # make the plot for a well
    well_img <- data %>% dplyr::filter(well == i) %>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = x, y = y, label = as.character(w_id)) +
      ggplot2::annotation_custom(grid::rasterGrob(img, width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc")), 0, w, 0, -h) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::geom_text(size = 2, nudge_y = 15) +
      ggplot2::geom_point(shape = 21, size = 1, alpha = 0.5, aes(color = model_select)) +
      ggplot2::annotate("path", x = w/2 + well_radius * cos(seq(0, 2 * pi, length.out = 100)),
                        y = h/2 + well_radius * sin(seq(0, 2 * pi, length.out = 100)),
                        color = "red", 
                        alpha = 0.25) +
      ggplot2::scale_x_continuous(expand=c(0,0),limits=c(0,w)) +
      ggplot2::scale_y_reverse(expand=c(0,0),limits=c(h,0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")
    
    # save the plot
    cowplot::ggsave2(well_img, filename = glue::glue("CellProfiler/20220427_dauerProtocol3/output/Analysis-20220505/imageJ_overlays/20220427-dauerProtocol3-p001-m2X_{i}.png"),
                     dpi = 300, width = 6.827, height = 6.827) # set to 2048 pixels
  }
}

# use it
d_overlay(data = dauer_df, wells = wells)

#==============================================#
# Test joining of truth set from imageJ
#==============================================#
# join the truth set - i'm seeing nearly 100% dauers in high dose
truth_join <- data.table::fread("CellProfiler/20220427_dauerProtocol3/dauerProtocol3_20220514_dauer_truthset.csv") %>%
  dplyr::rename(xt = x, yt = y, w_id = worm_id) %>%
  dplyr::left_join(., dauer_df) %>%
  dplyr::mutate(x_diff = xt-x,
                y_diff = yt-y,
                abs_diff = case_when(sign(x_diff) == 1 & sign(y_diff) == 1 ~ x_diff + y_diff,
                                     sign(x_diff) == 1 & sign(y_diff) == -1 ~ (x_diff + (-1*y_diff)),
                                     sign(x_diff) == -1 & sign(y_diff) == 1 ~ -1*x_diff + y_diff,
                                     sign(x_diff) == -1 & sign(y_diff) == -1 ~ -1*x_diff + -1*y_diff,
                                     TRUE ~ NA_real_))

# the positions look good for the truth set - I think they are joined properly
ggplot(truth_join) +
  geom_histogram(aes(x = abs_diff))

#=======================================#
# setup the classifier
#=======================================#
# explore histograms for dauer and non-dauer in truth set
hist_df <- truth_join %>%
  tidyr::pivot_longer(cols = names(dplyr::select(truth_join, sd_int_gut:area_um)))

# plot the histograms
all_trait_hist <- ggplot(hist_df %>% dplyr::filter(truth_class != "prune" & !(name %in% c("area", "length", "width")))) +
  geom_histogram(aes(x = value, y = stat(density), fill = truth_class), bins = 50, position = "identity", alpha = 0.5) +
  geom_density(aes(x = value, fill = truth_class), position = "identity", alpha = 0.5) +
  facet_wrap(~name, scales = "free", ncol = 3) +
  labs(title = "Worm traits by dauer class: dauer n=112, non-dauer n=119") +
  theme_bw() +
  theme(plot.title = element_text(size=12))
all_trait_hist

cowplot::ggsave2(all_trait_hist, filename = "plots/all_trait_distributions_truth_class.png", width = 7.5, height = 7.5)

#========================================+#
# OK now classifier
#=========================================#
# make the classifier data
dauer_classifier_df <- truth_join %>%
  dplyr::select(plate, well, w_id, Class = truth_class, sd_int_gut:area_um) %>%
  dplyr::select(-length, -width, -area) %>%
  dplyr::filter(Class != "prune") %>%
  dplyr::arrange(plate, well, w_id) %>%
  dplyr::select(-plate, -well, -w_id)

# setup the training and testing observations
set.seed(1)
inTraining <- caret::createDataPartition(dauer_classifier_df$Class, p = .5, list = FALSE) # take 50% of data randomly
training <- dauer_classifier_df[ inTraining,]
testing  <- dauer_classifier_df[-inTraining,]

# setup fit control for classifier
fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# make a boosted tree model (sbm) with gbm package - no reason for this other than it's in the example
gbmFit1 <- caret::train(Class ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)
gbmFit1

# make a new tuning grid
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:20)*10, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

# make a new model with the new tuning grid
set.seed(1)
gbmFit2 <- train(Class ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = gbmGrid)
gbmFit2

# see accuracy of models over tuning grid
trellis.par.set(caretTheme())
plot(gbmFit1)

trellis.par.set(caretTheme())
plot(gbmFit2)

# predict on testing set
predict(gbmFit2, newdata = testing, type = "prob")

#=============================================#
# caret example for classifier
#============================================#
# get data
library(mlbench)
data(Sonar)
str(Sonar[, 1:10])
# partition to training and testing
library(caret)
set.seed(998)
inTraining <- caret::createDataPartition(Sonar$Class, p = .75, list = FALSE) # take 75% of data randomly
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]
# setup fit control
fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

