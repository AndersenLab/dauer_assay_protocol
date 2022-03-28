library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Define a vector of your experiment directories. Can load more
dirs <- c("data/20170912_dauerGWA1A",
          "data/20170913_dauerGWA2A",
          "data/20170917_dauerGWA3A",
          "data/20170918_dauerGWA4A",
          "data/20171109_dauerGWA5A",
          "data/20170912_dauerGWA1B",
          "data/20170913_dauerGWA2B",
          "data/20170917_dauerGWA3B",
          "data/20170918_dauerGWA4B",
          "data/20171109_dauerGWA5B",
          "data/20170912_dauerGWA1C",
          "data/20170913_dauerGWA2C",
          "data/20170917_dauerGWA3C",
          "data/20170918_dauerGWA4C",
          "data/20171109_dauerGWA5C")

# Read in the data
raw <- easysorter::read_data(dirs)

`# Remove all data from the contaminated wells
raw_nocontam <- easysorter::remove_contamination(raw)

# Summarize the data to see whats up with strains
summedraw <- easysorter::sumplate(raw_nocontam, directories = TRUE, quantiles = TRUE)

# look at raw data
score <- raw_nocontam[[5]][[1]]
setup <- raw_nocontam[[5]][[2]]
# OK these are identical need, pick one and process with EM (y = GFP, x = TOF within well)

# EM with GWA 5
proc_dat <- score %>%
  dplyr::filter(strain == "JU310" | strain == "MY2693" | strain == "N2") %>%
  dplyr::group_by(round, assay, plate, row, col) %>%
  dplyr::mutate(well_n = n(),
                well_name = paste0(round, assay, plate, row, col)) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(well_n > 20 & well_n < 80) %>% 
  dplyr::filter(condition %in% c("ascr5.800nM","EtOH"))
  
# view bubbles?
table(proc_dat$call50, proc_dat$stage)

# well counts
table(proc_dat$well_name, proc_dat$strain)

# plot GFP by TOF
ggplot(proc_dat %>% dplyr::group_by(strain)) +
  aes(x = TOF, y = green) +
  geom_point() +
  facet_grid(~strain + condition)









