library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Define a vector of your experiment directories. Can load more
dirs <- c("data/raw/20220407_dauerProtocol1")

# Read in the data
raw <- easysorter::read_data(dirs)

# look at score data
score_proc <- as.data.frame(raw[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE")

# plot GFP by TOF
ggplot(score_proc %>% dplyr::group_by(condition)) +
  aes(x = TOF, y = green) +
  geom_point() +
  facet_grid(~condition)

