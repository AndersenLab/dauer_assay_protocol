d3 <- score_proc1 %>%
  dplyr::distinct(date, plate, row, col, strain, condition) %>%
  dplyr::mutate(well = paste0(row))

write.csv(d3, file = "~/Desktop/dauer_protocol3_design.csv")

# assay 4
dirs1 <- c("/Users/timcrombie/repos/dauer_assay_protocol/data/raw/20220504_dauerProtocol4") 

# Read in the data using easysorter: https://github.com/AndersenLab/easysorter
standard_raw <- easysorter::read_data(dirs1)

# look at score data for each color channel
score_proc1 <- as.data.frame(standard_raw[1]) %>%
  dplyr::filter(call50 != "bubble", contamination != "TRUE", !is.na(strain)) %>%
  dplyr::mutate(well = paste(plate,row,col, sep = ""))

d4 <- score_proc1 %>%
  dplyr::distinct(date, plate, row, col, strain, condition) %>%
  dplyr::mutate(well = paste0(row, col))

write.csv(d4, file = "~/Desktop/dauer_protocol4_design.csv")
