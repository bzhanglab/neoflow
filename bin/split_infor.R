library(tidyverse)

args=commandArgs(T)
path <- args[1]
experiment <- args[2]
sample <- args[3]

merge_data <- read.delim(paste(path, experiment, "_germline_somatic-varInfo.txt", sep = ""))
column_name <- as.name(paste("sample.", sample, sep = ""))

filter_infor <- merge_data %>%
  filter(str_detect(UQ(column_name), "1,"))

write.table(filter_infor, file = paste(path, sample, "_new_somatic_varInfo.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
