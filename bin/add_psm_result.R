library(tidyverse)

args <- commandArgs(TRUE)
peptide_fdr_file <- args[1]
psm_fdr_file <- args[4]
mgf_index_file <- args[2]
output_file <- args[3]

peptide_fdr_data <- read.delim(peptide_fdr_file)
mgf_index_data <- read.delim(mgf_index_file)
filter_psm_data <- peptide_fdr_data %>%
  filter(Qvalue <= 0.01, str_detect(protein, "NM"), !str_detect(protein, "cont_"))
filter_var_psm_data <- peptide_fdr_data %>%
  filter(Qvalue <= 0.01, !str_detect(protein, "NM"), !str_detect(protein, "cont_"))

psm_fdr_data <- read.delim(psm_fdr_file)
psm_fdr_pass_peptide <- psm_fdr_data %>%
        filter(peptide %in% filter_psm_data$peptide)

process_data <- psm_fdr_pass_peptide %>%
	select(index, peptide, mods, evalue) %>%
	mutate(index = as.numeric(index))
process_var_data <- filter_var_psm_data %>%
        select(index, peptide, mods, evalue) %>%
        mutate(index = as.numeric(index))

process_data_rt <- left_join(process_data, mgf_index_data, by = "index")
process_var_data_rt <- left_join(process_var_data, mgf_index_data, by = "index")
process_data_rt$rt <- process_data_rt$rt/60
process_var_data_rt$rt <- process_var_data_rt$rt/60

write.table(process_data_rt, paste(output_file, "_normal_psm.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(process_var_data_rt, paste(output_file, "_variant_psm.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")
