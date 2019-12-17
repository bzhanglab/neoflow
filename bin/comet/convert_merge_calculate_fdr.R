library(PGA)
library(tidyverse)
library(data.table)

args <- commandArgs(TRUE)
input_files <- dir(args[1], pattern="_rawResults.txt")
database <- args[2]
output <- args[1]
batchID <- args[3]

all_data <- data_frame()
for (file in input_files){
	sample=str_split(file, "_rawResults")[[1]][1]
	comet_result <- read.delim(file) %>%
		select(scan, "xcorr", charge, exp_neutral_mass, plain_peptide, protein, modifications)
	colnames(comet_result) <- c("index", "score", "charge", "mass", "peptide", "protein", "mods")
	comet_result$protein <- str_replace_all(comet_result$protein, ",", ";")
	comet_result$mods <- str_replace_all(comet_result$mods, ",", ";")
	comet_result <- comet_result %>% mutate(index = as.numeric(index))
	comet_result$index <- paste(sample, comet_result$index, sep=":")
	all_data <- bind_rows(all_data, comet_result)
}

write.table(all_data, paste(output, batchID, "-rawPSMs.txt", sep=""), row.names=F, quote=F, sep="\t")

calculateFDR(psmfile=paste(output, batchID, "-rawPSMs.txt", sep=""),
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=TRUE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output, "peptide_level/global_fdr", sep=""),
        xmx=20)
calculateFDR(psmfile=paste(output, batchID, "-rawPSMs.txt", sep=""),
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=F,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output, "psm_level/global_fdr", sep=""),
        xmx=20)
