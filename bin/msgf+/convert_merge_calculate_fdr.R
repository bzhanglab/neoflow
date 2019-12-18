library(PGA)
library(tidyverse)
library(data.table)

args <- commandArgs(TRUE)
input_files <- dir(args[1], pattern="mzid")
database <- args[2]
output <- args[1]
batchID <- args[3]

all_data <- data_frame()
for (file in input_files){
	sample=str_split(file, ".mzid")[[1]][1]
	parserGear(file=paste(args[1], file, sep=""),
        	db=database,
	        decoyPrefix="XXX_",
        	novelPrefix="VAR",
	        prefix=sample,
        	xmx=10,
	        thread=6,
        	alignment=0,
	        outdir=output)
	data <- fread(paste(output, sample, "-rawPSMs.txt", sep=""))
	data <- data %>% mutate(index = as.numeric(index) + 1)
	data$index <- paste(sample, data$index, sep=":")
	colnames(data)[2] <- "score"
	all_data <- bind_rows(all_data, data)
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
        xmx=10)
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
        xmx=10)
