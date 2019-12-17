library(tidyverse)

args <- commandArgs(T)
results_folder <- args[1]

pga_result <- read.delim(paste(results_folder, "/global_fdr/pga-peptideSummary.txt", sep=""))
colnames(pga_result)[1] <- "psm_id"
pga_result <- separate(pga_result, psm_id, into=c("fraction", "index"), sep=":", remove=F)

split_result <- split(pga_result, pga_result$fraction)

a <- lapply(names(split_result), function(x){
            write.table(split_result[[x]], paste(results_folder, "/global_fdr/", x, "_fraction.txt", sep = ""), row.names=F, quote=F, sep="\t")})
