library(PGA)

args <- commandArgs(T)

folder <- args[1]
experiment <- args[2]
reference_folder <- args[3]

buildTargetDecoyDB(db=paste(folder, experiment, "_germline_somatic-var.fasta.new.fasta", sep=""),
	decoyPrefix="XXX_",
	cont_file=paste(reference_folder, "contaminants.fasta", sep=""),
	output=paste(folder, experiment, "_target_decoy.fasta", sep=""))

