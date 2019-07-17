#!/usr/bin/env nextflow

/* Prints help when asked for and exits */


def helpMessage() {
    log.info"""
    =========================================
    neoflow => variant annotation and customized database construction
    =========================================
    Usage:
    nextflow run neoflow-db.nf
    Arguments:
      --customprodbj_file         CustomProdbj information result file
      --customprodbj_fasta		  CustomProdbj db file
      --hla_class_results		  Optitype hla class results file
      --out_dir                   Output folder, default is "./"
      --sample_ID                 Sample ID

    """.stripIndent()
}

params.help       = false
params.out_dir 	  = "./"
params.customprodbj_file = "./"
params.sample_ID  = "."
params.hla_class_results = "."
params.customprodbj_fasta = "./"

customprodbj_infor_results = file(params.customprodbj_file)
customprodbj_fasta = file(params.customprodbj_fasta)
hla_class_file = file(params.hla_class_results)
sample_id = params.sample_ID
out_dir = file(params.out_dir)


process Split_infor {
	tag "Separate customprodbj results"

	input:
	file customprodbj_infor_results

	output:
	file("${sample_id}_new_somatic_varInfo.txt") into sample_somatic_chanel

	script:
	"""
	#!/usr/local/bin/Rscript

	library(tidyverse)
	library(data.table)

	merge_infor <- fread("$customprodbj_infor_results")
	column_name <- as.name(paste("sample:", "$sample_id", sep = ""))
	filter_infor <- merge_infor %>%
  		filter(str_detect(UQ(column_name), "1,"))
  	write.table(filter_infor, paste("$sample_id", "_new_somatic_varInfo.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
	"""
}


process MHC_peptide_prediction {
	tag "netmhc pan prediction"

	publishDir "${out_dir}", mode: "copy", overwrite: true
	
	input:
	file hla_class_file
	file customprodbj_fasta
	file(somatic_file) from sample_somatic_chanel

	output:
	file("*.csv") into final_db

	script:
	"""
	python3 ${baseDir}/binding_prediction_debug.py \
		-id $sample_id \
		-a $hla_class_file \
		-fa $customprodbj_fasta \
		-txt $somatic_file \
		-o . \
		-netmhcpan /usr/local/netMHCpan-4.0/netMHCpan \
		-nt 6
	"""


}
