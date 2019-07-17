#!/usr/bin/env nextflow

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => variant annotation and customized database construction
    =========================================
    Usage:
    nextflow run neoflow-db.nf
    Mandatory arguments:
      --file                      The prefix for novel proteins in database, default is "VAR"
      --out_dir                   Output folder, default is "./"
      --cpu                       The number of CPUs

    Annovar arguments:
      --buildver                  Enzyme used for protein digestion. Default is trypsin
      --protocol                  The number of allowed missed cleavage sites. Default is 2
      --operation                 The error window on experimental peptide mass values, default is 10
      --anno_dir                  The unit of --tol, ppm or Da. Default is ppm

    """.stripIndent()
}


params.help       = false 
// Annovar parameters
params.buildver   = "hg38"
params.protocol   = "refGene"
params.operation  = "g"
params.anno_dir = false // annotation data folder for annovar 


params.cpu = 0
params.file = false // mapping file
params.out_dir = "./"
params.prefix = "test"

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


out_dir = file(params.out_dir)
cpus = params.cpu

// parameters for varaint annotation: Annovar
annovar_buildver = params.buildver // 
annovar_protocol = params.protocol //
annovar_operation = params.operation 
annovar_anno_dir  = file(params.anno_dir)

out_dir = file(params.out_dir)
out_prefix = params.prefix
// sample\texperiment\tfile\tfile_type
// a\te1\ta.vcf;b.vcf\tsomatic;germline
mapping_file = file(params.file)



process pre_processing {
	tag "${mapping_file}"

	input:
	file mapping_file

	output:
	file "*-mapping_file.tsv" into single_experiment_map_files

	script:
	"""
	#!/usr/bin/env /usr/local/bin/Rscript
	library(dplyr)
	library(readr)
	a = read.delim("${mapping_file}")
	experiment_names = unique(a[,"experiment"])
	for(f in experiment_names){
		dat = a %>% filter(experiment == f)
		out_file = paste(f,"-mapping_file.tsv",sep="")
		write_tsv(dat,out_file)
	}
	
	"""
}


process variant_annotation {
	tag "${single_experiment_map_file}"

	input:
	file single_experiment_map_file from single_experiment_map_files.flatten()
	file annovar_anno_dir

	output:
	set val(experiment_name),file("${ofile}") into anno_file
	file("*_multianno.txt") into v_multianno_file

	script:
	experiment_name = "${single_experiment_map_file}".replaceFirst(/-mapping_file.tsv/, "")
	ofile = experiment_name + "_anno.txt"
	"""
	#!/bin/bash
	#!/usr/bin/env /usr/local/bin/python
	variant_annotation.py -i ${single_experiment_map_file} \
		-d ${annovar_anno_dir} \
		-b ${annovar_buildver} \
		-c ${cpus} \
		-o ./ \
		-f ${ofile} \
		-a /usr/local/annovar/table_annovar.pl
	"""
}

process database_construction {
	tag "database_construction"

	//container "biocontainers/bwa:0.7.15"

	//afterScript "find $workflow.workDir -name ${fastq1} -delete; find $workflow.workDir -name ${fastq2} -delete"

	input:
	set val(experiment_name), file(anno_f) from anno_file
	file annovar_anno_dir
	file v_multianno_file

	output:
	file("*.fasta") into final_db

	script:
	mrna_fa   = "${annovar_anno_dir}/${annovar_buildver}_${annovar_protocol}Mrna.fa"
	gene_anno = "${annovar_anno_dir}/${annovar_buildver}_${annovar_protocol}.txt"
	"""
	CPJ=`which customprodbj.jar`
	java -jar \$CPJ \
		-f ${anno_f} \
		-d ${mrna_fa} \
		-r ${gene_anno} \
		-t \
		-o ./ \
		-p2 ${experiment_name} \
		-ref ref.fasta

	"""
}

