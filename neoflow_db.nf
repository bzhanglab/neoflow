#!/usr/bin/env nextflow

params.help       = false 
// Annovar parameters
params.ref_ver   = "hg19"
params.protocol   = "refGene"
//params.operation  = "g"
params.ref_dir = false // annotation data folder for annovar 
params.annovar_dir = false


params.cpu = 0
params.vcf_file = null // mapping file
params.out_dir = "./output"


def get_absolute_path(f,out_file){
	myFile = file(f)

	File new_File = new File(out_file)

	allLines  = myFile.readLines()
	new_File.write(allLines[0]+"\n")

	headline = allLines[0].split("\t")
	hmap = [:]
	for(int i=0;i<headline.size();i++){
		hmap[headline[i]] = i
	}
	allLines.remove(0)
	for(line in allLines ) {
    	d = line.split("\t")
    	v_files = d[hmap["file"]].split(";");
    	v_files_full_path = []
    	for(int i=0;i<v_files.size();i++){
    		v_files_full_path.add(new File(v_files[i]).absolutePath)
    	}
    	d[hmap["file"]] = v_files_full_path.join(";")
    	new_File << d.join("\t") + "\n"
	}
}

def format_fasta_db(db, new_db){
	newDB = File(new_db)
	newDB.write("")
	String regex = /^>/
	new File(db).eachLine { line ->
    	if(line ==~ regex){
    		d = line.split(" ")
    		newDB.append(d[0]+"\n")
    	}else{
    		newDB.append(line+"\n")
    	}
	}
}

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => variant annotation and customized database construction
    =========================================
    Usage:
    nextflow run neoflow_db.nf
    Arguments:
      --vcf_file              A txt file contains VCF file(s)
      --annovar_dir           ANNOVAR folder
      --protocol              The parameter of "protocol" for ANNOVAR, default is "refGene"
      --ref_dir               ANNOVAR annotation data folder
      --ref_ver               The genome version, hg19 or hg38, default is "hg19"
      --out_dir               Output folder, default is "./output"
      --cpu                   The number of CPUs
      --help                  Print help message

    """.stripIndent()
}


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


out_dir = file(params.out_dir)
cpus = params.cpu

// parameters for varaint annotation: Annovar
annovar_buildver = params.ref_ver // 
annovar_protocol = params.protocol //
annovar_anno_dir  = file(params.ref_dir)

out_dir = file(params.out_dir)

if(!out_dir.isDirectory()){
	out_dir_result = out_dir.mkdirs()
	println out_dir_result ? "Create folder: $out_dir!" : "Cannot create directory: $myDir!"
}

// sample\texperiment\tfile\tfile_type
// a\te1\ta.vcf;b.vcf\tsomatic;germline
mapping_vcf_file = file(params.vcf_file)

//
annovar_operation = "g" 

/*
 * validate input files
 */
annovar_tool = file(params.annovar_dir + "/table_annovar.pl")
if( !annovar_tool.exists() ) exit 1, "ANNOVAR perl is invalid: ${annovar_tool}"
annovar_dir = file(params.annovar_dir)


if( !mapping_vcf_file.exists() ) exit 1, "variant mapping file is invalid: ${mapping_vcf_file}"


println "Format input mapping file: ${mapping_vcf_file}!"
new_mapping_file = params.out_dir +"/" + mapping_vcf_file.getName()
println "New mapping file: ${new_mapping_file}!"
get_absolute_path(params.vcf_file,new_mapping_file)
mapping_f = file(new_mapping_file)


process pre_processing {
	tag "${mapping_file}"

	echo true

	container "proteomics/pga:latest"

	input:
	file mapping_file from mapping_f

	output:
	file "*-mapping_file.tsv" into single_experiment_map_files

	script:

	"""
	#!/usr/bin/env /usr/local/bin/Rscript

	library(dplyr)
	library(readr)
	
	a = read.delim("${mapping_f}")
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

	container "proteomics/neoflow:latest"

	publishDir "${out_dir}/variant_annotation/", mode: "copy", overwrite: true

	input:
	file single_experiment_map_file from single_experiment_map_files.flatten()
	file annovar_dir

	output:
	set val(experiment_name),file("${ofile}") into anno_file
	file("*_multianno.txt") into v_multianno_file

	script:
	experiment_name = "${single_experiment_map_file}".replaceFirst(/-mapping_file.tsv/, "")
	ofile = experiment_name + "_anno.txt"
	//annovar_pl = 
	"""
	#!/bin/bash
	#!/usr/bin/env /usr/local/bin/python
	python $baseDir/bin/variant_annotation.py -i ${single_experiment_map_file} \
		-d ${annovar_anno_dir} \
		-b ${annovar_buildver} \
		-c ${cpus} \
		-o ./ \
		-f ${ofile} \
		-a "${annovar_dir}/table_annovar.pl"
	"""
}

process database_construction {
	tag "$experiment_name"

	container "proteomics/neoflow:latest"

	publishDir "${out_dir}/customized_database/", mode: "copy", overwrite: true, pattern: '*{.fasta,-varInfo.txt}'

	input:
	set val(experiment_name), file(anno_f) from anno_file
	file annovar_anno_dir
	file v_multianno_file

	output:
	file("*-var.fasta")
	set val(experiment_name), file("${experiment_name}_anno-var.fasta") into target_customized_db_fa
	//set val(experiment_name), file("${experiment_name}_anno-varInfo.txt") into target_customized_db_info
	file("ref.fasta")
	file("*-varInfo.txt")
	//file("${experiment_name}_anno-var_format.fasta") into target_customized_db_format_fa

	script:
	mrna_fa   = "${annovar_anno_dir}/${annovar_buildver}_${annovar_protocol}Mrna.fa"
	gene_anno = "${annovar_anno_dir}/${annovar_buildver}_${annovar_protocol}.txt"
	"""
	java -jar /opt/customprodbj.jar \
		-f ${anno_f} \
		-d ${mrna_fa} \
		-r ${gene_anno} \
		-t \
		-o ./ \
		-p2 ${experiment_name} \
		-ref ref.fasta
	"""

}

process format_db {
	tag "$experiment_name"

	container "proteomics/neoflow:latest"

	publishDir "${out_dir}/customized_database/", mode: "copy", overwrite: true

	input:
	set val(experiment_name), file(target_customized_db) from target_customized_db_fa

	output:
	set val(experiment_name), file("*_format.fasta") into format_db_set

	script:
	"""
	python $baseDir/bin/format_db.py ${target_customized_db}
	"""

}

process generate_decoy_db{

	tag "$experiment_name"

	container "proteomics/pga:latest"

	publishDir "${out_dir}/customized_database/", mode: "copy", overwrite: true

	input:
	set val(experiment_name), file(format_db) from format_db_set

	output:
	file out_target_decoy_db 

	script:
	out_target_decoy_db = "${experiment_name}_target_decoy.fasta"
	cont_db = "${baseDir}/database/contaminants.fasta"
	"""
	#!/usr/bin/env /usr/local/bin/Rscript
	library(PGA)
	buildTargetDecoyDB(db="${format_db}",
		decoyPrefix="XXX_",
		cont_file="${cont_db}",
		output="${out_target_decoy_db}")

	"""
}







