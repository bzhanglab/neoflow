#!/usr/bin/env nextflow

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => MS/MS searching
    =========================================
    Usage:
    nextflow run neoflow-msms.nf
    Mandatory arguments:
      --novel_protein_prefix      The prefix for novel proteins in database, default is "VAR"
      --decoy_protein_prefix      The prefix for decoy proteins in database, default is "###REV###"
      --db                        The protein database
      --cont                      The contaminant protein database
      --ms                        MS/MS data
      --refdb                     Reference protein database
      --para_file                 Parameter file for MS/MS searching
      --out_dir                   Output folder, default is "./"
      --prefix                    The prefix of output files
      --cpu                       The number of CPUs
      --se                        The ID number of search engines

    PepQuery arguments:
      --enzyme                    Enzyme used for protein digestion. Default is trypsin
      --c                         The number of allowed missed cleavage sites. Default is 2
      --tol                       The error window on experimental peptide mass values, default is 10
      --tolu                      The unit of --tol, ppm or Da. Default is ppm
      --fixmod                    Fixed modification. The format is like : 1,2,3. Different modification is represented by different number
      --varmod                    Variable modification. The format is the same with --fixMod;

    """.stripIndent()
}


params.help = false 
params.novel_protein_prefix =  "VAR"
params.decoy_protein_prefix = "###REV###"

params.db = false
params.cont = false
params.ms = false
params.refdb = false
params.para_file = false
params.out_dir = "./"
params.prefix = "test"
params.cpu = 2
params.se = 1

// 
params.enzyme = 1 // enzyme used for PepQuery
params.c = 2 // The number of allowed missed cleavage sites. Default is 2;
params.tol = 10
params.tolu = "ppm"
params.itol = 0.05 // da
params.varmod = false
params.fixmod = false


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


// Input parameters
db = file(params.db)
cont_db = file(params.cont)
// a folder which contains multiple MGF files or a single MGF file
ms_data = file(params.ms)
msms_searching_para_file = file(params.para_file)
out_dir = file(params.out_dir)
out_prefix = params.prefix
novel_protein_prefix = params.novel_protein_prefix
cpus = params.cpu
// search engine: 1=>MS-GF+, 2=>X!Tandem, 3=>Comet
search_engine_ID = params.se

reference_protein_database = file(params.refdb) // reference protein database

// parameters for novel peptide validation (PepQuery)
pepquery_tol = params.tol
pepquery_tolu = params.tolu // ppm or da
pepquery_enzyme = params.enzyme 
pepquery_c = params.c
pepquery_varmod = params.varmod
pepquery_fixmod = params.fixmod
pepquery_itol = params.itol






process build_target_decoy_db {
	tag "build_target_decoy_db"

	//container "zhanglab18/pga"

	input:
	file db
	file cont_db

	output:
	file "target_decoy_db.fasta" into target_decoy_db

	script:
	"""
	#!/usr/bin/env /usr/local/bin/Rscript
	library(PGA)
	target_db = "${db}"
	cont_db = "${cont_db}"
	final_db = "target_decoy_db.fasta"
	buildTargetDecoyDB(target_db,cont_file=cont_db,decoyPrefix="###REV###",output=final_db)

	"""
}


if(ms_data.isFile()){
	println "Process single MS/MS file."
	ms_data_file = file(params.ms)
}else{
	println "Process multiple MS/MS files."
	ms_data_file = Channel.fromPath("${params.ms}/*.mgf")
}

process msms_searching{
	tag "msms_searching"

	//container "zhanglab18/neoflow:1.0.0"

	input:
	file ms_file from ms_data_file
	file para_file from msms_searching_para_file
	file target_decoy_db

	output:
	file "${ms_file.baseName}.mzid" into psm_raw_file

	script:
	if (search_engine_ID == 1){
		// MS-GF+
		"""
		java -jar /usr/local/msgf/MSGFPlus.jar \
			-conf ${para_file} \
			-s ${ms_file} \
			-d ${target_decoy_db} \
			-tda 0 \
			-o ${ms_file.baseName}.mzid

		"""
	}else if(search_engine_ID == 2){
		// X!Tandem
		"""
		tandem.exe ${para_file}
		## convert xml to mzid
		java -jar /usr/local/bin/mzidentml-lib-1.6.12-SNAPSHOT.jar Tandem2mzid  \
			tandem_result.xml \
			${ms_file.baseName}.mzid \
			-outputFragmentation false \
			-decoyRegex "${params.decoy_protein_prefix}" \
			-databaseFileFormatID "MS:1001348" \
			-massSpecFileFormatID MS:1001062 \
			-idsStartAtZero false \
			-compress false \
			-proteinCodeRegex "\\S+"
	
		"""
	}else if(search_engine_ID == 3){
		// Comet
		"""
		comet.exe -p ${para_file}
		## convert comet to tsv
		comet2tsv x.csv ${ms_file.baseName}.tsv
		"""
	}
}


process post_process{
	tag "post_process"

	//container "zhanglab18/neoflow:1.0.0"

	input:
	file psm_raw_file
	file target_decoy_db

	output:
	file "./sample-rawPSMs.txt" into psm_raw_tsv

	script:
	if (search_engine_ID == 1){
		// MS-GF+
		"""
		#!/usr/bin/env /usr/local/bin/Rscript
		library(PGA)
		parserGear(file="${psm_raw_file}",
			db="${target_decoy_db}",
			decoyPrefix="${params.decoy_protein_prefix}",
			novelPrefix="${params.novel_protein_prefix}",
			prefix="sample",
			xmx=60,
			thread=8,
			alignment=0,
			outdir="./")
		raw_psm <- read.delim("./sample-rawPSMs.txt")
		colnames(raw_psm)[2] <- "score"
		write.table(raw_psm,"./sample-rawPSMs.txt", row.names = FALSE, quote = FALSE, sep = "\t")

		"""
	}else if(search_engine_ID == 2){
		// X!Tandem
		"""
		echo "PASS"
		"""
	}else if(search_engine_ID == 3){
		// Comet
		"""
		echo "PASS"
		"""
	}
}




combine_psm_raw = psm_raw_tsv.collectFile(name:"combine_psm_raw.tsv",keepHeader:true)


process fdr_estimation {
	tag "fdr_estimation"

	//container "zhanglab18/neoflow:1.0.0"

	//afterScript "find $workflow.workDir -name ${fastq1} -delete; find $workflow.workDir -name ${fastq2} -delete"

	input:
	file combine_psm_raw
	file target_decoy_db

	output:
	file "*pga-peptideSummary.txt" into psm_filtered
	file "novel_peptide.txt" into novel_peptide_file

	script:
	"""
	#!/usr/bin/env /usr/local/bin/Rscript
	library(PGA)
	psm_raw_file = "${combine_psm_raw}" 
	db = "${target_decoy_db}"
	decoyPrefix = "${params.decoy_protein_prefix}"
	novelPrefix = "${params.novel_protein_prefix}"
	calculateFDR(psmfile=psm_raw_file,
			db=db,
			peptide_level=TRUE,
			fdr=0.01,
			decoyPrefix=decoyPrefix,
			novelPrefix=novelPrefix,
			better_score_lower=FALSE,
			remap=FALSE,
			out_dir="./",
			protein_inference=FALSE,
			xmx=4)
	psm_data = read.delim("./pga-peptideSummary.txt",stringsAsFactors=FALSE)
	library(dplyr)
	psm_qvalue = psm_data %>% filter(Qvalue <= 0.01,isSAP == "true", isdecoy == "false") %>% 
		select(peptide) %>% 
		distinct()
	write.table(psm_qvalue,file="novel_peptide.txt",sep="\t",col.names=FALSE,quote=FALSE,row.names=FALSE)

	"""
}

process novel_peptide_validation {
	tag "novel_peptide_validation"

	//container "zhanglab18/neoflow:1.0.0"

	input:
	file novel_peptide_file
	file ms_data_file

	output:
	file "validation/*" into validation_result

	script:
	"""
	java -jar /usr/local/bin/pepquery-1.2-jar-with-dependencies.jar \
		-c ${pepquery_c} \
		-cpu 0 \
		-e ${pepquery_enzyme} \
		-fixMod ${pepquery_fixmod} \
		-fragmentMethod 1 \
		-itol ${pepquery_itol} \
		-m 1 \
		-minScore 12 \
		-n 10000 \
		-tol ${pepquery_tol} \
		-um \
		-varMod ${pepquery_varmod} \
		-ms ${ms_data_file} \
		-pep ${novel_peptide_file} \
		-o validation/ \
		-db ${reference_protein_database}
	"""
}


