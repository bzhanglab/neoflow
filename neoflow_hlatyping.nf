#!/usr/bin/env nextflow

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => hla class typing by optitype
    =========================================
    Usage:
    nextflow run neoflow-db.nf
    Arguments:
      ### Input bam file or two fastq files
      --bam_file                  Input bam file
      --fq_file_r1				  Input pair end read 1 fastq file
      --fq_file_r2				  Input pair end read 2 fastq file
      --out_dir                   Output folder, default is "./"
      --hla_reference			  HLA reference file
      --cpu                       The number of CPUs

    """.stripIndent()
}


params.help       = false
params.cpu 		  = 0
params.out_dir    = "./"
params.prefix 	  = "test"
params.bam 		  = false
params.fq_file_r1 	  = false
params.fq_file_r2 	  = false

hla_reference 	  = file(params.hla_reference)
results_folder    = file(params.out_dir)
prefix 			  = params.prefix
threads 		  = params.cpu

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

if (params.bam != false) {
	bam_file = file(params.bam)

	process convert_bam {
		tag "$prefix"

		container "zhanglab18/samtools:1.9"

		input:
		file bam_file

		output:
		set file("${prefix}.pair_1.fastq"), file("${prefix}.pair_2.fastq") into pair_end_fastq

		script:
		"""
		samtools fastq -1 ${prefix}.pair_1.fastq -2 ${prefix}.pair_2.fastq ${bam_file}
		"""
	}
} else {
	fq_file_r1 = file(params.fq_file_r1)
	fq_file_r2 = file(params.fq_file_r2)
	pair_end_fastq = Channel.from(fq_file_r1, fq_file_r2)
}



process bwa_index {
	tag "$prefix"

	container "biocontainers/bwa:0.7.15"

	if (params.bam != false) {
		//afterScript "find $results_folder -name ${fastq1} -delete; find $results_folder -name ${fastq2} -delete"
	}

	input:
	set file(fastq1), file(fastq2) from pair_end_fastq
	file hla_reference

	output:
	set file("${prefix}.pair_1.sam"), file("${prefix}.pair_2.sam") into pair_end_sam

	script:
	"""

	bwa index ${hla_reference}

	bwa mem -t ${threads} -M ${hla_reference} ${fastq1} > ${prefix}.pair_1.sam
	bwa mem -t ${threads} -M ${hla_reference} ${fastq2} > ${prefix}.pair_2.sam
	"""
}

process samtools_view{
	tag "$prefix"

	container "zhanglab18/samtools:1.9"

	afterScript "find $results_folder -name ${sam1} -delete; find $results_folder -name ${sam2} -delete"

	input:
	set file(sam1), file(sam2) from pair_end_sam

	output:
	set file("${prefix}.pair_1.bam"), file("${prefix}.pair_2.bam") into pair_end_bam

	script:
	"""
	samtools view -bS ${sam1} > ${prefix}.pair_1.bam
	samtools view -bS ${sam2} > ${prefix}.pair_2.bam
	"""
}

process samtools_fastq {
	tag "$prefix"

	container "zhanglab18/samtools:1.9"

	afterScript "find $results_folder -name ${bam2} -delete; find $results_folder -name ${bam2} -delete"	

	input:
	set file(bam1), file(bam2) from	pair_end_bam

	output:
	set file("${prefix}.pair_1_input.fastq"), file("${prefix}.pair_2_input.fastq") into pair_end_input_fastq

	script:
	"""
	samtools fastq ${bam1} > ${prefix}.pair_1_input.fastq
	samtools fastq ${bam2} > ${prefix}.pair_2_input.fastq
	"""
}

process optitype {
	tag "$prefix"

	container "zhanglab18/optitype:1.3.1"

	afterScript "find $results_folder -name ${bam2} -delete; find $results_folder -name ${bam2} -delete"

	publishDir "${results_folder}", mode: "copy", overwrite: true

	input:
	set file(fastq1_input), file(fastq2_input) from pair_end_input_fastq

	output:
	file "*" into final_db

	script:
	"""
	python /opt/OptiType/OptiTypePipeline.py \
		--input ${fastq1_input} ${fastq1_input} \
		--prefix ${prefix}_optitype \
		--outdir . \
		--dna

	"""

}
