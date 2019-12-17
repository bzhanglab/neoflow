#!/usr/bin/env nextflow

params.help = null

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => HLA tuping
    =========================================
    Usage:
    nextflow run neoflow_hla_typing.nf
    """.stripIndent()
}


genome_bam_folder = file(params.hla.bam_folder)
hla_reference = file(params.hla.hla_reference)
results_folder = file(params.hla.optitype_output)
map_file = file(params.map_file)
threads = params.threads

process convert_bam {
	tag "params.experiment"

	container "zhanglab18/samtools:1.9"

	input:
	file genome_bam_folder

	output:
	file genome_bam_folder into pair_end_fastq

	script:
	"""
	#!/bin/sh
  	sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
    	mkdir -p ${results_folder}/\${sample}
		samtools fastq -1 ${results_folder}/\${sample}/\${sample}.pair_1.fastq -2 ${results_folder}/\${sample}/\${sample}.pair_2.fastq ${genome_bam_folder}/\${sample}.bam
	done
	"""
}

process bwa_index {
	tag "params.experiment"

	container "biocontainers/bwa:0.7.15"

	input:
	file genome_bam_folder from pair_end_fastq
	file hla_reference

	output:
	file genome_bam_folder into pair_end_sam

	script:
	"""
	#!/bin/sh
	bwa index ${hla_reference}/hla_reference_dna.fasta
	sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
		bwa mem -t ${threads} -M ${hla_reference}/hla_reference_dna.fasta ${results_folder}/\${sample}/\${sample}.pair_1.fastq > ${results_folder}/\${sample}/\${sample}.pair_1.sam
		bwa mem -t ${threads} -M ${hla_reference}/hla_reference_dna.fasta ${results_folder}/\${sample}/\${sample}.pair_2.fastq > ${results_folder}/\${sample}/\${sample}.pair_2.sam
		rm ${results_folder}/\${sample}/\${sample}.pair_1.fastq ${results_folder}/\${sample}/\${sample}.pair_2.fastq
	done
	"""
}

process samtools_view{
	tag "params.experiment"

	container "zhanglab18/samtools:1.9"

	input:
	file genome_bam_folder from pair_end_sam

	output:
	file genome_bam_folder into pair_end_bam

	script:
	"""
	#!/bin/sh
	sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
		samtools view -bS ${results_folder}/\${sample}/\${sample}.pair_1.sam > ${results_folder}/\${sample}/\${sample}.pair_1.bam
		samtools view -bS ${results_folder}/\${sample}/\${sample}.pair_2.sam > ${results_folder}/\${sample}/\${sample}.pair_2.bam
		rm ${results_folder}/\${sample}/\${sample}.pair_1.sam ${results_folder}/\${sample}/\${sample}.pair_2.sam
	done
	"""
}

process samtools_fastq {
	tag "params.experiment"

	container "zhanglab18/samtools:1.9"

	input:
	file genome_bam_folder from	pair_end_bam

	output:
	file genome_bam_folder into pair_end_input_fastq

	script:
	"""
	#!/bin/sh
	sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
		samtools fastq ${results_folder}/\${sample}/\${sample}.pair_1.bam > ${results_folder}/\${sample}/\${sample}.pair_1_input.fastq
		samtools fastq ${results_folder}/\${sample}/\${sample}.pair_2.bam > ${results_folder}/\${sample}/\${sample}.pair_2_input.fastq
		rm ${results_folder}/\${sample}/\${sample}.pair_1.bam ${results_folder}/\${sample}/\${sample}.pair_2.bam
	done
	"""
}

process optitype {
	tag "params.experiment"

	container "zhanglab18/optitype:1.3.1"

	afterScript "find $workflow.workDir -name ${bam2} -delete; find $workflow.workDir -name ${bam2} -delete"

	input:
	file genome_bam_folder from pair_end_input_fastq

	output:
	
	script:
	"""
	#!/bin/sh
	sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
		python /opt/OptiType/OptiTypePipeline.py \
			--input ${results_folder}/\${sample}/\${sample}.pair_1_input.fastq ${results_folder}/\${sample}/\${sample}.pair_2_input.fastq \
			--prefix \${sample}_optitype \
			--outdir ${results_folder}/\${sample} \
			--dna
		rm ${results_folder}/\${sample}/\${sample}.pair_1_input.fastq ${results_folder}/\${sample}/\${sample}.pair_2_input.fastq
	done
	"""

}
