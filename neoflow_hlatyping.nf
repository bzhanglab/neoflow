#!/usr/bin/env nextflow

params.help = false
params.reads = null
params.hla_ref_dir = "./"
params.seqtype = "dna"
params.singleEnd = false
params.out_dir = "./"
params.cpu = 6



/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => HLA typing
    =========================================
    Usage:
    nextflow run neoflow_hlatyping.nf
    Arguments:
      --reads                     Reads data in fastq.gz format. For example, "*_{1,2}.fq.gz"
      --hla_ref_dir               HLA reference folder
      --seqtype                   Read type, dna or rna. Default is dna.
      --singleEnd                 Single end or not, default is false (pair end reads)
      --cpu                       The number of CPUs, default is 6.
      --out_dir                   Output folder, default is "./"
      --help                      Print help message

    """.stripIndent()
}


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

out_dir = file(params.out_dir)

if(!out_dir.isDirectory()){
	out_dir_result = out_dir.mkdirs()
	println out_dir_result ? "Create folder: $out_dir!" : "Cannot create directory: $myDir!"
}



hla_reference_dir = file(params.hla_ref_dir)

Channel.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs" + "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { input_data }


process unzip {

	echo true
	
	input:
    set val(pattern), file(reads) from input_data

    output:
    set val(pattern), "unzipped_{1,2}.fastq" into raw_reads

    script:
    if(params.singleEnd == true)
    """
    zcat ${reads[0]} > unzipped_1.fastq
    """
    else
    """
    zcat ${reads[0]} > unzipped_1.fastq
    zcat ${reads[1]} > unzipped_2.fastq
    """
}


process reads_mapping {

	echo true

	tag "reads_mapping"

	container "biocontainers/bwa:0.7.15"

	input:
	file hla_reference_dir
	set val(pattern), file(reads) from raw_reads

	output:
	set val(pattern), "mapped_{1,2}.sam" into mapped_reads

	script:
	if(params.singleEnd == true)
	"""
	bwa mem -t ${params.cpu} -M ${hla_reference_dir}/hla_reference_dna.fasta ${reads} > mapped_1.sam
	"""
	else
	"""
	bwa mem -t ${params.cpu} -M ${hla_reference_dir}/hla_reference_dna.fasta ${reads[0]} > mapped_1.sam
	bwa mem -t ${params.cpu} -M ${hla_reference_dir}/hla_reference_dna.fasta ${reads[1]} > mapped_2.sam
	"""
}

process run_samtools{

	container "zhanglab18/samtools:1.9"

	input:
	file hla_reference_dir
	set val(pattern), file(reads) from mapped_reads

	output:
	set val(pattern), "mapped_{1,2}.fastq" into fished_reads

	script:
	if(params.singleEnd == true)
	"""
	samtools view -bS mapped_1.sam > mapped_1.bam
	samtools fastq mapped_1.bam > mapped_1.fastq
	"""
	else
	"""
	samtools view -bS mapped_1.sam > mapped_1.bam
	samtools view -bS mapped_2.sam > mapped_2.bam

	samtools fastq mapped_2.bam > mapped_2.fastq
	samtools fastq mapped_2.bam > mapped_2.fastq
	"""

}


process run_optitype {
	tag "run_optitype"

	container "zhanglab18/optitype:1.3.1"

	publishDir "${out_dir}/hla_type/", mode: "copy", overwrite: true

	input:
	set val(pattern), file(reads) from fished_reads

	output:
	file("${pattern}") into res

	script:
	"""
	python /opt/OptiType/OptiTypePipeline.py \
		--input ${reads} \
		--prefix ${pattern} \
		--outdir ${pattern} \
		--${params.seqtype}
	"""


}
