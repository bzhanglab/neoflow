#!/usr/bin/env nextflow

params.help = null

println """\

N E O F L O W
=========================================================================================================
Start      : $workflow.start

BWA        : BWA-0.7.17
Samtools   : Samtools-1.7
Picard     : Picard-2.3.0-0
GATK       : GenomeAnalysisTK-3.8-0

launchDir  : $workflow.launchDir
WorkDir    : $workflow.workDir
ResultsDir : ${params.results}/${params.sampleID}
"""
.stripIndent()


/* Parse the input parameters */

genome_fasta           = file(params.genome_fasta) 
genome_fasta_sa        = file(params.genome_fasta + ".sa")
genome_fasta_fai       = file(params.genome_fasta + ".fai")
genome_fasta_bwt       = file(params.genome_fasta + ".bwt")
genome_fasta_ann       = file(params.genome_fasta + ".ann")
genome_fasta_amb       = file(params.genome_fasta + ".amb")
genome_fasta_pac       = file(params.genome_fasta + ".pac")
genome_fasta_dict      = file("${genome_fasta.parent}/ucsc.hg19.dict")
cosmic_file            = file(params.cosmic_file)
mills_indel_file       = file(params.mills_indel_file)
known_indel_file       = file(params.known_indel_file)
known_dbsnps_file      = file(params.known_dbsnps_file)
known_dbsnps_1000_file = file(params.known_dbsnps_1000_file)
hla_reference_dna      = file(params.hla_reference_dna)
refseq_human_hg19_gtf  = file(params.refseq_human_hg19_gtf)
modification_file_txt  = file(params.modification_file_txt)
gene_annotation_txt    = file(params.annovar + "/humandb/hg19_refGene.txt")
mrna_database_fa       = file(params.annovar + "/humandb/hg19_refGeneMrna.fa")


/* Load DNA-Seq read pairs */

if (params.paired_dna_seq) {

	normal_ch_1 = Channel.fromPath("${params.normal_wxs_read_1}")
					.ifEmpty{error "Cannot find any reads matching: ${params.normal_wxs_read_1}"}
					.map{file -> tuple("normal", file)}
	
	normal_ch_2 = Channel.fromPath("${params.normal_wxs_read_2}")
					.ifEmpty{error "Cannot find any reads matching: ${params.normal_wxs_read_2}"}
					.map{file -> tuple("normal", file)}

	tumor_ch_1 = Channel.fromPath("${params.tumor_wxs_read_1}")
					.ifEmpty{error "Cannot find any reads matching: ${params.tumor_wxs_read_1}"}
					.map{file -> tuple("tumor", file)}
	
	tumor_ch_2 = Channel.fromPath("${params.tumor_wxs_read_2}")
					.ifEmpty{error "Cannot find any reads matching: ${params.tumor_wxs_read_2}"}
					.map{file -> tuple("tumor", file)}

	normal_ch_1.concat(normal_ch_2, tumor_ch_1, tumor_ch_2)
		.groupTuple()
		.into{grouped_raw_reads_wxs_ch1; grouped_raw_reads_wxs_ch2}

} else {

	normal_ch_1 = Channel.fromPath("${params.normal_wxs_read_1}")
					.ifEmpty{error "Cannot find any reads matching: ${params.normal_wxs_read_1}"}
					.map{file -> tuple("normal", file)}

	tumor_ch_1 = Channel.fromPath("${params.tumor_wxs_read_1}")
					.ifEmpty{error "Cannot find any reads matching: ${params.tumor_wxs_read_1}"}
					.map{file -> tuple("tumor", file)}

	normal_ch_1.concat(tumor_ch_1)
		.groupTuple()
		.into{grouped_raw_reads_wxs_ch1; grouped_raw_reads_wxs_ch2}

}


/* Load RNASeq read pairs */

if (params.paired_rna_seq) {

	rna_seq_ch_1 = Channel.fromPath("${params.rna_seq_read_1}")
					.ifEmpty{error "Cannot find any reads matching: ${params.rna_seq_read_1}"}
					.map{file -> tuple("rna-seq", file)}

	rna_seq_ch_2 = Channel.fromPath("${params.rna_seq_read_2}")
					.ifEmpty{error "Cannot find any reads matching: ${params.rna_seq_read_2}"}
					.map{file -> tuple("rna-seq", file)}

	rna_seq_ch_1.concat(rna_seq_ch_2)
		.groupTuple()
		.into{grouped_raw_reads_rna_seq_ch1; grouped_raw_reads_rna_seq_ch2; grouped_raw_reads_rna_seq_ch3}

} else {

	Channel
		.fromPath("${params.rna_seq_read_1}")
		.ifEmpty{error "Cannot find any reads matching: ${params.rna_seq_read_1}"}
		.map{file -> tuple("rna-seq", file)}
		.into{grouped_raw_reads_rna_seq_ch1; grouped_raw_reads_rna_seq_ch2; grouped_raw_reads_rna_seq_ch3}

}


/* 
 * Read alignment via Burrows-Wheeler Aligner - MEM algorithm and
 * use of Samtools to convert the SAM output from BWA to BAM format  
 */

 if(params.paired_dna_seq) {

	/* Runs for paired-end reads when paired_dna_seq == true */

	process bwa_mem { 
		
		tag "$sample"

		container "biocontainers/bwa:0.7.15"

		/*afterScript "find $workflow.workDir -name ${reads[0]} -delete; find $workflow.workDir -name ${reads[1]} -delete"*/

		input:
		file genome_fasta     
		file genome_fasta_sa  
		file genome_fasta_fai 
		file genome_fasta_bwt 
		file genome_fasta_ann 
		file genome_fasta_amb 
		file genome_fasta_pac
		
		set sample, file(reads) from grouped_raw_reads_wxs_ch1

		output: 
		set sample, file("${params.sampleID}.dna.${sample}.sam") into bwa_aligned_sam_wxs_ch

		script:
		"""
		bwa mem -t ${params.threads.mapping} -M ${genome_fasta} ${reads[0]} ${reads[1]} > "${params.sampleID}.dna.${sample}.sam"
		"""
	}
}

else {

	/* Runs for paired-end reads when paired_dna_seq == false */

	process bwa_mem { 
		
		tag "$sample"

		container "biocontainers/bwa:0.7.15"

		/*afterScript "find $workflow.workDir -name ${reads} -delete"*/

		input:
		file genome_fasta     
		file genome_fasta_sa  
		file genome_fasta_fai 
		file genome_fasta_bwt 
		file genome_fasta_ann 
		file genome_fasta_amb 
		file genome_fasta_pac
		
		set sample, file(reads) from grouped_raw_reads_wxs_ch1

		output: 
		set sample, file("${params.sampleID}.dna.${sample}.sam") into bwa_aligned_sam_wxs_ch

		script:
		"""
		bwa mem -t ${params.threads.mapping} -M ${genome_fasta} ${reads} > "${params.sampleID}.dna.${sample}.sam"
		"""
	}
}


/* Prints all alignments in the specified input alignment file into a BAM file */

process samtools_view {

	tag "$sample"

	container "biocontainers/samtools:1.3.1"

	/* afterScript "find $workflow.workDir -name ${sam_file_raw} -delete" */

	input:
	set sample, file(sam_file_raw) from bwa_aligned_sam_wxs_ch

	output: 
	set sample, file("${sam_file_raw.baseName}.aligned.bam") into bwa_aligned_bam_wxs_ch

	script:
	"""
	samtools view -bS ${sam_file_raw} > "${sam_file_raw.baseName}.aligned.bam"
	"""
}


/* Use Picard tools to add group to raw bam file */

process picard_add_or_replace_read_groups {
	
	tag "$sample"

	container "biocontainers/picard:2.3.0"

	/*afterScript "find $workflow.workDir -name ${bam_file_raw} -delete"*/
	
	input:
	set sample, file(bam_file_raw) from bwa_aligned_bam_wxs_ch
	
	output:
	set sample, file("${bam_file_raw.baseName}.grouped.sorted.bam") into picard_added_group_bam_wxs_ch1
	set sample, file("${bam_file_raw.baseName}.grouped.sorted.bam") into picard_added_group_bam_wxs_ch2
	
	script:
	"""
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar AddOrReplaceReadGroups \
		I=${bam_file_raw} \
		O="${bam_file_raw.baseName}.grouped.sorted.bam" \
		SO=coordinate \
		RGLB=${params.sampleID}.${sample} \
		RGPL=illumina \
		RGPU=${params.sampleID}.${sample} \
		RGSM=${params.sampleID}.${sample}
	"""
}


/* Create summary mapping using Samtools */

process samtools_flagstat {
	
	tag "$sample"

	container "biocontainers/samtools:1.3.1"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	publishDir "${params.results}/${params.sampleID}", mode: "move", overwrite: true
	
	input:
	set sample, file(bam_file_added_group) from picard_added_group_bam_wxs_ch1
	
	output:
	set sample, file("*.mapping.statistic.txt") into samtools_flagstat_wxs_ch
	
	script:
	"""	
	samtools flagstat ${bam_file_added_group} > "${bam_file_added_group.baseName}.mapping.statistic.txt"
	"""
}


/* Use Picard tools to mark duplicates */

process picard_mark_duplicates {
	
	tag "$sample"

	container "biocontainers/picard:2.3.0"

	/*afterScript "find $workflow.workDir -name ${bam_file_added_group} -delete"*/
	
	input:
	set sample, file(bam_file_added_group) from picard_added_group_bam_wxs_ch2
	
	output:
	set sample, file("${bam_file_added_group.baseName}.deduplicated.bam"), file("${bam_file_added_group.baseName}.deduplicated.bai") into picard_mark_duplicates_wxs_ch1
	set sample, file("${bam_file_added_group.baseName}.deduplicated.bam"), file("${bam_file_added_group.baseName}.deduplicated.bai") into picard_mark_duplicates_wxs_ch2
	
	script:
	"""	
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar MarkDuplicates \
		I=${bam_file_added_group} \
		O="${bam_file_added_group.baseName}.deduplicated.bam" \
		METRICS_FILE="${bam_file_added_group.baseName}.output.metrics" \
		VALIDATION_STRINGENCY=SILENT \
		CREATE_INDEX=true
	"""
}


/* GATK realigner target creator */

process gatk_realigner_target_creator {

	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	set sample, file(bam_file_marked_duplicate), file(bai_file_marked_duplicate) from picard_mark_duplicates_wxs_ch1
	
	output:
	set sample, file("${bam_file_marked_duplicate.baseName}.intervals.list") into gatk_target_creator_wxs_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-nt ${params.threads.target_creator} \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-o "${bam_file_marked_duplicate.baseName}.intervals.list" 		
	"""
}	


/* GATK perform local realignment of reads around indels */

process gatk_indel_realignment {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${interval_list} -delete; find $workflow.workDir -name ${bam_file_marked_duplicate} -delete; find $workflow.workDir -name ${bai_file_marked_duplicate} -delete" */
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	set sample, file(interval_list) from gatk_target_creator_wxs_ch
	set sample, file(bam_file_marked_duplicate), file(bai_file_marked_duplicate) from picard_mark_duplicates_wxs_ch2
	
	output:
	set sample, file("${bam_file_marked_duplicate.baseName}.realigned.bam"), file("${bam_file_marked_duplicate.baseName}.realigned.bai") into gatk_indel_realignment_wxs_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-targetIntervals ${interval_list} \
		-o "${bam_file_marked_duplicate.baseName}.realigned.bam"
	"""
}


/* Use GATK to detect systematic errors in base quality scores */

process gatk_base_recalibrator {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file known_indel_file
	file mills_indel_file
	file known_dbsnps_file
	file known_dbsnps_1000_file	
	
	set sample, file(indel_realigned_bam), file(indel_realigned_bai) from gatk_indel_realignment_wxs_ch
	
	output:
	set sample, file(indel_realigned_bam), file(indel_realigned_bai), file("${indel_realigned_bam.baseName}.data.table") into gatk_base_recalibrator_wxs_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-nct ${params.threads.base_recalibrator} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		--knownSites ${known_dbsnps_1000_file} \
		--knownSites ${known_dbsnps_file} \
		--knownSites ${mills_indel_file} \
		--knownSites ${known_indel_file} \
		-o "${indel_realigned_bam.baseName}.data.table"
	"""
}


/* Use GATK to write out sequence read data (for filtering, merging, subsetting etc) */

process gatk_print_reads {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${indel_realigned_bam} -delete; find $workflow.workDir -name ${indel_realigned_bai} -delete; find $workflow.workDir -name ${recalibration_table} -delete" */
		
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	set sample, file(indel_realigned_bam), file(indel_realigned_bai), file(recalibration_table) from gatk_base_recalibrator_wxs_ch
	
	output:
	file("${indel_realigned_bam.baseName}.processed.bam") into processed_bam_wxs_ch
	file("${indel_realigned_bam.baseName}.processed.bai") into processed_bai_wxs_ch

	file("${indel_realigned_bam.baseName}.processed.bam") into processed_bam_wxs_ch2
	file("${indel_realigned_bam.baseName}.processed.bai") into processed_bai_wxs_ch2
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T PrintReads \
		-nct ${params.threads.print_reads} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		-BQSR ${recalibration_table} \
		-o "${indel_realigned_bam.baseName}.processed.bam"
	"""
}


Channel
		.from("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY")
		.set{chromosome_wxs_ch}


/* Call somatic SNPs and indels via local re-assembly of haplotypes */
	
process mutect_2 {
	
	tag "$chromosome"

	container "broadinstitute/gatk3:3.8-0"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file cosmic_file
	file known_dbsnps_file
	
	file("*") from processed_bam_wxs_ch.toList()
	file("*") from processed_bai_wxs_ch.toList()
	
	val chromosome from chromosome_wxs_ch
	
	output:
	file("${params.sampleID}.${chromosome}.vcf") into mutect2_vcf_wxs_ch
		
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T MuTect2 \
		-R ${genome_fasta} \
		-I:tumor  ${params.sampleID}.dna.tumor.aligned.grouped.sorted.deduplicated.realigned.processed.bam \
		-I:normal ${params.sampleID}.dna.normal.aligned.grouped.sorted.deduplicated.realigned.processed.bam \
		-L ${chromosome} \
		--dbsnp ${known_dbsnps_file} \
		--cosmic ${cosmic_file} \
		-o "${params.sampleID}.${chromosome}.vcf" \
		-nt ${params.threads.mutect2}	
	"""	
}


/* Concatenate VCF files of non-overlapping genome intervals, all with the same set of samples */

process gatk_cat_variants {

	tag "wxs"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name '${params.sampleID}.chr*.vcf' -delete"*/
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict

	file("*") from mutect2_vcf_wxs_ch.toList()
	
	output:
	file("${params.sampleID}.mutect2.vcf") into mutect2_final_vcf_wxs_ch
	
	script:
	"""	
	java -Xmx6g -cp /usr/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
		-R ${genome_fasta} \
		-V ${params.sampleID}.chr1.vcf  \
		-V ${params.sampleID}.chr2.vcf  \
		-V ${params.sampleID}.chr3.vcf  \
		-V ${params.sampleID}.chr4.vcf  \
		-V ${params.sampleID}.chr5.vcf  \
		-V ${params.sampleID}.chr6.vcf  \
		-V ${params.sampleID}.chr7.vcf  \
		-V ${params.sampleID}.chr8.vcf  \
		-V ${params.sampleID}.chr9.vcf  \
		-V ${params.sampleID}.chr10.vcf \
		-V ${params.sampleID}.chr11.vcf \
		-V ${params.sampleID}.chr12.vcf \
		-V ${params.sampleID}.chr13.vcf \
		-V ${params.sampleID}.chr14.vcf \
		-V ${params.sampleID}.chr15.vcf \
		-V ${params.sampleID}.chr16.vcf \
		-V ${params.sampleID}.chr17.vcf \
		-V ${params.sampleID}.chr18.vcf \
		-V ${params.sampleID}.chr19.vcf \
		-V ${params.sampleID}.chr20.vcf \
		-V ${params.sampleID}.chr21.vcf \
		-V ${params.sampleID}.chr22.vcf \
		-V ${params.sampleID}.chrX.vcf  \
		-V ${params.sampleID}.chrY.vcf  \
		-V ${params.sampleID}.chrM.vcf  \
		-out "${params.sampleID}.mutect2.vcf" \
		-assumeSorted
	"""
}


/* Sort VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. */

process picard_sort_vcf {

	tag "wxs"

	container "biocontainers/picard:2.3.0"

	/*afterScript "find $workflow.workDir -name ${original_vcf} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file(original_vcf) from mutect2_final_vcf_wxs_ch
	
	output:
	file("${original_vcf.baseName}.sorted.vcf") into picard_sort_vcf_wxs_ch

	script:
	"""
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar SortVcf \
		I=${original_vcf} \
		O="${original_vcf.baseName}.sorted.vcf" \
		SEQUENCE_DICTIONARY=${genome_fasta_dict}
	"""
}


/* Filter variant calls based on INFO and FORMAT annotations. */

process gatk_variant_filtration {

	tag "wxs"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${sorted_vcf} -delete"*/

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict

	file(sorted_vcf) from picard_sort_vcf_wxs_ch

	output:
	file("${sorted_vcf.baseName}.filtered.vcf") into gatk_variant_filtration_vcf_wxs_ch

	script:
	"""
	java -jar /usr/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R ${genome_fasta} \
		-V ${sorted_vcf} \
		-window 35 \
		-cluster 3 \
		-filterName FS \
		-filter "FS > 30.0" \
		-filterName QD \
		-filter "QD < 2.0" \
		-o "${sorted_vcf.baseName}.filtered.vcf"
	"""
}


/* Removes all sites with a FILTER flag other than PASS. */

process vcftools_remove_filtered {

	tag "wxs"

	container "zhanglab18/vcftools:0.1.16"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	/*afterScript "find $workflow.workDir -name ${sorted_vcf} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file(sorted_vcf) from gatk_variant_filtration_vcf_wxs_ch

	output:
	file("${sorted_vcf.baseName}.recode.vcf") into vcftools_filtered_vcf_wxs_ch

	script:
	"""
	vcftools \
		--vcf ${sorted_vcf} \
		--remove-filtered-all \
		--recode \
		--out "${sorted_vcf.baseName}"
	"""
}


/* Functionally annotate genetic variants detected from diverse genomes using ANNOVAR. */

process annovar {

	tag "wxs"

	/*afterScript "find $workflow.workDir -name ${mutect2_vcf} -delete"*/

	input:
	file(mutect2_vcf) from vcftools_filtered_vcf_wxs_ch

	output:
	file("${params.sampleID}.wxs.hg19_multianno.txt") into annovar_wxs_ch

	script:
	"""
	perl ${params.annovar}/table_annovar.pl ${mutect2_vcf} ${params.annovar}/humandb \
		-buildver hg19 \
		-out "${params.sampleID}.wxs" \
		-protocol refGene \
		-operation g \
		-nastring . \
		-vcfinput \
		--thread ${params.threads.annovar} \
		--maxgenethread ${params.threads.annovar} \
		-polish \
		-otherinfo
	"""
}


/* Scan microsatellites from reference genome. */

process msisensor_scan {

	tag "msi"

	container "zhanglab18/msisensor:latest"

	input:
	file genome_fasta

	output:
	file("${params.sampleID}.microsatellites.list") into msisensor_scan_ch

	script:
	"""
	/opt/msisensor/binary/msisensor.linux scan \
		-d ${genome_fasta} \
		-o "${params.sampleID}.microsatellites.list"
	"""
}


/* MSI scoring. */

process msisensor_score {

	tag "msi"

	container "zhanglab18/msisensor:latest"

	input:
	file("*") from processed_bam_wxs_ch2.toList()
	file("*") from processed_bai_wxs_ch2.toList()

	file(microsatelittes_list) from msisensor_scan_ch

	output:
	set file("${params.sampleID}.msi.prefix_somatic"), file("${params.sampleID}.msi.prefix_dis") into msisensor_score_ch

	script:
	"""
	/opt/msisensor/binary/msisensor.linux msi \
		-d ${microsatelittes_list} \
		-n ${params.sampleID}.dna.normal.aligned.grouped.sorted.deduplicated.realigned.processed.bam \
		-t ${params.sampleID}.dna.tumor.aligned.grouped.sorted.deduplicated.realigned.processed.bam \
		-o ${params.sampleID}.msi.prefix \
		-b 40
	"""
}


/* Extract somatic MSI INDEL data. */

process convert_msi_for_annovar {

	tag "msi"

	container "zhanglab18/msisensor:latest"

	input:
	set file(prefix_somatic), file(prefix_dis_tab) from msisensor_score_ch

	output:
	file("${params.sampleID}.msi.annovar.input.txt") into msisensor_to_annovar_ch

	script:
	"""
	perl /opt/convert_msi_for_annovar.pl ${prefix_somatic} ${prefix_dis_tab} ${params.sampleID}.msi.annovar.input.txt
	"""
}


/* Functionally annotate genetic variants detected from diverse genomes using ANNOVAR. */

process annovar {

	tag "msi"

	input:
	file(msi_to_annovar) from msisensor_to_annovar_ch

	output:
	file("${params.sampleID}.msi.hg19_multianno.txt") into annovar_msi_ch

	script:
	"""
	perl ${params.annovar}/table_annovar.pl ${msi_to_annovar} ${params.annovar}/humandb \
		-buildver hg19 \
		-out "${params.sampleID}.msi" \
		-protocol refGene \
		-operation g \
		-nastring . \
		--thread ${params.threads.annovar} \
		--maxgenethread ${params.threads.annovar} \
		-polish \
		-otherinfo
	"""
}


/* Combine results from multiple MS/MS search engines with Customprodbj and build customized protein database. */

process protein_db_construction {

	tag "wxs"

	container "zhanglab18/customprodbj:1.1.0"

	/*afterScript "find $workflow.workDir -name ${hg19_multianno_msi} -delete; find $workflow.workDir -name ${hg19_multianno_wxs} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file mrna_database_fa
	file gene_annotation_txt

	file(hg19_multianno_wxs) from annovar_wxs_ch
	file(hg19_multianno_msi) from annovar_msi_ch

	output:
	file("${params.sampleID}.msgfplus.txt") into protein_db_construction_wxs_chx

	file("${params.sampleID}.merge-var.fasta") into protein_db_construction_wxs_ch1
	file("${params.sampleID}.merge-var.fasta") into protein_db_construction_wxs_ch2
	file("${params.sampleID}.pro-ref.fasta") into protein_db_reference
	
	set file("${params.sampleID}.merge-var.fasta"), file("${params.sampleID}.merge-varInfo.txt") into protein_db_construction_wxs_ch3

	script:
	"""
	python /opt/concat.py --file-one ${hg19_multianno_wxs} --file-two ${hg19_multianno_msi} --prefix ${params.sampleID}

	java -jar /opt/customprodbj.jar \
		-i ${params.sampleID}.msgfplus.txt \
		-d ${mrna_database_fa} \
		-r ${gene_annotation_txt} \
		-o . \
		-ref ${params.sampleID}.pro-ref.fasta \
		-t

	mv merge-var.fasta ${params.sampleID}.merge-var.fasta
	mv merge-varInfo.txt ${params.sampleID}.merge-varInfo.txt
	"""
}


/* 
 * Use MSGF+ to perform peptide identification by scoring MS/MS 
 * spectra against peptides derived from a protein sequence database. 
 */

Channel.fromPath("${params.mzml_files_dir}/*.mzML")
	.into{mzml_files_ch1; mzml_files_ch2}

process msgfplus {

	tag "proteomics"

	container "zhanglab18/msgfplus:2018.07.17"

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file modification_file_txt

	file(database_file_fasta) from protein_db_construction_wxs_ch1
	file(mzml_files) from mzml_files_ch1.collect()

	output:
	file("*.mzid") into msgfplus_mzid_wxs_ch
	file("${params.sampleID}.msgfplus.tsv") into msgfplus_wxs_ch

	script:
	"""
	for mzml in ${mzml_files}
	do
		basename=`basename \$mzml .mzML`

		java -Xmx6g -jar /opt/MSGFPlus.jar \
			-thread ${params.threads.msgfplus} \
			-s \$mzml \
			-d ${database_file_fasta} \
			-mod ${modification_file_txt} \
			-o \$basename.mzid \
			-m ${params.msgfplus.m} \
			-maxLength ${params.msgfplus.maxLength} \
			-n ${params.msgfplus.n} \
			-ti ${params.msgfplus.ti} \
			-minLength ${params.msgfplus.minLength} \
			-addFeatures ${params.msgfplus.addFeatures} \
			-t ${params.msgfplus.t} \
			-minCharge ${params.msgfplus.minCharge} \
			-maxCharge ${params.msgfplus.maxCharge} \
			-tda ${params.msgfplus.tda} \
			-protocol ${params.msgfplus.protocol} \
			-ntt ${params.msgfplus.ntt} \
			-inst ${params.msgfplus.inst} \
			-e ${params.msgfplus.e} \
			-maxMissedCleavages ${params.msgfplus.maxMissedCleavages}

		java -Xmx3500M -cp /opt/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv \
			-i \$basename.mzid \
			-o \$basename.msgfplus.tsv
	done

	python /opt/concat.py --path . --prefix ${params.sampleID}
	"""
}


/* 
 * Use moFF to extract apex MS1 intensity using a set of identified MS2 peptides. Uses a Go library to directly 
 * extract data from Thermo Raw spectrum files, eliminating the need for conversions from other formats.
 */

process deepflq {

	tag "wxs"

	container "zhanglab18/deeplfq:1.0.0"

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file("*") from mzml_files_ch2.collect()
	file("*") from msgfplus_mzid_wxs_ch.collect()
	file(database_file_fasta) from protein_db_construction_wxs_ch2

	output:
	file("${params.sampleID}.ibaq.txt") into deepflq_wxs_ch

	script:
	"""
	java -jar /opt/DeepLFQ-1.0.0.jar \
		-d ${database_file_fasta} \
		-i . \
		-s . \
		-o . \
		-moff /opt/moFF-1.2 \
		-py /usr/local/bin/python \
		-m 2 \
		-tol 10 \
		-rt_p 1 \
		-rt_w 3 

	mv iBAQ.txt ${params.sampleID}.ibaq.txt
	"""
}


/* 
 * Use BWA to index the HLA reference file that comes with the optitype 
 * and map the reads against this reference sequence.
 */

process bwa_index {

	tag "wxs"

	container "biocontainers/bwa:0.7.15"

	input:
	file hla_reference_dna
	set sample, file(reads) from grouped_raw_reads_wxs_ch2

	output:
	set file("${params.sampleID}.optitype.pair_1.sam"), file("${params.sampleID}.optitype.pair_2.sam") into bwa_index_and_mem_sam_wxs_ch

	when:
	"$sample" == "tumor"

	script:
	"""
	bwa index ${hla_reference_dna}

	bwa mem -t ${params.threads.mapping} -M ${hla_reference_dna} ${reads[0]} > "${params.sampleID}.optitype.pair_1.sam"
	bwa mem -t ${params.threads.mapping} -M ${hla_reference_dna} ${reads[1]} > "${params.sampleID}.optitype.pair_2.sam"
	"""
}


/* Prints all alignments in the specified input alignment file into a BAM file */

process samtools_view {

	tag "wxs"

	container "biocontainers/samtools:1.3.1"

	/*afterScript "find $workflow.workDir -name ${sam_file_raw_1} -delete; find $workflow.workDir -name ${sam_file_raw_2} -delete"*/

	input:
	set file(sam_file_raw_1), file(sam_file_raw_2) from bwa_index_and_mem_sam_wxs_ch

	output: 
	set file("${sam_file_raw_1.baseName}.bam"), file("${sam_file_raw_2.baseName}.bam") into optitype_bam_wxs_ch

	script:
	"""
	samtools view -bS ${sam_file_raw_1} > "${sam_file_raw_1.baseName}.bam"
	samtools view -bS ${sam_file_raw_2} > "${sam_file_raw_2.baseName}.bam"
	"""
}


process samtools_fastq  {

	tag "wxs"

	container "biocontainers/samtools:1.3.1"

	/*afterScript "find $workflow.workDir -name ${bam_1} -delete; find $workflow.workDir -name ${bam_2} -delete"*/

	input:
	set file(bam_1), file(bam_2) from optitype_bam_wxs_ch

	output:
	set file("${bam_1.baseName}.fastq"), file("${bam_2.baseName}.fastq") into samtools_filtered_wxs_ch

	script:
	"""
	samtools fastq ${bam_1} > ${bam_1.baseName}.fastq
	samtools fastq ${bam_2} > ${bam_2.baseName}.fastq
	"""
}


/* HLA genotyping using OptiType. */

process optitype {

	tag "wxs"

	container "zhanglab18/optitype:1.3.1"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	/*afterScript "find $workflow.workDir -name ${filtered_fastq_1} -delete; find $workflow.workDir -name ${filtered_fastq_2} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	set file(filtered_fastq_1), file(filtered_fastq_2) from samtools_filtered_wxs_ch

	output:
	set file("*.tsv"), file("*.pdf") into optitype_final_wxs_ch

	script:
	"""
	python /opt/OptiType/OptiTypePipeline.py \
		--input ${filtered_fastq_1} ${filtered_fastq_2} \
		--prefix ${params.sampleID}_optitype \
		--outdir . \
		--dna 
	"""
}


/* HLA-peptide binding prediction and matrix construction using netMHCpan. */

process peptide_binding_prediction {

	tag "wxs"

	// afterScript "find $workflow.workDir -name ${optitype_result_tsv} -delete; find $workflow.workDir -name ${optitype_coverage_plot} -delete"

	input:
	set file(database_file_fasta), file(var_info_txt) from protein_db_construction_wxs_ch3
	set file(optitype_result_tsv), file(optitype_coverage_plot) from optitype_final_wxs_ch

	output:
	file("${params.sampleID}_somatic_mutation_binding_result.csv") into peptide_binding_prediction_wxs_ch

	script:
	"""
	python3 ${baseDir}/bin/binding_prediction.py \
		-id ${params.sampleID} \
		-a ${optitype_result_tsv} \
		-fa ${database_file_fasta} \
		-txt ${var_info_txt} \
		-o . \
		-netmhcpan ${params.netmhcpan} \
		-nt ${params.threads.mapping}
	"""
}

/* 
 * Generate STAR genome index files 
 */

 if(params.paired_rna_seq) {

	/* Runs for paired-end reads when paired_rna_seq == true */

	process star_first_pass {

		tag "$sample"

		container "quay.io/biocontainers/star:2.6.0b--0"

		input:
		file genome_fasta

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch1

		output:
		file("${params.sampleID}.SJ.out.tab") into star_first_pass_rna_seq_ch

		script:
		"""
		STAR --runMode genomeGenerate \
			--genomeDir . \
			--genomeFastaFiles ${genome_fasta} \
			--runThreadN ${params.threads.mapping} 

		STAR --genomeDir . \
			--runThreadN ${params.threads.mapping} \
			--outSAMtype BAM Unsorted \
			--readFilesCommand zcat \
			--readFilesIn ${reads[0]} ${reads[1]} \
			--outFileNamePrefix "${params.sampleID}."
		"""
	}
}

else {

	/* Runs for single-end reads when paired_rna_seq == false */

	process star_first_pass {

		tag "$sample"

		container "quay.io/biocontainers/star:2.6.0b--0"

		input:
		file genome_fasta

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch1

		output:
		file("${params.sampleID}.SJ.out.tab") into star_first_pass_rna_seq_ch


		script:
		"""
		STAR --runMode genomeGenerate \
			--genomeDir . \
			--genomeFastaFiles ${genome_fasta} \
			--runThreadN ${params.threads.mapping} 

		STAR --genomeDir . \
			--runThreadN ${params.threads.mapping} \
			--outSAMtype BAM Unsorted \
			--readFilesCommand zcat \
			--readFilesIn ${reads} \
			--outFileNamePrefix "${params.sampleID}."
		"""
	}
}


/* 
 * STAR second pass.
 * A new index is then created using splice junction information contained in the file SJ.out.tab from the first pass.
 */

if(params.paired_rna_seq) {

	/* Runs for paired-end reads when paired_rna_seq == true */

	 process star_second_pass {

		tag "$sample"

		container "quay.io/biocontainers/star:2.6.0b--0"

		/*afterScript "find $workflow.workDir -name ${reads[0]} -delete; find $workflow.workDir -name ${reads[1]} -delete; find $workflow.workDir -name ${params.sampleID}.SJ.out.tab -delete"*/

		input:
		file genome_fasta
		file("*") from star_first_pass_rna_seq_ch

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch2

		output:
		file("${params.sampleID}.rna.aligned.bam") into star_second_pass_rna_seq_ch
		val 'STAR mapping for variant calling is done' into star_mapping_variant_call

		script:
		"""
		STAR --runMode genomeGenerate \
			--genomeDir . \
			--genomeFastaFiles ${genome_fasta} \
			--sjdbFileChrStartEnd ${params.sampleID}.SJ.out.tab \
			--sjdbOverhang 75 \
			--runThreadN ${params.threads.mapping}

		STAR --genomeDir . \
			--runThreadN ${params.threads.mapping} \
			--outSAMtype BAM Unsorted \
			--readFilesCommand zcat \
			--readFilesIn ${reads[0]} ${reads[1]} \
			--outFileNamePrefix "${params.sampleID}.rna."

		mv ${params.sampleID}.rna.Aligned.out.bam ${params.sampleID}.rna.aligned.bam
		"""
	}
}

else {

	/* Runs for paired-end reads when paired_rna_seq == false */

	process star_second_pass {

		tag "$sample"

		container "quay.io/biocontainers/star:2.6.0b--0"

		/*afterScript "find $workflow.workDir -name ${reads} -delete; find $workflow.workDir -name ${params.sampleID}.SJ.out.tab -delete"*/

		input:
		file genome_fasta
		file("*") from star_first_pass_rna_seq_ch

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch2

		output:
		file("${params.sampleID}.rna.aligned.bam") into star_second_pass_rna_seq_ch
		val 'STAR mapping for variant calling is done' into star_mapping_variant_call

		script:
		"""
		STAR --runMode genomeGenerate \
			--genomeDir . \
			--genomeFastaFiles ${genome_fasta} \
			--sjdbFileChrStartEnd ${params.sampleID}.SJ.out.tab \
			--sjdbOverhang 75 \
			--runThreadN ${params.threads.mapping}

		STAR --genomeDir . \
			--runThreadN ${params.threads.mapping} \
			--outSAMtype BAM Unsorted \
			--readFilesCommand zcat \
			--readFilesIn ${reads} \
			--outFileNamePrefix "${params.sampleID}.rna."

		mv ${params.sampleID}.rna.Aligned.out.bam ${params.sampleID}.rna.aligned.bam
		"""
	}
}


/* Use Picard tools to add group to raw sam file */

process picard_add_or_replace_read_groups {
	
	tag "rna-seq"

	container "biocontainers/picard:2.3.0"

	/*afterScript "find $workflow.workDir -name ${bam_file_raw} -delete"*/
	
	input:
	file(bam_file_raw) from star_second_pass_rna_seq_ch
	
	output:
	file("${bam_file_raw.baseName}.grouped.sorted.bam") into picard_added_group_bam_rna_seq_ch
	
	script:
	"""
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar AddOrReplaceReadGroups \
		I=${bam_file_raw} \
		O="${bam_file_raw.baseName}.grouped.sorted.bam" \
		SO=coordinate \
		RGLB=${params.sampleID}.rna-seq \
		RGPL=illumina \
		RGPU=${params.sampleID}.rna-seq \
		RGSM=${params.sampleID}.rna-seq
	"""
}


/* Use Picard tools to mark duplicates */

process picard_mark_duplicates {
	
	tag "rna-seq"

	container "biocontainers/picard:2.3.0"
	
	input:
	file(bam_file_added_group) from picard_added_group_bam_rna_seq_ch
	
	output:
	file(bam_file_added_group) into picard_mark_duplicates_rna_seq_ch1
	set file("${bam_file_added_group.baseName}.deduplicated.bam"), file("${bam_file_added_group.baseName}.deduplicated.bai") into picard_mark_duplicates_rna_seq_ch2
	
	script:
	"""	
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar MarkDuplicates \
		I=${bam_file_added_group} \
		O="${bam_file_added_group.baseName}.deduplicated.bam" \
		METRICS_FILE="${bam_file_added_group.baseName}.output.metrics" \
		VALIDATION_STRINGENCY=SILENT \
		CREATE_INDEX=true
	"""
}


/* Create summary mapping using Samtools */

process samtools_flagstat {
	
	tag "rna-seq"

	container "biocontainers/samtools:1.3.1"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	/*afterScript "find $workflow.workDir -name ${bam_file_added_group} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "move", overwrite: true
	
	input:
	file(bam_file_added_group) from picard_mark_duplicates_rna_seq_ch1
	
	output:
	file("*.mapping.statistic.txt") into samtools_flagstat_rna_seq_ch
	
	script:
	"""	
	samtools flagstat ${bam_file_added_group} > "${bam_file_added_group.baseName}.mapping.statistic.txt"
	"""
}


/* 
 * Use GATK SplitNCigarReads to split exon segments (getting rid of Ns but maintaining grouping information) 
 * and hard­clip any sequences overhanging into the intronic regions 
 */

process gatk_split_n_cigar_reads {
	
	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${bam_file_marked_duplicate} -delete; find $workflow.workDir -name ${bai_file_marked_duplicate} -delete"*/
		
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	set file(bam_file_marked_duplicate), file(bai_file_marked_duplicate) from picard_mark_duplicates_rna_seq_ch2
	
	output:
	set file("${bam_file_marked_duplicate.baseName}.splitted.bam"), file("${bam_file_marked_duplicate.baseName}.splitted.bai") into split_n_cigar_rna_seq_ch1
	set file("${bam_file_marked_duplicate.baseName}.splitted.bam"), file("${bam_file_marked_duplicate.baseName}.splitted.bai") into split_n_cigar_rna_seq_ch2
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T SplitNCigarReads \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate} \
		-o "${bam_file_marked_duplicate.baseName}.splitted.bam" \
		-RMQF 255 \
		-RMQT 60 \
		-rf ReassignOneMappingQuality \
		-U ALLOW_N_CIGAR_READS
	"""
}


/* GATK realigner target creator */

process gatk_realigner_target_creator {

	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	set file(bam_file_marked_duplicate_splitted), file(bai_file_marked_duplicate_splitted) from split_n_cigar_rna_seq_ch1
	
	output:
	file("${bam_file_marked_duplicate_splitted.baseName}.intervals.list") into gatk_target_creator_rna_seq_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-nt ${params.threads.target_creator} \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate_splitted} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-o "${bam_file_marked_duplicate_splitted.baseName}.intervals.list" 		
	"""
}	


/* GATK perform local realignment of reads around indels */

process gatk_indel_realignment {
	
	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${interval_list} -delete; find $workflow.workDir -name ${bam_file_marked_duplicate_splitted} -delete; find $workflow.workDir -name ${bai_file_marked_duplicate_splitted} -delete"*/
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	file(interval_list) from gatk_target_creator_rna_seq_ch
	set file(bam_file_marked_duplicate_splitted), file(bai_file_marked_duplicate_splitted) from split_n_cigar_rna_seq_ch2
	
	output:
	set file("${bam_file_marked_duplicate_splitted.baseName}.realigned.bam"), file("${bam_file_marked_duplicate_splitted.baseName}.realigned.bai") into gatk_indel_realignment_rna_seq_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate_splitted} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-targetIntervals ${interval_list} \
		-o "${bam_file_marked_duplicate_splitted.baseName}.realigned.bam"
	"""
}


/* Use GATK to detect systematic errors in base quality scores */

process gatk_base_recalibrator {
	
	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file known_indel_file
	file mills_indel_file
	file known_dbsnps_file
	file known_dbsnps_1000_file	
	
	set file(indel_realigned_bam), file(indel_realigned_bai) from gatk_indel_realignment_rna_seq_ch
	
	output:
	set file(indel_realigned_bam), file(indel_realigned_bai), file("${indel_realigned_bam.baseName}.data.table") into gatk_base_recalibrator_rna_seq_ch
	
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-nct ${params.threads.base_recalibrator} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		--knownSites ${known_dbsnps_1000_file} \
		--knownSites ${known_dbsnps_file} \
		--knownSites ${mills_indel_file} \
		--knownSites ${known_indel_file} \
		-o "${indel_realigned_bam.baseName}.data.table"
	"""
}


/* Use GATK to write out sequence read data (for filtering, merging, subsetting etc) */

process gatk_print_reads {
	
	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${indel_realigned_bam} -delete; find $workflow.workDir -name ${indel_realigned_bai} -delete; find $workflow.workDir -name ${recalibration_table} -delete"*/
		
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	set file(indel_realigned_bam), file(indel_realigned_bai), file(recalibration_table) from gatk_base_recalibrator_rna_seq_ch
	
	output:
	set file("${indel_realigned_bam.baseName}.processed.bam"), file("${indel_realigned_bam.baseName}.processed.bai") into processed_bam_bai_rna_seq_ch
	
	script:
	"""
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T PrintReads \
		-nct ${params.threads.print_reads} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		-BQSR ${recalibration_table} \
		-o "${indel_realigned_bam.baseName}.processed.bam"
	"""
}


/* Call somatic SNPs and indels via local re-assembly of haplotypes */
	
process haplotype_caller {
	
	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${processed_bam} -delete; find $workflow.workDir -name ${processed_bai} -delete"*/
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file known_dbsnps_file
	
	set file(processed_bam), file(processed_bai) from processed_bam_bai_rna_seq_ch
		
	output:
	file("${params.sampleID}.haplotypecaller.vcf") into haplotype_caller_vcf_rna_seq_ch
		
	script:
	"""	
	java -Xmx6g -jar /usr/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ${genome_fasta} \
		-I ${processed_bam} \
		--dbsnp ${known_dbsnps_file} \
		-dontUseSoftClippedBases \
		-stand_call_conf 20.0 \
		-o "${params.sampleID}.haplotypecaller.vcf"
	"""	
}


/* Sort VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. */

process picard_sort_vcf {

	tag "rna-seq"

	container "biocontainers/picard:2.3.0"

	/*afterScript "find $workflow.workDir -name ${original_vcf} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file(original_vcf) from haplotype_caller_vcf_rna_seq_ch
	
	output:
	file("${original_vcf.baseName}.sorted.vcf") into picard_sort_vcf_rna_seq_ch

	script:
	"""
	java -Xmx6g -jar /opt/conda/share/picard-2.3.0-0/picard.jar SortVcf \
		I=${original_vcf} \
		O="${original_vcf.baseName}.sorted.vcf" \
		SEQUENCE_DICTIONARY=${genome_fasta_dict}
	"""
}


/* Filter variant calls based on INFO and FORMAT annotations. */

process gatk_variant_filtration {

	tag "rna-seq"

	container "broadinstitute/gatk3:3.8-0"

	/*afterScript "find $workflow.workDir -name ${sorted_vcf} -delete"*/

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict

	file(sorted_vcf) from picard_sort_vcf_rna_seq_ch

	output:
	file("${sorted_vcf.baseName}.filtered.vcf") into gatk_variant_filtration_vcf_rna_seq_ch

	script:
	"""
	java -jar /usr/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R ${genome_fasta} \
		-V ${sorted_vcf} \
		-window 35 \
		-cluster 3 \
		-filterName FS \
		-filter "FS > 30.0" \
		-filterName QD \
		-filter "QD < 2.0" \
		-o "${sorted_vcf.baseName}.filtered.vcf"
	"""
}


/* Removes all sites with a FILTER flag other than PASS. */

process vcftools_remove_filtered {

	tag "rna-seq"

	container "zhanglab18/vcftools:0.1.16"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	/*afterScript "find $workflow.workDir -name ${sorted_vcf} -delete"*/

	publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

	input:
	file(sorted_vcf) from gatk_variant_filtration_vcf_rna_seq_ch

	output:
	file("${sorted_vcf.baseName}.recode.vcf") into haplotype_caller_final_vcf_rna_seq_ch

	script:
	"""
	vcftools \
		--vcf ${sorted_vcf} \
		--remove-filtered-all \
		--recode \
		--out "${sorted_vcf.baseName}"
	"""
}


/* Functionally annotate genetic variants detected from diverse genomes using ANNOVAR. */

process annovar {

	tag "rna-seq"

	/*afterScript "find $workflow.workDir -name ${haplotype_caller_vcf} -delete"*/

	input:
	file(haplotype_caller_vcf) from haplotype_caller_final_vcf_rna_seq_ch

	output:
	file("${haplotype_caller_vcf.baseName}.hg19_multianno.txt") into annovar_rna_seq_ch

	script:
	"""
	perl ${params.annovar}/table_annovar.pl ${haplotype_caller_vcf} ${params.annovar}/humandb \
		-buildver hg19 \
		-out "${haplotype_caller_vcf.baseName}" \
		-protocol refGene \
		-operation g \
		-nastring . \
		-vcfinput \
		--thread ${params.threads.annovar} \
		--maxgenethread ${params.threads.annovar} \
		-polish
	"""
}


/* Combine results from multiple MS/MS search engines with Customprodbj and build customized protein database. */

process protein_db_construction {

	tag "rna-seq"

	container "zhanglab18/customprodbj:1.1.0"

	/*afterScript "find $workflow.workDir -name ${hg19_multianno} -delete"*/

	input:
	file mrna_database_fa
	file gene_annotation_txt

	file(hg19_multianno) from annovar_rna_seq_ch

	output:
	file("merge-varInfo.txt") into protein_db_construction_rna_seq_ch

	script:
	"""
	java -jar /opt/customprodbj.jar \
		-i ${hg19_multianno} \
		-d ${mrna_database_fa} \
		-r ${gene_annotation_txt} \
		-o . \
		-t
	"""
}


/* 
 * Use RSEM to estimate gene and isoform expression levels from RNA-Seq data. 
 */

if(params.paired_rna_seq) {

	/* Runs for paired-end reads when paired_rna_seq == true */

	process rsem {

		tag "$sample"

		container "zhanglab18/rsem:1.2.26"

		beforeScript "mkdir -p ${params.results}/${params.sampleID}"

		publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

		input:
		file genome_fasta
		file refseq_human_hg19_gtf
		val x from star_mapping_variant_call

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch3

		output:
		file("${params.sampleID}.genes.results") into rsem_rna_seq_ch

		script:
		"""
		/opt/RSEM-1.2.26/rsem-prepare-reference \
			--gtf ${refseq_human_hg19_gtf} \
			--star \
			--star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
			-p ${params.threads.mapping} \
			${genome_fasta} \
			/tmp/${params.sampleID}

		/opt/RSEM-1.2.26/rsem-calculate-expression \
			--paired-end \
			--star \
			--star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
			-p ${params.threads.mapping} \
			--gzipped-read-file \
			${reads[0]} \
			${reads[1]} \
			/tmp/${params.sampleID} \
			${params.sampleID}
		"""
	}
}


else {

	/* Runs for paired-end reads when paired_rna_seq == false */

	process rsem {

		tag "$sample"

		container "zhanglab18/rsem:1.2.26"

		beforeScript "mkdir -p ${params.results}/${params.sampleID}"

		publishDir "${params.results}/${params.sampleID}", mode: "copy", overwrite: true

		input:
		file genome_fasta
		file refseq_human_hg19_gtf
		val x from star_mapping_variant_call

		set sample, file(reads) from grouped_raw_reads_rna_seq_ch3

		output:
		file("${params.sampleID}.genes.results") into rsem_rna_seq_ch

		script:
		"""
		/opt/RSEM-1.2.26/rsem-prepare-reference \
			--gtf ${refseq_human_hg19_gtf} \
			--star \
			--star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
			-p ${params.threads.mapping} \
			${genome_fasta} \
			/tmp/${params.sampleID}

		/opt/RSEM-1.2.26/rsem-calculate-expression \
			--star \
			--star-path /opt/STAR-2.6.0a/bin/Linux_x86_64 \
			-p ${params.threads.mapping} \
			--gzipped-read-file \
			${reads} \
			/tmp/${params.sampleID} \
			${params.sampleID}
		"""
	}
}


process evidence_and_expression_addition {

	container "zhanglab18/neoflow:latest"

	beforeScript "mkdir -p ${params.results}/${params.sampleID}"

	//afterScript "find $workflow.workDir -name ${snv_data} -delete; find $workflow.workDir -name ${reference_matrix} -delete; find $workflow.workDir -name ${gene_expression_data} -delete; find $workflow.workDir -name ${variant_peptide_data} -delete; find $workflow.workDir -name ${ibaq_txt} -delete"

	publishDir "${params.results}/${params.sampleID}", mode: "move", overwrite: true

	input:
	file(ibaq_txt) from deepflq_wxs_ch
	file(gene_expression_data) from rsem_rna_seq_ch
	file(variant_peptide_data) from msgfplus_wxs_ch
	file(snv_data) from protein_db_construction_rna_seq_ch
	file(reference_matrix) from peptide_binding_prediction_wxs_ch
	file(pro_ref_db) from protein_db_reference

	output:
	file("*") into evidence_and_expression_addition_ch

	script:
	"""
	Rscript /opt/combine_multiple_evidences.R \
	    -i ${params.sampleID} \
	    --pro_exp ${ibaq_txt} \
	    --pro_var ${variant_peptide_data} \
	    --rna_var ${snv_data} \
	    --rna_exp ${gene_expression_data} \
	    -m ${reference_matrix} \
	    --db ${pro_ref_db} \
	    -f ${params.mhc_affinity_cutoff} \
	    -o .
	"""
}


/* Prints help when asked for and exits */

if (params.help) {
	log.info ''
	log.info '---------------------------------------------------------------------------------------------------------------------'
	log.info '                                               N E O F L O W                                                         '
	log.info '---------------------------------------------------------------------------------------------------------------------'
	log.info ''
	log.info 'Usage:'
	log.info ''
	log.info '	nextflow run neoflow.nf --rna_seq_read_1 STRING --rna_seq_read_2 STRING --normal_wxs_read_1 STRING --normal_wxs_read_2 STRING --tumor_wxs_read_1 STRING --tumor_wxs_read_2 STRING --paired_dna_seq BOOLEAN --paired_rna_seq BOOLEAN --sampleID STRING  --results STRING --genome_fasta FILE --known_dbsnps_file FILE --cosmic_file FILE --known_indel_file FILE --known_dbsnps_1000_file FILE --mills_indel_file FILE --hla_reference_dna FILE --refseq_human_hg19_gtf FILE'
	log.info ''
	log.info '	Mandatory arguments:'
	log.info ''
	log.info '	--sampleID                STRING     a unique string to identify the sample source or ID.'
	log.info '	--paired_dna_seq          BOOLEAN    <true> for paired-end DNA-SEQ reads. <false> otherwise.'
	log.info '	--paired_rna_seq          BOOLEAN    <true> for paired-end RNA-SEQ reads. <false> otherwise.'
	log.info '	--genome_fasta            FILE       indexed hg19 Human reference genome.'
	log.info '	--known_dbsnps_file       FILE       indexed dbSNP file hg19 VCF.'
	log.info '	--cosmic_file             FILE       indexed VCF file of COSMIC sites.'
	log.info '	--known_indel_file        FILE       indexed input VCF file(s) with known indels.'
	log.info '	--known_dbsnps_1000_file  FILE       indexed dbSNP file high confidence hg19 sites VCF.'
	log.info '	--mills_indel_file        FILE       indexed mills indel reference from the GATK resource bundle.'
	log.info '	--hla_reference_dna       FILE       the HLA reference fasta files. Available at /path/to/OptiType/data/'
	log.info '	--refseq_human_hg19_gtf   FILE       RefSeq genes data from the UCSC table browser for genome version hg19.'
	log.info '	--mzml_files_dir          DIR        directory containing mzML files (mass spectrometer output files)'
	log.info ''
	log.info '	Optional arguments:'
	log.info ''
	log.info '	--help                               show this message and exit.'
	log.info '	--results                 STRING     output results directory. Defaults to <results> in the current working directory.'
	log.info ' 	--rna_seq_read_1          STRING     file path to gzipped RNA-SEQ sample read 1'
	log.info ' 	--rna_seq_read_2          STRING     file path to gzipped RNA-SEQ sample read 2'
	log.info ' 	--tumor_wxs_read_1        STRING     file path to gzipped tumor sample read 1.'
	log.info ' 	--tumor_wxs_read_2        STRING     file path to gzipped tumor sample read 2.'
	log.info ' 	--normal_wxs_read_1       STRING     file path to gzipped normal sample read 1.'
	log.info ' 	--normal_wxs_read_2       STRING     file path to gzipped normal sample read 2.'	
	log.info '	--msgfplus.t              STRING     precursor Mass Tolerance (e.g. 2.5Da, 20ppm or 0.5Da,2.5Da, Default: 20ppm).' 																							
	log.info '	--msgfplus.tda            BINARY     if 0: don nott search decoy database (default), 1: search decoy database.' 																									
	log.info '	--msgfplus.m              INT        fragment Method ID (0: as written in the spectrum or CID if no info (default), 1: CID, 2: ETD, 3: HCD).' 																	
	log.info '	--msgfplus.inst           INT        instrument ID (0: low-res LCQ/LTQ (default), 1: High-res LTQ, 2: TOF, 3: Q-Exactive).' 																					
	log.info '	--msgfplus.e              INT        enzyme ID (0: unspecific cleavage, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: glutamyl endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage).' 
	log.info '	--msgfplus.protocol       INT        protocol ID (0: NoProtocol (default), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho).' 																					
	log.info '	--msgfplus.minLength      INT        minimum peptide length to consider, default: 6.' 																																
	log.info '	--msgfplus.maxLength      INT        maximum peptide length to consider, default: 40.' 																																
	log.info '	--msgfplus.n              INT        number of matches per spectrum to be reported, default: 1.' 																													
	log.info '	--msgfplus.ti             STRING     isotope Error Range. Range of allowed isotope peak errors, default: 0,1.' 																								
	log.info '	--msgfplus.addFeatures    BINARY     0: output basic scores only (default), 1: output additional features.' 																									
	log.info '	--msgfplus.minCharge      INT        minimum precursor charge to consider if charges are not specified in the spectrum file, default: 2.' 																		
	log.info '	--msgfplus.maxCharge      INT        maximum precursor charge to consider if charges are not specified in the spectrum file, default: 3.' 																		
	log.info '	--msgfplus.ntt            INT        number of Tolerable Termini, default: 2.'
	log.info '	--annovar                 FILE       file path to the ANNOVAR directory. i.e. path/to/annovar/.'
	log.info '	--netmhcpan               FILE       file path to the netmhcpan executable. Should be found at .../netMHCpan-4.0/netMHCpan.'
	log.info '	--mhc_affinity_cutoff     DOUBLE     The cutoff of peptide MHC-binding affinity. Default: 150.'
	log.info ''
	exit 0
}
