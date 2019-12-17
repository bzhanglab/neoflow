#!/usr/bin/env nextflow

params.help = null

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => variant annotation and customized database construction
    =========================================
    Usage:
    nextflow run neoflow_db.nf
    """.stripIndent()
}

// parameters for varaint annotation: Annovar
annovar_buildver = params.annovar.buildver // 
annovar_protocol = params.annovar.protocol //
annovar_operation = params.annovar.operation 
annovar_reference = params.annovar.reference
annovar_anno_dir  = file(params.annovar.anno_dir)
annovar_soft = params.annovar.annovar
experiment_name = params.experiment


map_file = file(params.map_file)

out_dir = file(params.cpj_out_dir)

/* Pre-process vcf files for annovar input */

process run_annovar {

	tag "$experiment_name"


	input:
	file map_file

	output:
	file(map_file) into annovar_results_ch

	script:
	"""
	#!/bin/sh
	mkdir -p ${annovar_anno_dir}/somatic
	mkdir -p ${annovar_anno_dir}/germline
	input_file=${map_file}_germline_somatic.txt
	echo "sample	somatic	germline" > \$input_file

	sed "s|base_dir|$baseDir|g" ${map_file} > ${map_file}_used

	sed 1d ${map_file}_used | while read -r experiment sample germline_vcf somatic_vcf
	do
		echo "\${sample}	${annovar_anno_dir}/somatic/\${sample}.${annovar_buildver}_multianno.txt	${annovar_anno_dir}/germline/\${sample}.${annovar_buildver}_multianno.txt" >> \$input_file
		perl ${annovar_soft}/convert2annovar.pl \
			-format vcf4 \
			\$germline_vcf > \${germline_vcf}.avinput
		perl ${annovar_soft}/table_annovar.pl \${germline_vcf}.avinput ${annovar_reference} \
			-buildver ${annovar_buildver} \
			-out ${annovar_anno_dir}/germline/\${sample} \
			-protocol ${annovar_protocol} \
			-operation ${annovar_operation} \
			-nastring . \
            --thread 8 \
            --maxgenethread 8 \
            -polish \
            -otherinfo

        perl ${annovar_soft}/convert2annovar.pl \
			-format vcf4 \
			\$somatic_vcf > \${somatic_vcf}.avinput
		perl ${annovar_soft}/table_annovar.pl \${somatic_vcf}.avinput ${annovar_reference} \
			-buildver ${annovar_buildver} \
			-out ${annovar_anno_dir}/somatic/\${sample} \
			-protocol ${annovar_protocol} \
			-operation ${annovar_operation} \
			-nastring . \
            --thread 8 \
            --maxgenethread 8 \
            -polish 
            -otherinfo
    done

    cp \$input_file ${annovar_anno_dir}/${experiment_name}_germline_somatic.txt
	"""
}

process database_construction {

	tag "$experiment_name"

	container "zhanglab18/customprodbj:1.1.0"

	input:
	file(annovar_result_file) from annovar_results_ch

	output:
	file(annovar_result_file) into db_cons_ch

	script:
	
	"""
	#!/bin/sh
	java -jar /opt/customprodbj.jar \
		-f ${annovar_anno_dir}/${experiment_name}_germline_somatic.txt \
		-d ${annovar_reference}/${annovar_buildver}_${annovar_protocol}Mrna.fa \
		-r ${annovar_reference}/${annovar_buildver}_${annovar_protocol}.txt \
		-ref ${out_dir}/protein.pro-ref.fasta \
		-t \
		-o ${out_dir}/
	"""
}

process generate_new_db{
	
	tag "$experiment_name"

	input:
	file(annovar_result_file) from db_cons_ch

	output:
	file(annovar_result_file) into new_db_cons_ch

	script:
	"""
	#!/bin/sh
	java -cp ${baseDir}/bin/ thanks ${out_dir}/${experiment_name}_germline_somatic-var.fasta
	"""
}

process generate_decoy_db{
	
	tag "$experiment_name"

	input:
	file(annovar_result_file) from new_db_cons_ch

	output:

	script:
	"""
	Rscript ${baseDir}/bin/pga_generate_decoy.R ${out_dir}/ ${experiment_name} ${annovar_reference}/../
	"""
}



