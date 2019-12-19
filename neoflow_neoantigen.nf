#!/usr/bin/env nextflow

params.help = null

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => Neoantigen prediction
    =========================================
    Usage:
    nextflow run neoflow_neoantigen.nf
    """.stripIndent()
}

map_file = file(params.map_file)
out_dir = file(params.neoantigen.out_dir)
hla_class_dir = file(params.neoantigen.hla_class_results_dir)
cpj_result_dir  = file(params.cpj_out_dir)
netMHC_path = params.neoantigen.netMHC_path


process split_infor {
  tag "Separate customprodbj results"

  container "proteomics/pga:latest"

  input:
  file map_file

  output:
  file map_file into new_info_ch

  script:
  """
  #!/bin/sh
  sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
        Rscript ${baseDir}/bin/split_infor.R ${cpj_result_dir}/ \${experiment} \$sample
    done
  """
}

process MHC_peptide_prediction {
    tag "netmhc pan prediction"

    //publishDir "${out_dir}", mode: "copy", overwrite: true

    container "bzhanglab/python:3.6.8"
    
    input:
    file map_file from new_info_ch

    output:

    script:
    """
    #!/bin/sh
    sed 1d $map_file | while read -r experiment sample germline_vcf somatic_vcf
    do
        python3 ${baseDir}/bin/binding_prediction.py \
            -id \$sample \
            -a ${hla_class_dir}/\${sample}/\${sample}_optitype_result.tsv \
            -fa ${cpj_result_dir}/\${experiment}_germline_somatic-var.fasta.new.fasta \
            -txt ${cpj_result_dir}/\${sample}_new_somatic_varInfo.txt \
            -o $out_dir \
            -netmhcpan $netMHC_path \
            -nt 6
    done
    """

}
