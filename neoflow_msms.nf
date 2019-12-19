#!/usr/bin/env nextflow

params.help = null

/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    neoflow => Variant peptide identification
    =========================================
    Usage:
    nextflow run neoflow_msms.nf
    """.stripIndent()
}

experiment = params.experiment
mgf_dir = file(params.mgf_dir)
cpj_result_dir = file(params.cpj_out_dir)
software = params.search_engine
pepquery_result_dir = file(params.pepquery.out_dir)
autort_input_folder = file(params.autoRt.input)
autort_output_folder = file(params.autoRt.output)
autoRT_path = params.autoRt.autoRT_path
autoRT = params.autoRt.run

out_dir = file(params.search_result_dir)

process msms_search{

    tag "$experiment"

    container "bzhanglab/neoflow:1.0"

    publishDir "$out_dir", mode: "copy", overwrite: true

    input:
    file(mgf_dir)

    output:
    set file("*") into msgf_results_ch

    script:
    if (software == "msgf") {
        """
        #!/bin/sh

        for mgf_file in ${mgf_dir}/*.mgf
        do
            basename=`basename \${mgf_file} .mgf`
            java -Xmx10g -jar /opt/MSGFPlus.jar \
                -s \$mgf_file \
                -d ${cpj_result_dir}/${experiment}_target_decoy.fasta \
                -conf ${baseDir}/bin/msgf+/msgf_conf.txt \
                -o \${basename}.mzid
        done
        """
    } else if (software == "comet") {
        """
        #!/bin/sh
        database="${cpj_result_dir}/${experiment}_target_decoy.fasta"
        sed "s|database_path|\$database|g" ${baseDir}/bin/comet/decoy_model.params > ./${experiment}.params
        for mgf_file in ${mgf_dir}/*.mgf
        do
            basename=`basename \${mgf_file} .mgf`
            /opt/comet.2018014.linux.exe -P./${experiment}.params -N./\${basename}_rawResults \${mgf_file}
            sed -i '1d' \${basename}_rawResults.txt
            sed -i '1 s/\$/\tna/' \${basename}_rawResults.txt
        done
        """
        
    } else if (software == "xtandem") {
        """
        #!/bin/sh
        database="${cpj_result_dir}/${experiment}_target_decoy.fasta"
        sed "s|database_path|\$database|g" ${baseDir}/bin/xtandem/decoy_db_model.xml > ${experiment}_db.xml
        for mgf_file in ${mgf_dir}/*.mgf
        do
            basename=`basename \${mgf_file} .mgf`
            db_par="${experiment}_db.xml"
            result="\${basename}.xml"
            sed "s|fraction_mgf|\$mgf_file|g; s|db_parameter|\$db_par|g; s|result_path|\$result|g" ${baseDir}/bin/xtandem/decoy_parameter_model.xml > \${basename}_parameter.xml
        
            /opt/tandem-linux-17-02-01-4/bin/./tandem.exe \${basename}_parameter.xml

            java -Xmx10g -jar /opt/mzidlib-1.7/mzidlib-1.7.jar Tandem2mzid \
                \${basename}.xml \
                \${basename}.mzid \
                -outputFragmentation false \
                -decoyRegex XXX_ \
                -databaseFileFormatID MS:1001348 \
                -massSpecFileFormatID MS:1001062 \
                -idsStartAtZero false \
                -compress false \
                -proteinCodeRegex "\\S+"
        done
        """       
    }
    
}

process convert_mzid_calculata_fdr{

    tag "$experiment"

    container "proteomics/pga:latest"

    input:
    file(mzid_files) from msgf_results_ch.collect()
    file(mgf_dir)

    output:
    file(mgf_dir) into pga_result_ch

    script:
    if (software == "msgf") {
        """
        mkdir -p ${out_dir}/peptide_level/
        mkdir -p ${out_dir}/psm_level/
        mkdir -p ${out_dir}/peptide_level/global_fdr
        mkdir -p ${out_dir}/psm_level/global_fdr
        Rscript ${baseDir}/bin/msgf+/convert_merge_calculate_fdr.R ${out_dir}/ ${cpj_result_dir}/${experiment}_germline_somatic-var.fasta.new.fasta $experiment
        """
    } else if (software == "comet") {
        """
        mkdir -p ${out_dir}/peptide_level/
        mkdir -p ${out_dir}/psm_level/
        mkdir -p ${out_dir}/peptide_level/global_fdr
        mkdir -p ${out_dir}/psm_level/global_fdr
        Rscript ${baseDir}/bin/comet/convert_merge_calculate_fdr.R ${out_dir}/ ${cpj_result_dir}/${experiment}_germline_somatic-var.fasta.new.fasta $experiment
        """
    
    } else if (software == "xtandem") {
        """
        mkdir -p ${out_dir}/peptide_level/
        mkdir -p ${out_dir}/psm_level/
        mkdir -p ${out_dir}/peptide_level/global_fdr
        mkdir -p ${out_dir}/psm_level/global_fdr
        Rscript ${baseDir}/bin/xtandem/convert_merge_calculate_fdr.R ${out_dir}/ ${cpj_result_dir}/${experiment}_germline_somatic-var.fasta.new.fasta $experiment
        """
    }
}

process generate_mgf_index{

    tag "$experiment"

    container "bzhanglab/neoflow:1.0"

    input:
    file(mgf_dir) from pga_result_ch

    output:
    file(mgf_dir) into generate_index_ch

    script:
    """
    #!/bin/sh

    java -cp ${baseDir}/bin/PDV-1.5.4_generate_index/PDV-1.5.4_generate_index.jar PDVGUI.generate_index $mgf_dir $mgf_dir
    """
}

process prepare_pepquery_input{

  tag "$experiment"

  container "proteomics/pga:latest"

  input:
  file(mgf_dir) from generate_index_ch

  output:
  file (mgf_dir) into pre_process_ch

  script:
  """
  #!/bin/sh
  rm ${out_dir}/peptide_level/global_fdr/*_fraction.txt_var_pep.txt
  Rscript ${baseDir}/bin/generate_pepquery_input.R ${out_dir}/peptide_level/
  """
}

process run_pepquery{

  tag "$experiment"

  container "bzhanglab/neoflow:1.0"

  input:
  file (mgf_dir) from pre_process_ch

  output:
  file (mgf_dir) into pepquery_ch

  script:
  """
  #!/bin/sh

  for file in ${mgf_dir}/*mgf
  do
    fraction=`basename \${file} .mgf`
    mkdir -p ${pepquery_result_dir}/\${fraction}
    java -Xmx10g -jar /opt/PepQuery_v1.3.0/pepquery-1.3.jar \
      -pep ${out_dir}/peptide_level/global_fdr/\${fraction}_fraction.txt_var_pep.txt \
      -db ${cpj_result_dir}/protein.pro-ref.fasta \
      -ms \${file} \
      -fixMod 6 \
      -varMod 107 \
      -cpu 8 \
      -minScore 12 \
      -tol 20 \
      -itol 0.05 \
      -n 10000 \
      -um \
      -m 1 \
      -prefix \$fraction \
      -o ${pepquery_result_dir}/\${fraction}/
  done
  """
}

if (autoRT == "Yes") {
    
    process psm_fdr_fraction{

    tag "$experiment"

    container "proteomics/pga:latest"

    input:
    file(mgf_dir) from pepquery_ch

    output:
    file(mgf_dir) into psm_fdr_ch

    script:
    """
    #!/bin/sh
    Rscript ${baseDir}/bin/psm_level_fraction.R ${out_dir}/psm_level/
    """
    }

    process combine_psm_input{

    tag "$experiment"

    container "proteomics/pga:latest"

    input:
    file(mgf_dir) from psm_fdr_ch

    output:
    file(mgf_dir) into combine_fdr_ch

    script:
    """
    #!/bin/sh
    mkdir -p ${autort_input_folder}/intermediate_files/
    for file in ${out_dir}/peptide_level/global_fdr/*fraction.txt
    do
    basename=`basename \${file}`
    fraction=`echo "\${basename/_fraction\\.txt/}"`
    psm_file=${out_dir}/psm_level/global_fdr/\${basename}
    Rscript ${baseDir}/bin/add_psm_result.R \$file ${mgf_dir}/\${fraction}.mgf_index.txt ${autort_input_folder}/intermediate_files/\${fraction} \$psm_file
    done
    """
    }

    process prepare_autoRT_input{

    tag "$experiment"

    container "bzhanglab/python:3.6.8"

    input:
    file(mgf_dir) from combine_fdr_ch  

    output:
    file(mgf_dir) into prepare_autoRT_ch

    script:
    """
    #!/bin/sh
    for file in ${autort_input_folder}/intermediate_files/*normal_psm.txt
    do
    basename=`basename \${file}`
    fraction=`echo "\${basename/_normal_psm\\.txt/}"`
    python3 ${baseDir}/bin/prepare_train_data_${software}.py \$fraction \$file ${autort_input_folder}/ "train"
    done

    for file in ${autort_input_folder}/intermediate_files/*variant_psm.txt
    do
    basename=`basename \${file}`
    fraction=`echo "\${basename/_variant_psm\\.txt/}"`
    python3 ${baseDir}/bin/prepare_train_data_${software}.py \$fraction \$file ${autort_input_folder}/ "prediction"
    done

    """
    }

    process run_autoRT{
    
    tag "$experiment"

    container "bzhanglab/python:3.6.8"

    input:
    file(mgf_dir) from prepare_autoRT_ch

    output:


    """
    #!/bin/sh
    for file in ${autort_input_folder}/*train.txt
    do
    basename=`basename \${file} _train.txt`
    mkdir -p ${autort_output_folder}/\${basename}
    python3 $autoRT_path/autort.py train -i \$file -o ${autort_output_folder}/\${basename}/ -e 40 -b 64 -u m -m ${baseDir}/bin/autoRT/models/base_models_PXD006109/model.json -rlr -n 10
    python3 $autoRT_path/autort.py predict -t ${autort_input_folder}/\${basename}_prediction.txt -s ${autort_output_folder}/\${basename}/model.json -o ${autort_output_folder}/\${basename}/ -p \${basename}
    done
    """

    }
}




