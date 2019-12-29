# NeoFlow

## Overview

**NeoFlow: a proteogenomics pipeline for neoantigen discovery**

NeoFlow includes four modules:

1. Variant annotation and customized database construction: neoflow_db.nf;
2. Variant peptide identification: neoflow_msms.nf;

   * MS/MS searching. Three search engines are available: [MS-GF+](https://github.com/MSGFPlus/msgfplus), [X!Tandem](https://www.thegpm.org/tandem/) and [Comet](http://comet-ms.sourceforge.net/);
   * FDR estimation: global FDR estimation;
   * Novel peptide validation by [PepQuery](http://pepquery.org/);
   * RT based validation for novel peptide identifications using [AutoRT](https://github.com/bzhanglab/AutoRT): optional (GPU required).

3. HLA typing: neoflow_hlatyping.nf;
4. Neoantigen prediction: neoflow_neoantigen.nf


NeoFlow supports both label free and iTRAQ/TMT data.

## Installation

1. Download neoflow:

```sh
git clone https://github.com/bzhanglab/neoflow
```

2. Install [Docker](https://docs.docker.com/install/) (>=19.03).

3. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

4. Install **ANNOVAR** by following the instruction at [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/).

5. Install **netMHCpan 4.0** by following the instruction at [http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme). 

6. Install [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) (>=2.2.2) for [**AutoRT**](https://github.com/bzhanglab/AutoRT/) by following the instruction at [https://github.com/NVIDIA/nvidia-docker](https://github.com/NVIDIA/nvidia-docker). This is optional and it is only required when users want to use the RT based validation for novel peptide identifications using AutoRT.

All other tools used by NeoFlow have been dockerized and will be automatically installed when NeoFlow is run in the first time on a computer.

## Usage

### 1. Variant annotation and customized database construction

```sh
 $ nextflow run neoflow_db.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `neoflow_db.nf` [irreverent_faggin] - revision: 741bf1a931
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
```

The txt file for parameter `--vcf_file` contains the path of VCF file(s) and its format is shown below:

| experiment | sample | file | file_type |
|---|---|---|---|
| TMT01 | T1 | T1_somatic.vcf;T1_rna.vcf | somatic;rna |
| TMT01 | T2 | T2_somatic.vcf;T2_rna.vcf | somatic;rna |
| TMT02 | T3 | T3_somatic.vcf;T3_rna.vcf | somatic;rna |
| TMT02 | T4 | T4_somatic.vcf;T4_rna.vcf | somatic;rna |

The column of `experiment` is TMT or iTRAQ experiment name and the column of `sample` is sample name. If it's iTRAQ or TMT data, the samples from the same iTRAQ or TMT experiment should have the same `experiment` name.  All variant files (for example, somatic variant vcf file and variant calling result vcf file based on RNA-Seq data) for the same sample should be in the same row (column `file`) and different files should be separated by ";". The column of `file_type` indicates the corresponding variant types for the vcf files in column `file`. 

The ANNOVAR annotation data (`--annovar_dir`) can be downloaded following the instruction at [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). 

The output files of `neoflow_db.nf` include customized protein databases in FASTA format for each experiment, variant annotation result files for each sample.

#### Example

```sh
nextflow run neoflow_db.nf --ref_dir /data/tools/annovar/humandb_hg19/ \
                           --vcf_file example_data/test_vcf_files.tsv \
                           --annovar_dir /data/tools/annovar/ \
                           --ref_ver hg19 \
                           --out_dir output
```
Please update  inputs for parameters `--ref_dir`  and `--annovar_dir` before run the above example. The input file for `--vcf_file` can be downloaded from the [example data](http://pdv.zhang-lab.org/data/download/neoflow_example_data/example_data.tar.gz) prepared for testing. After the example data is downloaded to users' computer, unzip the data and all the testing data are available in the **example_data** folder.

The running time of above example is less than 5 minutes on a Linux server with 40 cores.

### 2. Variant peptide identification
Please note that the customized database generated in the first step will be used in this step. 
```sh
 $ ./nextflow run neoflow_msms.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `neoflow_msms.nf` [drunk_nobel] - revision: 6d58fb19bd
=========================================
neoflow => Variant peptide identification
=========================================
Usage:
nextflow run neoflow-msms.nf
MS/MS searching arguments:
  --db                        The customized protein database (target + decoy sequences) in FASTA format which is generated by neoflow_db.nf
  --ms                        MS/MS data in MGF format
  --msms_para_file            Parameter file for MS/MS searching
  --out_dir                   Output folder, default is "./"
  --prefix                    The prefix of output files
  --search_engine             The search engine used for MS/MS searching, comet=Comet, msgf=MS-GF+ or xtandem=X!Tandem

PepQuery arguments:
  --pv_enzyme                 Enzyme used for protein digestion. 0:Non enzyme, 1:Trypsin (default), 2:Trypsin (no P rule), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C
  --pv_c                      The max missed cleavages, default is 2
  --pv_tol                    Precursor ion m/z tolerance, default is 10
  --pv_tolu                   The unit of --tol, ppm or Da. Default is ppm
  --pv_itol                   The error window for fragment ion, default is 0.5
  --pv_fixmod                 Fixed modification. The format is like : 1,2,3. Different modification is represented by different number
  --pv_varmod                 Variable modification. The format is the same with --fixMod;
  --pv_refdb                  Reference protein database

AutoRT parameters:
  --rt_validation             Perform RT based validation
  
  --help                      Print help message
```

The output files of `neoflow_msms.nf` include MS/MS searching raw identification files, FDR estimation result files at both PSM and peptide levels, PepQuery validation result files. 

#### Example

```sh
nextflow run neoflow_msms.nf --ms example_data/mgf/ \
               --msms_para_file example_data/comet_parameter.txt \
               --search_engine comet \
               --db output/customized_database/neoflow_crc_target_decoy.fasta \
               --out_dir output \
               --pv_refdb output/customized_database/ref.fasta \
               --pv_tol 20 \
               --pv_itol 0.05
```
The variant peptide identification result is in this file `output/novel_peptide_identification/novel_peptides_psm_pepquery.tsv`.

The running time of above example is less than 15 minutes on a Linux server with 40 cores.

### 3. HLA typing
```sh
 $ ./nextflow run neoflow_hlatyping.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `neoflow_hlatyping.nf` [spontaneous_hawking] - revision: 5fd970e701
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
```
The  output of `neoflow_hlatyping.nf` is a txt format file containing HLA alleles for a sample. This file is generated by [OptiType](https://github.com/FRED-2/OptiType).

#### Example
```sh
nextflow run neoflow_hlatyping.nf --hla_ref_dir example_data/hla_reference \
                  --reads "example_data/dna/*_{1,2}.fastq.gz" \
                  --out_dir output/ \
                  --cpu 40
```
The HLA typing result is in this file `output/hla_type/sample1/sample1_result.tsv`.

The running time of above example is less than 10 minutes on a Linux server with 40 cores.


### 4. Neoantigen prediction
Please note that the results generated in step 1-3 will be used in this step. 
```sh
 $ ./nextflow run neoflow_neoantigen.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `neoflow_neoantigen.nf` [mighty_roentgen] - revision: e4261baca3
=========================================
neoflow => Neoantigen prediction
=========================================
Usage:
nextflow run neoflow_neoantigen.nf
Arguments:
  --var_db                  Variant (somatic) database in fasta format generated by neoflow_db.nf
  --var_info_file           Variant (somatic) information in txt format generated by neoflow_db.nf
  --ref_db                  Reference (known) protein database
  --hla_type                HLA typing result in txt format generated by Optitype
  --netmhcpan_dir           NetMHCpan 4.0 folder
  --var_pep_file            Variant peptide identification result generated by neoflow_msms.nf, optional.
  --var_pep_info            Variant information in txt format for customized database used for variant peptide identification
  --prefix                  The prefix of output files
  --out_dir                 Output directory
  --cpu                     The number of CPUs
  --help                    Print help message
```

The output of `neoflow_neoantigen.nf` is a tsv format file containing neoantigen prediction result as shown below:

Variant\_ID|Chr|Start|End|Ref|Alt|Variant\_Type|Variant\_Function|Gene|mRNA|Neoepitope|Variant\_Start|Variant\_End|AA\_before|AA\_after|HLA\_type|netMHCpan\_binding\_affinity\_nM|netMHCpan\_precentail\_rank|protein\_var\_evidence\_pep
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
VAR\|NM\_002536\|10054|chrX|48418659|48418659|G|A|nonsynonymous SNV|protein-altering|TBC1D25|NM\_002536|TGFGGHRG|1|1|A|T|HLA-A*01:01|44216.6|88.5537|-
VAR\|NM\_002536\|10054|chrX|48418659|48418659|G|A|nonsynonymous SNV|protein-altering|TBC1D25|NM\_002536|TGFGGHRG|1|1|A|T|HLA-C*07:01|43330|73.7774|-
VAR\|NM\_002536\|10054|chrX|48418659|48418659|G|A|nonsynonymous SNV|protein-altering|TBC1D25|NM\_002536|TGFGGHRG|1|1|A|T|HLA-B*08:01|35925.8|70.8561|-
VAR\|NM\_001348265\|10055|chrX|48418659|48418659|G|A|nonsynonymous SNV|protein-altering|TBC1D25|NM\_001348265|TGFGGHRG|1|1|A|T|HLA-A*01:01|44216.6|88.5537|-
VAR\|NM\_001348265\|10055|chrX|48418659|48418659|G|A|nonsynonymous SNV|protein-altering|TBC1D25|NM\_001348265|TGFGGHRG|1|1|A|T|HLA-C*07:01|43330|73.7774|-

#### Example
```sh
nextflow run neoflow_neoantigen.nf --prefix sample1 \
                   --hla_type output/hla_type/sample1/sample1_result.tsv \
                   --var_db output/customized_database/sample1-somatic-var.fasta \
                   --var_info_file output/customized_database/sample1-somatic-varInfo.txt \
                   --out_dir output/ \
                   --netmhcpan_dir /data/tools/netMHCpan-4.0/ \
                   --cpu 40 \
                   --ref_db output/customized_database/ref.fasta \
                   --var_pep_file output/novel_peptide_identification/novel_peptides_psm_pepquery.tsv \
                   --var_pep_info output/customized_database/neoflow_crc_anno-varInfo.txt
```
Please update  input for parameter `--netmhcpan_dir` before run the above example. 

The neoantigen prediction result is in this file `output/novel_peptide_identification/novel_peptides_psm_pepquery.tsv`.

The running time of above example is less than 30 minutes on a Linux server with 40 cores.

##  Example data

The test data used for above examples can be downloaded by clicking [test data ](http://pdv.zhang-lab.org/data/download/neoflow_example_data/example_data.tar.gz). 
