# NeoFlow: a proteogenomics pipeline for neoantigen discovery

## Description

NeoFlow is a standalone tool developed based on Nextflow and Docker and is developed for tumor neoantigens identification using DNA-Seq, RNA-Seq and MS/MS data.

## Installation

1. Install [Docker](https://docs.docker.com/install/) engine 1.10.x (or higher).

2. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). Note that Nextflow requires **BASH** and **Java 8** or higher to be installed or run. More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

3. Install [Python 3.0](https://www.python.org/downloads/) (or higher). The following packages are required:

	- os
	- Bio
	- csv
	- glob
	- argparse
	- subprocess
	- numpy
	- pandas
	- multiprocessing

4. Install **ANNOVAR** by following the instructions provide at [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/).

5. Install **netMHCpan** by following the instructions provide at [http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme). 

6. It's recommended to have the directory structure of the working directory as shown below.

7. Edit the [*nextflow.config*](https://github.com/bingzhang16/neoflow/blob/master/nextflow.config) as required.

## Usage

```sh
nextflow run neoflow.nf --help
```

The input parameters could be passed as command line parameters, however it's suggested that users pass these parameters by editing the [*nextflow.config*](https://github.com/bingzhang16/neoflow/blob/master/nextflow.config).

| Parameters 				| Type 			| Brief description 																																							|
|:--------------------------|:--------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|--help 					| 				| show help message and exit. 																																					|
|--sampleID 				| STRING 		| a unique string to identify the sample source or ID. 																															|
|--genome_fasta 			| FILE 			| indexed hg19 Human reference genome. 																																			|
|--known_dbsnps_file 		| FILE 			| indexed dbSNP file hg19 VCF. 																																					|
|--cosmic_file 				| FILE 			| indexed VCF file of COSMIC sites. 																																			|
|--known_indel_file 		| FILE 			| indexed input VCF file(s) with known indels. 																																	|
|--known_dbsnps_1000_file 	| FILE 			| indexed dbSNP file high confidence hg19 sites VCF. 																															|
|--mills_indel_file 		| FILE 			| indexed mills indel reference from the GATK resource bundle. 																													|
|--hla_reference_dna 		| FILE 			| the HLA reference fasta files. Available at `/path/to/OptiType/data/`																											|
|--refseq_human_hg19_gtf 	| FILE 			| RefSeq genes data from the UCSC table browser for genome version hg19. 																										|
|--results 					| STRING 		| output results directory. Defaults to `<results>` in the current working directory. 																							| 
|--paired_dna_seq           | BOOLEAN       | `<true>` for paired-end DNA-SEQ reads. `<false>` otherwise. 																													|
|--paired_rna_seq			| BOOLEAN	    | `<true>` for paired-end RNA-SEQ reads. `<false>` otherwise. 																													|
|--rna_seq_read_1 			| FILE 			| file path to gzipped RNA-SEQ sample read 1. 																																	|
|--rna_seq_read_2 			| FILE 			| file path to gzipped RNA-SEQ sample read 2. 																																	|
|--tumor_wxs_read_1 		| FILE 			| file path to gzipped tumor sample read 1. 																																	|
|--tumor_wxs_read_2 		| FILE 			| file path to gzipped tumor sample read 2. 																																	|
|--normal_wxs_read_1 		| FILE 			| file path to gzipped normal sample read 1.																																	|
|--normal_wxs_read_2 		| FILE 			| file path to gzipped normal sample read 2.																																	|
|--msgfplus.t 				| STRING		| Precursor Mass Tolerance (e.g. 2.5Da, 20ppm or 0.5Da,2.5Da, Default: 20ppm). 																									|
|--msgfplus.tda 			| BINARY 		| if 0: don't search decoy database (default), 1: search decoy database. 																										|
|--msgfplus.m 				| INT 			| Fragment Method ID (0: as written in the spectrum or CID if no info (default), 1: CID, 2: ETD, 3: HCD). 																		|
|--msgfplus.inst 			| INT 			| Instrument ID (0: low-res LCQ/LTQ (default), 1: High-res LTQ, 2: TOF, 3: Q-Exactive). 																						|
|--msgfplus.e 				| INT 			| Enzyme ID (0: unspecific cleavage, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: glutamyl endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage). 	|
|--msgfplus.protocol 		| INT 			| Protocol ID (0: NoProtocol (default), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho). 																						|
|--msgfplus.minLength 		| INT 			| minimum peptide length to consider, default: 6. 																																|
|--msgfplus.maxLength 		| INT 			| maximum peptide length to consider, default: 40. 																																|
|--msgfplus.n 				| INT 			| number of matches per spectrum to be reported, default: 1. 																													|
|--msgfplus.ti 				| STRING 		| Isotope Error Range. Range of allowed isotope peak errors, default: 0,1. 																										|
|--msgfplus.addFeatures 	| BINARY 		| 0: output basic scores only (default), 1: output additional features. 																										|
|--msgfplus.minCharge 		| INT 			| minimum precursor charge to consider if charges are not specified in the spectrum file, default: 2. 																			|
|--msgfplus.maxCharge 		| INT 			| maximum precursor charge to consider if charges are not specified in the spectrum file, default: 3. 																			|
|--msgfplus.ntt 			| INT 			| number of Tolerable Termini, default: 2. 																																		|
|--annovar					| FILE			| file path to the ANNOVAR directory. i.e. `path/to/annovar/`. 																													|
|--netmhcpan 				| FILE 			| file path to the netmhcpan executable. Should be found at `.../netMHCpan-4.0/netMHCpan`.|	
|--mhc_affinity_cutoff		| DOUBLE 		| The cutoff of peptide MHC-binding affinity. Default: 150.|

#### The directory structure of the working directory.

```
├── neoflow.nf
├── nextflow.config
├── reads/
|   ├── sampleID
│   │   ├── RNA.fastq
│   │   ├── tumor_1.fastq
│   │   ├── tumor_2.fastq
│   │   ├── normal_1.fastq
│   │   ├── normal_2.fastq
├── data/
|   ├── gatk-resource/
│   │   ├── 1000G_phase1.indels.hg19.sites.vcf
│   │   ├── 1000G_phase1.indels.hg19.sites.vcf.idx
│   │   ├── 1000G_phase1.snps.high_confidence.hg19.sites.vcf
│   │   ├── 1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx
│   │   ├── b37_cosmic_v54_120711.vcf
│   │   ├── b37_cosmic_v54_120711.vcf.idx
│   │   ├── dbsnp_138.hg19.vcf
│   │   ├── dbsnp_138.hg19.vcf.idx
│   │   ├── hg19_cosmic_v54_120711.vcf
│   │   ├── hg19_cosmic_v54_120711.vcf.idx
│   │   ├── Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
│   │   ├── Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
│   │   ├── ucsc.hg19.dict
│   │   ├── ucsc.hg19.fasta
|   │   ├── ucsc.hg19.fasta.sa
│   │   ├── ucsc.hg19.fasta.amb
│   │   ├── ucsc.hg19.fasta.ann
│   │   ├── ucsc.hg19.fasta.bwt
│   │   ├── ucsc.hg19.fasta.fai
│   │   ├── ucsc.hg19.fasta.pac
|   ├── optitype/
│   │   ├── hla_reference_dna.fasta
│   │   ├── hla_reference_dna.fasta.sa
│   │   ├── hla_reference_dna.fasta.amb
│   │   ├── hla_reference_dna.fasta.ann
│   │   ├── hla_reference_dna.fasta.bwt
│   │   ├── hla_reference_dna.fasta.pac
|   ├── msgfplus/
│   │   ├── test.mgf
│   │   ├── test.mods.txt
│   ├── rsem/
│   │   ├── refseq_human_hg19_refGene_20170329.gtf
│   ├── annovar/
│   │   ├── annotate_variation.pl
│   │   ├── coding_change.pl
│   │   ├── convert2annovar.pl
│   │   ├── retrieve_seq_from_fasta.pl
│   │   ├── table_annovar.pl
│   │   ├── variants_reduction.pl
│   │   ├── myanno.hg19_multianno.csv
│   │   ├── humandb/
│   │   │   ├── GRCh37_MT_ensGeneMrna.fa
│   │   │   ├── GRCh37_MT_ensGene.txt
│   │   │   ├── hg18_refGeneMrna.fa
│   │   │   ├── hg18_refGene.txt
│   │   │   ├── hg18_refGeneVersion.txt
│   │   │   ├── hg19_avsnp147.txt
│   │   │   ├── hg19_avsnp147.txt.idx
│   │   │   ├── hg19_cytoBand.txt
│   │   │   ├── hg19_dbnsfp30a.txt
│   │   │   ├── hg19_dbnsfp30a.txt.idx
│   │   │   ├── hg19_exac03.txt
│   │   │   ├── hg19_exac03.txt.idx
│   │   │   ├── hg19_example_db_generic.txt
│   │   │   ├── hg19_example_db_gff3.txt
│   │   │   ├── hg19_MT_ensGeneMrna.fa
│   │   │   ├── hg19_MT_ensGene.txt
│   │   │   ├── hg19_refGeneMrna.fa
│   │   │   ├── hg19_refGene.txt
│   │   │   ├── hg19_refGeneVersion.txt
│   │   │   ├── hg19_refGeneWithVerMrna.fa
│   │   │   ├── hg19_refGeneWithVer.txt
│   │   ├── netMHCpan-4.0
│   │   │   ├── data  
│   │   │   ├── Linux_x86_64  
│   │   │   ├── netMHCpan  
│   │   │   ├── netMHCpan.1  
│   │   │   ├── netMHCpan-4.0.readme  
│   │   │   ├── scratch  
│   │   │   ├── test
```

#### The Configuration File.

The [*nextflow.config*](https://github.com/bingzhang16/neoflow/blob/master/nextflow.config) file should be updated with the following information:

```
params {
	sampleID               = "XXXX"
	results                = "path/to/results"
	paired_dna_seq         = BOOLEAN
	paired_rna_seq         = BOOLEAN
	rna_seq_read_1         = "path/to/rna_seq_read_1"
	rna_seq_read_2         = "path/to/rna_seq_read_2"
	tumor_wxs_read_1       = "path/to/tumor_wxs_read_1"
	tumor_wxs_read_2       = "path/to/tumor_wxs_read_2"
	normal_wxs_read_1      = "path/to/normal_wxs_read_1"
	normal_wxs_read_2      = "path/to/normal_wxs_read_1"
	genome_fasta           = "data/gatk-resource/ucsc.hg19.fasta"
	known_dbsnps_file      = "data/gatk-resource/dbsnp_138.hg19.vcf"
	cosmic_file            = "data/gatk-resource/hg19_cosmic_v54_120711.vcf"
	known_indel_file       = "data/gatk-resource/1000G_phase1.indels.hg19.sites.vcf"
	known_dbsnps_1000_file = "data/gatk-resource/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
	mills_indel_file       = "data/gatk-resource/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"	
	spectrum_file_mgf      = "data/msgfplus/test.mgf"
	modification_file_txt  = "data/msgfplus/test.mods.txt"
	hla_reference_dna      = "data/optitype/hla_reference_dna.fasta"
	refseq_human_hg19_gtf  = "data/rsem/refseq_human_hg19_refGene_20170329.gtf"
	annovar                = "path/to/annovar"
	netmhcpan              = "path/to/netMHCpan"
	mhc_affinity_cutoff    = 150

	msgfplus {
		t           = "20.0ppm"
		tda         = 1
		m           = 1
		inst        = 1
		e           = 1
		protocol    = 5
		minLength   = 7
		maxLength   = 45
		n           = 1
		ti          = "0,0"
		addFeatures = 1		
		minCharge   = 2
		maxCharge   = 3
		ntt         = 2
	}

	threads {
		mutect2           = 1
		annovar           = 5
		mapping           = 5
		download          = 5
		print_reads       = 5
		target_creator    = 5
		base_recalibrator = 5
		msgfplus          = 5
    }
}
```



#### Running the pipeline from a Unix/Linux host.

Provides support for single-end and paired-end raw reads.

```{bash}
$ nextflow run neoflow.nf --paired_dna_seq true --paired_rna_seq true
```

```{bash}
$ nextflow run neoflow.nf --paired_dna_seq false --paired_rna_seq false
```

#### Output
The final neoantigen prediction result is included in file "\*-filtered_neoantigens.tsv".

## Example

Please download a test data using this link: [https://doi.org/10.6084/m9.figshare.7022462.v1](https://doi.org/10.6084/m9.figshare.7022462.v1)

## Citation
To cite NeoFlow in publications, please use:

## Contribution
Contributions to the package are more than welcome.

