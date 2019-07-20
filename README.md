# NeoFlow

## Description

NeoFlow: a proteogenomics pipeline for neoantigen discovery

NeoFlow includes four modules:


1. HLA typing: neoflow_hlatyping.nf
2. Variant annotation and customized database construction: neoflow_db.nf
3. Variant peptide identification: neoflow_msms.nf
4. Neoantigen prediction: neoflow_neoantigen.nf


NeoFlow supports both label free and iTRAQ/TMT data.

## Installation

1. Download neoflow:

```sh
git clone https://github.com/bingzhang16/neoflow
```

2. Install [Docker](https://docs.docker.com/install/) engine 1.10.x (or higher).

3. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). Note that Nextflow requires **BASH** and **Java 8** or higher to be installed or run. More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

4. Install [Python 3.0](https://www.python.org/downloads/) (or higher). The following packages are required:

	- os
	- Bio
	- csv
	- glob
	- argparse
	- subprocess
	- numpy
	- pandas
	- multiprocessing

5. Install **ANNOVAR** by following the instructions provide at [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/).

6. Install **netMHCpan** by following the instructions provide at [http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme). 


## Usage

```sh
nextflow run neoflow-*.nf --help
```
