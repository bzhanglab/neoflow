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
git clone https://github.com/bzhanglab/neoflow
```

2. Install [Docker](https://docs.docker.com/install/) engine 1.10.x (or higher).

3. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). Note that Nextflow requires **BASH** and **Java 8** or higher to be installed or run. More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

4. Install **AutoRt** by following the instructions provide at [https://github.com/bzhanglab/AutoRT](https://github.com/bzhanglab/AutoRT).

5. Install **ANNOVAR** by following the instructions provide at [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/).

6. Install **netMHCpan** by following the instructions provide at [http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme). 

7. Add **AutoRt**, **ANNOVAR** and **netMHCpan** path into configure file _nextflow.config_.

8. Download test data from [http://pdv.zhang-lab.org/data/download/neoflow_example_data/test_data.zip](http://pdv.zhang-lab.org/data/download/neoflow_example_data/test_data.zip) and uncompress it. Then put whole _test_data_ folder with neoflow scripts. 


## Usage

At first, users need build customized database.
```sh
nextflow run neoflow_db.nf
```
And then, run msms search.
```sh
nextflow run neoflow_msms.nf
```
Users need run hlatyping before neoantigen prediction.
```sh
nextflow run neoflow_hlatyping.nf
```
```sh
nextflow run neoflow_neoantigen.nf
```
Most of parameters and files' path can be modified in _nextflow.config_. For msms search engines' parameters, there are different parameter files in _bin/msgf+_,  _bin/comet_ and  _bin/xtandem_.
