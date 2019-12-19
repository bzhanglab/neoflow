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

2. Install [Docker](https://docs.docker.com/install/) engine 1.10.x (or higher).

3. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). Note that Nextflow requires **BASH** and **Java 8** or higher to be installed or run. More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

4. Install **ANNOVAR** by following the instruction at [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/).

5. Install **netMHCpan 4** by following the instruction at [http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme). 

6. Install **AutoRT** by following the instruction at [https://github.com/bzhanglab/AutoRT](https://github.com/bzhanglab/AutoRT). This is optional and it is only required when users want to use the RT based validation for novel peptide identifications using AutoRT.

After **ANNOVAR**, **netMHCpan**  and **AutoRT** are installed, users must update the path for each of these tools in the configure file **nextflow.config** as shown below before users run NeoFlow.

```
annovar{
   ...
   annovar = "/data/tools/annovar/" // update this accordingly
   ...
}

neoantigen{
   ...
   netMHC_path = "/data/tools/netMHCpan-4.0/netMHCpan" // update this accordingly
   ...
}

autoRt{
    ...
    autoRT_path = "/data/tools/autoRT/"
    ...
}
```

All other tools used by NeoFlow have been dockerized and will be automatically installed when NeoFlow is run in the first time on a computer.

## Usage

### 1. Variant annotation and customized database construction

```sh
nextflow run neoflow_db.nf
```
### 2. Variant peptide identification
Please note that the customized database generated in the first step will be used in this step. 
```sh
nextflow run neoflow_msms.nf
```
### 3. HLA typing
```sh
nextflow run neoflow_hlatyping.nf
```
### 4. Neoantigen prediction
Please note that the results generated in step 1-3 will be used in this step. 
```sh
nextflow run neoflow_neoantigen.nf
```
For each of these modules, all the parameters are contained in the **nextflow.config** file.

##  Example data

An example data can be downloaded by clicking [test data ](http://pdv.zhang-lab.org/data/download/neoflow_example_data/test_data.zip). 
