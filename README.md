# DNAMetabarcoding

Pipeline for Amplicon Sequence Variant (ASV) resolution from microbial community fastq-formatted data. 

This Pipeline expects fastq data and outputs a comprhensive .csv file with taxonomic and abundance information for each ASV, as well as various Qiime objects needed for downstream processing in the Qiime platform.

This pipeline is based on the DNAmetabarcoding workflow developed by Jessica Parks and Jonathan Fleming and can be found here:
https://github.com/jessicaparks/DNAmetabarcoding.git

## Table of Contents

* [Installation](#installation)
* [Databases](#databases)
  * [DADA2 taxonomy databases](#dada2-taxonomy-databases)
  * [NCBI nt BLAST database](#ncbi-nt-blast-database)
  * [Taxize NCBI database](#taxize-ncbi-database)
* [DNAmetabarcoding Pipeline](#dnametabarcoding-pipeline)
  * [Pipeline Structure](#pipeline-structure)
  * [Primer trimming: cutadapt](#primer-trimming-cutadapt)
  * [Primer trimming: DADA2](#primer-trimming-cutadapt)
  * [ASV identification: DADA2](#asv-identification-dada2)
  * [Taxonomy: DADA2](#taxonomy-dada2)
  * [Taxonomy: BLAST and taxizedb](#taxonomy-blast-and-taxizedb)
 
## Installation
To install the pipeline on your local computer or server use the command listed here provided that git is installed. 
```bash
git clone 
```
This will download the code to a directory named `DNAmetabarcoding`. After downloading the code run the following command from within the folder. This will make a conda environment with the appropriate packages. 
```bash
conda create -n dnametabarcoding -f environment.yml
```

## Databases
After all dependencies have been downloaded please make sure to download all necessary databases for your chosen taxonomic identification method. 
The BLAST taxonomic method utilizes both the NCBI nt BLAST database and the taxize and the Taxize NCBI database.

All methods utilize the taxize NCBI database, and the pipeline will then merge this information with either the taxonomic information from the NCBI nt BLAST database or those results from DADA2 taxonomy functionality. 

### Taxize NCBI database
The taxize NCBI database from [taxizedb](https://ropensci.github.io/taxizedb/) is used to assign taxonomy to the taxids identified with BLAST. To download and then update this database run, `taxizedb_download.sh`. This script will download the most recent version of the database to the specified location.
```bash
cd DNAmetabarcoding
conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime
Rscript taxizedb_download.R
```


## DNAmetabarcoding pipeline
The DNAmetabarcoding pipeline has three major steps those being primer identification and trimming with CUTADAPT, resolution of amplicon seqeunce variants with DADA2, and taxonomic identification with either DADA2 or BLAST. There is also an option to skip the primer trimming step with CUTADAPT and instead remove the primers with the in-house DADA2 functions. 
WARNING: to use this method you will need to know the location of the primers within your reads. 

### Pipeline Structure
The main processing script for the pipeline is either `main.py` or `main_cutadapt.py`. Both scripts are dependent on functions from the scripts `dada2.R`, `dada2_taxonomy.R`, and `taxizedb.R`. The only difference between these two scirpts is that one utilizes cutadapt for primer trimming and the other utilizes DADA2 processes. These main processing scripts will analyze a single fastq file. A job submission script `submit_main.csh`, allows to user to run `main.py` on all files within a directory, and submit those jobs to a cluster running a IBM Spectrum LSF system. 



