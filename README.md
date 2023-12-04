# DNAMetabarcoding

Pipeline for Amplicon Sequence Variant (ASV) resolution from microbial community fastq-formatted data. 

This Pipeline expects fastq data and outputs a comprhensive .csv file with taxonomic and abundance information for each ASV.

This pipeline is based on the DNAmetabarcoding workflow developed by Jessica Parks and Jonathan Fleming and can be found here:
https://github.com/jessicaparks/DNAmetabarcoding.git

## Table of Contents

* [Installation](#installation)
 
## Installation
To install the pipeline on your local computer or server use the command listed here provided that git is installed. 
```bash
git clone 
```
This will download the code to a directory named `DNAmetabarcoding`. After downloading the code run the following command from within the folder. This will make a conda environment with the appropriate packages. 
```bash
conda create -n dnametabarcoding -f environment.yml
```
