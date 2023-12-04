#!/usr/bin/env Rscript

# load DADA2 library
library(dada2, quietly=TRUE)

# get input file path, output path, and reference file path from arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]
ref <- args[3]

# read input file into a data frame
seqtab <- read.csv(file)

# assign the taxonomy
taxa <- assignTaxonomy(
    seqtab,
    ref,
    multithread = TRUE,
    tryRC = TRUE
)

# write taxonomy to file
dir.create(dirname(output), showWarnings=FALSE, recursive=TRUE)
write.csv(taxa, output)
