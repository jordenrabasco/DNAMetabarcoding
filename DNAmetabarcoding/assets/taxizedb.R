#!/usr/bin/env Rscript

# set HOME environment variable to user directory under the
# group directory because the taxize database was downloaded
# to this location
#user <- Sys.getenv("USER")
user <- 'tmtiedge'
Sys.setenv(HOME = paste("/usr/local/usrapps/trnL_blast/", user, sep=""))

# load libraries
library(taxizedb, quietly=TRUE)
library(data.table)

# get the input and output file paths
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]

# load the ncbi taxize db connection
src <- src_ncbi()

# read input data (list of taxids)
data <- read.csv(file, header=FALSE, col.names="taxid")

# taxonomic classification for each taxid
results <- vector('list', length(data))
i <- 1
for (taxid in data$taxid){
  result <- classification(taxid)[[1]]
  result$taxid <- taxid
  results[[i]] <- result
  i <- i + 1
}
results <- rbindlist(results, fill=TRUE)

# write results to output file
write.csv(results, output, row.names=FALSE)
