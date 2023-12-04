#!/usr/bin/env Rscript

# load the DADA2 library
library(dada2, quietly=TRUE)

# get the input file path and output path from the arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]
feat_table_output <- args[3]
rep_seqs_output <- args[4]
trim_length <- args[5]
tempdir <- args[6]
# get the sample name from the input file path
samplename <- gsub(pattern='.fastq','',basename(file))

# filter and trim the reads (DADA2 requires no Ns)
filteredreads <- paste0(tempdir, '/', samplename, '_filtered.fastq')
filterAndTrim(
    fwd = file,
    filt = filteredreads,
    truncLen = as.integer(trim_length),
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
)

if(file.exists(filteredreads) == TRUE){

# learn the error rates
err <- learnErrors(
    fls = filteredreads,
    multithread = TRUE
)

# dereplicate identical reads
derep <- derepFastq(fls = filteredreads)

# apply sample inference
dadaresult <- dada(
    derep = derep,
    err = err,
    multithread = TRUE
)

# create sequence table and remove chimeras
samplelist <- list(dadaresult)
names(samplelist) <- c(samplename)
seqtab <- makeSequenceTable(
    samples = samplelist
)
seqtab.nochim <- removeBimeraDenovo(
    unqs = seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
)

# transpose sequence table and write to file
# data is now in columns named 'sequence' and 'abundance'
df <- t(seqtab.nochim)
df <- cbind(rownames(df), data.frame(df, row.names=NULL))
names(df) <- c('sequence', 'abundance')
dir.create(dirname(output), showWarnings=FALSE, recursive=TRUE)
write.csv(df, output)
write.table(t(seqtab.nochim), feat_table_output, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout=rep_seqs_output, ids=colnames(seqtab.nochim))
}
