#!/usr/bin/env Rscript

# set HOME environment variable to user directory under the
# group directory because the taxize database will be
# downloaded to this location
user <- Sys.getenv("USER")
Sys.setenv(HOME = paste("/usr/local/usrapps/trnL_blast/", user, sep=""))

# load taxizedb library
library(taxizedb, quietly=TRUE)

# download the NCBI taxize database
path <- db_download_ncbi(overwrite=TRUE, verbose=TRUE)
cat("NCBI taxize db downloaded to", path, "\n")
