#!/usr/bin/env Rscript
library(phyloseq, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(dplyr, quietly=TRUE, warn.conflicts = FALSE)

#read in arguments (passed from visualizations.py)
args <- commandArgs(trailingOnly = TRUE)
ASVfile <- args[1]
Taxafile <- args[2]
taxalevel <- args[3]
outputFile <- args[4]

#read in abundance and taxa file
ASV_mat <- read.csv(ASVfile)
taxa_mat <- read.csv(Taxafile)
#samples_df <- read.csv(sample file location and name here) #(uncomment for nmds plots if there are variables for sample data)


#name abundance rows
row.names(ASV_mat) <- ASV_mat$sequence
ASV_mat <- ASV_mat %>% select (-sequence)


#name taxa rows
row.names(taxa_mat) <- taxa_mat$sequence
taxa_mat <- taxa_mat %>% select (-sequence) 

#name sample rows (if you have a sample file displaying where each sample came from and other sample info, with the id column named sample)
#row.names(samples_df) <- samples_df$sample   #uncomment for sample file
#samples_df <- samples_df %>% select (-sample)    #uncomment for sample file

#transform into matrix
ASV_mat <- as.matrix(ASV_mat)
taxa_mat <- as.matrix(taxa_mat)

#prep tables to create a phyloseq object
ASV = otu_table(ASV_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
#samples = sample_data(samples_df) uncomment if you have a sample file

#create phyloseq object
carbom <- phyloseq(ASV, TAX) #or use carbom <- phyloseq(ASV, TAX, samples) if you have a sample file

#normalize the data
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)

#make abundance plot based on taxa level chosen
pdf(NULL)
if(taxalevel == "superkingdom"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=superkingdom, fill=superkingdom), stat="identity", position="stack")
} else if(taxalevel == "kingdom"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=kingdom, fill=kingdom), stat="identity", position="stack")
} else if(taxalevel == "phylum"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
} else if(taxalevel == "class"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=class, fill=class), stat="identity", position="stack")
} else if(taxalevel == "order"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=order, fill=order), stat="identity", position="stack")
} else if(taxalevel == "family"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=family, fill=family), stat="identity", position="stack")
} else if(taxalevel == "genus"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
} else if(taxalevel == "species"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=species, fill=species), stat="identity", position="stack")
}
ggsave(outputFile, height=8, width=12)

#make nmds plots
#uncomment for nmds plots
#carbom.ord <- ordinate(carbom, "NMDS", "bray")
#plot_ordination(carbom, carbom.ord, type="samples", color="fraction", shape="level", title="Samples") + geom_point(size=3) #uncomment for nmds plots
    #Note: fill in variable to color by for "fraction", fill in variable to shape by for "level", remove shape and/or color if not desired and you just want circular black samples
#ggsave(location for nmds plot, height=8, width=12)
