

cat("\f")
 
rm(list = ls())

library("ggplot2")
library(tidyr)
library(mia)
library(miaViz)
library(tidyverse)
library(scater)
library(ape)
library(phyloseq)
library(microbiome)
library(mia)


#Read the data 

samples_df <- read.csv(file ="metadata.csv", 
                       header = TRUE)

otu_mat <- read.csv(file ="genus.csv", 
                    header = TRUE)

phyloseq <- readRDS("physeq.rds")

phyloseq


sample_names(phyloseq)

ntaxa(phyloseq)

nsamples(phyloseq)

rank_names(phyloseq)

#https://microbiome.github.io/tutorials/Preprocessing.html

sample_variables(phyloseq)

metadata <- meta(phyloseq)
taxonomy <- tax_table(phyloseq)
# Absolute abundances
otu.absolute <- abundances(phyloseq)
# Relative abundances
otu.relative <- abundances(phyloseq, "compositional")


#total read counts

reads_sample <- readcount(phyloseq)



# convert phyloseq to TSE
tse <- makeTreeSummarizedExperimentFromPhyloseq(phyloseq) 
tse


assays(tse)


se <- relAbundanceCounts(tse)

rowData(se)
colData(se)
rowTree(se)
rowLinks(se)


dim(se) ##

# inspect possible values for SampleType
unique(se$Country)


library(mia)
crc <- meltAssay(se,
                       add_row_data = TRUE,
                       add_col_data = TRUE,
                       abund_values = "relabundance")
crc
