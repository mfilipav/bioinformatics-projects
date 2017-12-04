####################################################
## Analysis of immune repertoires with R
## Modestas Filipavicius
## ETH Zurich
## mfilipav@ethz.ch
####################################################

## This file walks you through analysis of NGS of immune repertoires.


###################  PART I  #######################
####################################################
##### Quick start in R and immune repertoires ######
####################################################

## Set working directory where you keep the folder
## of the analysis. For example

wd = "~/immune_reps"
setwd(wd)

## Install R packages
source("http://www.bioconductor.org/biocLite.R")
biocLite()
install.packages(c("pcaMethods", "reshape", "phylotools", "VennDiagram", 
                   "fastmatch", "HDMD","pcaMethods","ConsensusClusterPlus", 
                   "corrgram", "igraph", "NMF"))

library(pcaMethods)
library(ggplot2)
library(knitr)
library(seqinr)
library(RColorBrewer)
library(xtable)
library(plyr)
library(ggplot2)
library(ShortRead)
library(grid)
library(reshape)
library(ape)
library(phylotools)
library(stringr)
library(gridExtra)
library(hexbin)
library(data.table)
library(VennDiagram)
library(scales)
library(fastmatch)
library(HDMD)
library(Biobase)
library(pcaMethods)
library(Biostrings)
library(stringdist)
library(ConsensusClusterPlus)
library(Hmisc)
library(gplots)
library(corrgram)
library(igraph)
library(NMF)

################## READ IMGT FILES ################
###################################################
### Read-in IMGT CSV output. Set "eval=TRUE" to run the code.
# Insert path to 1_Summary.txt files output from IMGT/High-VQUEST
# Files to read in are CSV IMGT output files from dataset2 and dataset4

files1 <- c("code/IMGT/dataset2_pandaseq_aa/1_Summary.txt")

# Read in selected columns: apply the function read.table
# to the list of files save in the variable "files1".
# "character" instructs the read in,
# "NULL" skips the column.
ds_files1 <- lapply(1:length(files1), function(x) read.table(files1[x],
    colClasses = c("NULL", "NULL", rep("character", 2), "NULL",
    "character", "NULL", "NULL", "NULL", "character",
    "NULL", "NULL", "NULL", "character", rep("NULL", 4),
    rep("NULL", 2), rep("NULL", 9)),
    sep="\t", header=TRUE, strip.white=TRUE,
    na.strings=TRUE, fill=TRUE, nrows=500001))
length(ds_files1)

#ds_files1[[1]][2]

# Output of the columns that have been read-in
colnames(ds_files1[[1]])
# [1] "Functionality", "V.GENE.and.allele", "V.REGION.identity..",
# "J.GENE.and.allele", "D.GENE.and.allele"
# Label the donor for each sequence in each CSV row of the read-in files
hiv_iavi17 <- rep(c("HIV-1 IAVI donor 17"),
                  times = c(nrow(ds_files1[[1]]) + nrow(ds_files1[[2]]) +
                                nrow(ds_files1[[3]]) + nrow(ds_files1[[4]]) +
                                nrow(ds_files1[[5]]) + nrow(ds_files1[[6]]) +
                                nrow(ds_files1[[7]])))
healthy_donor1 <- rep(c("Healthy donor 1"),
                      times = c(nrow(ds_files1[[8]]) + nrow(ds_files1[[9]]) +
                                    nrow(ds_files1[[10]]) + nrow(ds_files1[[11]]) +
                                    nrow(ds_files1[[12]]) + nrow(ds_files1[[13]]) +
                                    nrow(ds_files1[[14]])))

Donor <- c(hiv_iavi17, healthy_donor1)
# Bind the all the rows with "rbind"
summary_file <- do.call("rbind", ds_files1)
# Change the name to the column "V.REGION.identity.."
colnames(summary_file)[3] <- c("V.REGION.identity")


files3 <- c("code/IMGT/dataset2_pandaseq_aa/3_Nt-sequences.txt")
ds_files3 <- lapply(1:length(files3), function(x) read.csv(files3[x],
                                                           colClasses = c(rep("NULL", 6), "character", "NULL",
                                                                          rep("character", 8), rep("NULL", 7),
                                                                          "character", rep("NULL", 16),
                                                                          rep("character", 2), rep("NULL", 72)),
                                                           sep="\t", header=TRUE, strip.white=TRUE,
                                                           na.strings=TRUE, fill=TRUE, nrows=500001))


nt_sequences_file <- do.call("rbind", ds_files3)

colnames(nt_sequences_file)
# [1] "V.D.J.REGION" "V.REGION" "FR1.IMGT"
# "CDR1.IMGT" "FR2.IMGT" "CDR2.IMGT" "FR3.IMGT"
# [8] "CDR3.IMGT" "JUNCTION" "D.REGION" "J.REGION"
# "FR4.IMGT"
colnames(nt_sequences_file)[1:12] <- c("V.D.J.REGION.nt", "V.REGION.nt",
                                       "FR1.IMGT.nt", "CDR1.IMGT.nt", "FR2.IMGT.nt", "CDR2.IMGT.nt", "FR3.IMGT.nt",
                                       "CDR3.IMGT.nt", "JUNCTION.nt", "D.REGION.nt", "J.REGION.nt", "FR4.IMGT.nt")
files5 <- c("code/IMGT/dataset2_pandaseq_aa/5_AA-sequences.txt")

ds_files5 <- lapply(1:length(files5), function(x) read.table(files5[x],
                                                             colClasses = c(rep("NULL", 6), "character", "NULL",
                                                                            rep("character", 10)),
                                                             sep="\t", header=TRUE, strip.white=TRUE,
                                                             na.strings=TRUE, fill=TRUE, nrows=101))
aa_sequence_file <- do.call("rbind", ds_files5)
colnames(aa_sequence_file)
# [1] "V.D.J.REGION" "V.REGION" "FR1.IMGT" "CDR1.IMGT"
# "FR2.IMGT" "CDR2.IMGT" "FR3.IMGT"
# [8] "CDR3.IMGT" "JUNCTION" "J.REGION" "FR4.IMGT" "Donor"
files8 <- c("code/IMGT/dataset2_pandaseq_aa/8_V-REGION-nt-mutation-statistics.txt")

ds_files8 <- lapply(1:length(files8), function(x) read.table(files8[x],
                                                             colClasses = c(rep("NULL", 7), rep("character",3), rep("NULL", 15),
                                                                            rep("character", 3), rep("NULL", 15), rep("character", 3),
                                                                            rep("NULL", 15), rep("character", 3), rep("NULL", 15),
                                                                            rep("character", 3), rep("NULL", 15), rep("character", 3),
                                                                            rep("NULL", 15), rep("character",3), rep("NULL", 12)),
                                                             sep="\t", header=TRUE, strip.white=TRUE,
                                                             na.strings=TRUE, fill=TRUE, nrows=500001))
nt_mutation_file <- do.call("rbind", ds_files8)
colnames(nt_mutation_file)
# [1] "V.REGION.Nb.of.mutations" "V.REGION.Nb.of.silent.mutations"
# "V.REGION.Nb.of.nonsilent.mutations" [4] "FR1.IMGT.Nb.of.mutations"
# "FR1.IMGT.Nb.of.silent.mutations" "FR1.IMGT.Nb.of.nonsilent.mutations"
# [7] "CDR1.IMGT.Nb.of.mutations" "CDR1.IMGT.Nb.of.silent.mutations"
# "CDR1.IMGT.Nb.of.nonsilent.mutations" [10] "FR2.IMGT.Nb.of.mutations"
# "FR2.IMGT.Nb.of.silent.mutations" "FR2.IMGT.Nb.of.nonsilent.mutations"
# [13] "CDR2.IMGT.Nb.of.mutations" "CDR2.IMGT.Nb.of.silent.mutations"
# "CDR2.IMGT.Nb.of.nonsilent.mutations" [16] "FR3.IMGT.Nb.of.mutations"
# "FR3.IMGT.Nb.of.silent.mutations" "FR3.IMGT.Nb.of.nonsilent.mutations"
# [19] "CDR3.IMGT.Nb.of.mutations" "CDR3.IMGT.Nb.of.silent.mutations"
# "CDR3.IMGT.Nb.of.nonsilent.mutations"
files9 <- c("code/IMGT/dataset2_pandaseq_aa/9_V-REGION-AA-change-statistics.txt")

ds_files9 <- lapply(1:length(files9), function(x) read.table(files9[x],
                                                             colClasses = c(rep("NULL", 5), "character", "NULL", "character",
                                                                            rep("NULL", 12), "character", "NULL", "character", rep("NULL", 12),
                                                                            "character", "NULL", "character", rep("NULL", 12), "character",
                                                                            "NULL", "character", rep("NULL", 12), "character", "NULL",
                                                                            "character", rep("NULL", 12), "character", "NULL", "character",
                                                                            rep("NULL", 12), "character", "NULL", "character", rep("NULL", 11)),
                                                             sep="\t", header=TRUE, strip.white=TRUE,
                                                             na.strings=TRUE, fill=TRUE, nrows=101))

aa_change_file <- do.call("rbind", ds_files9)

colnames(aa_change_file)
# [1] "V.REGION.Nb.of.AA" "V.REGION.Nb.of.AA.changes" "FR1.IMGT.Nb.of.AA"
# [4] "FR1.IMGT.Nb.of.AA.changes" "CDR1.IMGT.Nb.of.AA" "CDR1.IMGT.Nb.of.AA.changes"
# [7] "FR2.IMGT.Nb.of.AA" "FR2.IMGT.Nb.of.AA.changes" "CDR2.IMGT.Nb.of.AA"
# [10] "CDR2.IMGT.Nb.of.AA.changes" "FR3.IMGT.Nb.of.AA" "FR3.IMGT.Nb.of.AA.changes"
# [13] "CDR3.IMGT.Nb.of.AA" "CDR3.IMGT.Nb.of.AA.changes"
### Extracted parameters from IMGT files
# nt_sequences_files
# aa_sequence_files
# nt_mutation_files
# aa_change_files
# Bind the the aa files (all rows together) to the "Donor" column
# that indicates the donor for each sequence-row
imgt_files <- cbind(summary_file, aa_sequence_file, aa_change_file, Donor)


### Save your data for faster and access without re-runing read-in
save(imgt_files, file = "imgt_files_mf.RData")

#!



#!
################### FILTERING #####################
###################################################

# Start timing of the code
start <- proc.time()

# dataset has 3,35 and 3,10 mln healthy and hiv sequences
load("code/imgt_files.RData")

# Match samples. Check the function by running: unique(imgt_files$Donor)
# and levels(factor(imgt_files$Donor))
match_stage <- match(unique(imgt_files$Donor), levels(factor(imgt_files$Donor)))
# Count annotated total sequences per sample
annotated_seq_initial <- data.frame(summary(imgt_files$Donor)[match_stage])

# Select columns
selection <- imgt_files[,c("Functionality", "V.GENE.and.allele",
                           "D.GENE.and.allele", "J.GENE.and.allele", "V.D.J.REGION",
                           "V.REGION", "JUNCTION", "J.REGION", "FR1.IMGT", "FR2.IMGT",
                           "FR3.IMGT", "FR4.IMGT", "CDR1.IMGT", "FR1.IMGT","CDR2.IMGT",
                           "CDR3.IMGT", "V.REGION.Nb.of.AA", "V.REGION.Nb.of.AA.changes",
                           "FR1.IMGT.Nb.of.AA.changes", "FR2.IMGT.Nb.of.AA.changes",
                           "FR3.IMGT.Nb.of.AA.changes", "CDR1.IMGT.Nb.of.AA.changes",
                           "CDR2.IMGT.Nb.of.AA.changes", "CDR3.IMGT.Nb.of.AA.changes",
                           "Donor")]

# How many NA in data?
sum(is.na(selection))
# [1] 0

# How many sequences in total?
nrow(selection)
# [1] 6455246

selection_df <- data.frame(summary(selection$Donor)[match_stage])
# Select for productive reads according to IMGT functionality
selection_productive <- selection$Functionality=="productive"
productive <- selection[(selection_productive),] # nrow(productive) [1] 195242
sum(is.na(productive))
# [1] 0
nrow(productive)
#[1] 441540
productive_df <- data.frame(summary(productive$Donor)[match_stage])

# Add lengths of reagions as additional columns of the data frame
productive$V.REGION.Length <- nchar(productive[,"V.REGION"])
productive$V.D.J.REGION.Length <- nchar(productive[,"V.D.J.REGION"])
productive$JUNCTION.Length <- nchar(productive[,"JUNCTION"])
productive$J.REGION.Length <- nchar(productive[,"J.REGION"])
productive$FR1.Length <- nchar(productive[,"FR1.IMGT"])
productive$FR2.Length <- nchar(productive[,"FR2.IMGT"])
productive$FR3.Length <- nchar(productive[,"FR3.IMGT"])
productive$FR4.Length <- nchar(productive[,"FR4.IMGT"])
productive$CDR1.Length <- nchar(productive[,"CDR1.IMGT"])
productive$CDR2.Length <- nchar(productive[,"CDR2.IMGT"])
productive$CDR3.Length <- nchar(productive[,"CDR3.IMGT"])


# Extract number-s of AA and AA changed (mutations) from IMGT columns
productive$V.REGION.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"V.REGION.Nb.of.AA.changes"], "[0-9]{1,3}"))
productive$FR1.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"FR1.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
productive$FR2.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"FR2.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
    
productive$FR3.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"FR3.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
productive$CDR1.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"CDR1.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
productive$CDR2.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"CDR2.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
productive$CDR3.IMGT.Nb.of.AA.changes <- as.numeric(
    str_extract(productive[,"CDR3.IMGT.Nb.of.AA.changes"], "[0-9]{1,3}"))
# Extract V(D)J genes
V.Gene.and.allele <- str_extract(productive$V.GENE.and.allele,
                                 pattern = "[[:alnum:]]{1,6}-(.*?)\\*[0-9]{1,2}")
productive$V.Gene.and.allele <- V.Gene.and.allele
# Check output of extraction
sample(table(V.Gene.and.allele),3)
# V.Gene.and.allele
# IGHV3-30*10 IGHV7-4-1*02 IGLV1-51*02
# 6 228 5021
V.Gene.Subgroup <- str_extract(productive$V.GENE.and.allele,
                               pattern = "[[:alnum:]]{1,6}-(.*?)[[:alnum:]]{1,2}")
productive$V.Gene.Subgroup <- V.Gene.Subgroup
V.Gene <- str_extract(productive$V.GENE.and.allele,
                      pattern = "IG[A-Z]V[0-9]{1,2}")
productive$V.Gene <- V.Gene
J.Gene.and.allele <- str_extract(productive$J.GENE.and.allele,
                                 pattern = "[[:alnum:]]{1,6}\\*[0-9]{1,2}")
productive$J.Gene.and.allele <- J.Gene.and.allele
J.Gene <- str_extract(productive$J.GENE.and.allele,
                      pattern = "IG[A-Z]J[0-9]{1,2}")
productive$J.Gene <- J.Gene
productive2 <- productive

# From productive2 variable, extract only heavy chains (productive)
V.Gene.and.allele <- str_extract(productive$V.GENE.and.allele, "IG[H]V")
productive$HV.Gene.and.allele <- V.Gene.and.allele

# From productive2 variable, extract only light chains (productive2)
V.Gene.and.allele <- str_extract(productive2$V.GENE.and.allele, "IG[K-L]V")
productive2$V.Gene.and.allele <- V.Gene.and.allele
productive2$HV.Gene.and.allele <- NULL
row.has.na2 <- apply(productive2, 1, function(x){any(is.na(x))})

sum(row.has.na2)
# [1] 291131
# Take out the rows with NA (from VH assignment)
productive_lv <- productive2[!row.has.na2,]
nrow(productive_lv)
# [1] 150409
nrow(productive)
# [1] 441540
row.has.na <- apply(productive, 1, function(x){any(is.na(x))})
productive_hv <- productive[!row.has.na,]
nrow(productive_hv) # Number of productive heavy chains
# [1] 73474
# Extract D gene
D.Gene.and.allele <- str_extract(productive_hv$D.GENE.and.allele,
                                 pattern = "[[:alnum:]]{1,6}-(.*?)\\*[0-9]{1,2}")
productive_hv$D.Gene.and.allele <- D.Gene.and.allele
D.Gene.Subgroup <- str_extract(productive_hv$D.GENE.and.allele,
                               pattern = "[[:alnum:]]{1,6}-(.*?)[[:alnum:]]{1,2}")
productive_hv$D.Gene.Subgroup <- D.Gene.Subgroup
D.Gene <- str_extract(productive_hv$D.GENE.and.allele,
                      pattern = "IG[A-Z]D[0-9]{1,2}")
productive_hv$D.Gene <- D.Gene
# Select columns after V(D)J extraction
extracted_hc <- productive_hv[,c("Functionality", "V.Gene.and.allele",
                                 "V.Gene.Subgroup", "V.Gene", "D.Gene.and.allele",
                                 "D.Gene.Subgroup", "D.Gene", "J.Gene.and.allele",
                                 "J.Gene", "V.REGION", "V.D.J.REGION", "JUNCTION",
                                 "J.REGION", "FR1.IMGT", "FR2.IMGT", "FR3.IMGT",
                                 "FR4.IMGT", "CDR1.IMGT", "CDR2.IMGT", "CDR3.IMGT",
                                 "V.REGION.Nb.of.AA.changes", "FR1.IMGT.Nb.of.AA.changes",
                                 "FR2.IMGT.Nb.of.AA.changes", "FR3.IMGT.Nb.of.AA.changes",
                                 "CDR1.IMGT.Nb.of.AA.changes", "CDR2.IMGT.Nb.of.AA.changes",
                                 "CDR3.IMGT.Nb.of.AA.changes", "Donor")]
extracted_lc <- productive_lv[,c("Functionality", "V.Gene.and.allele",
                                 "V.Gene.Subgroup", "V.Gene", "J.Gene.and.allele", "J.Gene",
                                 "V.REGION", "V.D.J.REGION", "JUNCTION", "J.REGION", "FR1.IMGT",
                                 "FR2.IMGT", "FR3.IMGT", "FR4.IMGT", "CDR1.IMGT", "CDR2.IMGT",
                                 "CDR3.IMGT", "V.REGION.Nb.of.AA.changes", "FR1.IMGT.Nb.of.AA.changes",
                                 "FR2.IMGT.Nb.of.AA.changes", "FR3.IMGT.Nb.of.AA.changes",
                                 "CDR1.IMGT.Nb.of.AA.changes", "CDR2.IMGT.Nb.of.AA.changes",
                                 "CDR3.IMGT.Nb.of.AA.changes", "Donor")]

# Find rows with NA
row.has.na_hc <- apply(extracted_hc, 1, function(x){any(is.na(x))})
# How many rows have NA?
sum(row.has.na_hc)
# [1] 487
# To check the NAs in the data:
head(extracted_hc[row.has.na_hc,])
# Note that there was no D gene assigned, where the NA are.
# Exclude (!) rows with NA in the heavy chains data frame
nona_hc <- extracted_hc[!row.has.na_hc,]
nrow(nona_hc) # [1] 72987
# Exclude (!) rows with NA in the light chains data frame
row.has.na_lc <- apply(extracted_lc, 1, function(x){any(is.na(x))})
nona_lc <- extracted_lc[!row.has.na_lc,]
nrow(nona_lc) # [1] 150409
### Some data present "*" in annotated sequences, which are stop codons.
### To filter out VDJ with '*'
# s <- nona_hc$V.D.J.REGION
# my_pattern='\\*'
# strings_without_asterisk <- s[!grep(my_pattern, s)]
# with_asterisk <- grep(my_pattern, s)
# selection <- initial_filtered[-c(with_asterisk),]
### Some data present "X" in an annotated sequences, which indicates
### undetermined aminoacid.
### Exclusion of sequences that have a "X" in the V region or in the CDR3 sequence
# nona_c3X <- str_detect(nona_hc$CDR3.IMGT, pattern = c("X")) #sum(nona_c3X) # [1] 0
# initial <- nona_hc[!nona_c3X,]
initial_hc <- nona_hc
initial_lc <- nona_lc

# Additional CDR3.Length for heavy chains
initial_hc$FR1.Length <- nchar(initial_hc$FR1.IMGT)
initial_hc$FR2.Length <- nchar(initial_hc$FR2.IMGT)
initial_hc$FR3.Length <- nchar(initial_hc$FR3.IMGT)
initial_hc$FR4.Length <- nchar(initial_hc$FR4.IMGT)
initial_hc$CDR1.Length <- nchar(initial_hc$CDR3.IMGT)
initial_hc$CDR2.Length <- nchar(initial_hc$CDR3.IMGT)
initial_hc$CDR3.Length <- nchar(initial_hc$CDR3.IMGT)

# Additional CDR3.Length for light chains
initial_lc$FR1.Length <- nchar(initial_lc$FR1.IMGT)
initial_lc$FR2.Length <- nchar(initial_lc$FR2.IMGT)
initial_lc$FR3.Length <- nchar(initial_lc$FR3.IMGT)
initial_lc$FR4.Length <- nchar(initial_lc$FR4.IMGT)


initial_lc$CDR1.Length <- nchar(initial_lc$CDR3.IMGT)
initial_lc$CDR2.Length <- nchar(initial_lc$CDR3.IMGT)
initial_lc$CDR3.Length <- nchar(initial_lc$CDR3.IMGT)
# Numbers of preprocessed sequences at the different stages
selection_df <- data.frame(summary(selection$Donor)[match_stage])
productive_hc_df <- data.frame(summary(productive_hv$Donor)[match_stage])
productive_lc_df <- data.frame(summary(productive_lv$Donor)[match_stage])
initial_hc_df <- data.frame(summary(initial_hc$Donor)[match_stage])
initial_lc_df <- data.frame(summary(initial_lc$Donor)[match_stage])

# Counts of CDR3 per donor. Apply table funtion, to CDR3 sequences according donor
# Calculates the frequency of unique CDR3s
cdr3_count_pre_hc <- lapply(tapply(initial_hc$CDR3.IMGT,
                                   initial_hc$Donor, function(x)
                                       sort(table(x), decreasing = TRUE)), function(x) x)
cdr3_count_pre_lc <- lapply(tapply(initial_lc$CDR3.IMGT, initial_lc$Donor,
                                   function(x) sort(table(x), decreasing = TRUE)), function(x) x)

# How many unique CDR3 are in list 1?
length(cdr3_count_pre_hc[[1]])
# [1] 4264

# How many unique CDR3 are each donor?
sapply(cdr3_count_pre_hc, length)[match_stage]
# HIV-1 IAVI donor 17 Healthy donor 1
# 3161 4264
sapply(cdr3_count_pre_lc, length)[match_stage]
# HIV-1 IAVI donor 17 Healthy donor 1
# 4203 8093

## Add for table
cdr3_count_pr_hc <- sapply(cdr3_count_pre_hc, length)
cdr3_count_pr_lc <- sapply(cdr3_count_pre_lc, length)
cdr3_count_pr_hc_df <- data.frame(cdr3_count_pr_hc[match_stage])
cdr3_count_pr_lc_df <- data.frame(cdr3_count_pr_lc[match_stage])

# Include only CDR3 that appear more than once in data (exclude singletons)
cdr3_exclusion_hc <- lapply(1:length(cdr3_count_pre_hc), function(x)
    names(which(cdr3_count_pre_hc[[x]]>1))) #

# How many CDR3 are there after exclusion of singletons in heavy chains?
sapply(cdr3_exclusion_hc, length)
# [1] 1862 1557
cdr3_exclusion_lc <- lapply(1:length(cdr3_count_pre_lc), function(x)
    names(which(cdr3_count_pre_lc[[x]]>1)))

sapply(cdr3_exclusion_lc, length)
# [1] 3773 2377
# Filter-in from the "initial" stage only the CDR3 that have passed checks
cdr3_exclusion_final_hc <- unlist(lapply(1:length(levels(initial_hc$Donor)), function(x)
    initial_hc$CDR3.IMGT[initial_hc$Donor==levels(initial_hc$Donor)[x]] %in%
        cdr3_exclusion_hc[[x]])[match_stage])
length(cdr3_exclusion_final_hc)
# [1] 72987
cdr3_exclusion_final_lc <- unlist(lapply(1:length(levels(initial_lc$Donor)), function(x)
    initial_lc$CDR3.IMGT[initial_lc$Donor==levels(initial_lc$Donor)[x]] %in%
        cdr3_exclusion_lc[[x]])[match_stage])

length(cdr3_exclusion_final_hc)
# [1] 72987

cdr3_included_hc <- initial_hc[cdr3_exclusion_final_hc,]
length(cdr3_included_hc[[1]])
# [1] 68981

tapply(cdr3_included_hc$CDR3.IMGT, cdr3_included_hc$Donor, length)
# Healthy donor 1 HIV-1 IAVI donor 17
# 31800 37181

cdr3_included_lc <- initial_lc[cdr3_exclusion_final_lc,]
tapply(cdr3_included_lc$CDR3.IMGT, cdr3_included_lc$Donor, length)
# Healthy donor 1 HIV-1 IAVI donor 17
# 60094 84169

## Add for table, exclude singletons (CDR3 that appear once in data)
cdr3_included_hc_df <- data.frame(summary(cdr3_included_hc$Donor)[match_stage])
cdr3_included_lc_df <- data.frame(summary(cdr3_included_lc$Donor)[match_stage])
# Filter for CDR3 that are longer than 2 a.a.
cdr3_length_hc <- nchar(str_trim(cdr3_included_hc$CDR3.IMGT))>2
length(cdr3_length_hc)
#[1] 68981

cdr3_length_lc <- nchar(str_trim(cdr3_included_lc$CDR3.IMGT))>2
length(cdr3_length_lc)
# [1] 144263

# All CDR3s that have passed filters but not unique
cdr3_hc <- cdr3_included_hc[cdr3_length_hc,]
tapply(cdr3_hc$CDR3.IMGT, cdr3_hc$Donor, length)

# Healthy donor 1 HIV-1 IAVI donor 17
# 31800 37181
cdr3_lc <- cdr3_included_lc[cdr3_length_lc,]

#including only cdr3>2aa
cdr3_hc_df <- data.frame(summary(cdr3_hc$Donor)[match_stage])
cdr3_lc_df <- data.frame(summary(cdr3_lc$Donor)[match_stage])

# Preprocess in unique CDR3s
cdr3_processed <- lapply(tapply(cdr3_hc$CDR3.IMGT, cdr3_hc$Donor,
                                function(x) sort(table(x), decreasing = TRUE)),
                         function(y) y)
# How many unique CDR3s in each dataset?
sapply(cdr3_processed, length)

# Healthy donor 1 HIV-1 IAVI donor 17
# 1862 1557
cdr3_hc <- cdr3_hc[,c("Functionality", "V.Gene.and.allele",
                      "V.Gene.Subgroup", "V.Gene", "J.Gene.and.allele",
                      "J.Gene", "V.REGION", "V.D.J.REGION", "JUNCTION",
                      "J.REGION", "FR1.IMGT","FR2.IMGT", "FR3.IMGT",
                      "FR4.IMGT", "CDR1.IMGT", "CDR2.IMGT", "CDR3.IMGT",
                      "V.REGION.Nb.of.AA.changes", "FR1.IMGT.Nb.of.AA.changes",
                      "FR2.IMGT.Nb.of.AA.changes", "FR3.IMGT.Nb.of.AA.changes",
                      "CDR1.IMGT.Nb.of.AA.changes", "CDR2.IMGT.Nb.of.AA.changes",
                      "CDR3.IMGT.Nb.of.AA.changes", "FR1.Length", "FR2.Length",
                      "FR3.Length", "FR4.Length", "CDR1.Length", "CDR2.Length",
                      "CDR3.Length", "Donor")]
nrow(cdr3_hc)
# [1] 68981
## Add for table
cdr3_df <- data.frame(summary(cdr3_hc$Donor)[match_stage])
cdr3_processed_df <- data.frame(summary(cdr3_processed)[match_stage])
productive_hc_df
# Build table of numbers of sequences through preprocessing steps
processing_seq <- cbind(selection_df, productive_df, productive_hc_df,
                        initial_hc_df, cdr3_included_hc_df,
                        cdr3_df, cdr3_processed_df)
# Set names of columns
colnames(processing_seq) <- c("All Paired", "All Productive", "Productive HC",
                              "NO Missing, X", "NO Singletons", "CDR3 > 2a.a.",
                              "Unique")
### Write a CSV file of the preprocessed results
write.csv(processing_seq, "processing_seq.csv")
### Make a plot
seqs_df <- data.frame(rownames(processing_seq), processing_seq)

# Eliminate rownames
rownames(seqs_df) <- NULL

# Select columns to show in plot
seqs_df <- seqs_df[,c(1,2,4,8)]
colnames(seqs_df)[1] <- c("Sample")

# Melt data
melt_seq <- melt(seqs_df, id.vars="Sample")
print(melt_seq)
# Sample variable value
# 1 HIV-1 IAVI donor 17 All.Paired 3104454
# 2 Healthy donor 1 All.Paired 3350792
# 3 HIV-1 IAVI donor 17 Productive.HC 38994
# 4 Healthy donor 1 Productive.HC 34480
# 5 HIV-1 IAVI donor 17 Unique 1557
# 6 Healthy donor 1 Unique 1862

plot_seq <- ggplot(melt_seq, aes(x=variable, y=as.numeric(value)))
plot_seq <- plot_seq + geom_bar(stat="identity", aes(fill=variable), position="dodge") +
    facet_wrap(~Sample, nrow=3) +
    theme_bw() +
    labs(x= "Processing stages", y="Sequence count (log10)") +
    theme(plot.title = element_text(face="bold", size=rel(1), hjust=0),
          plot.margin = unit(c(2, 2, 2, 2),"points"),
          strip.text = element_text(size = rel(0.75))) +
    theme(strip.background = element_rect(fill = scales:::alpha("blue", 0.3))) +
    scale_fill_discrete(name="Processing stage") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
# pdf("figure/processed_seqs.pdf", width=7, heigh=7)
plot_seq
# dev.off()
# Prints a table that can be imported in LaTeX
print(xtable(processing_seq, caption = "Processed sequences",
             label="tab:process-seqs", digits = 2),
      table.placement="h", include.rownames=TRUE, size = "tiny")
# Stop the clock (in seconds). Code takes ~ 2 minutes filtering.
proc.time() - start
# user system elapsed
# 105.194 7.846 119.172

#!
################# 6.2 ROCESSING ###############
###################################################
# Unique CDR3 are matched to mutations, germline genes (V, D, J) and regions (FRs, CDRs, V region).
# Results are saved in .RData files that can be loaded and continue other analysis and visualization.
process_time <- proc.time()

# Preprocess heavy chains
# Repeat for light chains with: cdr3 <- cdr3_lc
cdr3 <- cdr3_hc

# Mutations by CDR3
mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.REGION.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

# Generate CDR3 table in order to match with CDR3 mutation-table:
# mutations_by_cdr3_pre
# Match by CDR3
cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    # order in descending order by CDR3 abundance
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

mutations_by_cdr3 <- lapply(1:length(mutations_by_cdr3_pre), function(x) {
    mutations <- match(cdr3_table[[x]], names(mutations_by_cdr3_pre[[x]]))
    mutations_by_cdr3_pre[[x]][mutations]})

mutations_by_cdr3_stage <- mutations_by_cdr3[match_stage]

# FR1 mutations
fr1_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function (y)
    tapply(cdr3$FR1.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

fr1_mutations_by_cdr3 <- lapply(1:length(fr1_mutations_by_cdr3_pre), function(x) {
    fr1_mutations <- match(cdr3_table[[x]], names(fr1_mutations_by_cdr3_pre[[x]]))
    fr1_mutations_by_cdr3_pre[[x]][fr1_mutations]})

fr1_mutations_by_cdr3_stage <- fr1_mutations_by_cdr3[match_stage]

# FR2 mutations
fr2_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$FR2.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

fr2_mutations_by_cdr3 <- lapply(1:length(fr2_mutations_by_cdr3_pre), function(x) {
    fr2_mutations <- match(cdr3_table[[x]], names(fr2_mutations_by_cdr3_pre[[x]]))
    fr2_mutations_by_cdr3_pre[[x]][fr2_mutations]})

fr2_mutations_by_cdr3_stage <- fr2_mutations_by_cdr3[match_stage]

# FR3 mutations
fr3_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$FR3.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

fr3_mutations_by_cdr3 <- lapply(1:length(fr3_mutations_by_cdr3_pre), function(x) {
    fr3_mutations <- match(cdr3_table[[x]], names(fr3_mutations_by_cdr3_pre[[x]]))
    fr3_mutations_by_cdr3_pre[[x]][fr3_mutations]})

fr3_mutations_by_cdr3_stage <- fr3_mutations_by_cdr3[match_stage]

# CDR1 mutations
cdr1_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$CDR1.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

cdr1_mutations_by_cdr3 <- lapply(1:length(cdr1_mutations_by_cdr3_pre), function(x) {
    cdr1_mutations <- match(cdr3_table[[x]], names(cdr1_mutations_by_cdr3_pre[[x]]))
    cdr1_mutations_by_cdr3_pre[[x]][cdr1_mutations]})

cdr1_mutations_by_cdr3_stage <- cdr1_mutations_by_cdr3[match_stage]

# CDR2 mutations
cdr2_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$CDR2.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

cdr2_mutations_by_cdr3 <- lapply(1:length(cdr2_mutations_by_cdr3_pre), function(x) {
    cdr2_mutations <- match(cdr3_table[[x]], names(cdr2_mutations_by_cdr3_pre[[x]]))
    cdr2_mutations_by_cdr3_pre[[x]][cdr2_mutations]})

cdr2_mutations_by_cdr3_stage <- cdr2_mutations_by_cdr3[match_stage]

# CDR3 MUTATIONS
cdr3_mutations_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$CDR3.IMGT.Nb.of.AA.changes[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               median_mutations <- median(as.numeric(x), na.rm = TRUE)
               median_mutations
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

cdr3_mutations_by_cdr3 <- lapply(1:length(cdr3_mutations_by_cdr3_pre), function(x) {
    cdr3_mutations <- match(cdr3_table[[x]], names(cdr3_mutations_by_cdr3_pre[[x]]))
    cdr3_mutations_by_cdr3_pre[[x]][cdr3_mutations]})

cdr3_mutations_by_cdr3_stage <- cdr3_mutations_by_cdr3[match_stage]

# CDR3 abundance
cdr3_abundance <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    abundance <- table(x)
    # order in descending order by CDR3 abundance
    abundance <- abundance[order(abundance, decreasing = TRUE)]
    abundance
})
cdr3_abundance_stage <- cdr3_abundance[match_stage]

# CDR3 frequency
cdr3_frequency <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    abundance <- table(x)
    abundance <- abundance[order(abundance, decreasing = TRUE)]
    frequency <- abundance/sum(abundance) *100
    frequency
})

cdr3_frequency_stage <- cdr3_frequency[match_stage]

# V region by CDR3 and donor
vregion_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.REGION[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x) {
               # cat(.) # Uncomment this to see a dot printed for each list processed
               v.region <- table(x)
               names(v.region)[which.max(v.region)]
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

vregion_by_cdr3 <- lapply(1:length(vregion_by_cdr3_pre), function(x) {
    vregion <- match(cdr3_table[[x]], names(vregion_by_cdr3_pre[[x]]))
    vregion_by_cdr3_pre[[x]][vregion]})

vregion_by_cdr3_stage <- vregion_by_cdr3[match_stage]

# VDJ region by CDR3 and donor
vdjregion_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.D.J.REGION[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               # cat(.)
               vdj.region <- table(x)
               names(vdj.region)[which.max(vdj.region)]
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

vdjregion_by_cdr3 <- lapply(1:length(vdjregion_by_cdr3_pre), function(x) {
    vdjregion <- match(cdr3_table[[x]], names(vdjregion_by_cdr3_pre[[x]]))
    vdjregion_by_cdr3_pre[[x]][vdjregion]})

vdjregion_by_cdr3_stage <- vdjregion_by_cdr3[match_stage]

# FR1 by CDR3 and donor
fr1_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$FR1.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               fr1.region <- table(x)
               names(fr1.region)[which.max(fr1.region)]
           }))

cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    # order in descending order by CDR3 abundance
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})

fr1_by_cdr3 <- lapply(1:length(fr1_by_cdr3_pre), function(x) {
    fr1 <- match(cdr3_table[[x]], names(fr1_by_cdr3_pre[[x]]))
    fr1_by_cdr3_pre[[x]][fr1]})
fr1_by_cdr3_stage <- fr1_by_cdr3[match_stage]
# FR2 by CDR3 and donor
fr2_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$FR2.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               fr2.region <- table(x)
               names(fr2.region)[which.max(fr2.region)]
           }))
cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})
fr2_by_cdr3 <- lapply(1:length(fr2_by_cdr3_pre), function(x) {
    fr2 <- match(cdr3_table[[x]], names(fr2_by_cdr3_pre[[x]]))
    fr2_by_cdr3_pre[[x]][fr2]})
fr2_by_cdr3_stage <- fr2_by_cdr3[match_stage]
# FR3 BY CDR3 AND DONOR
fr3_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$FR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               # cat(.)
               fr3.region <- table(x)
               names(fr3.region)[which.max(fr3.region)]
           }))
cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})
fr3_by_cdr3 <- lapply(1:length(fr3_by_cdr3_pre), function(x) {
    fr3 <- match(cdr3_table[[x]], names(fr3_by_cdr3_pre[[x]]))
    fr3_by_cdr3_pre[[x]][fr3]})
fr3_by_cdr3_stage <- fr3_by_cdr3[match_stage]
# CDR1 by CDR3 and donor
cdr1_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$CDR1.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               cdr1.region <- table(x)
               names(cdr1.region)[which.max(cdr1.region)]
           }))
cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})
cdr1_by_cdr3 <- lapply(1:length(cdr1_by_cdr3_pre), function(x) {
    cdr1 <- match(cdr3_table[[x]], names(cdr1_by_cdr3_pre[[x]]))
    cdr1_by_cdr3_pre[[x]][cdr1]})
cdr1_by_cdr3_stage <- cdr1_by_cdr3[match_stage]
# CDR2 by CDR3 and donor
cdr2_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$CDR2.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               cdr2.region <- table(x)
               names(cdr2.region)[which.max(cdr2.region)]
           }))
cdr3_table <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    names(cdr3.table)})
cdr2_by_cdr3 <- lapply(1:length(cdr2_by_cdr3_pre), function(x) {
    cdr2 <- match(cdr3_table[[x]], names(cdr2_by_cdr3_pre[[x]]))
    cdr2_by_cdr3_pre[[x]][cdr2]})
cdr2_by_cdr3_stage <- cdr2_by_cdr3[match_stage]
# V gene assignment to unique CDR3 analysed
vgene_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.Gene[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               v.gene <- table(x)
               names(v.gene)[which.max(v.gene)]
           }))
vgene_by_cdr3 <- lapply(1:length(vgene_by_cdr3_pre), function(x) {
    table.vgene <- match(cdr3_table[[x]], names(vgene_by_cdr3_pre[[x]]))
    vgene_by_cdr3_pre[[x]][table.vgene]
})
vgene_by_cdr3_stage <- vgene_by_cdr3[match_stage]
# V gene subgroup assignment to unique CDR3 analysed
vgene_subgr_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.Gene.Subgroup[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               v.gene.sub <- table(x)
               names(v.gene.sub)[which.max(v.gene.sub)]
           }))
vgene_subgr_by_cdr3 <- lapply(1:length(vgene_subgr_by_cdr3_pre), function(x) {
    match_by_cdr3 <- match(cdr3_table[[x]], names(vgene_subgr_by_cdr3_pre[[x]]))
    vgene_subgr_by_cdr3_pre[[x]][match_by_cdr3]
})
vgene_subgr_by_cdr3_stage <- vgene_subgr_by_cdr3[match_stage]
# V gene and allele assignment to unique CDR3 analysed
vgene_allele_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$V.Gene.and.allele[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               v.gene.allele <- table(x)
               names(v.gene.allele)[which.max(v.gene.allele)]
           }))
vgene_allele_by_cdr3 <- lapply(1:length(vgene_allele_by_cdr3_pre), function(x) {
    match_by_cdr3 <- match(cdr3_table[[x]], names(vgene_allele_by_cdr3_pre[[x]]))
    vgene_allele_by_cdr3_pre[[x]][match_by_cdr3]
})
vgene_allele_by_cdr3_stage <- vgene_allele_by_cdr3[match_stage]
# J gene assignment to unique cdr3s analysed
jgene_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$J.Gene[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               v.gene <- table(x)
               names(v.gene)[which.max(v.gene)]
           }))
jgene_by_cdr3 <- lapply(1:length(jgene_by_cdr3_pre), function(x) {
    jgene <- match(cdr3_table[[x]], names(jgene_by_cdr3_pre[[x]]))
    jgene_by_cdr3_pre[[x]][jgene]
})
jgene_by_cdr3_stage <- jgene_by_cdr3[match_stage]


# J gene and allele assignment to unique cdr3s analysed
jgene_allele_by_cdr3_pre <- lapply(1:length(levels(cdr3$Donor)), function(y)
    tapply(cdr3$J.Gene.and.allele[cdr3$Donor==levels(cdr3$Donor)[y]],
           cdr3$CDR3.IMGT[cdr3$Donor==levels(cdr3$Donor)[y]], function(x){
               j.gene.allele <- table(x)
               names(j.gene.allele)[which.max(j.gene.allele)]
           }))

jgene_allele_by_cdr3 <- lapply(1:length(jgene_allele_by_cdr3_pre), function(x) {
    j.gene.allele <- match(cdr3_table[[x]], names(jgene_allele_by_cdr3_pre[[x]]))
    jgene_allele_by_cdr3_pre[[x]][j.gene.allele]
})

jgene_allele_by_cdr3_stage <- jgene_allele_by_cdr3[match_stage]

# CDR3 rank
cdr3_table_rank <- tapply(cdr3$CDR3.IMGT, cdr3$Donor, function(x) {
    cdr3.table <- table(x)
    # Order in descending order by CDR3 abundance
    cdr3.table <- cdr3.table[order(cdr3.table, decreasing = TRUE)]
    cdr3.table/sum(cdr3.table) *100
})

cdr3_table_rank_stage <- cdr3_table_rank[match_stage]

save(cdr3_table_rank_stage, file="cdr3_table_rank_stage.RData")



###############################################################
#################### PARAMETERS BY CDR3 AA ####################
# Put all the values into a list
working_list <- list()
for(k in 1:length(cdr3_table_rank_stage)){
    working_list[[k]] <- data.frame(cdr3 = names(cdr3_table_rank_stage[[k]]),
                                    cdr3_rank = c(1:length(cdr3_table_rank_stage[[k]])),
                                    cdr3_frequency = as.numeric(cdr3_frequency_stage[[k]]),
                                    cdr3_abundance = as.numeric(cdr3_abundance_stage[[k]]),
                                    vregion_mutations = mutations_by_cdr3_stage[[k]],
                                    fr1_mutations = fr1_mutations_by_cdr3_stage[[k]],
                                    fr2_mutations = fr2_mutations_by_cdr3_stage[[k]],
                                    fr3_mutations = fr3_mutations_by_cdr3_stage[[k]],
                                    cdr1_mutations = cdr1_mutations_by_cdr3_stage[[k]],
                                    cdr2_mutations = cdr2_mutations_by_cdr3_stage[[k]],
                                    cdr3_mutations = cdr3_mutations_by_cdr3_stage[[k]],
                                    vgeneallele = vgene_allele_by_cdr3_stage[[k]],
                                    vgene = vgene_by_cdr3_stage[[k]],
                                    vgene_subgroup = vgene_subgr_by_cdr3_stage[[k]],
                                    jgeneallele = jgene_allele_by_cdr3_stage[[k]],
                                    jgene = jgene_by_cdr3_stage[[k]],
                                    vregion_by_cdr3 = vregion_by_cdr3_stage[[k]],
                                    fr1_by_cdr3 = fr1_by_cdr3_stage[[k]],
                                    fr2_by_cdr3 = fr2_by_cdr3_stage[[k]],
                                    fr3_by_cdr3 = fr3_by_cdr3_stage[[k]],
                                    cdr1_by_cdr3 = cdr1_by_cdr3_stage[[k]],
                                    cdr2_by_cdr3 = cdr2_by_cdr3_stage[[k]],
                                    cdr1_length = nchar(cdr1_by_cdr3_stage[[k]]),
                                    cdr2_length = nchar(cdr2_by_cdr3_stage[[k]]),
                                    cdr3_length = nchar(names(cdr3_table_rank_stage[[k]])),
                                    fr1_length = nchar(fr1_by_cdr3_stage[[k]]),
                                    fr2_length = nchar(fr2_by_cdr3_stage[[k]]),
                                    fr3_length = nchar(fr3_by_cdr3_stage[[k]]))}

names(working_list) <- names(cdr3_table_rank_stage)

# Stop the clock (in seconds). Code takes less than 2 minutes for processing.
proc.time() - process_time
# user system elapsed
# 17.364 0.621 27.921
# save(working_list, file = "working_list.RData")

    
#!
################# HIV EXAMPLE MIXCR ###############
###################################################

colclasses <- rep("NULL", 33)
colclasses[c(1,5,6,7,8,9,10,11,32)] <- rep("character", length(c(1,5,6,7,8,9,10,11,32)))
cdr3_clones <- read.csv2("examples/clones_dataset2.txt", row.names=NULL, sep = "\t", colClasses = colclasses)




#!

###################  PART II #######################
####################################################
### Immune repertoire analysis and visualization ###
####################################################



## ----diversity (pg 46)
# Load in data calculated previously
load("cdr3_table_rank_stage.RData")
### Plot frequency and cumulative frequency
require(scales)
# Get titles for plot
rank_plots <- as.list(0)
titles <- c(levels(cdr3_hc$Donor)[match_stage])
for (i in 1:length(cdr3_table_rank_stage)){
    cdr3_freq <- data.frame(as.numeric(cdr3_table_rank_stage[[i]]))
    cdr3_freq$index <- 1:nrow(cdr3_freq)
    colnames(cdr3_freq)[1] <- "V1"
    cdr3_freq <- rbind(cdr3_freq, data.frame(V1=cumsum(cdr3_freq$V1),
                                             index=1:nrow(cdr3_freq)))
    cdr3_freq$facet <- factor(rep(factor(c("Frequency","Cumulative frequency")),
                                  each=length(cdr3s_table_rank_stage[[i]])),
                              levels = c("Frequency", "Cumulative frequency"))
    rank_plots[[i]] <- ggplot(cdr3_freq, aes(x=index, y= V1))
    rank_plots[[i]] <- rank_plots[[i]] + geom_point(size=1, alpha = 0.8) +
        facet_grid(facet ~., scales = "free") +
        labs(x= "CDR3 index", y="CDR3 frequency (%)", title = titles[i]) +
        theme_bw() +
        theme(plot.title = element_text(face="bold", size=rel(1), hjust=0),
              plot.margin = unit(c(2, 2, 2, 2),"points"),
              strip.text = element_text(size = rel(0.75)),
              legend.position = "none") +
        theme(strip.background = element_rect(fill = scales:::alpha("blue", 0.3)))
}














## ----SHMs (pg 48)

load("code/clon90")
mean_vmutations <- tapply(clon90$mean_mutations, clon90$L1, function(x)
    mean(x, na.rm=TRUE))
median_vmutations <- tapply(clon90$mean_mutations, clon90$L1, function(x)
    median(x, na.rm=TRUE))

# pdf("figure/density_shm_clon90.pdf") # Uncomment to produce plot
ggplot(clon90) +
    geom_freqpoly(binwidth=1, aes(x=mean_mutations, y=..density.., colour=L1), size=1.2) +
    scale_color_manual("Donor", values=c("red",
                                         "#56B4E9", "darkgreen", "yellow","#999999", "#E69F00"),
                       labels=c("HIV-1 IAVI donor 17 HC", "HIV-1 IAVI donor 17 LC",
                                "Uninfected donor 1 HC", "Uninfected donor 1 LC",
                                "Uninfected donor 2 HC", "Uninfected donor 2 LC")) +
    theme(axis.text.x = element_text(angle=360,size=5)) +
    theme_bw(base_size=12, base_family = "Helvetica") +
    labs (title = "", y="Density", x="Mutations clones (mean)") +
    theme(strip.background = element_rect(fill="white",
                                          linetype="solid", color="white"),
          strip.text=element_text(face="bold", size=rel(1.2), hjust=0.5, vjust=1)) +
    theme(strip.background=element_rect(fill = "white")) +
    theme(strip.background=element_rect(fill = "white")) +
    theme(plot.title=element_text(family="Helvetica",
                                  face="bold", size=rel(1.5), hjust=-0.035, vjust=3.5)) +
    theme(axis.text=element_text(family="Helvetica", size=rel(1.2)),
          axis.title=element_text(family="Helvetica", size=rel(1.2)),
          legend.text=element_text(family="Helvetica", size=rel(1.1)),
          legend.title=element_text(family="Helvetica", face="bold", size=rel(1.2))) +
    theme(axis.title.x = element_text(vjust=-0.5, size=rel(1.2)),
          axis.title.y = element_text(vjust=1, size=rel(1.2))) +
    theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
          axis.line=element_line(size = 0.7, color = "black"),
          text = element_text(size = 14)) +
    theme(plot.margin = unit(c(0, 0.2, 0.5, 0.2),"cm")) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +
    theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
          axis.line = element_line(size = 0.7, color = "black"),
          text = element_text(size = 12))+ theme(legend.position=c(0.805,0.85)) +
    geom_vline(xintercept=mean_vmutations, linetype="dashed", size=1,
               color=c("red", "#56B4E9", "darkgreen", "yellow","#999999", "#E69F00"))
# dev.off()










## ----cdr3_length_distribution (pg 50)

load("/code/clon90")

mean_cdr3_lengths <- tapply(clon90$cdr3_length, clon90$L1, function(x)
  mean(x, na.rm=TRUE))

median_cdr3_lengths <- tapply(clon90$cdr3_length, clon90$L1, function(x)
  median(x, na.rm=TRUE))

pdf("figure/geom_density_cdr3_lengths_clon90.pdf")
ggplot(clon90, aes(x=cdr3_length, color=L1)) + geom_density(size=1.2) +
  scale_color_manual("Donor", values=c("red", "#56B4E9", "darkgreen",
                     "yellow","#999999", "#E69F00"),
                     labels=c("HIV-1 IAVI donor 17 HC", "HIV-1 IAVI donor 17 LC",
                              "Uninfected donor 1 HC", "Uninfected donor 1 LC",
                              "Uninfected donor 2 HC", "Uninfected donor 2 LC")) +
  theme(axis.text.x = element_text(angle=360, size=5)) +
  theme_bw(base_size=12, base_family = "Helvetica") +
  labs (title = "", y="Density", x="CDR3 lengths") +
  theme(strip.background = element_rect(fill="white", linetype="solid", color="white"),
        strip.text=element_text(face="bold", size=rel(1.2), hjust=0.5, vjust=1)) +
  theme(strip.background=element_rect(fill = "white")) +
  theme(strip.background=element_rect(fill = "white")) +
  theme(plot.title=element_text(family="Helvetica", face="bold", size=rel(1.5),
                                hjust=-0.035, vjust=3.5)) +
  theme(axis.text=element_text(family="Helvetica", size=rel(1.2)),
        axis.title=element_text(family="Helvetica", size=rel(1.2)),
        legend.text=element_text(family="Helvetica", size=rel(1.1)),
        legend.title=element_text(family="Helvetica", face="bold", size=rel(1.2))) +
  theme(axis.title.x = element_text(vjust=-0.5, size=rel(1.2)),
        axis.title.y = element_text(vjust=1, size=rel(1.2))) +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line=element_line(size = 0.7, color = "black"),
        text = element_text(size = 14)) +
  theme(plot.margin = unit(c(0, 0.2, 0.5, 0.2),"cm")) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.7, color = "black"),
        text = element_text(size = 12))+ theme(legend.position=c(0.805,0.85)) +
  geom_vline(xintercept=mean_cdr3_lengths, linetype="dashed", size=1,
             color=c("red", "#56B4E9", "darkgreen", "yellow","#999999", "#E69F00"))
dev.off()




##----vgene_frequency_plot

load("working_hiv1_iavi17.RData")

# Germline genes - V gene frequency analysis

vgene_freq <- prop.table(table(working_hiv1_iavi17$vgene_subgroup))

vgene_plot <- ggplot(data.frame(vgene_freq), aes(x=Var1, y=Freq*100)) +
  theme_bw() + geom_bar(stat="identity", width=0.5, fill="dodgerblue4") +
  labs(x="", y="Frequency (%)") +
  theme(axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1)) +
  theme(axis.text.x=element_text(angle=90, size=20, hjust=1, vjust=0.5)) +
  theme(axis.title.y = element_text(colour="black", size=25, vjust=1.5, hjust=0.4),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=20)) +
  theme(plot.margin = unit(c(12, 12, 12, 12),"points"))

pdf("figure/vgene_frequency.pdf", width=12)
vgene_plot
dev.off()




## ----extraction and subsetting techniques
# The following code shows how to select for top frequency CDR3 
# in a dataset, and how to subset for different parameters

# Load data
load("code/working_hiv1_iavi17.RData")

# As we have ordered our dataframes by CDR3 ranks,
# we can subset the 100 most frequent clones by
top_frequency_cdr3 <- as.character(head(working_hiv1_iavi17[["cdr3"]], 100))

# In order to select all the columns of the top 100 most frequent
# CDR3s in the original dataset:
all_top_cdr3 <- working_hiv1_iavi17[top_frequency_cdr3,]

# Calculate the median mutations in this subset of CDR3s
median(all_top_cdr3$vregion_mutations)
# [1] 21

# FR2 for CDR3 that have V gene IGHV1-2 in the top 100 most frequent CDR3s
data.frame(fr2=all_top_cdr3[all_top_cdr3$vgene_subgroup=="IGHV1-2",]$fr2)

# A faster way is by using the subset function.
# We subset the fr3 in the original dataframe and V genes
# by the clones that have CDR3 longer than 27 a.a.
subset(working_hiv1_iavi17, cdr3_length > 27, select = c(fr3,vgene_subgroup))

# Select for the "CDR3" with all its characteristics (V gene, CDR3 frequency, etc)
working_hiv1_iavi17[working_hiv1_iavi17$cdr3 == "VRDGAYGCSGASCYFGALGNFVYYYYMDV",]

# Find sequence in dataset
which(working_hiv1_iavi17$cdr3 == "VRDGAYGCSGASCYFGALGNFVYYYYMDV")
# [1] 210

# Find similar sequences
cdr3_distance <- stringDist(as.character(working_hiv1_iavi17$cdr3),
                            method = "levenshtein")

# Assign 1 to all similar sequences (max. 1 a.a. different)
# and zero to all sequences that differ by more than 1 a.a.
cdr3_matrix <- ifelse(as.matrix(cdr3_distance) == 1, 1, 0)
colnames(cdr3_matrix) <- as.character(working_hiv1_iavi17$cdr3)
rownames(cdr3_matrix) <- as.character(working_hiv1_iavi17$cdr3)

which(cdr3_matrix[rownames(cdr3_matrix)=="ARRGIAGPDYYSYHGLDV"]==1)
# [1]   52  106  122 1186 1436 1437 1438 1552


## ----save_output

# Example loading presaved data

load("working_hiv1_iavi17.RData")

# Example of saving results as RData

subset_fr3 <- subset(working_hiv1_iavi17, cdr3_length > 27, select = c(fr3,vgene_subgroup))

save(subset_fr3, file="subset_fr3.RData")

# Example of writing results as CSV

top_cdr3 <- as.character(head(working_hiv1_iavi17[["cdr3"]], 100))

write.csv(top_cdr3, "top_cdr3.csv")

# Example of loading and reading-in previous results

load("subset_fr3.RData")

# Print the first row of the data

print(subset_fr3[1,])

nrow(subset_fr3)









## ----visualization (pg 58)

# Load clones based on 90 % CDR3 a.a. identity
load("clon90")
install.packages("e1071")
library(e1071)

skewness(clon90[clon90$L1=="HIV-1 IAVI donor 17 HC",]$mean_mutations)
# [1] 1.301486

#placeholder plot - prints nothing at all
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
     theme(
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank()
     )

#scatterplot of x and y variables
scatter <- ggplot(clon90,aes(cdr3_length, mean_mutations)) +
geom_point(aes(color=L1)) +
scale_color_manual(values = c("red", "#56B4E9", "darkgreen",
                              "yellow","#999999", "#E69F00")) +
theme(legend.position=c(1,1),legend.justification=c(1,1))
  # scatter <- hc_abundance_mutations_plot

#marginal density of x - plot on top
plot_top <- ggplot(clon90,aes(cdr3_length, fill=L1)) +
  geom_density(alpha=.5) +
  scale_fill_manual(values = c("red", "#56B4E9", "darkgreen",
                               "yellow","#999999", "#E69F00")) +
  theme(legend.position = "none")

#marginal density of y - plot on the right
plot_right <- ggplot(clon90,aes(mean_mutations, fill=L1)) +
  geom_density(alpha=.5) +
  coord_flip() +
  scale_fill_manual(values = c("red", "#56B4E9", "darkgreen",
                               "yellow","#999999", "#E69F00")) +
  theme(legend.position = "none")

# Arrange the plots together, with appropriate
# height and width for each row and column
pdf("figure/cdr3_length_mean_mutations.pdf")
grid.arrange(plot_top, empty, scatter, plot_right,
             ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()



## ----networks

load("code/working_hiv1_iavi17.RData")

library(igraph)

cdr3 <- as.character(working_hiv1_iavi17$cdr3)
cdr3_dist <- stringDist(cdr3, method = "levenshtein")
cdr3_mat <- as.matrix(cdr3_dist)
cdr3_bol <- cdr3_mat
cdr3_bol[cdr3_bol<=1] <- 1
cdr3_bol[cdr3_bol>1] <- 0
colnames(cdr3_bol) <- cdr3
rownames(cdr3_bol) <- cdr3

# CDR3 corresponding to the top aligned VDJ with the bNAb database
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3)

# The entire row
working_hiv1_iavi17[working_hiv1_iavi17$vdj_region==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSRVV
ISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]

cdr3_graph <- igraph::simplify(graph.adjacency(cdr3_bol, weighted=T, mode = "undirected"))

# Name of the first node of the graph
V(cdr3_graph)$name[1]

nodesize <- ifelse(V(cdr3_graph)$name ==
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3),
3, working_hiv1_iavi17$vregion_mutations/11)

# Mark in red the top hit (best aligned) with bNAb database sequences
nodecolor <- ifelse(V(cdr3_graph)$name ==
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3),
"red", "grey")

pdf("figure/cdr3_network.pdf")
set.seed(11)
plot(cdr3_graph, vertex.frame.color=NA,
     layout=layout_with_fr(cdr3_graph, grid="nogrid", niter=200),
     vertex.label=NA, vertex.color=nodecolor, edge.width=1,
     vertex.size=nodesize, edge.color = "black")
dev.off()

