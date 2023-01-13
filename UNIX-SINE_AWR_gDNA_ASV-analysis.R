# ---- Description ----
# Analysis of ASV output

# ---- Set Working Directory ----
setwd("/Users/samanthabeal/Documents/MSc/Bioinformatics/")

# ---- Load Packages ----
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

##########################################################################
#do this the first time, but don't need to reinstall everything each time
#for ranacapa
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("multtest")
#BiocManager::install("phyloseq")

#packageVersion("phyloseq")

#if (!requireNamespace("devtools", quietly = TRUE))
  #install.packages('devtools')
#devtools::install_github("gauravsk/ranacapa")
##########################################################################

# ---- Load Rest of Packages ----
library(multtest)
library(phyloseq)
library(devtools)
library(ranacapa)

# ---- Check File Number ----
# Create list with each .tsv as dataframe with parameter settings specified

# Get path and filenames of all .fastq.gz files in output
#allFiles(filepath)
allFiles = Sys.glob("AWF_gDNA/input/*.fastq.gz")
# Keep only file name
allFiles2 = sub(".*/", "", allFiles)
# Keep only location, site number, site visit, and replicate
allFiles3 = sub('^([^-]+-[^-]+-[^-]+).*', "\\1", allFiles2)
# Replace all "-" with "_" to match eDNA Master
allFiles4 = gsub('\\-', '_', allFiles3)

# Create Dataframe
files = ldply(allFiles4, data.frame)
colnames(files) = "ID"

#check
files

# Count number of occurances
# Note: total occurrences for each ID must be divided by 2 because paired-end reads (R1 & R2)
fileNumber = as.data.frame(table(files$ID))
colnames(fileNumber) = c("sample_ID", "PCR_reps")

# Divide PCR_reps by 2 because paired-end reads
fileNumber$PCR_reps = fileNumber$PCR_reps / 2
dim(fileNumber)
#50,2

##########################################################################
#since only looking at Aquatron samples, not sure if i need to do this all 
# Cross-reference with eDNA-MASTER (i.e., number of sites)
# Import Data
master = read.csv("eDNA_metadata_MASTER_220802.csv",
                  header = TRUE)
# Data Check
head(master)

sites = master[c("sample_location_ID", "site_number", "sample_ID")]

# Remove duplicate sample_IDs
sites2 = sites[!duplicated(sites$sample_ID), ]

# If needed, cross-reference with plate layout (i.e., number of expected PCR replicates)

# Merge sites2 and fileNumber
# Use all = TRUE then filter by sample_location_ID to check each site
rep_check = merge(sites2, fileNumber, by = "sample_ID", all = TRUE)

# Check a specific lake, site or sample series
ASP = rep_check[which(rep_check$sample_location_ID == 'ASP'), ]
#& rep_check$sampling_event <= 3), ]

# save output
write.csv(ASP, "ASP-rep_check-220811.csv",
          row.names = FALSE)


##########################################################################


# ---- Import Data ----
# file with number of reads for each ASV per sample 
# .tsv
# the names in this table aren't shortened. figure out how to fix this
asv_table <- read.delim("AWF_gDNA/output-SINE-anchored/ASV/ASV_AWFgDNA.tsv",
                        header = TRUE)

# file with blast output (most likely species) 
# .txt
blast <- read.table("AWF_gDNA/output-SINE-anchored/ASV/ASV_AWFgDNA.txt")

# file with taxonomy
# use to cross-reference with blast output
taxon <- read.delim("Ref_Databases/SINE_taxonomy.txt", header = FALSE)

# --- Data Check ---
head(asv_table)
head(blast)
head(taxon)

# name columns
names(asv_table)[1] <- "ASV"
names(taxon)[1] <- "ID"
names(taxon)[2] <- "taxonomy"
names(blast)[1] <- "ASV"
names(blast)[2] <- "ID"
names(blast)[3] <- "perc_ident_matches"
names(blast)[4] <- "length"

# the smaller the evalue, the better the match
names(blast)[11] <- "evalue" 

# bit-score adjusted to the sequence database size
# the higher the bit-score, the better the sequence similarity
names(blast)[12] <- "bitscore"  

# split column ASV in blast
# create columns for size and ASV
blast2 <- separate(data = blast, col = ASV, into = c("ASV", "size"), sep = "\\;")

# remove 'size=' from column size
blast2$size <- gsub("size=([0-9]+).*", "\\1", blast2$size)


# ---- Index Species ID ----
# add taxonomy to blast
# cross-reference ID; taxon V1 with blast V2
# add taxon V2 to blast
blastID <- merge(taxon, blast2, by.x = "ID")

# --- Data Check ---
head(blastID)
View(blastID)



# ---- Index Taxonomy ----
# add taxonomy to asv_table
# cross-reference ASV; asv_table V1 with blastID V3
# first hit will be imported to asv_table
asv_out <- merge(blastID, asv_table, by.x = "ASV")

# --- Data Check ---
head(asv_out)


# --- Sort ---
# sort asv_out by bitscore
# the higher the bit-score, the better the sequence similarity
asv_final <- asv_out[order(asv_out$bitscore,
                           decreasing = TRUE), ]
# sort blastID by perc_ident_matches 
# highest value first (i.e. descending)
asv_final <- asv_final[order(asv_final$perc_ident_matches,
                             decreasing = TRUE), ]
# sort by ASV
asv_final <- asv_final[order(asv_final$ASV,
                             decreasing = FALSE), ]

# --- Data Check ---
head(asv_final)
View(asv_final)


# ---- Remove Duplicate ASVs ----
# keep only the first (i.e. best match; highest percent identity match and bitscore)
# .keep_all = TRUE will keep all columns (if FALSE keeps only 'ASV')
top_asv <- distinct(asv_final, ASV, .keep_all = TRUE)


# --- Data Check ---
head(top_asv)
View(top_asv)

dim(top_asv)
#387,64

# check minimum perc_ident_matches
# should be > 97 (cut-off used in pipeline)
min(top_asv[ ,5]) 
#97.015

# compare number of rows to asv_table
dim(asv_table)
#2328, 51

# if different
# find which ASVs are missing compared to asv_table
anti_join(asv_table, top_asv)

##########################################################################
# ---- WORKING: Identify Species Missing from Database ----
# find which ASVs are missing from database
# cross reference with .fasta
# open .fasta and search blast manually 
# add species to database if needed
##########################################################################

# ---- Format Output ----
# split column taxonomy 
# Kingdom;Phylum;Class;Order;Family;Genus;Species
top_asv2 <- separate(data = top_asv,
                     col = taxonomy,
                     into = c("Kingdom",
                              "Phylum",
                              "Class",
                              "Order",
                              "Family",
                              "Genus",
                              "Species"), 
                     sep = "\\;")

# save output
write.csv(top_asv2, "AWF_gDNA/output-SINE-anchored/ASV/AWFgDNA.csv",
          row.names = FALSE)





# ---- ranacapa ----
# see github: https://github.com/gauravsk/ranacapa
# input file must contain one column for each sample, full taxonomic path, and sequence number (optional)

# modify top_asv
# remove V2, V4 - V14
rana_in <- subset(top_asv, select = -c(ID, size, perc_ident_matches, length, V5, V6, V7, V8, V9, V10, evalue, bitscore))
# change name of ASV (V1) to 12S_seq_number
#names(rana_in)[1] <- "12S_seq_number"
# move taxonomy (V3) to end
rana_in <- rana_in %>% relocate(taxonomy, .after = last_col())
names(rana_in)[50] <- "sum.taxonomy"


# ---- Import Data ----
# site metadata


















