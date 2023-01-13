#!/bin/bash

################
# Primer removal 
################

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX

# list files in folder
# extract part of name before second underscore and find unique
ls input-AWF_gDNA | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable and check contents
samples=$(ls input-AWF_gDNA | cut -d_ -f1,2 | sort | uniq)
echo $samples

# count number of sequences across all files in folder
cd input-AWF_gDNA
gzcat *.fastq.gz | grep -c "^@M00" 
#2217624

################################################################################################
#SINE-----
cd ..
mkdir output-SINE-anchored

# do not anchor (^) primers - anchoring tells cutadapt that the primer is at the beginning of ther read
# set minimum overlap parameter (-O) to 18 (i.e. must find, in order at least 18 bases of primer seq)
# add 'X' to indicate nothing (i.e. sequence should start at first bp of primer)

# SINE (Smal cor-II) - only R1 sequenced
for s in $samples;
do
	cutadapt -g ^TAGCTCAGCTGGTAGAGCAC \
	-o output-SINE-anchored/${s}_L001_R1.fastq.gz --discard-untrimmed \
	input/${s}_L001_R1_001.fastq.gz;
done

# count number of sequences across all files in folder
cd output-SINE-anchored
gzcat *.fastq.gz | grep -c "^@M00" 
#unanchored = 1484656
#anchored 1385527
#what about -g XT... ? 

#####################
# Sequence Quality
#####################

# IF fusing rather than merging do this after quality filtering
# make sure in folder 'trimmed'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

# Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/AWF_gDNA/output-SINE-anchored/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/AWF_gDNA/output-SINE-anchored/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/AWF_gDNA/output-SINE-anchored/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/AWF_gDNA/output-SINE-anchored/qc

cd qc
multiqc .

#####################
# Trim Sequences
#####################

cd .. #(back to output)
mkdir trimmed

# left oligo removed in previous step
# right oligo not sequenced...run with 151 cycles 
# expected fragment length 90bp
# trim to be conservative and have consistent length
# left oligo 20bp
# qc indicated drop-off at ~130bp
# 130 - 20bp = 110bp
#trimmed to 100

ls ./*_R1.fastq.gz | sed 's/_R1.fastq.gz$//' | parallel 'java -jar /Users/samanthabeal/Trimmomatic-0.39/trimmomatic-0.39.jar SE {}_R1.fastq.gz /Users/samanthabeal/Documents/MSc/Bioinformatics/AWF_gDNA/output-SINE-anchored/trimmed/{}_R1.fastq.gz CROP:100 MINLEN:100'


#####################
# Concatenate Samples
#####################

mkdir ASV #(in output)

cd trimmed
gunzip *.fastq.gz

# add sample name to sequences
# combine all sequences in one file
for f in *;
do                         
        sed -e "s/\(^@M00.*\) .*$/\1;sample=${f%.*};/" $f \
        >> ../ASV/concatenated_AWFgDNA.fastq;
done


###################
# Quality Filtering
###################

# move up in folder structure
cd ..
cd ASV

# quality filtering
# ONLY do this if merged rather than fused
# if fused will discard all sequences
vsearch --fastx_filter concatenated_AWFgDNA.fastq --fastq_maxee 1 --fastaout concatenated_AWFgDNA.fasta
#1,295,321 sequences kept (of which 0 truncated), 90206 sequences discarded.

################
# Dereplication
################

# header of each unique sequence records number of copies in dataset
vsearch --derep_fulllength concatenated_AWFgDNA.fasta --output derep_AWFgDNA.fasta --sizeout --relabel uniq
#280335 unique sequences, avg cluster 4.6, median 1, max 134070
#count number of unique sequences
grep -c "^>" derep_AWFgDNA.fasta
#280335

############
# Denoising
############
vsearch --cluster_unoise derep_AWFgDNA.fasta --minsize 8 --unoise_alpha 2 --centroids denoised_AWFgDNA.fasta
#1621800 nt in 16218 seqs, min 100, max 100, avg 100
#minsize 8: 264117 sequences discarded
#Masking 100% 
#Sorting by abundance 100%
#Counting k-mers 100% 
#Clustering 100%  
#Sorting clusters 100%
#Writing clusters 100% 
#Clusters: 2356 Size min 1, max 11823, avg 6.9
#Singletons: 2053, 12.7% of seqs, 87.1% of clusters


#####################
# Chimera Filtering
#####################
vsearch --uchime3_denovo denoised_AWFgDNA.fasta --nonchimeras ASV_AWFgDNA.fasta
#Found 28 (1.2%) chimeras, 2328 (98.8%) non-chimeras
#and 0 (0.0%) borderline sequences in 2356 unique sequences
#Taking abundance information into account, this corresponds to
#832 (0.3%) chimeras, 242975 (99.7%) non-chimeras,
#and 0 (0.0%) borderline sequences in 243807 total sequences

#######################
# Mapping Reads to ASVs
#######################

vsearch --search_exact concatenated_AWFgDNA.fasta -db ASV_AWFgDNA.fasta -otutabout ASV_AWFgDNA.tsv
#232800 nt in 2328 seqs, min 100, max 100, avg 100
#Masking 100% 
#Hashing database sequences 100%  
#Searching 100%  
#Matching query sequences: 242975 of 1295321 (18.76%)
#Writing OTU table (classic) 100%  

##############################
# Identifying ASVs with Blast
##############################

# -outmt 6 is standard blast tabulated output
# set strict evalue and perc_identity since expect many closely related species
# -evalue 
# -perc_identity

cd ../..

blastn -query output-SINE-anchored/ASV/ASV_AWFgDNA.fasta -subject SINE.fasta -outfmt 6 -out output-SINE-anchored/ASV/ASV_AWFgDNA.txt \
-num_threads 1 -evalue 0.001 -perc_identity 97










################################################################################################
################################################################################################

