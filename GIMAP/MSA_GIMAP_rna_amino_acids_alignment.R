#MSA_GIMAP_amino_acids_RNAseq_alignments.R

#load packages

source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")
library(msa)
library(Biostrings)
system.file("tex", "texshade.sty", package="msa") #see where required LATEX file is located
install.packages("ape")
library(ape)
install.packages("seqinr")
library(seqinr)

## Upload sequences
# all found GIMAP sequences were used for this 
GIMAP_amino_acid_sequence <- readAAStringSet('/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline/Cvir_apoptosis_manual_annotation/GIMAP/GIMAP_XP_all_amino_acid_sequence_shortest_header.fasta')
GIMAP_rna_sequence <- readAAStringSet('/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline/Cvir_apoptosis_manual_annotation/GIMAP/GIMAP_XP_all_rna_sequence_fixed_shortest_headers.fasta')
# Create alignments with clustalw in MSA
GIMAP_amino_acid_alignment <- msa(GIMAP_amino_acid_sequence)
GIMAP_rna_alignment <- msa(GIMAP_rna_sequence)


## Print the alignments in a pretty way and trim regions at ends with low consensus 
msaPrettyPrint(GIMAP_amino_acid_alignment, y=c(200,600), output="pdf", showLogo="none", consensusColor="ColdHot")

#this aligment is too big for LaTex to print
#msaPrettyPrint(GIMAP_rna_alignment, output="tex", showConsensus= "none", askForOverwrite=TRUE, verbose=FALSE)
#library(tools)
#texi2pdf('GIMAP_rna_alignment.tex', clean=TRUE)

## Add Gene names onto the appropriate rows of the alignment using Paint

# Create phylogenetic trees using MSA 


