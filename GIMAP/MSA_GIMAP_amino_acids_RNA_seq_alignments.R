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
GIMAP_amino_acid_sequence <- readAAStringSet('/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline/Cvir_apoptosis_manual_annotation/GIMAP/GIMAP_XP_all_amino_acid_sequence_shortened_header.fasta')
GIMAP_rna_sequence <- readAAStringSet('/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline/Cvir_apoptosis_manual_annotation/GIMAP/GIMAP_XP_all_rna_sequence_fixed_shortened_headers.fasta')
# Create alignments with clustalw in MSA
GIMAP_amino_acid_alignment <- msa(GIMAP_amino_acid_sequence)
GIMAP_rna_alignment <- msa(GIMAP_rna_sequence)

## Upload alignments into MSA format using 

## Print the alignments in a pretty way 
msaPrettyPrint(GIMAP_amino_acid_alignment, output="pdf", showLogo="none", askForOverwrite = FALSE)
msaPrettyPrint(GIMAP_rna_alignment, output="pdf", showLogo="none", askForOverwrite = FALSE)

## Add Gene names onto the appropriate rows of the alignment using 
