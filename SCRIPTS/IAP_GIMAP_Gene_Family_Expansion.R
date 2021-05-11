#### R script to identify and test gene family expansion
# Erin Roberts, 2020
# PhD Candidate University of Rhode Island 

#### Load packages ####
library(ape)
library(Biostrings)
library(ggrepel)
library(phylotools)
library(treeio)
library(tidytree)
library(ggimage)
library(tidyverse)
library(tidytext)
library(rtracklayer)
library(data.table)
library(chopper)
library(alakazam)
library(phylotools)
library(viridis)
library(ggpubr)
library(forcats) 
library(cowplot)
#install.packages("remotes") 
#remotes::install_github("YuLab-SMU/ggtree")
library(ggtree) # install the dev version to get the get.tree function
library(aplot)
library(RColorBrewer)
library(gtable)
library(gridExtra)
library(egg)
library(grid)
library(ggmsa)
library(ggtext)
library(ggplotify)
library(UpSetR)
library(gt)
library(GenomicFeatures)
library(GenomicRanges)
library(RCircos)
setwd("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics")

#### IMPORT GENOMES AND ANNOTATIONS #####
load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
# Load the gff file for all mollusc genomes
All_molluscs_CDS_gff <- readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/All_molluscs_CDS.gff")
All_molluscs_CDS_gff <- as.data.frame(All_molluscs_CDS_gff)
All_mollusc_gene_gff <- readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/All_mollusc_gene.gff")
All_mollusc_gene_gff <- as.data.frame(All_mollusc_gene_gff)
All_mollusc_exon_gff <- readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/All_mollusc_exon.gff")
All_mollusc_exon_gff <- as.data.frame(All_mollusc_exon_gff)

# remove from memory when not using
rm(All_molluscs_CDS_gff,All_mollusc_gene_gff,All_mollusc_exon_gff)

### MANUAL SEARCH OF IAPS AND GIMAPS IN REFERENCE ANNOTATION ####
Cvir_gff_IAP_family <- C_vir_rtracklayer[grepl("inhibitor of apoptosis", C_vir_rtracklayer$product, ignore.case=TRUE) | grepl("XIAP", C_vir_rtracklayer$product, ignore.case=TRUE) |
                                            grepl("baculoviral", C_vir_rtracklayer$product, ignore.case=TRUE),]
Cvir_gff_IAP_family_XP <- Cvir_gff_IAP_family[!duplicated(Cvir_gff_IAP_family$protein_id),]
length(unique(Cvir_gff_IAP_family_XP$protein_id)) # 136
length(unique(Cvir_gff_IAP_family_XP$gene)) # 67

Cvir_gff_IAP_family_XP <- Cvir_gff_IAP_family_XP[!is.na(Cvir_gff_IAP_family_XP$protein_id),]
Cgig_gff_IAP_family <- C_gig_rtracklayer[grepl("inhibitor of apoptosis", C_gig_rtracklayer$product, ignore.case=TRUE) | grepl("XIAP", C_gig_rtracklayer$product, ignore.case=TRUE)|
                                           grepl("baculoviral", C_gig_rtracklayer$product, ignore.case=TRUE),]
Cgig_gff_IAP_family_XP <- Cgig_gff_IAP_family[!duplicated(Cgig_gff_IAP_family$protein_id),]
Cgig_gff_IAP_family_XP <- Cgig_gff_IAP_family_XP[!is.na(Cgig_gff_IAP_family_XP$protein_id),]
length(unique(Cgig_gff_IAP_family_XP$protein_id)) #64
length(unique(Cgig_gff_IAP_family_XP$gene)) # 35

Cvir_gff_GIMAP_family <- C_vir_rtracklayer[grepl("IMAP", C_vir_rtracklayer$product, ignore.case=TRUE) | grepl("immune-associated nucleotide-binding protein", C_vir_rtracklayer$product, ignore.case=TRUE),]
Cvir_gff_GIMAP_family_XP <- Cvir_gff_GIMAP_family[!duplicated(Cvir_gff_GIMAP_family$protein_id),]
Cvir_gff_GIMAP_family_XP <- Cvir_gff_GIMAP_family_XP[!is.na(Cvir_gff_GIMAP_family_XP$protein_id),]
length(unique(Cvir_gff_GIMAP_family_XP$protein_id)) #121
length(unique(Cvir_gff_GIMAP_family_XP$gene)) # 65

Cgig_gff_GIMAP_family <- C_gig_rtracklayer[grepl("IMAP", C_gig_rtracklayer$product, ignore.case=TRUE) | grepl("immune-associated nucleotide-binding protein", C_gig_rtracklayer$product, ignore.case=TRUE),]
Cgig_gff_GIMAP_family_XP <- Cgig_gff_GIMAP_family[!duplicated(Cgig_gff_GIMAP_family$protein_id),]
Cgig_gff_GIMAP_family_XP <- Cgig_gff_GIMAP_family_XP[!is.na(Cgig_gff_GIMAP_family_XP$protein_id),]
length(unique(Cgig_gff_GIMAP_family_XP$protein_id)) # 39
length(unique(Cgig_gff_GIMAP_family_XP$gene)) # 35

CV_CG_IAP <- c(Cvir_gff_IAP_family_XP$protein_id, Cgig_gff_IAP_family_XP$protein_id)
CV_CG_GIMAP <- c(Cvir_gff_GIMAP_family_XP$protein_id, Cgig_gff_GIMAP_family_XP$protein_id)

#### IAP and GIMAP HMMER and INTERPROSCAN ANALYSIS ####
## LOAD HMM and Interproscan XPS ##
BIR_XP_gff <- readGFF(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_hmmsearch_eval3_XP_list.fa.gff3")
BIR_XP_gff <- as.data.frame(BIR_XP_gff)
length(unique(BIR_XP_gff$seqid)) # 800

AIG1_XP_ALL_gff <- readGFF(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_hmmsearch_eval3_XP_list.fa.gff3")
AIG1_XP_ALL_gff <- as.data.frame(AIG1_XP_ALL_gff)
length(unique(AIG1_XP_ALL_gff$seqid)) # 1410 unique proteins 

## Filter out IAP proteins that don't have a BIR repeat domain
BIR_XP_gff <- BIR_XP_gff %>% group_by(seqid) %>% filter(any(grepl("BIR", signature_desc)))
length(unique(BIR_XP_gff$seqid)) # 674 if you only include BIR from the CDD search, but 791 if you include all BIR repeat profiles from any source
colnames(BIR_XP_gff)[1] <- "protein_id"

## Filter out AIG Interproscan results for proteins that have a line with coil and a line with AIG1 from CDD and from Coils
AIG1_XP_ALL_gff_GIMAP <- AIG1_XP_ALL_gff %>% group_by(seqid) %>% filter(any(source == "CDD" & signature_desc == "AIG1") | any(source == "Pfam" & signature_desc == "AIG1 family"))
AIG1_XP_ALL_gff_GIMAP <- AIG1_XP_ALL_gff_GIMAP %>% group_by(seqid) %>% filter(any(source == "Coils" & Name == "Coil"))
length(unique(AIG1_XP_ALL_gff_GIMAP$seqid)) # 403

# Check annotations
colnames(AIG1_XP_ALL_gff_GIMAP)[1] <- "protein_id"
class(All_molluscs_CDS_gff$protein_id) # character
AIG1_XP_ALL_gff_GIMAP$protein_id <- as.character(AIG1_XP_ALL_gff_GIMAP$protein_id)
AIG1_XP_ALL_gff_GIMAP_annot <- left_join(AIG1_XP_ALL_gff_GIMAP[,c("protein_id")], All_molluscs_CDS_gff[,c("protein_id","product","gene")])
AIG1_XP_ALL_gff_GIMAP_annot <- AIG1_XP_ALL_gff_GIMAP_annot[!duplicated(AIG1_XP_ALL_gff_GIMAP_annot$protein_id),]

class(BIR_XP_gff$protein_id) # factor
BIR_XP_gff$protein_id <- as.character(BIR_XP_gff$protein_id)
BIR_XP_gff_annot <- left_join(BIR_XP_gff[,c("protein_id")], All_molluscs_CDS_gff[,c("protein_id","product","gene")])
BIR_XP_gff_annot <- BIR_XP_gff_annot[!duplicated(BIR_XP_gff_annot$protein_id),]

## Get the species name for each sequence
# Get species names using seqid (assembly ID)
BIR_XP_gff_seqid <- left_join(BIR_XP_gff, unique(All_molluscs_CDS_gff[,c("protein_id","seqid")]))
AIG1_XP_ALL_gff_GIMAP_seqid <- left_join(AIG1_XP_ALL_gff_GIMAP, unique(All_molluscs_CDS_gff[,c("protein_id","seqid")]))
length(unique(BIR_XP_gff_seqid$seqid)) # 194
length(unique(AIG1_XP_ALL_gff_GIMAP_seqid$seqid)) # 119

write.table(unique(BIR_XP_gff_seqid$seqid), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_seqid.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)
write.table(unique(AIG1_XP_ALL_gff_GIMAP_seqid$seqid), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_XP_ALL_gff_GIMAP_seqid.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)
# In batch entrez load file and search Nucleotide database, export summary to file
# edit
#$ sed '/^$/d' BIR_species_list.txt > BIR_species_list_blanks_rm.txt
#$ sed '/^$/d' GIMAP_species_list.txt > GIMAP_species_list_blanks_rm.txt
## paste every three lines together
#$ paste -d ", "  - - - < GIMAP_species_list_blanks_rm.txt > GIMAP_species_list_blanks_rm_paste.txt
#$ paste -d ", "  - - - < BIR_species_list_blanks_rm.txt > BIR_species_list_blanks_rm_paste.txt
# edited in excel

# Load into R
GIMAP_seqid <- read.csv(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_species_list_Acc.csv")
IAP_seqid <- read.csv(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_species_list_Acc.csv")

# Combine lists
Species_seqid <- rbind(GIMAP_seqid,IAP_seqid)

# combine with genome to add on species ID
All_molluscs_CDS_gff <- left_join(All_molluscs_CDS_gff, unique(Species_seqid))

# Join species ID with proteins and gene info
BIR_XP_gff_species <- left_join(BIR_XP_gff, unique(All_molluscs_CDS_gff[,c("protein_id","Species","product","gene","locus_tag")]))
length(unique(BIR_XP_gff_species$protein_id)) # 791
BIR_XP_gff_species %>% filter(is.na(Species)) #0
AIG1_XP_ALL_gff_GIMAP_species <- left_join(AIG1_XP_ALL_gff_GIMAP, unique(All_molluscs_CDS_gff[,c("protein_id","Species","product","gene","locus_tag")]))
length(unique(AIG1_XP_ALL_gff_GIMAP_species$protein_id)) # 403
AIG1_XP_ALL_gff_GIMAP_species %>% filter(is.na(Species)) #0

## Were any proteins added or missing with HMMER/Interproscan compare to genome? 
BIR_XP_gff_CG <- BIR_XP_gff_species %>% filter(Species=="Crassostrea_gigas")
CG_IAP_Hmmer_added <- BIR_XP_gff_CG[!(BIR_XP_gff_CG$protein_id %in% Cgig_gff_IAP_family_XP$protein_id),] # 15 proteins added by HMMER/interproscan that were not annotated in genome
CG_IAP_Hmmer_added[!duplicated(CG_IAP_Hmmer_added$protein_id),] # 15
CG_IAP_Hmmer_missing <- Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$protein_id %in% BIR_XP_gff_CG$protein_id),] # 5 proteins from original list not found
nrow(CG_IAP_Hmmer_missing[!duplicated(CG_IAP_Hmmer_missing$protein_id),]) # 5

BIR_XP_gff_CV <- BIR_XP_gff_species %>% filter(Species=="Crassostrea_virginica")
CV_IAP_Hmmer_added <- BIR_XP_gff_CV[!(BIR_XP_gff_CV$protein_id %in% Cvir_gff_IAP_family_XP$protein_id),] # 39 proteins added by HMMER/interproscan that were not annotated in genome
CV_IAP_Hmmer_added[!duplicated(CV_IAP_Hmmer_added$protein_id),] # 39
CV_IAP_Hmmer_missing <- Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$protein_id %in% BIR_XP_gff_CV$protein_id),] # 11 proteins from original list not found
nrow(CV_IAP_Hmmer_missing[!duplicated(CV_IAP_Hmmer_missing$protein_id),]) # 11

AIG1_XP_ALL_gff_GIMAP_CG <- AIG1_XP_ALL_gff_GIMAP_species %>% filter(Species=="Crassostrea_gigas")
CG_GIMAP_Hmmer_added <- AIG1_XP_ALL_gff_GIMAP_CG[!(AIG1_XP_ALL_gff_GIMAP_CG$protein_id %in% Cgig_gff_GIMAP_family_XP$protein_id),] # 4 proteins added by HMMER/interproscan that were not annotated in genome
CG_GIMAP_Hmmer_added[!duplicated(CG_GIMAP_Hmmer_added$protein_id),] # 4
CG_GIMAP_Hmmer_missing <- Cgig_gff_GIMAP_family_XP[!(Cgig_gff_GIMAP_family_XP$protein_id %in% AIG1_XP_ALL_gff_GIMAP_CG$protein_id),] #  proteins from original list not found
nrow(CG_GIMAP_Hmmer_missing[!duplicated(CG_GIMAP_Hmmer_missing$protein_id),]) # 12 proteins not found from original list

AIG1_XP_ALL_gff_GIMAP_CV <- AIG1_XP_ALL_gff_GIMAP_species %>% filter(Species=="Crassostrea_virginica")
CV_GIMAP_Hmmer_added <- AIG1_XP_ALL_gff_GIMAP_CV[!(AIG1_XP_ALL_gff_GIMAP_CV$protein_id %in% Cvir_gff_GIMAP_family_XP$protein_id),] # 23 proteins added by HMMER/interproscan that were not annotated in genome
CV_GIMAP_Hmmer_added[!duplicated(CV_GIMAP_Hmmer_added$protein_id),] # 23
CV_GIMAP_Hmmer_missing <- Cvir_gff_GIMAP_family_XP[!(Cvir_gff_GIMAP_family_XP$protein_id %in% AIG1_XP_ALL_gff_GIMAP_CV$protein_id),] # 35 proteins from original list not found
nrow(CV_GIMAP_Hmmer_missing[!duplicated(CV_GIMAP_Hmmer_missing$gene),]) # 35 proteins not found from original list, 18 genes
  # checking these genes in NCBI to look for the presence of the AIG1 domain, they seem to just have the Ploop NTPase superfamily domain and not AIG1

# Export protein lists for use with making Trees with MAFFT and RAxmL
write.table(unique(BIR_XP_gff_species$protein_id), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(unique(AIG1_XP_ALL_gff_GIMAP_species$protein_id), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

### Export XP and XM lists for CV and CG only to use to look up transcripts in the DESeq dataset ###
BIR_XP_gff_CG_uniq_XP <- unique(BIR_XP_gff_CG$protein_id)
BIR_XP_gff_CV_uniq_XP <- unique(BIR_XP_gff_CV$protein_id)
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP <- unique(AIG1_XP_ALL_gff_GIMAP_CG$protein_id)
AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP <- unique(AIG1_XP_ALL_gff_GIMAP_CV$protein_id)

#save(BIR_XP_gff_CG_uniq_XP, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_protein_list_CG.Rdata")
#save(BIR_XP_gff_CV_uniq_XP, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_protein_list_CV.Rdata")
#save(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_XP_ALL_gff_GIMAP_protein_list_CG.Rdata")
#save(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_XP_ALL_gff_GIMAP_protein_list_CV.Rdata")

# Join with XM information 
BIR_XP_gff_CG_uniq_XP <- as.data.frame(BIR_XP_gff_CG_uniq_XP)
BIR_XP_gff_CV_uniq_XP <- as.data.frame(BIR_XP_gff_CV_uniq_XP)
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP)
AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP)

colnames(BIR_XP_gff_CG_uniq_XP)[1] <- "protein_id" 
colnames(BIR_XP_gff_CV_uniq_XP)[1] <- "protein_id" 
colnames(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP)[1] <- "protein_id" 
colnames(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP)[1] <- "protein_id" 

BIR_XP_gff_CG_uniq_XP_XM <- left_join(BIR_XP_gff_CG_uniq_XP, C_gig_rtracklayer)
BIR_XP_gff_CV_uniq_XP_XM <- left_join(BIR_XP_gff_CV_uniq_XP, C_vir_rtracklayer)
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM <- left_join(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP, C_gig_rtracklayer)
AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM <- left_join(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP, C_vir_rtracklayer)

BIR_XP_gff_CG_uniq_XP_XM <- BIR_XP_gff_CG_uniq_XP_XM %>% distinct(protein_id, Parent, .keep_all = TRUE)
BIR_XP_gff_CV_uniq_XP_XM <- BIR_XP_gff_CV_uniq_XP_XM %>% distinct(protein_id, Parent, .keep_all = TRUE)
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM <- AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM %>% distinct(protein_id, Parent, .keep_all = TRUE)
AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM <- AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM %>% distinct(protein_id, Parent, .keep_all = TRUE)

# Remove "rna-" from Parent column for CG 
BIR_XP_gff_CG_uniq_XP_XM$Parent <- str_remove(BIR_XP_gff_CG_uniq_XP_XM$Parent, "rna-")
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM$Parent <- str_remove(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM $Parent, "rna-")

colnames(BIR_XP_gff_CG_uniq_XP_XM)[20] <- "transcript_id"
BIR_XP_gff_CG_uniq_XP_XM <- BIR_XP_gff_CG_uniq_XP_XM[,-23]

BIR_XP_gff_CV_uniq_XP_XM <- BIR_XP_gff_CV_uniq_XP_XM[,-10]
colnames(BIR_XP_gff_CV_uniq_XP_XM)[23] <- "ID"

colnames(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM)[20] <- "transcript_id"
AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM <- AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM[,-23]

AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM <- AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM[,-10]
colnames(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM)[23] <- "ID"

# EXPORT XM XP lists - use this to join IAP transcripts to apoptosis list to serch for in DESeq2 and WGCNA
save(BIR_XP_gff_CG_uniq_XP_XM, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_CG_uniq_XP_XM.Rdata")
save(BIR_XP_gff_CV_uniq_XP_XM, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_CV_uniq_XP_XM.Rdata")
save(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP_XM.Rdata")
save(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP_XM.Rdata")

# How many IAP proteins in each 
length(BIR_XP_gff_CG_uniq_XP) # 74
length(BIR_XP_gff_CV_uniq_XP) # 164 (is 158 later after haplotigs are collapsed)
length(AIG1_XP_ALL_gff_GIMAP_CG_uniq_XP) # 31
length(AIG1_XP_ALL_gff_GIMAP_CV_uniq_XP) # 109

# how many uncharacterized?
BIR_XP_gff_CG %>% distinct(protein_id, .keep_all = TRUE) %>% filter(grepl("uncharacterized", product)) # 15
BIR_XP_gff_CV %>% distinct(protein_id, .keep_all = TRUE) %>% filter(grepl("uncharacterized", product)) # 39

#Count IAP genes across species
BIR_XP_gff_species_gene_count <- BIR_XP_gff_species %>% group_by(Species) %>% filter(is.na(locus_tag)) %>% distinct(gene)  %>% summarise(gene_count = n())
BIR_XP_gff_species_locus_tag_count <- BIR_XP_gff_species %>% group_by(Species) %>% filter(is.na(gene)) %>%  distinct(locus_tag) %>% summarise(locus_tag_count = n())
colnames(BIR_XP_gff_species_locus_tag_count)[2] <- "gene_count"
BIR_XP_gff_species_gene_locus_tag_count <- rbind(BIR_XP_gff_species_gene_count, BIR_XP_gff_species_locus_tag_count)

# total gene count across species 
sum(BIR_XP_gff_species_gene_locus_tag_count$gene_count) # 380


#Count GIMAP and IAN genes across species to compare with Lu et al. 2020 paper
AIG1_XP_ALL_gff_GIMAP_species_gene_count <- AIG1_XP_ALL_gff_GIMAP_species %>% group_by(Species) %>% filter(is.na(locus_tag)) %>% distinct(gene)  %>% summarise(gene_count = n())
AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count <- AIG1_XP_ALL_gff_GIMAP_species %>% group_by(Species) %>% filter(is.na(gene)) %>%  distinct(locus_tag) %>% summarise(locus_tag_count = n())
colnames(AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count)[2] <- "gene_count"
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_count <- rbind(AIG1_XP_ALL_gff_GIMAP_species_gene_count, AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count)
# All are between five and 1 over

## EXPORT GENE LISTS PER SPECIES TO EXAMINE POTENTIAL ARTIFACTS
BIR_XP_gff_species_genes <- BIR_XP_gff_species %>% filter(is.na(locus_tag)) %>% distinct(gene, Species)
BIR_XP_gff_species_locus_tag <- BIR_XP_gff_species %>% filter(is.na(gene)) %>%  distinct(locus_tag, Species) 
colnames(BIR_XP_gff_species_locus_tag)[3] <- "gene"
BIR_XP_gff_species_gene_locus_tag <- rbind(BIR_XP_gff_species_genes, BIR_XP_gff_species_locus_tag)

length(unique(BIR_XP_gff_species_gene_locus_tag$gene)) # 380
write.table(unique(BIR_XP_gff_species_gene_locus_tag$gene), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_genes_HMMER_Interpro_BIR.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

AIG1_XP_ALL_gff_GIMAP_species_gene <- AIG1_XP_ALL_gff_GIMAP_species %>% filter(is.na(locus_tag)) %>% distinct(gene, Species)
AIG1_XP_ALL_gff_GIMAP_species_locus_tag <- AIG1_XP_ALL_gff_GIMAP_species  %>% filter(is.na(gene)) %>%  distinct(locus_tag, Species)
colnames(AIG1_XP_ALL_gff_GIMAP_species_locus_tag )[3] <- "gene"
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag <- rbind(AIG1_XP_ALL_gff_GIMAP_species_gene, AIG1_XP_ALL_gff_GIMAP_species_locus_tag)
length(unique(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$gene)) #252

write.table(unique(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$gene), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_genes_HMMER_Interpro_AIG.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Convert locus tag and gene LOC to Extrez ID to get sequences with Batch join by Name column 
colnames(BIR_XP_gff_species_gene_locus_tag)[3] <- "Name"
BIR_XP_gff_species_gene_locus_tag_convert <- left_join(BIR_XP_gff_species_gene_locus_tag, All_mollusc_gene_gff[,c("Name","Dbxref","start","end")])
BIR_XP_gff_species_gene_locus_tag_convert$Dbxref <- str_remove(BIR_XP_gff_species_gene_locus_tag_convert$Dbxref,"GeneID:")
BIR_XP_gff_species_gene_locus_tag_convert <- BIR_XP_gff_species_gene_locus_tag_convert %>% filter(Dbxref != "character(0)")

colnames(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag)[3] <- "Name"
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert <-  left_join(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag, All_mollusc_gene_gff[,c("Name","Dbxref",'start',"end")])
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert$Dbxref <- str_remove(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert$Dbxref,"GeneID:")
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert <- AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert %>% filter(Dbxref != "character(0)")

# Export All Entrez Gene IDs from all species except Elychlor
write.table(BIR_XP_gff_species_gene_locus_tag_convert$Dbxref, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_genes_HMMER_Interpro_BIR_except_Elychlor.txt",
            quote = FALSE,col.names = FALSE, row.names=FALSE)
write.table(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert$Dbxref, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_genes_HMMER_Interpro_AIG_except_Elychlor.txt",
            quote = FALSE,col.names = FALSE, row.names=FALSE)

# Export elysia chlorotica sequences start and end
BIR_XP_gff_species_gene_locus_tag_convert_Elchlor <- BIR_XP_gff_species_gene_locus_tag_convert %>% filter(Species == "Elysia_chlorotica") %>% select(Name,start,end)
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert_Elchlor <-AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert %>% filter(Species == "Elysia_chlorotica") %>% select(Name,start,end)

write.table(BIR_XP_gff_species_gene_locus_tag_convert_Elchlor[,c("Name","start","end")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_genes_HMMER_Interpro_BIR_Elychlor.txt",
            quote = FALSE,col.names = FALSE, row.names=FALSE)
write.table(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_convert_Elchlor[,c("Name","start","end")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_genes_HMMER_Interpro_AIG_Elychlor.txt",
            quote = FALSE,col.names = FALSE, row.names=FALSE)

# Export gene lists by species 
#by(BIR_XP_gff_species_gene_locus_tag, BIR_XP_gff_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
#paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/BIR_Gene_list_", i$Species[1], ".txt"), 
#quote = FALSE,col.names = FALSE, row.names=FALSE))

#by(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag, AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
#paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/AIG_GIMAP_Gene_list_", i$Species[1], ".txt"), 
#quote = FALSE,col.names = FALSE, row.names=FALSE))

## EXPORT ONLY C. VIRGINICA GENE LIST AS BED FILE WITH THE START AND END COORDINATES TO LOOK AT MAPPING COVERAGE AND COMPARE IDENTITY
BIR_XP_gff_species_gene_locus_tag_Cvir <- BIR_XP_gff_species_gene_locus_tag %>% filter(Species =="Crassostrea_virginica")
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir <- AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag %>% filter(Species =="Crassostrea_virginica")
colnames(BIR_XP_gff_species_gene_locus_tag_Cvir)[3] <- "gene"
colnames(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir)[3] <- "gene"

BIR_XP_gff_species_gene_locus_tag_Cvir_BED <- left_join(BIR_XP_gff_species_gene_locus_tag_Cvir, unique(All_mollusc_gene_gff[,c("gene","seqid","start","end")]))
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir_BED <- left_join(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir, unique(All_mollusc_gene_gff[,c("gene","seqid","start","end")]))

IAP_BED <- BIR_XP_gff_species_gene_locus_tag_Cvir_BED[,c(3:5)]
GIMAP_BED <- AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir_BED[,c(3:5)]

# Write out table
write.table(IAP_BED, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_IAP_EMR.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE)
write.table(GIMAP_BED, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_GIMAP_EMR.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE)

IAP_BED_name  <- BIR_XP_gff_species_gene_locus_tag_Cvir_BED[,c(3:5,1)]
GIMAP_BED_name <- AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir_BED[,c(3:5,1)]

# Write out table
write.table(IAP_BED_name, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_IAP_EMR_Name.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(GIMAP_BED_name, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_GIMAP_EMR_Name.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")

save(IAP_BED_name, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_IAP_EMR_Name.Rdata")
save(GIMAP_BED_name, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_GIMAP_EMR_Name.Rdata")


## Review Matches
BIR_XP_gff_species 
AIG1_XP_ALL_gff_GIMAP_species
View(unique(BIR_XP_gff_species$product)) # has one phosphatase and actin regulator 4-B-like but does have the IAP repeats 
View(unique(AIG1_XP_ALL_gff_GIMAP_species$product))

# export this list 
save(BIR_XP_gff_species, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species.Rdata")

##################### ORTHOGROUP ANALYSIS ###############################

### Load May15th Orthogroup Analysis of 10 Mollusc species from Orthogroup.tsv file ###
# Load tsv
#Orthogroups <- read_tsv("/Volumes/My Passport for Mac/OrthoFinder_3_25_2020_Bluewaves_Backup/Results_May15/Orthogroups/Orthogroups.tsv",
 #                       col_names = c("Orthogroup","Elysia_chlorotica", "Aplysia_californica", "Crassostrea_gigas", "Lottia_gigantea", 
#                                      "Biomphalaria_glabrata", "Octopus_bimaculoides",
#                                      "C_virginica", "Mizuhopecten_yessoensis",	"Pomacea_canaliculata",	"Octopus_sinensis"))

#### ASSESS ORTHOGROUPS USING ONLY THE CV AND CG REFERENCE ANNOTATIONS ####
CV_CG_IAP_list <- as.list(CV_CG_IAP)
CV_CG_GIMAP_list <- as.list(CV_CG_GIMAP)

#CV_CG_IAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_IAP_list, collapse="|"), i))),]
#length(CV_CG_IAP_list_lookup$Orthogroup)
# 27 orthogroups

#CV_CG_GIMAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_GIMAP_list, collapse="|"), i))),]
#length(CV_CG_GIMAP_list_lookup$Orthogroup)
#9

#### USE FULL IAP AND GIMAP LISTS TO PULL OUT ALL MOLLUSC ORTHOGROUPS ####
BIR_XP_gff_species_list <- as.list(unique(BIR_XP_gff_species$protein_id))
#BIR_XP_gff_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_species_list, collapse="|"), i))),]
#length(BIR_XP_gff_species_list_lookup$Orthogroup)
# 35 orthogroups (got one extra orthogroup when setting hmmer to eval -3)

#Compare to list from original orthogroup search using only proteins in annotation
#setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_species_list_lookup$Orthogroup) # 4 missed "OG0003807" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
#setdiff(BIR_XP_gff_species_list_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # 13 added "OG0001642" "OG0002611" "OG0004344" "OG0007118" "OG0011865" "OG0011926" "OG0012919" "OG0013878" "OG0014276" "OG0016483" "OG0017158" "OG0017983" "OG0018491"

AIG1_XP_ALL_gff_GIMAP_species_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_species$protein_id))
length(AIG1_XP_ALL_gff_GIMAP_species_list)
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_XP_ALL_gff_GIMAP_species_list, collapse="|"), i))),]
#length(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup)
#12

#setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup) # 2 not found "OG0013109" "OG0016155"
#setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 5 added

#### USE CG AND CV SPECIFIC IAP AND GIMAP HMMER/INTERPROSCAN LISTS TO PULL OUT ORTHOGROUPS ###
BIR_XP_gff_CG_CV <- rbind(BIR_XP_gff_CG, BIR_XP_gff_CV)
AIG1_XP_ALL_gff_GIMAP_CG_CV <- rbind(AIG1_XP_ALL_gff_GIMAP_CG, AIG1_XP_ALL_gff_GIMAP_CV)

# Use CV and CG only lists to pull out orthogroups
BIR_XP_gff_CG_CV_list <- as.list(unique(BIR_XP_gff_CG_CV$protein_id))
#BIR_XP_gff_CG_CV_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_CG_CV_list, collapse="|"), i))),]
#length(BIR_XP_gff_CG_CV_lookup$Orthogroup)
# 26 orthogroups, 11 are added when you include all the mollusc proteins

#Compare to list from original orthogroup search using only proteins in annotation
#setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_CG_CV_lookup$Orthogroup) #  "OG0003807" "OG0006831" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
#setdiff(BIR_XP_gff_CG_CV_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # "OG0004344" "OG0011865" "OG0012919" "OG0013878" "OG0017983"

AIG1_CDD_GIMAP_only_CV_CG_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_CG_CV$protein_id))
length(AIG1_CDD_GIMAP_only_CV_CG_list) # 140
#AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_CDD_GIMAP_only_CV_CG_list, collapse="|"), i))),]
#length(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup)
# 9 orthogroups,

#Compare to list from original orthogroup search using only proteins in annotation
#setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup) # "OG0013109" "OG0016155"
#setdiff(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 0"OG0007435" "OG0014793"

#### GET ALL MOLLUSCS IAP ORTHOGROUP HITS #####

## IAP genes 
# Get full list of proteins for each species by transposing and uniting
# Transpose the rows and column 
#BIR_XP_gff_species_list_lookup_transpose <- t(BIR_XP_gff_species_list_lookup)
#class(BIR_XP_gff_species_list_lookup_transpose) # matrix
#BIR_XP_gff_species_list_lookup_transpose <- as.data.frame(BIR_XP_gff_species_list_lookup_transpose)
## unite all columns into one column 
#BIR_XP_gff_species_list_lookup_transpose_united <- unite(BIR_XP_gff_species_list_lookup_transpose, full_protein_list, sep=",")
## remove NAs
#BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub("NA,", "",BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
#BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub(",NA", "",BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
#BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub("NA", "", BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
## Put all into single vector for annot and export to make tree
## Concatenate each into single vector
#BIR_XP_gff_species_list_lookup_transpose_united_all <- BIR_XP_gff_species_list_lookup_transpose_united %>% summarise(combined =paste(full_protein_list, collapse=","))
#BIR_XP_gff_species_list_lookup_transpose_united_all_col <- data.frame(protein_id = unlist(strsplit(as.character(BIR_XP_gff_species_list_lookup_transpose_united_all$combined), ",")))
## trimws and remove orthogroups
#BIR_XP_gff_species_list_lookup_united_all_col <- BIR_XP_gff_species_list_lookup_transpose_united_all_col[-c(1:35),1]
#BIR_XP_gff_species_list_lookup_united_all_col <- trimws(BIR_XP_gff_species_list_lookup_united_all_col, which="left")

# How many XPs identified?
#length(BIR_XP_gff_species_list_lookup_united_all_col) # 601
#length(BIR_XP_gff_species_list) # 791

# Are there any duplicated proteins found in orthogroups?
#BIR_XP_gff_species_list_lookup_united_all_col[duplicated(BIR_XP_gff_species_list_lookup_united_all_col)] # 0 duplicated

# Were all the original proteins found? - NO
#setdiff(BIR_XP_gff_species_list,  BIR_XP_gff_species_list_lookup_united_all_col) # 272 proteins are present in original HMMER list that are not in the Orthogroup list
#length(setdiff(BIR_XP_gff_species_list_lookup_united_all_col, BIR_XP_gff_species_list)) # 82 are added in that are not in the original orthogroup list

## Find genes in original HMMER list for each protein NOTE THAT Elysia chlorotica, Lottia gigantea ONLY HAVE LOCUS TAGS AND NOT GENES
BIR_XP_gff_species_join <- left_join(unique(BIR_XP_gff_species[,c("protein_id","product","Species")]), All_molluscs_CDS_gff[,c("protein_id","gene")])
BIR_XP_gff_species_join <- left_join(BIR_XP_gff_species_join[,c("protein_id","product","gene","Species")], All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene")])
BIR_XP_gff_species_join <- unique(BIR_XP_gff_species_join)
  # are any still unfilled?
#BIR_XP_gff_species_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
# remove duplicates
BIR_XP_gff_species_genes <- BIR_XP_gff_species_join[!duplicated(BIR_XP_gff_species_join[,c("gene","locus_tag")]),]

# How many total genes or locus tags?
nrow(BIR_XP_gff_species_genes) # 380

## Find genes in full Orthogroup list and compare
#BIR_XP_gff_species_list_lookup_united_all_col_df <- as.data.frame(BIR_XP_gff_species_list_lookup_united_all_col)
#colnames(BIR_XP_gff_species_list_lookup_united_all_col_df)[1]<-"protein_id"
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- left_join(BIR_XP_gff_species_list_lookup_united_all_col_df, All_molluscs_CDS_gff[,c("protein_id","gene","product")])
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- left_join(BIR_XP_gff_species_list_lookup_united_all_col_df_genes, All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene","product")])
## are any still unfilled?
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
## remove duplicates
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes[!duplicated(BIR_XP_gff_species_list_lookup_united_all_col_df_genes[,c("gene","locus_tag")]),]
## how many genes
#nrow(BIR_XP_gff_species_list_lookup_united_all_col_df_genes ) # 308

## Write out to table the gene list to use for gathering sequences for MAFFT 
# Collapse the locus tag and gene columns into one
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes %>% 
#        unite(gene_locus_tag, c("gene","locus_tag"))
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "NA_")
#BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "_NA")
#
#write.table(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_mollusc_orthogroups_gene_locus_tag_list.txt",
#            row.names = FALSE, col.names = FALSE, quote=FALSE)

### GIMAP genes
# Get full list of proteins for each species by transposing and uniting
# Transpose the rows and column 
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose <- t(AIG1_XP_ALL_gff_GIMAP_species_list_lookup)
#class(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose) # matrix
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose)
## unite all columns into one column 
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united <- unite(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose, full_protein_list, sep=",")
## remove NAs
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub("NA,", "",AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub(",NA", "",AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub("NA", "", AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
## Put all into single vector for annot and export to make tree
## Concatenate each into single vector
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united %>% summarise(combined =paste(full_protein_list, collapse=","))
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all_col <- data.frame(protein_id = unlist(strsplit(as.character(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all$combined), ",")))
## trimws and remove orthogroups
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all_col[-c(1:12),1]
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col <- trimws(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col, which="left")
#
# How many XPs identified?
#length(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col) # 438
#length(AIG1_XP_ALL_gff_GIMAP_species_list ) # 403

# Are there any duplicated proteins found in orthogroups?
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col[duplicated(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col)] #0 duplicated

# Were all the original proteins found? - NO
#setdiff(AIG1_XP_ALL_gff_GIMAP_species_list,  AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col) # 132 proteins are present in original HMMER list that are not in the Orthogroup list
#length(setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col, AIG1_XP_ALL_gff_GIMAP_species_list)) # 167 are added in that are not in the original orthogroup list

## Find genes in original HMMER list for each protein NOTE THAT Elysia chlorotica, Lottia gigantea ONLY HAVE LOCUS TAGS AND NOT GENES
AIG1_XP_ALL_gff_GIMAP_species_join <- left_join(unique(AIG1_XP_ALL_gff_GIMAP_species[,c("protein_id","product","Species")]), All_molluscs_CDS_gff[,c("protein_id","gene")])
AIG1_XP_ALL_gff_GIMAP_species_join <- left_join(AIG1_XP_ALL_gff_GIMAP_species_join[,c("protein_id","product","gene","Species")], All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene")])
AIG1_XP_ALL_gff_GIMAP_species_join <- unique(AIG1_XP_ALL_gff_GIMAP_species_join)
# remove duplicates
AIG1_XP_ALL_gff_GIMAP_species_genes <- AIG1_XP_ALL_gff_GIMAP_species_join[!duplicated(AIG1_XP_ALL_gff_GIMAP_species_join[,c("gene","locus_tag")]),]

# How many total genes or locus tags?
nrow(AIG1_XP_ALL_gff_GIMAP_species_genes) # 252

## Find genes in full Orthogroup list and compare
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col)
#colnames(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df)[1]<-"protein_id"
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- left_join(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df, All_molluscs_CDS_gff[,c("protein_id","gene","product")])
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- left_join(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes, All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene","product")])
## are any still unfilled?
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
## remove duplicates
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[!duplicated(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[,c("gene","locus_tag")]),]
## how many genes
#nrow(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes ) # 308

## Write out to table the gene list to use for gathering sequences for MAFFT 
# Collapse the locus tag and gene columns into one
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes %>% 
#  unite(gene_locus_tag, c("gene","locus_tag"))
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "NA_")
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "_NA")
#
#write.table(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_mollusc_orthogroups_gene_locus_tag_list.txt",
 #           row.names = FALSE, col.names = FALSE, quote=FALSE)

###### COMPARATIVE STATISTICS BETWEEN HMMER, ANNOTATION, AND ORTHOFINDER #####

### IAP ###
## Were all the original HMMER/Interproscan genes found in the Orthofinder results?
#IAP_Orthogroup_missing_genes <- BIR_XP_gff_species_genes[!(BIR_XP_gff_species_genes$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]
#IAP_Orthogroup_missing_genes <- IAP_Orthogroup_missing_genes %>% filter(!is.na(gene)) %>% distinct(gene)
#length(setdiff(BIR_XP_gff_species_genes$locus_tag,  BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag)) # 2  locus tag added
#nrow(IAP_Orthogroup_missing_genes) #120
# 120 genes + 2 locus tags is 122 genes
  # Mostly missing PREDICTED genes, partial genes, and uncharacterized protein genes. All missing are Lottia gigantea and Elysia chlorotica genes

## Were any genes added in Orthogroups that weren't in HMMER?
#IAP_Orthogroup_added_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes[!(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag %in% BIR_XP_gff_species_genes$gene),]
## need to add in line to take into account locus tag vs gene here 
#setdiff(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$locus_tag,  BIR_XP_gff_species_genes$locus_tag) # 0 new locus tags added

## Are all the CV and CG IAP genes found from the genome in the Orthogroup search ?
#View(Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]) # 7 CG genes were missed by Orthofinder that were annotated in genome
#View(Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]) # 12 CV genes were missed by Orthofinder that were annotated in genome

## Were any genes added by HMMER that were not in the genome? 
BIR_XP_gff_species_genes_CG <- BIR_XP_gff_species_genes %>% filter(Species=="Crassostrea_gigas")
BIR_XP_gff_species_genes_CG[!(BIR_XP_gff_species_genes_CG$gene %in% Cgig_gff_IAP_family_XP$gene),] # 9 uncharacterized Loci genes were added by HMMER that were not annotated in genome

BIR_XP_gff_species_genes_CV <- BIR_XP_gff_species_genes %>% filter(Species=="Crassostrea_virginica")
BIR_XP_gff_species_genes_CV[!(BIR_XP_gff_species_genes_CV$gene %in% Cvir_gff_IAP_family_XP$gene),] # 14 uncharacterized Loci genes were added by HMMER that were not annotated in genome

## Are all the CV and CG IAP genes found from the genome in the HMMER/Interproscan search ?
Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_genes$gene),] # 4 CG genes were not in HMMER that were annotated in genome
    #LOC105333301 # only the zinc finger domain
    #LOC105336740 # only the RingUbox domain 
    #LOC105338773 # only has the UBCc no BIR
    #LOC105338774 # only has UBCc domains no BIR
Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_genes$gene),] # 6 CV genes were not in by HMMER that were annotated in genome
    #LOC111136287 # only the pfam zinc ring finger domain        
    #LOC111101682 # only the RingUbox domain no BIR repeat           
    #LOC111100802 # has the RING-HC_BIRC2_3_7; RING finger, HC subclass, found in apoptosis protein c-IAP1, c-IAP2, livin, and similar proteins but not the BIR repeat        
    #LOC111104430 # no BIR repeat only RingUbox in the NCBI domains            
    #LOC111109770 # no BIR repeat only RingUbox in the NCBI domains            
    #LOC111106726 # only has the RINGUbox           

## GIMAP ###
## Were all the original HMMER/Interproscan genes found in the Orthofinder results?
#GIMAP_Orthogroup_missing_genes <- AIG1_XP_ALL_gff_GIMAP_species_genes[!(AIG1_XP_ALL_gff_GIMAP_species_genes$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]
#GIMAP_Orthogroup_missing_genes <- GIMAP_Orthogroup_missing_genes %>% filter(!is.na(gene)) %>% distinct(gene)
#length(setdiff(AIG1_XP_ALL_gff_GIMAP_species_genes$locus_tag,  AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag)) # 3  locus tag added
#nrow(GIMAP_Orthogroup_missing_genes) #65
# 65 genes + 3 locus tags is 68 genes

## Were any genes added in Orthogroups that weren't in HMMER?
#GIMAP_Orthogroup_added_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[!(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag %in% AIG1_XP_ALL_gff_GIMAP_species_genes$gene),]
# need to add in line to take into account locus tag vs gene here 
#setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$locus_tag,  AIG1_XP_ALL_gff_GIMAP_species_genes$locus_tag) # 0 new locus tags added

## Are all the CV and CG IAP genes found from the genome in the Orthogroup search ?
#GIMAP_missing_CG<- Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),] #  CG genes were missed by Orthofinder that were annotated in genome
#length(unique(GIMAP_missing_CG$gene)) #35
#GIMAP_missing_CV <- Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),] #  CV genes were missed by Orthofinder that were annotated in genome
#length(unique(GIMAP_missing_CV$gene)) #67

## Were any genes added by HMMER that were not in the genome? 
AIG1_XP_ALL_gff_GIMAP_species_genes_CG <- AIG1_XP_ALL_gff_GIMAP_species_genes %>% filter(Species=="Crassostrea_gigas")
GIMAP_HMMER_v_annot_CG <- AIG1_XP_ALL_gff_GIMAP_species_genes_CG[!(AIG1_XP_ALL_gff_GIMAP_species_genes_CG$gene %in% Cgig_gff_IAP_family_XP$gene),] #  uncharacterized Loci genes were added by HMMER that were not annotated in genome
length(unique(GIMAP_HMMER_v_annot_CG$gene)) # 27

AIG1_XP_ALL_gff_GIMAP_species_genes_CV <- AIG1_XP_ALL_gff_GIMAP_species_genes %>% filter(Species=="Crassostrea_virginica")
GIMAP_HMMER_v_annot_CV <- AIG1_XP_ALL_gff_GIMAP_species_genes_CV[!(AIG1_XP_ALL_gff_GIMAP_species_genes_CV$gene %in% Cvir_gff_IAP_family_XP$gene),] # uncharacterized Loci genes were added by HMMER that were not annotated in genome
length(unique(GIMAP_HMMER_v_annot_CV$gene))

## Are all the CV and CG IAP genes found from the genome in the HMMER/Interproscan search ?
GIMAP_annot_v_HMMER_CG <- Cgig_gff_GIMAP_family_XP[!(Cgig_gff_GIMAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_genes$gene),] # 4 CG genes were not in HMMER that were annotated in genome
length(unique(GIMAP_annot_v_HMMER_CG$gene))

GIMAP_annot_v_HMMER_CV <- Cvir_gff_GIMAP_family_XP[!(Cvir_gff_GIMAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_genes$gene),] # 6 CV genes were not in by HMMER that were annotated in genome
length(unique(GIMAP_annot_v_HMMER_CV$gene))

### RUN PROTEIN SEQUENCES IN MAFFT AND RAXML ###
 # RAxML and HMMER full results used

#### IDENTICAL PROTEINS REMOVED BY CD-HIT ####
# Remember from investigation of CD-HITs behavior that sequences that are shorter but exactly the same are still removed as well. This could be explaining the discrepancy in the numbers below in the software 
# Load removed duplicates file
# phylotools reads the sequences in as a dataframe, treeio reads in as vectors
AIG_seq_rm_dup <- treeio::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa")
BIR_seq_rm_dup <- treeio:::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup.fa")
AIG_seq_rm_dup_phylo <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa")
BIR_seq_rm_dup_phylo <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup.fa")

length(names(AIG_seq_rm_dup)) # 313
length(names(BIR_seq_rm_dup)) # 499

# get protein IDs
AIG_seq_rm_dup_headers <- labels(AIG_seq_rm_dup)
AIG_seq_rm_dup_headers <- str_split_fixed(AIG_seq_rm_dup_headers, " ", 2)
AIG_seq_rm_dup_headers <- as.data.frame(AIG_seq_rm_dup_headers)
colnames(AIG_seq_rm_dup_headers)[1] <- "protein_id"
BIR_seq_rm_dup_headers <- labels(BIR_seq_rm_dup)
BIR_seq_rm_dup_headers <- str_split_fixed(BIR_seq_rm_dup_headers, " ", 2)
BIR_seq_rm_dup_headers <- as.data.frame(BIR_seq_rm_dup_headers)
colnames(BIR_seq_rm_dup_headers)[1] <- "protein_id"

# Compare protein IDs
AIG1_dup_seq_rm <- AIG1_XP_ALL_gff_GIMAP_species_join[!(AIG1_XP_ALL_gff_GIMAP_species_join$protein_id %in% AIG_seq_rm_dup_headers$protein_id),] # 90 proteins removed for being duplicated
BIR_dup_seq_rm <- BIR_XP_gff_species_join[!(BIR_XP_gff_species_join$protein_id %in% BIR_seq_rm_dup_headers$protein_id),] # 20 removed

# Protein IDs that were kept 
AIG1_dup_seq_rm_kept <- AIG1_XP_ALL_gff_GIMAP_species_join[AIG1_XP_ALL_gff_GIMAP_species_join$protein_id %in% AIG_seq_rm_dup_headers$protein_id,] # 90 proteins removed for being duplicated
BIR_dup_seq_rm_kept <- BIR_XP_gff_species_join[BIR_XP_gff_species_join$protein_id %in% BIR_seq_rm_dup_headers$protein_id,] # 20 removed

## Parse the CD-HIT cluster file
# Followed code from this site: https://rpubs.com/rmurdoch/cdhit_to_mapping_file
AIG_seq_rm_dup_clustering <- phylotools:::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa.clstr")
BIR_seq_rm_dup_clustering <- phylotools:::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup.fa.clstr")

AIG_seq_rm_dup_clstr <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr)
BIR_seq_rm_dup_clstr <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr)

# parse
AIG_seq_rm_dup_clstr2 <- AIG_seq_rm_dup_clstr
n = nrow(AIG_seq_rm_dup_clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(AIG_seq_rm_dup_clstr2[row,1]) == TRUE) {
    AIG_seq_rm_dup_clstr2[row,1] <- x}
  else {NULL}
  x <- AIG_seq_rm_dup_clstr2[row,1]
}
head(AIG_seq_rm_dup_clstr2 )

BIR_seq_rm_dup_clstr2 <- BIR_seq_rm_dup_clstr
n = nrow(BIR_seq_rm_dup_clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(BIR_seq_rm_dup_clstr2[row,1]) == TRUE) {
    BIR_seq_rm_dup_clstr2[row,1] <- x}
  else {NULL}
  x <- BIR_seq_rm_dup_clstr2[row,1]
}
head(BIR_seq_rm_dup_clstr2)

# remove rows with empty V2
AIG_seq_rm_dup_clstr4 <- AIG_seq_rm_dup_clstr2[-which(AIG_seq_rm_dup_clstr2$V2 == ""), ]
BIR_seq_rm_dup_clstr4 <- BIR_seq_rm_dup_clstr2[-which(BIR_seq_rm_dup_clstr2$V2 == ""), ]

# remove “>” * split V2 by “aa,” and by “… at” and by “%,”
AIG_seq_rm_dup_clstr5 <- AIG_seq_rm_dup_clstr4
AIG_seq_rm_dup_clstr5[] <- lapply(AIG_seq_rm_dup_clstr5, gsub, pattern='>', replacement='')
AIG_seq_rm_dup_clstr5.2 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5$V2, "aa, ", 2))
AIG_seq_rm_dup_clstr5.3 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5.2$X2, "... ", 2))
AIG_seq_rm_dup_clstr6 <- cbind(AIG_seq_rm_dup_clstr5[1],AIG_seq_rm_dup_clstr5.2[1],AIG_seq_rm_dup_clstr5.3[1:2])
colnames(AIG_seq_rm_dup_clstr6) <- c("cluster","aa","protein_id","stat")
head(AIG_seq_rm_dup_clstr6)

BIR_seq_rm_dup_clstr5 <- BIR_seq_rm_dup_clstr4
BIR_seq_rm_dup_clstr5[] <- lapply(BIR_seq_rm_dup_clstr5, gsub, pattern='>', replacement='')
BIR_seq_rm_dup_clstr5.2 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5$V2, "aa, ", 2))
BIR_seq_rm_dup_clstr5.3 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5.2$X2, "... ", 2))
BIR_seq_rm_dup_clstr6 <- cbind(BIR_seq_rm_dup_clstr5[1],BIR_seq_rm_dup_clstr5.2[1],BIR_seq_rm_dup_clstr5.3[1:2])
colnames(BIR_seq_rm_dup_clstr6) <- c("cluster","aa","protein_id","stat")
head(BIR_seq_rm_dup_clstr6)

# The representative sequence is indicated by a "*"

# Join with information about product and gene
AIG_seq_rm_dup_clstr6 <- left_join(AIG_seq_rm_dup_clstr6, AIG1_XP_ALL_gff_GIMAP_species_join[,c("protein_id","gene","product","locus_tag","Species")], by = "protein_id")
BIR_seq_rm_dup_clstr6 <- left_join(BIR_seq_rm_dup_clstr6, BIR_XP_gff_species_join[,c("protein_id","gene","product","locus_tag","Species")], by = "protein_id")

# Check if any proteins collapsed came from different genes 
# look at duplicated clusters and keep all rows 
AIG_seq_rm_dup_clstr6_dup <- AIG_seq_rm_dup_clstr6 %>% group_by(cluster) %>% filter(n() > 1)
BIR_seq_rm_dup_clstr6_dup <- BIR_seq_rm_dup_clstr6 %>% group_by(cluster) %>% filter(n() > 1)

# Look at distinct clusters where there may be two different genes in a cluster
AIG_seq_rm_dup_clstr6_dup_diff_gene <-AIG_seq_rm_dup_clstr6_dup %>% distinct(cluster,gene,locus_tag) %>% group_by(cluster) %>% filter(n() > 1) 
    # 0 diff genes in same cluster for AIG1
BIR_seq_rm_dup_clstr6_dup_diff_gene <- BIR_seq_rm_dup_clstr6_dup %>% distinct(cluster,gene,locus_tag) %>% group_by(cluster) %>% filter(n() > 1)
BIR_seq_rm_dup_clstr6_dup_diff_gene_lookup <- unique(BIR_seq_rm_dup_clstr6_dup_diff_gene$cluster)
BIR_seq_rm_dup_clstr6_dup_diff_gene_product <- BIR_seq_rm_dup_clstr6[BIR_seq_rm_dup_clstr6$cluster %in% BIR_seq_rm_dup_clstr6_dup_diff_gene_lookup,]

# subset for C virginica
BIR_seq_rm_dup_clstr6_dup_diff_gene_product_Cvir <- BIR_seq_rm_dup_clstr6_dup_diff_gene_product %>% filter(Species == "Crassostrea_virginica")

# Join start and end coordinates to investigate if they are tandem duplicates
BIR_seq_rm_dup_clstr6_dup_diff_gene_product_Cvir <- left_join(BIR_seq_rm_dup_clstr6_dup_diff_gene_product_Cvir, unique(All_mollusc_gene_gff[,c("gene","seqid","start","end")]))

## all are the same protein names, going to keep the tree annotation the same!

#### IDENTIFY 95 % IDENTICAL PROTIENS CD-HIT ####
## Parse the CD-HIT cluster file
# Followed code from this site: https://rpubs.com/rmurdoch/cdhit_to_mapping_file
AIG_seq_rm_dup_clstr_95 <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup_95.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr_95)
BIR_seq_rm_dup_clstr_95 <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup_95.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr_95)

# parse
AIG_seq_rm_dup_clstr2_95 <- AIG_seq_rm_dup_clstr_95
n = nrow(AIG_seq_rm_dup_clstr_95)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(AIG_seq_rm_dup_clstr2_95[row,1]) == TRUE) {
    AIG_seq_rm_dup_clstr2_95[row,1] <- x}
  else {NULL}
  x <- AIG_seq_rm_dup_clstr2_95[row,1]
}

BIR_seq_rm_dup_clstr2_95 <- BIR_seq_rm_dup_clstr_95
n = nrow(BIR_seq_rm_dup_clstr_95)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(BIR_seq_rm_dup_clstr2_95[row,1]) == TRUE) {
    BIR_seq_rm_dup_clstr2_95[row,1] <- x}
  else {NULL}
  x <- BIR_seq_rm_dup_clstr2_95[row,1]
}

# remove rows with empty V2
AIG_seq_rm_dup_clstr4_95 <- AIG_seq_rm_dup_clstr2_95[-which(AIG_seq_rm_dup_clstr2_95$V2 == ""), ]
BIR_seq_rm_dup_clstr4_95 <- BIR_seq_rm_dup_clstr2_95[-which(BIR_seq_rm_dup_clstr2_95$V2 == ""), ]

# remove “>” * split V2 by “aa,” and by “… at” and by “%,”
AIG_seq_rm_dup_clstr5_95 <- AIG_seq_rm_dup_clstr4_95
AIG_seq_rm_dup_clstr5_95[] <- lapply(AIG_seq_rm_dup_clstr5_95, gsub, pattern='>', replacement='')
AIG_seq_rm_dup_clstr5.2_95 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5_95$V2, "aa, ", 2))
AIG_seq_rm_dup_clstr5.3_95 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5.2_95$X2, "... ", 2))
AIG_seq_rm_dup_clstr6_95 <- cbind(AIG_seq_rm_dup_clstr5_95[1],AIG_seq_rm_dup_clstr5.2_95[1],AIG_seq_rm_dup_clstr5.3_95[1:2])
colnames(AIG_seq_rm_dup_clstr6_95) <- c("cluster","aa","protein_id","stat")
head(AIG_seq_rm_dup_clstr6_95)

BIR_seq_rm_dup_clstr5_95 <- BIR_seq_rm_dup_clstr4_95
BIR_seq_rm_dup_clstr5_95[] <- lapply(BIR_seq_rm_dup_clstr5_95, gsub, pattern='>', replacement='')
BIR_seq_rm_dup_clstr5.2_95 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5_95$V2, "aa, ", 2))
BIR_seq_rm_dup_clstr5.3_95 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5.2_95$X2, "... ", 2))
BIR_seq_rm_dup_clstr6_95 <- cbind(BIR_seq_rm_dup_clstr5_95[1],BIR_seq_rm_dup_clstr5.2_95[1],BIR_seq_rm_dup_clstr5.3_95[1:2])
colnames(BIR_seq_rm_dup_clstr6_95) <- c("cluster","aa","protein_id","stat")
head(BIR_seq_rm_dup_clstr6_95)

# The representative sequence is indicated by a "*"

# Join with information about product and gene
AIG_seq_rm_dup_clstr6_95 <- left_join(AIG_seq_rm_dup_clstr6_95, AIG1_XP_ALL_gff_GIMAP_species_join[,c("protein_id","gene","product","locus_tag","Species")], by = "protein_id")
BIR_seq_rm_dup_clstr6_95 <- left_join(BIR_seq_rm_dup_clstr6_95, BIR_XP_gff_species_join[,c("protein_id","gene","product","locus_tag","Species")], by = "protein_id")

# Check if any proteins collapsed came from different genes 
# look at duplicated clusters and keep all rows 
AIG_seq_rm_dup_clstr6_dup_95 <- AIG_seq_rm_dup_clstr6_95 %>% group_by(cluster) %>% filter(n() > 1)
BIR_seq_rm_dup_clstr6_dup_95 <- BIR_seq_rm_dup_clstr6_95 %>% group_by(cluster) %>% filter(n() > 1)

# Look at distinct clusters where there may be two different genes in a cluster
AIG_seq_rm_dup_clstr6_dup_diff_gene_95 <- AIG_seq_rm_dup_clstr6_dup_95 %>% distinct(cluster,gene,locus_tag) %>% group_by(cluster) %>% filter(n() > 1) 
AIG_seq_rm_dup_clstr6_dup_diff_gene_95_lookup_95 <- unique(AIG_seq_rm_dup_clstr6_dup_diff_gene_95 $cluster)
AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95 <- AIG_seq_rm_dup_clstr6_95[AIG_seq_rm_dup_clstr6_95$cluster %in% AIG_seq_rm_dup_clstr6_dup_diff_gene_95_lookup_95,]

BIR_seq_rm_dup_clstr6_dup_diff_gene_95 <- BIR_seq_rm_dup_clstr6_dup_95 %>% distinct(cluster,gene,locus_tag) %>% group_by(cluster) %>% filter(n() > 1)
BIR_seq_rm_dup_clstr6_dup_diff_gene_lookup_95 <- unique(BIR_seq_rm_dup_clstr6_dup_diff_gene_95$cluster)
BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95 <- BIR_seq_rm_dup_clstr6_95[BIR_seq_rm_dup_clstr6_95$cluster %in% BIR_seq_rm_dup_clstr6_dup_diff_gene_lookup_95,]

# export to file for use in markdown for jon and marta
save(AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95.Rdata")
save(BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95.Rdata")


#### IDENTIFY 95-100% IDENTICAL NUCLEOTIDE GENE SEQUENCES CD-HIT ####
## Parse the CD-HIT cluster file
# Followed code from this site: https://rpubs.com/rmurdoch/cdhit_to_mapping_file
AIG_seq_rm_dup_clstr_NUC_95 <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_GIMAP_EMR_Name_CD_Hit_95.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr_NUC_95)
BIR_seq_rm_dup_clstr_NUC_95 <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Cvir_IAP_EMR_Name_CD_Hit_95.fa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
head(AIG_seq_rm_dup_clstr_NUC_95)

# parse
AIG_seq_rm_dup_clstr2_NUC_95 <- AIG_seq_rm_dup_clstr_NUC_95
n = nrow(AIG_seq_rm_dup_clstr_NUC_95)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(AIG_seq_rm_dup_clstr2_NUC_95[row,1]) == TRUE) {
    AIG_seq_rm_dup_clstr2_NUC_95[row,1] <- x}
  else {NULL}
  x <- AIG_seq_rm_dup_clstr2_NUC_95[row,1]
}

BIR_seq_rm_dup_clstr2_NUC_95 <- BIR_seq_rm_dup_clstr_NUC_95
n = nrow(BIR_seq_rm_dup_clstr_NUC_95)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(BIR_seq_rm_dup_clstr2_NUC_95[row,1]) == TRUE) {
    BIR_seq_rm_dup_clstr2_NUC_95[row,1] <- x}
  else {NULL}
  x <- BIR_seq_rm_dup_clstr2_NUC_95[row,1]
}

# remove rows with empty V2
AIG_seq_rm_dup_clstr4_NUC_95 <- AIG_seq_rm_dup_clstr2_NUC_95[-which(AIG_seq_rm_dup_clstr2_NUC_95$V2 == ""), ]
BIR_seq_rm_dup_clstr4_NUC_95 <- BIR_seq_rm_dup_clstr2_NUC_95[-which(BIR_seq_rm_dup_clstr2_NUC_95$V2 == ""), ]

# remove “>” * split V2 by “aa,” and by “… at” and by “%,”
AIG_seq_rm_dup_clstr5_NUC_95 <- AIG_seq_rm_dup_clstr4_NUC_95
AIG_seq_rm_dup_clstr5_NUC_95[] <- lapply(AIG_seq_rm_dup_clstr5_NUC_95, gsub, pattern='>', replacement='')
AIG_seq_rm_dup_clstr5.2_NUC_95 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5_NUC_95$V2, "aa, ", 2))
AIG_seq_rm_dup_clstr5.3_NUC_95 <- data.frame(str_split_fixed(AIG_seq_rm_dup_clstr5.2_NUC_95$X2, "... ", 2))
AIG_seq_rm_dup_clstr6_NUC_95 <- cbind(AIG_seq_rm_dup_clstr5_NUC_95[1],AIG_seq_rm_dup_clstr5.2_NUC_95[1],AIG_seq_rm_dup_clstr5.3_NUC_95[1:2])
colnames(AIG_seq_rm_dup_clstr6_NUC_95) <- c("cluster","aa","gene","stat")
head(AIG_seq_rm_dup_clstr6_NUC_95)

BIR_seq_rm_dup_clstr5_NUC_95 <- BIR_seq_rm_dup_clstr4_NUC_95
BIR_seq_rm_dup_clstr5_NUC_95[] <- lapply(BIR_seq_rm_dup_clstr5_NUC_95, gsub, pattern='>', replacement='')
BIR_seq_rm_dup_clstr5.2_NUC_95 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5_NUC_95$V2, "aa, ", 2))
BIR_seq_rm_dup_clstr5.3_NUC_95 <- data.frame(str_split_fixed(BIR_seq_rm_dup_clstr5.2_NUC_95$X2, "... ", 2))
BIR_seq_rm_dup_clstr6_NUC_95 <- cbind(BIR_seq_rm_dup_clstr5_NUC_95[1],BIR_seq_rm_dup_clstr5.2_NUC_95[1],BIR_seq_rm_dup_clstr5.3_NUC_95[1:2])
colnames(BIR_seq_rm_dup_clstr6_NUC_95) <- c("cluster","aa","gene","stat")
head(BIR_seq_rm_dup_clstr6_NUC_95)

# The representative sequence is indicated by a "*"

# Join with information about gene
AIG_seq_rm_dup_clstr6_NUC_95 <- left_join(AIG_seq_rm_dup_clstr6_NUC_95, unique(AIG1_XP_ALL_gff_GIMAP_species_join[,c("gene","locus_tag","Species")]), by = "gene")
BIR_seq_rm_dup_clstr6_NUC_95 <- left_join(BIR_seq_rm_dup_clstr6_NUC_95, unique(BIR_XP_gff_species_join[,c("gene","locus_tag","Species")]), by = "gene")

# save to Rdata
save(AIG_seq_rm_dup_clstr6_NUC_95, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_seq_rm_dup_clstr6_NUC_95.Rdata")
save(BIR_seq_rm_dup_clstr6_NUC_95, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_seq_rm_dup_clstr6_NUC_95.Rdata")

# Check if any proteins collapsed came from different genes 
# look at duplicated clusters and keep all rows 
AIG_seq_rm_dup_clstr6_dup_NUC_95 <- AIG_seq_rm_dup_clstr6_NUC_95 %>% group_by(cluster) %>% filter(n() > 1)
BIR_seq_rm_dup_clstr6_dup_NUC_95 <- BIR_seq_rm_dup_clstr6_NUC_95 %>% group_by(cluster) %>% filter(n() > 1)

# Join coordinates using BED file used to grab the sequences
AIG_seq_rm_dup_clstr6_dup_NUC_95_coord <- left_join(AIG_seq_rm_dup_clstr6_dup_NUC_95, GIMAP_BED_name, by = "gene"  )
BIR_seq_rm_dup_clstr6_dup_NUC_95_coord <- left_join(BIR_seq_rm_dup_clstr6_dup_NUC_95, IAP_BED_name, by="gene")

# Compare Nuc 95% match file to the Protein 95% file (before subsetting for only clusters from different genes)
colnames(AIG_seq_rm_dup_clstr6_dup_NUC_95_coord)[4] <-"gene_identity_stat"
colnames(BIR_seq_rm_dup_clstr6_dup_NUC_95_coord)[4] <- "gene_identity_stat"
colnames(AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95)[4] <- "prot_identity_stat"
colnames(BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95 )[4] <- "prot_identity_stat"

AIG_NUC_95_prot_95 <- left_join(AIG_seq_rm_dup_clstr6_dup_NUC_95_coord,AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95, by=c("gene","Species"))
BIR_NUC_95_prot_95 <- left_join(BIR_seq_rm_dup_clstr6_dup_NUC_95_coord,unique(BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95), by=c("gene","Species"))
BIR_NUC_95_prot_95[unique(BIR_NUC_95_prot_95$gene),]

# why is there not total overlap between protein and nucleotide file..because the clustering is different. I guess because I used the global identity to identity similarity 
AIG_seq_rm_dup_clstr6_95[grepl("LOC111105333", AIG_seq_rm_dup_clstr6_95$gene) | grepl("LOC111105339", AIG_seq_rm_dup_clstr6_95$gene),]


#### INVESTIGATE POTENTIAL GENE ARTIFACTS ####

# LOAD BED FILES AND HAPLOTIG FILES FROM JON PURITZ
## Notes: I created BED files for the locations of all Cvirginica IAP and GIMAP genes. Jon used those coordinates to pull out 
# the mean coverage values at each gene location. He also provided the file he created by running HaploMerger to identify haplotigs across the genome. I am going to compare 
# my results with his in a short report after this. 
# Exported data frames used: IAP_BED_nams, GIMAP_BED_name

# Load BED files 
Cvir_GIMAP_meanCov <- read.table(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Jon_Puritz_6_4_2020_regenefamilybedfiles/Cvir_GIMAP.meanCov.bed",
                                 sep="\t", col.names = c("seqid","start","end","gene","meanCov"))

Cvir_IAP_meanCov <- read.table(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Jon_Puritz_6_4_2020_regenefamilybedfiles/Cvir_IAP.meanCov.bed",
                               sep="\t", col.names = c("seqid","start","end","gene","meanCov"))

Cvir_haplotigs <- read.table(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Jon_Puritz_6_4_2020_regenefamilybedfiles/haplotigs.bed",
                             sep="\t", skip= 1, col.names = c("seqid","start","end","counts","dataset"))

# Join BED files coverage to the protein sequences file where more than two genes are in a cluster
AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95 
BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95 

Cvir_GIMAP_meanCov_CD_Hit_95 <- left_join(AIG_seq_rm_dup_clstr6_dup_diff_gene_95_product_95 , Cvir_GIMAP_meanCov) %>% filter(Species =="Crassostrea_virginica")
Cvir_IAP_meanCov_CD_Hit_95 <- left_join(BIR_seq_rm_dup_clstr6_dup_diff_gene_product_95 , Cvir_IAP_meanCov) %>% filter(Species =="Crassostrea_virginica")

# Join full gene length since Ximing mentioned that haplotigs are pretty long 
Cvir_GIMAP_meanCov_CD_Hit_95 <- left_join(Cvir_GIMAP_meanCov_CD_Hit_95, GIMAP_BED_name)
Cvir_IAP_meanCov_CD_Hit_95 <- left_join(Cvir_IAP_meanCov_CD_Hit_95, IAP_BED_name)

GIMAP_gene_length_aa <- AIG_seq_rm_dup_clstr6_NUC_95[,c("aa","gene")]
IAP_gene_length_aa   <- BIR_seq_rm_dup_clstr6_NUC_95[,c("aa","gene")]
colnames(GIMAP_gene_length_aa)[1] <- "gene_length"
colnames(IAP_gene_length_aa  )[1] <- "gene_length"

Cvir_GIMAP_meanCov_CD_Hit_95_length <- left_join(Cvir_GIMAP_meanCov_CD_Hit_95, GIMAP_gene_length_aa)
Cvir_IAP_meanCov_CD_Hit_95_length <- left_join(Cvir_IAP_meanCov_CD_Hit_95 , IAP_gene_length_aa)
Cvir_IAP_meanCov_CD_Hit_95_length <- unique(Cvir_IAP_meanCov_CD_Hit_95_length)

# Make unique for each gene
Cvir_GIMAP_meanCov_CD_Hit_95_length_unique <- Cvir_GIMAP_meanCov_CD_Hit_95_length %>% distinct(gene, .keep_all = TRUE)
Cvir_IAP_meanCov_CD_Hit_95_length_unique <-  Cvir_IAP_meanCov_CD_Hit_95_length %>% distinct(gene, .keep_all = TRUE)

# Is their overlap with the haplotigs file?
Cvir_GIMAP_haplomerger_haplotigs <- Cvir_GIMAP_meanCov_CD_Hit_95_length[Cvir_GIMAP_meanCov_CD_Hit_95_length$start %in% Cvir_haplotigs,]
# 0 in the overlap
Cvir_IAP_haplomerger_haplotigs <- Cvir_IAP_meanCov_CD_Hit_95_length[Cvir_IAP_meanCov_CD_Hit_95_length$start %in% Cvir_haplotigs,]
# 0 exact gene overlaps, need to check if my genes hit to any of these ranges 

# No exact overlaps with gene coordinates, check if it is within range
Cvir_GIMAP_meanCov_CD_Hit_95_length_unique$HM_found_start <- ifelse(sapply(Cvir_GIMAP_meanCov_CD_Hit_95_length_unique$start, function(p) 
  any(Cvir_haplotigs$start <= p & Cvir_haplotigs$end >= p)),"YES", NA)
Cvir_GIMAP_meanCov_CD_Hit_95_length_unique$HM_found_end <- ifelse(sapply(Cvir_GIMAP_meanCov_CD_Hit_95_length_unique$end, function(p) 
  any(Cvir_haplotigs$start <= p & Cvir_haplotigs$end >= p)),"YES", NA)

Cvir_IAP_meanCov_CD_Hit_95_length_unique$HM_found_start <- ifelse(sapply(Cvir_IAP_meanCov_CD_Hit_95_length_unique$start, function(p) 
  any(Cvir_haplotigs$start <= p & Cvir_haplotigs$end >= p)),"YES", NA)
Cvir_IAP_meanCov_CD_Hit_95_length_unique$HM_found_end <- ifelse(sapply(Cvir_IAP_meanCov_CD_Hit_95_length_unique$end, function(p) 
  any(Cvir_haplotigs$start <= p & Cvir_haplotigs$end >= p)),"YES", NA)

## Calculate mean coverage within gene clusters
Cvir_GIMAP_meanCov_CD_Hit_95_length_unique_mean <- Cvir_GIMAP_meanCov_CD_Hit_95_length_unique %>% group_by(cluster) %>% mutate(mean_Cov_clstr = mean(meanCov))
Cvir_IAP_meanCov_CD_Hit_95_length_unique_mean <- Cvir_IAP_meanCov_CD_Hit_95_length_unique %>% group_by(cluster) %>% mutate(mean_Cov_clstr = mean(meanCov)) 

#GIMAP results:
#  - Cluster 291: Mean coverage of 436. Look at the nucleotide sequences of these genes LOC111110115, LOC111106081

#IAP results:
#  - Cluster 62: mean coverage of 444. Includes LOC111100470 and LOC111101689
#- Cluster 328: mean coverage across cluster of 280. Includes LOC111132301 LOC111114013, LOC111103682, LOC111132489, LOC111132589, LOC111102106, LOC111114070
#- Cluster 344: mean coverage across cluster 484. Includes LOC111117856, LOC111116826, LOC111111659

# Make files with list of IDs to export to bluewaves to extract nucleotide sequences 
GIMAP_cluster_219_gene <- Cvir_GIMAP_meanCov_CD_Hit_95_length_unique_mean  %>% filter(cluster == "Cluster 219")
IAP_cluster_62_gene <- Cvir_IAP_meanCov_CD_Hit_95_length_unique_mean  %>% filter(cluster == "Cluster 62") 
IAP_cluster_328_gene <- Cvir_IAP_meanCov_CD_Hit_95_length_unique_mean  %>% filter(cluster == "Cluster 328")
IAP_cluster_344_gene <- Cvir_IAP_meanCov_CD_Hit_95_length_unique_mean  %>% filter(cluster == "Cluster 344") 

write.table(GIMAP_cluster_219_gene$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_cluster_219_gene.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)
write.table(IAP_cluster_62_gene$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_62_gene.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)
write.table(IAP_cluster_328_gene$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_328_gene.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)
write.table(IAP_cluster_344_gene$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_344_gene.txt",
            quote=FALSE, row.names=FALSE, col.names = FALSE)

## Export the coordinates of clusters of interest as BED files so I can show in my IGV track session
write.table(GIMAP_cluster_219_gene[,c("seqid","start","end","cluster")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_cluster_219_gene.bed",
            quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
write.table(IAP_cluster_62_gene[,c("seqid","start","end","cluster")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_62_gene.bed",
            quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
write.table(IAP_cluster_328_gene[,c("seqid","start","end","cluster")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_328_gene.bed",
            quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
write.table(IAP_cluster_344_gene[,c("seqid","start","end","cluster")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_cluster_344_gene.bed",
            quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")

## Extract nucleotide sequences of interesting clusters on bluewaves (see code in notes file)
## Align nucleotide sequences in bluewaves and view in IGV tracks 
## Put results into report compiled by Jon 

#### COLLAPSE GIMAP AND IAP GENES AND PROTEINS BASED ON HAPLOTIG RESULTS ####

## Conclusions and Collapsing the GIMAP Cluster analysis

#1. GIMAP Cluster 219: two genes need to be collapsed into 1 1. gene LOC111106081 is shorter and has lower coverage, this one should be removed 
# LOC111106081 protein XP_022296317.1 should be combined with XP_022302183.1
# Does LOC111106081 have other proteins?
AIG1_XP_ALL_gff_GIMAP_species_join %>% filter(gene=="LOC111106081") # No just this one! 

# First create for reference a data frame with all species and the haplotig collapsed
AIG1_XP_ALL_gff_GIMAP_species_join_haplotig_collapsed <- AIG1_XP_ALL_gff_GIMAP_species_join %>% filter(gene !="LOC111106081")

# Create collapsed data frame from the original MAFFT list where all the exact identical protein sequences from genes were removed (when CD-Hit was run only removing 100% identical sequences)
AIG1_dup_seq_rm_kept_haplotig_collapsed <- AIG1_dup_seq_rm_kept %>% filter(gene !="LOC111106081")

# Subset for only Mizuhopecten yessoensis (outgroup) and Crassostrea virginica and Crassostrea gigas sequences will be used to generate smaller tree for just these species 
AIG1_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG <- AIG1_dup_seq_rm_kept_haplotig_collapsed %>% filter(Species == "Mizuhopecten_yessoensis" | Species == "Crassostrea_virginica" |
                                                                                                         Species == "Crassostrea_gigas"  )

# Export MY, CV, CG Sequence list for MAFFT and RAxML 
write.table(unique(AIG1_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG$protein_id), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

## Conclusions and Collapsing the IAP Cluster analysis
# Remember that for the IAP analysis, the original tree had collapsed different genes with exact protein sequences. Which genes were collapsed? 

#1. IAP Cluster 62 should NOT be collapsed

#2. IAP Cluster 328 The two genes on top of the alignment (LOC111132489 and LOC111114013) form a cluster together and are most similar and should be collapsed into one gene. 
# gene LOC111132489 has much higher coverage and is likely the "real" gene. Going to remove gene LOC111114013 and it's protein XP_022308010.1 should be collapsed with 
# gene LOC111132489 protein XP_022336007.1 
# Does LOC111114013 have multiple proteins - NO
BIR_XP_gff_species_join %>% filter(gene=="LOC111114013")  # Has only 1 protein, XP_022308010.1

# IAP Cluster 328 other five genes should be collapsed together: LOC111132301 has the highest relative coverage, meaning LOC111103682, LOC111132589,LOC111102106,LOC111114070 should be collapsed into this one gene
# Do any of these gene have multiple proteins? 
BIR_XP_gff_species_join %>% filter(gene=="LOC111103682") # only 1 protein XP_022292821.1
BIR_XP_gff_species_join %>% filter(gene=="LOC111132589") # only 1 protein XP_022336127.1
BIR_XP_gff_species_join %>% filter(gene=="LOC111102106") # only 1 protein XP_022290466.1
BIR_XP_gff_species_join %>% filter(gene=="LOC111114070") # only 1 protein XP_022308067.1

#3. IAP Cluster 344: The two sequences with the greatest similarity in gene sequence, LOC111116826 and LOC111111659, should be collapsed. LOC111111659 has extremely low coverage and its protein XP_022304464.1
# should be collapsed into LOC111116826
# Does LOC111111659 have more proteins 
BIR_XP_gff_species_join %>% filter(gene=="LOC111111659")  # Has only 1 protein, XP_022304464.1

# First create for reference a data frame with all species and the 6 haplotigs collapsed
IAP_gene_remove <- c("LOC111111659", 
                     "LOC111114013" ,
                     "LOC111103682" ,
                     "LOC111132589" ,
                     "LOC111102106" ,
                     "LOC111114070")
BIR_XP_gff_species_join_haplotig_collapsed <- BIR_XP_gff_species_join[(!BIR_XP_gff_species_join$gene %in%  IAP_gene_remove),]

#### Get Gene and protein list after collapsing IAP genes ####
All_mollusc_IAP_gene_list_after_haplotig_collapsed <- BIR_XP_gff_species_join_haplotig_collapsed %>%
  ungroup() %>%
  mutate(gene_locus_tag = case_when(
    is.na(.$gene) ~ locus_tag, 
    TRUE ~ gene)) %>% 
  distinct(gene_locus_tag, Species) %>% 
  count(Species) %>%
  arrange(desc(n))
# Use this table for paper

All_mollusc_IAP_gene_list_after_haplotig_collapsed %>% gt::gt(groupname_col = "Species") %>%
  col_label(n = md("**Number of Genes**"))

# How many total proteins and total uncharacterized for C. virginica
BIR_XP_gff_species_join_haplotig_collapsed %>%
  ungroup() %>%
  mutate(gene_locus_tag = case_when(
    is.na(.$gene) ~ locus_tag, 
    TRUE ~ gene)) %>% 
  filter(Species == "Crassostrea_virginica") %>% 
  #grepl("uncharacterized", product)) %>%
count() %>% 
  arrange(desc(n)) # 39 uncharacterized

# Create collapsed data frame from the original MAFFT list where all the exact identical protein sequences from genes were removed (when CD-Hit was run only removing 100% identical sequences)
BIR_dup_seq_rm_kept_haplotig_collapsed <- BIR_dup_seq_rm_kept[(! BIR_dup_seq_rm_kept$gene %in%  IAP_gene_remove),]

# Subset for only Mizuhopecten yessoensis (outgroup) and Crassostrea virginica and Crassostrea gigas sequences will be used to generate smaller tree for just these species 
BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG <- BIR_dup_seq_rm_kept_haplotig_collapsed%>% filter(Species == "Mizuhopecten_yessoensis" | Species == "Crassostrea_virginica" |
                                                                                                      Species == "Crassostrea_gigas"  )

# Export MY, CV, CG Sequence list for MAFFT and RAxML 
write.table(unique(BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG$protein_id), file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#### EXPORT HAPLOTIG COLLAPSED IAP GENE SEQUENCES FOR GENE FUNCTIONAL DIVERSITY TREE ####

### EXPORT CDS SEQUENCES - didnt end up needing the CDS because interproscan does this internally 
# Use the full list of haplotig collapsed genes above to get exon coordinates to run MAFFT and RAxML in order to assess the potential full domain diversity not just what is expressed
    # in the transcript variants
# check for correct gene number
BIR_XP_gff_species_join_haplotig_collapsed %>% filter(Species == "Crassostrea_gigas" | Species == "Crassostrea_virginica") %>% ungroup() %>% distinct(gene, Species) %>% count(Species)

BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY <- BIR_XP_gff_species_join_haplotig_collapsed %>% filter(Species == "Crassostrea_gigas" | Species == "Crassostrea_virginica" |
                                                             Species == "Mizuhopecten_yessoensis") %>% ungroup() %>%  distinct(gene, Species)
  
# Find the CDS exon coordinates for these genes in the C. gig, C. vir and M. yes gff annotations
C_vir_rtracklayer %>% filter(type =="CDS") %>% head()
# Load the MY annotation
MY_gff <-  readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/GCF_002113885.1_ASM211388v2_genomic.gff")
MY_gff <- as.data.frame(MY_gff)

# get the CDS lines for each C. virginica gene 
BIR_XP_gff_species_join_haplotig_collapsed_CV_CDS <- C_vir_rtracklayer[C_vir_rtracklayer$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "CDS") %>%
  # sort in order of start position
  arrange(gene, start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)

# get the CDS lines for each C. gigas gene 
BIR_XP_gff_species_join_haplotig_collapsed_CG_CDS <- C_gig_rtracklayer[C_gig_rtracklayer$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "CDS") %>%
  # sort in order of start position
  arrange(gene,start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)

# get the CDS lines for each M. yessoensis gene 
BIR_XP_gff_species_join_haplotig_collapsed_MY_CDS <- MY_gff[MY_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "CDS") %>%
  # sort in order of start position
  arrange(gene,start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)

# Export CDS position lists and load into bluewaves for extraction of sequences from genomes, concatenation, alignment with MAFFT, and phylogenetic inference with RAxML
write.table(BIR_XP_gff_species_join_haplotig_collapsed_CV_CDS, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_CV_CDS.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_CG_CDS, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_CG_CDS.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_MY_CDS, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_MY_CDS.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")

### EXPORT FULL GENE SEQUENCE

# get the full gene line for each C. virginica gene 
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene <- C_vir_rtracklayer[C_vir_rtracklayer$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene) # 69

# get the full gene line for each C. gigas gene 
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene <- C_gig_rtracklayer[C_gig_rtracklayer$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene) # 40

# get the full gene line for each M. yessoensis gene 
BIR_XP_gff_species_join_haplotig_collapsed_MY_gene <- MY_gff[MY_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,] %>% filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_MY_gene) # 68

# Export gene position lists and load into bluewaves for extraction of sequences from genomes, concatenation, alignment with MAFFT, and phylogenetic inference with RAxML
write.table(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_CV_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_CG_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_MY_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_MY_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")

# unlist MY annotation to save space in workspace 
rm(MY_gff)

# make list with full gene and species
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_species <- BIR_XP_gff_species_join_haplotig_collapsed_CV_gene %>% mutate(Species = "Crassostrea_virginica") %>% dplyr::select(gene, Species)
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_species <- BIR_XP_gff_species_join_haplotig_collapsed_CG_gene %>% mutate(Species = "Crassostrea_gigas") %>% dplyr::select(gene, Species)
BIR_XP_gff_species_join_haplotig_collapsed_MY_gene_species <- BIR_XP_gff_species_join_haplotig_collapsed_MY_gene %>% mutate(Species = "Mizuhopecten_yessoensis") %>% dplyr::select(gene, Species)

CV_CG_MY_gene_species <- rbind(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_species,
                               BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_species,
                               BIR_XP_gff_species_join_haplotig_collapsed_MY_gene_species)

# total genes
nrow(CV_CG_MY_gene_species) # 177 because some were collapsed

### REPEAT TO EXPORT IAP GENE COORDINATES FROM OTHER 7 MOLLUSCS SO I CAN I MAKE A GENE TREE WITH ALL THE SPECIES ###

# which species do I need to add to tree
BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% distinct(Species)
    #2 Biomphalaria_glabrata  
    #4 Octopus_bimaculoides   
    #5 Pomacea_canaliculata   
    #6 Lottia_gigantea        
    #8 Aplysia_californica    
    #9 Octopus_vulgaris       
    #10 Elysia_chlorotica   

BIR_XP_gff_species_join_haplotig_collapsed_Octopus_bimaculoides <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Octopus_bimaculoides") %>% distinct(gene)  %>% mutate(Species = "Octopus_bimaculoides")
BIR_XP_gff_species_join_haplotig_collapsed_Pomacea_canaliculata <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Pomacea_canaliculata") %>% distinct(gene) %>% mutate(Species = "Pomacea_canaliculata")
BIR_XP_gff_species_join_haplotig_collapsed_Lottia_gigantea      <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Lottia_gigantea") %>% distinct(locus_tag) %>% rename(gene = locus_tag) %>% mutate(Species = "Lottia_gigantea")
BIR_XP_gff_species_join_haplotig_collapsed_Biomphalaria_glabrata  <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Biomphalaria_glabrata") %>% distinct(gene) %>% mutate(Species = "Biomphalaria_glabrata")
BIR_XP_gff_species_join_haplotig_collapsed_Aplysia_californica  <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Aplysia_californica") %>% distinct(gene) %>% mutate(Species = "Aplysia_californica")
BIR_XP_gff_species_join_haplotig_collapsed_Octopus_vulgaris     <- BIR_XP_gff_species_join_haplotig_collapsed %>% ungroup() %>% filter(Species== "Octopus_vulgaris") %>% distinct(gene) %>% mutate(Species = "Octopus_vulgaris")
BIR_XP_gff_species_join_haplotig_collapsed_Elysia_chlorotica   <- BIR_XP_gff_species_join_haplotig_collapsed %>%  ungroup() %>% filter(Species==  "Elysia_chlorotica") %>% distinct(locus_tag)  %>% rename(gene = locus_tag) %>% mutate(Species = "Elysia_chlorotica")

# get the full gene line for each Octopus_bimaculoides gene 
BIR_XP_gff_species_join_haplotig_collapsed_OB_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_Octopus_bimaculoides$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_OB_gene) # 11

# get the full gene line for each Pomacea_canaliculata gene 
BIR_XP_gff_species_join_haplotig_collapsed_PC_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_Pomacea_canaliculata$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_PC_gene) # 27 

# get the full gene line for each Lottia_gigantea gene 
BIR_XP_gff_species_join_haplotig_collapsed_LG_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$locus_tag %in% BIR_XP_gff_species_join_haplotig_collapsed_Lottia_gigantea$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, locus_tag) 
nrow(BIR_XP_gff_species_join_haplotig_collapsed_LG_gene) # 23

# get the full gene line for each Biomphalaria_glabrata gene 
BIR_XP_gff_species_join_haplotig_collapsed_BG_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_Biomphalaria_glabrata$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_BG_gene) # 88

# get the full gene line for each Aplysia_californica gene 
BIR_XP_gff_species_join_haplotig_collapsed_AC_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_Aplysia_californica$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_AC_gene) # 23

# get the full gene line for each Octopus_vulgaris gene 
BIR_XP_gff_species_join_haplotig_collapsed_OV_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_Octopus_vulgaris $gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, gene)
nrow(BIR_XP_gff_species_join_haplotig_collapsed_OV_gene) # 10

# get the full gene line for each Elysia_chlorotica gene 
BIR_XP_gff_species_join_haplotig_collapsed_EC_gene <- All_mollusc_gene_gff[All_mollusc_gene_gff$locus_tag %in% BIR_XP_gff_species_join_haplotig_collapsed_Elysia_chlorotica$gene,] %>% 
  filter(type == "gene") %>%
  # sort in order of start position
  arrange(start) %>%
  # extract in bed format for getfasta: seqid, start, stop, gene 
  distinct(seqid, start, end, locus_tag) 
nrow(BIR_XP_gff_species_join_haplotig_collapsed_EC_gene) # 15

# Export gene position lists and load into bluewaves for extraction of sequences from genomes, concatenation, alignment with MAFFT, and phylogenetic inference with RAxML
write.table(BIR_XP_gff_species_join_haplotig_collapsed_OB_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_OB_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_PC_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_PC_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_LG_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_LG_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_BG_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_BG_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_AC_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_AC_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_OV_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_OV_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")
write.table(BIR_XP_gff_species_join_haplotig_collapsed_EC_gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_EC_Gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")

# make list with full gene and species

All_gene_species <- rbind(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_species,
                          BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_species,
                          BIR_XP_gff_species_join_haplotig_collapsed_MY_gene_species,
                          BIR_XP_gff_species_join_haplotig_collapsed_Octopus_bimaculoides [,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Pomacea_canaliculata [,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Lottia_gigantea[,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Biomphalaria_glabrata[,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Aplysia_californica[,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Octopus_vulgaris[,c("gene","Species")],
                          BIR_XP_gff_species_join_haplotig_collapsed_Elysia_chlorotica[,c("gene","Species")])

# total genes
nrow(All_gene_species) # 374 because some were collapsed

#### PLOT ALL MOLLLUSC GENE TREE ####
### Load IAP Gene RAxML tree data 
IAP_GENE_all_species_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_XP_gff_species_join_haplotig_collapsed_all_species_Gene_MSA_RAxML")
IAP_GENE_all_species_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
IAP_GENE_all_species_raxml_tibble <- as_tibble(IAP_GENE_all_species_raxml)
colnames(IAP_GENE_all_species_raxml_tibble)[4] <- "gene"
# add species data 
IAP_GENE_all_species_raxml_tibble_join <- left_join(IAP_GENE_all_species_raxml_tibble, All_gene_species)
colnames(IAP_GENE_all_species_raxml_tibble_join)[4] <- "label"
# Convert to treedata object to store tree plus outside data
IAP_GENE_all_species_raxml_treedata <- as.treedata(IAP_GENE_all_species_raxml_tibble_join)

# check total gene number in tree
IAP_GENE_all_species_raxml_tibble_join %>% filter(!is.na(Species)) %>% nrow() # 373

## Plot IAP gene sequence tree as circular first 
IAP_GENE_all_species_raxml_treedata_circular_gene <- ggtree(IAP_GENE_all_species_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=label,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32",  "#a68340",
                                                 "#a3c763", "#c257b0", "#c083d0","#59a1cf","#c2134a","#ead76b"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis", 
                                                                                                                                             "Elysia_chlorotica","Lottia_gigantea", "Octopus_bimaculoides", "Octopus_vulgaris", "Pomacea_canaliculata", "Biomphalaria_glabrata","Aplysia_californica"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis", 
                                 "Elysia chlorotica","Lottia gigantea", "Octopus bimaculoides", "Octopus vulgaris", "Pomacea canaliculata", "Biomphalaria glabrata", "Aplysia californica")) 

# Export plot to file to put together with gene tree for paper
ggsave(filename = "IAP_all_species_GENE_circular_tree.tiff", plot= IAP_GENE_all_species_raxml_treedata_circular_gene, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/",
       width = 10 ,
       height = 10,
       units = "in",
       dpi=300)


#### PLOT MY, CV, CG GENE TREE WITH DOMAIN INFO ####

### Load IAP Gene RAxML tree data 
IAP_GENE_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene_MSA_RAxML")
IAP_GENE_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
IAP_GENE_raxml_tibble <- as_tibble(IAP_GENE_raxml)
colnames(IAP_GENE_raxml_tibble)[4] <- "gene"
# add species data 
IAP_GENE_raxml_tibble_join <- left_join(IAP_GENE_raxml_tibble, CV_CG_MY_gene_species)
colnames(IAP_GENE_raxml_tibble_join)[4] <- "label"
# Convert to treedata object to store tree plus outside data
IAP_GENE_raxml_treedata <- as.treedata(IAP_GENE_raxml_tibble_join)

# check total gene number in tree
IAP_GENE_raxml_tibble_join %>% filter(!is.na(Species)) %>% nrow() # 177

## Plot IAP gene sequence tree as circular first 
IAP_GENE_raxml_treedata_circular_gene <- ggtree(IAP_GENE_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=label,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis")) +
  guides(col = guide_legend(ncol =1, title.position = "top", override.aes = aes(label = "")) ) # need to override aes to get rid of "a"

# Export plot to file 
ggsave(filename = "IAP_CV_CG_MY_GENE_circular_tree.tiff", plot=IAP_GENE_raxml_treedata_circular_gene, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/",
       width = 10 ,
       height = 10,
       units = "in",
       dpi=300)

## Plot vertical tree for making IAP protein domain tree
IAP_GENE_raxml_treedata_vertical <- 
  ggtree(IAP_GENE_raxml_treedata, aes(color=Species, fill=Species),  branch.length = "none") + 
  geom_tiplab(aes(label=label), fontface="bold", size =3.5, offset=0) + # geom_tiplab2 flips the labels correctly
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=2.0) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=2.0) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=2.0) +
   ## Add clade labels for the 21 domain groups domain groups using the internal node number
 #geom_cladelabel(261, label="1",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') + # get node order from below 
 #geom_cladelabel(254, label="2",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(233, label="3",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(228, label="4",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(217, label="5",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(211, label="6",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(207, label="7",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(198, label="8",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(201, label="9",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(311, label="10", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(308, label="11", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(291, label="12", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(304, label="13", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(329, label="14", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(340, label="15", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(343, label="16", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(347, label="17", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(353, label="18", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(359, label="19", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(188, label="20", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
 #geom_cladelabel(366, label="21", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black', extend = 0.5) +
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=14, family="sans"),
        legend.title = element_text(size=16, family="sans")) +
  #geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 2.0, fontface="bold") + # allows for subset
  xlim(-70,31.8) + #change scaling so branch lengths are smaller and all alias labels are showing
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis")) +
  guides(col = guide_legend(ncol =1, title.position = "top", override.aes = aes(label = "")) ) # need to override aes to get rid of "a"

### Load Interproscan data for genes
IAP_gene_gff <- readGFF(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene.fa.gff3")
IAP_gene_gff <- as.data.frame(IAP_gene_gff)
length(unique(IAP_gene_gff$seqid)) # 892

# total number of genes in IAP gff
IAP_gene_gff %>% distinct(gene) # 176... 1 was taken out.not sure why  

# Join species information, first split seqid into gene and tag
IAP_gene_gff <- IAP_gene_gff %>% separate(seqid, into = c("gene","Interpro_tag"), sep = "_")
IAP_gene_gff_join <- left_join(IAP_gene_gff, CV_CG_MY_gene_species) 

# how many Interproscan domains were found?
unique(as.character(IAP_gene_gff$Dbxref)) # 55 unique Interproscan domain types found
IAP_gene_gff %>% distinct(Name) %>% filter(!grepl("LOC",Name)) %>% View()
unique(as.character(IAP_gene_gff$signature_desc)) # 57 unique signature descriptions 

## Were any domain descriptions found here that were not found for the IAP protein Interproscan results?
  # remember setdiff shows elements in X not in Y
# descriptions in transcript domain annotation not in gene annotation 
setdiff(BIR_XP_gff_Interpro_Domains_all$signature_desc, IAP_gene_gff$signature_desc)
  # [[1]] "Ubiquitin-associated domain (UBA) profile.", "Ubiquitin-conjugating enzymes family profile.", "UBCc", "Ubiquitin-conjugating enzyme","Trp-Asp (WD) repeats signature.", "UBA_IAPs"

# descriptions in transcribed genes not in the transcript level annotations 
domain_gene_not_transcript <- as.character(setdiff(IAP_gene_gff$signature_desc, BIR_XP_gff_Interpro_Domains_all$signature_desc)) # 45 domains

## compare Dbxref ids
# Dbxrefs in gene and not in transcript level
setdiff(as.character(IAP_gene_gff$Dbxref), as.character(BIR_XP_gff_Interpro_Domains_all$Dbxref))
setdiff(as.character(BIR_XP_gff_Interpro_Domains_all$Dbxref), as.character(IAP_gene_gff$Dbxref))


# It's clear that a variety of Interproscan domains are at the gene level and not at the protein level
# need to take a closer look at the domains for comparison of which types are wholly absent or just have a slightly different subtype betwene the two methods 
BIR_XP_gff_Interpro_Domains_all %>% filter(grepl("Ig", as.character(signature_desc), ignore.case = TRUE))
BIR_XP_gff_Interpro_Domains_all %>% filter(grepl("Im", as.character(signature_desc),ignore.case = TRUE))
BIR_XP_gff_Interpro_Domains_all %>% filter(grepl("AIG", as.character(signature_desc),ignore.case = TRUE))
BIR_XP_gff_Interpro_Domains_all %>% filter(grepl("transposase", as.character(signature_desc),ignore.case = TRUE))

IAP_gene_gff %>% filter(grepl("Immunoglobulin I-set domain", as.character(signature_desc), ignore.case = TRUE)) # only in 1 gene LOC111100407
IAP_gene_gff %>% filter(gene == "LOC111100407") %>%  View() 
IAP_gene_gff %>% filter(grepl("AIG1", as.character(signature_desc), ignore.case = TRUE)) # one gene has a GIMAP AIG1 domain LOC111100017
IAP_gene_gff %>% filter(gene == "LOC111100017") %>%  View() 
IAP_gene_gff %>% filter(grepl("reverse", as.character(signature_desc), ignore.case = TRUE)) %>% distinct(gene) # 12 genes 

## Plot CDD and Interproscan domains 
# Use combination of geom_segment and geom_rect and combine plot with vertical tree using ggarrange from ggpubr
IAP_GENE_raxml_tibble_join_na <- IAP_GENE_raxml_tibble_join %>% filter(!is.na(label)) # remove rows with just bootstrap information
colnames(IAP_GENE_raxml_tibble_join_na)[4] <- "gene"
IAP_GENE_Interpro_Domains <-  left_join(IAP_GENE_raxml_tibble_join_na[,c("gene","node","Species")], IAP_gene_gff)
IAP_GENE_Interpro_Domains_only <- IAP_GENE_Interpro_Domains %>% 
  filter(source =="CDD" | grepl("InterPro:IPR", Dbxref) | grepl("SSF57924",Name) | 
           grepl("G3DSA:1.10.1170.10", Name) | grepl("G3DSA:1.10.533.10", Name) |
           grepl("SSF57850",Name) | grepl("G3DSA:4.10.60.10", Name) | grepl("PF13920", Name)) 
          # keep Interproscan domain lines, CDD NCBI lines, and add in SUPERFAMILY IAP entry, Gene3D Death domain, RING/Ubox, Zinc finger
         
IAP_GENE_Interpro_Domains_all <- IAP_GENE_Interpro_Domains

# Which domains were removed
IAP_GENE_Interpro_Domains_Name <- IAP_GENE_Interpro_Domains %>% distinct(Name) %>% filter(!grepl("LOC", Name))
IAP_GENE_Interpro_Domains_only_Name <- IAP_GENE_Interpro_Domains_only %>% distinct(Name) 

IAP_GENE_Interpro_Domains_Name[!(IAP_GENE_Interpro_Domains_Name$Name %in% IAP_GENE_Interpro_Domains_only_Name$Name),] %>% View() # 22 were removed 


## Compare with list from the protein Interproscan
BIR_XP_gff_Interpro_Domains_all_IP_cdd <- BIR_XP_gff_Interpro_Domains_all %>% ungroup() %>% 
  # make Dbxref a character
  mutate(Dbxref = as.character(Dbxref)) %>% filter(grepl("cd", Dbxref) | grepl("InterPro",Dbxref)) %>% distinct(Dbxref) %>% mutate(Type = "protein")
BIR_XP_gff_Interpro_Domains_all_IP_cdd$Dbxref <- gsub('[\"]', '', BIR_XP_gff_Interpro_Domains_all_IP_cdd$Dbxref)

IAP_GENE_Interpro_Domains_all_IP_cdd <- IAP_GENE_Interpro_Domains_all %>% ungroup() %>% 
  # make Dbxref a character
  mutate(Dbxref = as.character(Dbxref)) %>%
  filter(grepl("cd", Dbxref) | grepl("InterPro",Dbxref)) %>% distinct(Dbxref) %>% mutate(Type = "gene")
IAP_GENE_Interpro_Domains_all_IP_cdd$Dbxref <- gsub('[\"]', '', IAP_GENE_Interpro_Domains_all_IP_cdd$Dbxref)

# Compare lists
gene_prot_shared <- merge(BIR_XP_gff_Interpro_Domains_all_IP_cdd[1],IAP_GENE_Interpro_Domains_all_IP_cdd[1]) # use merge to find what they both share
prot_unique <- BIR_XP_gff_Interpro_Domains_all_IP_cdd[!(BIR_XP_gff_Interpro_Domains_all_IP_cdd$Dbxref) %in% IAP_GENE_Interpro_Domains_all_IP_cdd$Dbxref,]
gene_unique <- IAP_GENE_Interpro_Domains_all_IP_cdd[!(IAP_GENE_Interpro_Domains_all_IP_cdd$Dbxref) %in% BIR_XP_gff_Interpro_Domains_all_IP_cdd$Dbxref,]

## Prioritize Gene Interproscan domains to plot: which domain terms are most common across genes?
IAP_GENE_Interpro_Domains_only  %>% ungroup() %>% 
  # make Dbxref a character
  mutate(Dbxref = as.character(Dbxref)) %>%
  # filter for cd and Interpro
  mutate(Dbxref = case_when(Dbxref == "character(0)"~ Name,
                            TRUE ~ Dbxref)) %>%
  # get distinct gene and term
  distinct(gene,Dbxref) %>% 
  # count frequency of each term across genes
  count(Dbxref) %>% arrange(desc(n)) %>% View()

    # I ORGANIZED THIS DATA AND ADDED DESCRIPTIONS FOR EACH TERM IN THE SPREADSHEET "Interproscan_term_Gene_Protein_level.xlsx" 
    # - came up with a list of 44 terms based on high frequency in genes, the ones most shared with the protein level, and the ones containing interesting transposon related domains 

## Get full length of gene from IAP_GENE_Interpro_Domains_all
IAP_GENE_Interpro_Domains_all_full_prot <- IAP_GENE_Interpro_Domains_all %>% filter(is.na(Interpro_tag)) %>% distinct()
nrow(IAP_GENE_Interpro_Domains_all_full_prot) # 177

# Fill in the CDD rows that have NULL for DBxref with the Name column
IAP_GENE_Interpro_Domains_only$Dbxref[IAP_GENE_Interpro_Domains_only$Dbxref == "character(0)" ] <- "CDD"
IAP_GENE_Interpro_Domains_all $Dbxref[IAP_GENE_Interpro_Domains_all$Dbxref == "character(0)" ] <- "CDD"

# unlist
IAP_GENE_Interpro_Domains_only <- IAP_GENE_Interpro_Domains_only %>% unnest(Dbxref)
IAP_GENE_Interpro_Domains_all <-  IAP_GENE_Interpro_Domains_all  %>% unnest(Dbxref)

# Change CDD rows to be the Name column
IAP_GENE_Interpro_Domains_only <- IAP_GENE_Interpro_Domains_only %>% mutate(Dbxref = ifelse(Dbxref == "CDD", Name, Dbxref))
IAP_GENE_Interpro_Domains_all <-  IAP_GENE_Interpro_Domains_all  %>% mutate(Dbxref = ifelse(Dbxref == "CDD", Name, Dbxref))

# how many unique non cd and Interpro terms now in the full domain list?
IAP_GENE_Interpro_Domains_all_SF_PF_G3 <- IAP_GENE_Interpro_Domains_all %>% ungroup() %>% distinct(Dbxref) %>% filter(!grepl("LOC", Dbxref)) %>% filter(!grepl("InterPro", Dbxref)) %>%  filter(!grepl("cd", Dbxref)) 


# Get the node order from IAP gene tree
IAP_GENE_raxml_treedata_tip  <- fortify(IAP_GENE_raxml_treedata) # not changing code from here down
IAP_GENE_raxml_treedata_tip <- subset(IAP_GENE_raxml_treedata_tip, isTip)
IAP_GENE_raxml_treedata_tip_order <- IAP_GENE_raxml_treedata_tip$label[order(IAP_GENE_raxml_treedata_tip$y, decreasing=TRUE)]
IAP_GENE_raxml_treedata_tip_order <- as.data.frame(IAP_GENE_raxml_treedata_tip_order)
colnames(IAP_GENE_raxml_treedata_tip_order)[1] <- "gene"

# Reorder the protein and polygon
IAP_GENE_Interpro_Domains_only_reorder <- full_join(IAP_GENE_raxml_treedata_tip_order, IAP_GENE_Interpro_Domains_only)
IAP_GENE_Interpro_Domains_all_reorder <- full_join(IAP_GENE_raxml_treedata_tip_order, IAP_GENE_Interpro_Domains_all)
IAP_GENE_Interpro_Domains_all_full_prot <- full_join(IAP_GENE_raxml_treedata_tip_order,IAP_GENE_Interpro_Domains_all_full_prot)

# Add polygon height
IAP_GENE_Interpro_Domains_only_reorder_ID  <- IAP_GENE_Interpro_Domains_only_reorder  %>% distinct(gene) 
IAP_GENE_Interpro_Domains_only_reorder_ID <- IAP_GENE_Interpro_Domains_only_reorder_ID %>% 
  mutate(height_start = rev(as.numeric(row.names(IAP_GENE_Interpro_Domains_only_reorder_ID )) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(IAP_GENE_Interpro_Domains_only_reorder_ID)) + .5))

IAP_GENE_Interpro_Domains_all_reorder_ID  <- IAP_GENE_Interpro_Domains_all_reorder  %>% distinct(gene) 
IAP_GENE_Interpro_Domains_all_reorder_ID <- IAP_GENE_Interpro_Domains_all_reorder_ID %>% 
  mutate(height_start = rev(as.numeric(row.names(IAP_GENE_Interpro_Domains_all_reorder_ID )) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(IAP_GENE_Interpro_Domains_all_reorder_ID)) + .5))

# Join back in height
IAP_GENE_Interpro_Domains_only_reorder <- left_join(IAP_GENE_Interpro_Domains_only_reorder , IAP_GENE_Interpro_Domains_only_reorder_ID )
IAP_GENE_Interpro_Domains_all_reorder <-  left_join(IAP_GENE_Interpro_Domains_all_reorder ,  IAP_GENE_Interpro_Domains_all_reorder_ID )

# How many unique domains in MY, CV, CG
IAP_GENE_Interpro_Domains_only_reorder %>% unnest(Dbxref) %>% distinct(Dbxref)

## Some genes are extremely long, shorten those that are over 1500 and add in railroad tracks 
# shorten any genes that are longer than 2500 base pairs
IAP_GENE_Interpro_Domains_all_full_prot_shortened <- IAP_GENE_Interpro_Domains_all_full_prot %>%
  mutate(end = case_when(
    end >= 1500 ~ as.numeric(1500),
    TRUE ~ as.numeric(end)
  ))

# which genes are longer than 1500?
IAP_GENE_Interpro_Domains_all_full_prot_shortened_list <- IAP_GENE_Interpro_Domains_all_full_prot_shortened %>% filter(end == 1500) %>% distinct(gene)

# what is the height start and end for each gene
IAP_GENE_Interpro_Domains_all_full_prot_shortened_list_coord <- left_join(IAP_GENE_Interpro_Domains_all_full_prot_shortened_list, IAP_GENE_Interpro_Domains_only_reorder) %>% distinct(gene, height_start, height_end)
          
# Add in set of railroad tracks for genes over 1500
IAP_GENE_Interpro_Domains_only_reorder_railroad <-  IAP_GENE_Interpro_Domains_all_full_prot_shortened_list_coord %>% mutate(Dbxref = "gene_shortened", start = 1490, end = 1495)

# Add in new rows 
IAP_GENE_Interpro_Domains_only_reorder <- plyr::rbind.fill(IAP_GENE_Interpro_Domains_only_reorder, IAP_GENE_Interpro_Domains_only_reorder_railroad)
class(IAP_GENE_Interpro_Domains_only_reorder$Dbxref)

# Subset for particular domains of interest
IAP_Gene_Interproscan_cdd_list <- c("InterPro:IPR006703",
                                    "InterPro:IPR022103",
                                    "InterPro:IPR001370","SSF57924","G3DSA:1.10.1170.10","InterPro:IPR016187","InterPro:IPR001304","InterPro:IPR016186",
                                    "G3DSA:1.10.533.10","InterPro:IPR011010","InterPro:IPR008593","cd00397","InterPro:IPR043502","InterPro:IPR036691",
                                    "InterPro:IPR000305","InterPro:IPR035901","InterPro:IPR011335","InterPro:IPR011604","InterPro:IPR009057","InterPro:IPR013098",
                                    "InterPro:IPR007110","InterPro:IPR036179","InterPro:IPR013783","InterPro:IPR041588","InterPro:IPR013762","InterPro:IPR001584",
                                    "InterPro:IPR002104","InterPro:IPR010998","InterPro:IPR043128","InterPro:IPR000477","InterPro:IPR026960","InterPro:IPR002156",
                                    "InterPro:IPR036397","InterPro:IPR012337","cd16713","SSF57850","cd09275","cd09274","cd09276","cd03714","cd01647","cd01650",
                                    "InterPro:IPR038717","InterPro:IPR027805","InterPro:IPR002492","InterPro:IPR019080","InterPro:IPR001841","InterPro:IPR013083",
                                    "G3DSA:4.10.60.10","PF13920","gene_shortened")

# subset for the Interproscan and cdd domains of interest 
IAP_GENE_Interpro_Domains_only_reorder$Dbxref <- gsub('[\"]', '', IAP_GENE_Interpro_Domains_only_reorder$Dbxref)

IAP_GENE_Interpro_Domains_only_reorder_subset <- IAP_GENE_Interpro_Domains_only_reorder %>% filter(Dbxref %in% IAP_Gene_Interproscan_cdd_list)

# Set factor level order of the nodes set levels in reverse order
IAP_GENE_Interpro_Domains_only_reorder_subset$node <- factor(IAP_GENE_Interpro_Domains_only_reorder_subset$node, levels = rev(unique(IAP_GENE_Interpro_Domains_only_reorder_subset$node)))
IAP_GENE_Interpro_Domains_only_reorder_subset$Dbxref <- factor(IAP_GENE_Interpro_Domains_only_reorder_subset$Dbxref, levels = unique(IAP_GENE_Interpro_Domains_only_reorder_subset$Dbxref))

IAP_GENE_Interpro_Domains_all_reorder$node <-   factor(IAP_GENE_Interpro_Domains_all_reorder$node, levels = unique(IAP_GENE_Interpro_Domains_all_reorder$node))
#IAP_GENE_Interpro_Domains_all_reorder$Dbxref <- factor(IAP_GENE_Interpro_Domains_all_reorder$Dbxref, levels = unique(IAP_GENE_Interpro_Domains_all_reorder$Dbxref))

IAP_GENE_Interpro_Domains_all_full_prot_shortened$node <- factor(IAP_GENE_Interpro_Domains_all_full_prot_shortened$node, levels = rev(unique(IAP_GENE_Interpro_Domains_all_full_prot_shortened$node)))

# Do all the genes have Interproscan BIR domains or IAP superfamily classification?
IAP_GENE_Interpro_Domains_only_reorder_subset_confirm <- IAP_GENE_Interpro_Domains_only_reorder_subset %>% filter(grepl("IPR001370",Dbxref) | grepl("SSF57924",Dbxref) | 
                                                            grepl("IPR022103", Dbxref) | grepl("G3DSA:1.10.1170.10", Dbxref)) %>% 
  distinct(gene) # 168 genes have a BIR repeat!
IAP_GENE_Interpro_Domains_all_reorder_gene <- unique(IAP_GENE_Interpro_Domains_all_reorder$gene)
# which genes don't have these BIR/IAP signatures? Do they have different ones? 
setdiff(IAP_GENE_Interpro_Domains_only_reorder_subset_confirm$gene, IAP_GENE_Interpro_Domains_all_reorder_gene)
IAP_GENE_not_found <- setdiff(IAP_GENE_Interpro_Domains_all_reorder_gene, IAP_GENE_Interpro_Domains_only_reorder_subset_confirm$gene) 

IAP_GENE_not_found_domain <- IAP_GENE_Interpro_Domains_all_reorder[IAP_GENE_Interpro_Domains_all_reorder$gene %in% IAP_GENE_not_found,]
IAP_GENE_not_found_domain %>% filter(!grepl("Interpro", as.character(Dbxref)) | !grepl("cd", as.character(Dbxref))) %>% View()

### Stats regarding domain usage in different categories 
# get distinct gene and Dbxref combinations
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct <- IAP_GENE_Interpro_Domains_only_reorder_subset %>% distinct(gene, Dbxref)

# Join with manually curated list of description for each term
Interproscan_term_Gene_Protein_level <- readxl::read_xlsx("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Interproscan_term_Gene_Protein_level.xlsx",
                                                          sheet = 2, col_names = c("Dbxref","Dbxref_list", "title","color","old_color"))
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title <- left_join(IAP_GENE_Interpro_Domains_only_reorder_subset_distinct, Interproscan_term_Gene_Protein_level[,c("Dbxref","title")])

## Search for particular terms
# transposase 
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("transposase", ignore.case = TRUE, title)) # 3 genes with transposase
    # 1 LOC110460644 InterPro:IPR002492                Transposase, Tc1-like
    # 2 LOC110460644 InterPro:IPR038717     Tc1-like transposase, DDE domain
    # 3 LOC110452306 InterPro:IPR002492                Transposase, Tc1-like
    # 4 LOC110465395 InterPro:IPR027805 Transposase, Helix-turn-helix domain
# Ty3, RNase_HI_RT_DIRS1, or ribonuclease
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("Ty3", ignore.case = TRUE, title) | 
                                    grepl("RNase_HI_RT_DIRS1", ignore.case = TRUE, title) | grepl("Ribonuclease", ignore.case = TRUE, title)) 
    #gene             Dbxref                                                          title
    #1 LOC111103270            cd09275                                              RNase_HI_RT_DIRS1 # C. vir
    #2 LOC110462612 InterPro:IPR002156                                          Ribonuclease H domain # mizuhopecten 
    #3 LOC110462612 InterPro:IPR036397                                     Ribonuclease H superfamily # mizuhopecten
    #4 LOC110462612 InterPro:IPR012337                                Ribonuclease H-like superfamily # mizuhopecten
    #5 LOC111112532 InterPro:IPR036397                                     Ribonuclease H superfamily # C. vir
    #6 LOC111112532 InterPro:IPR012337                                Ribonuclease H-like superfamily # C. vir
    #7 LOC111100017 InterPro:IPR036397                                     Ribonuclease H superfamily # C. vir
    #8 LOC111100017 InterPro:IPR012337                                Ribonuclease H-like superfamily # C. vir
    #9 LOC111100017            cd09274 Ty3/Gypsy family of RNase HI in long-term repeat retroelements # C. vir
  
#Immunoglobulin    
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("immunoglobulin", ignore.case = TRUE, title)) 
    #1 LOC111100407 InterPro:IPR036179 Immunoglobulin-like domain superfamily
    #2 LOC111100407 InterPro:IPR007110             Immunoglobulin-like domain
    #3 LOC111100407 InterPro:IPR013783               Immunoglobulin-like fold
    #4 LOC111100407 InterPro:IPR013098                   Immunoglobulin I-set
#C-lectin     
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("C-type lectin", ignore.case = TRUE, title)) 
    #gene             Dbxref                                      title
    #1 LOC111103155 InterPro:IPR016187                         C-type lectin fold
    #2 LOC111103155 InterPro:IPR016186 C-type lectin-like/link domain superfamily
    #3 LOC111103155 InterPro:IPR001304                         C-type lectin-like

# endonuclease
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("endonuclease", ignore.case = TRUE, title)) 
    #gene             Dbxref                                            title
    #1 LOC110457934 InterPro:IPR036691 Endonuclease/exonuclease/phosphatase superfamily
    #2 LOC110456394 InterPro:IPR036691 Endonuclease/exonuclease/phosphatase superfamily
    #3 LOC111100394 InterPro:IPR036691 Endonuclease/exonuclease/phosphatase superfamily
    #4 LOC111100394 InterPro:IPR035901                 GIY-YIG endonuclease superfamily
    #5 LOC111100394 InterPro:IPR000305                             GIY-YIG endonuclease
    #6 LOC111100394 InterPro:IPR011335            Restriction endonuclease type II-like
    #7 LOC110457936 InterPro:IPR036691 Endonuclease/exonuclease/phosphatase superfamily
# Integrase 
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("integrase", ignore.case = TRUE, title)) 
    #gene             Dbxref                                        title
    #1 LOC111103270 InterPro:IPR002104                  Integrase, catalytic domain
    #2 LOC111103270 InterPro:IPR013762 Integrase-like, catalytic domain superfamily
    #3 LOC111103270 InterPro:IPR010998            Integrase/recombinase, N-terminal
    #4 LOC111100017 InterPro:IPR001584                    Integrase, catalytic core
    #5 LOC111100017 InterPro:IPR041588                Integrase zinc-binding domain

# RT, Integrase, and RnaseH needed for retrotransposon
IAP_GENE_Interpro_Domains_only_reorder_subset_distinct_title %>% filter(grepl("integrase", ignore.case = TRUE, title) | grepl("RT", ignore.case = TRUE, title) | grepl("RNase", ignore.case = TRUE, title) |
                                                                          grepl("Ribonuclease", ignore.case = TRUE, title)) %>% 
  group_by(gene) %>% View()
  # LOC111103270 has Integrase, RNase_HI_RT_DIRS1, RT_nLTR_like, RT_DIRS1 
  # LOC111112532 has Ribonuclase H and RT_nLTR_like
  # LOC111100017 has Integrase, Ribonuclease H, RT_nLTR_like, Ty3/Gypsy family of RNase HI in long-term repeat retroelements,RT_LTR: Reverse transcriptases (RTs) from retrotransposons and retroviruses


### Create main plot ###
IAP_GENE_Interproscan_domain_plot_domain_subset <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data = IAP_GENE_Interpro_Domains_all_full_prot_shortened,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add protein domain boxes with geom_rect 
  geom_rect(data=IAP_GENE_Interpro_Domains_only_reorder_subset,inherit.aes = FALSE,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Dbxref)) +
  #add labels
  labs(y = NULL, x = "Nucleotide position") +
  # add theme
  theme_bw() + 
  # plot theme
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=14, family="sans"),
        legend.title = element_text(size=16, family="sans"),
        axis.text.x=element_text(size=18, family="sans"),
        axis.title.x.bottom = element_text(size=18, family="sans")) +
  # Change y axis ticks
  scale_x_continuous(breaks=c(0,500,1000,1500,1800), expand = c(0,0)) + 
  # Change domain labels 
  scale_fill_manual(values=c(                "#b2dbfe",
                                             "#5ba6a6",
                                             "#524ed4",
                                             "#524ed4",
                                             "#524ed4",
                                             "#b6eee2",
                                             "#b6eee2",
                                             "#b6eee2",
                                             "brown", 
                                             "#53e1a1",
                                             "#d393d4",
                                             "#b77853",
                                             "#57e5d6",
                                             "#e1c2aa",
                                             "#d04d8f",
                                             "#d04d8f",
                                             "#d04d8f",
                                             "#495f3a",
                                             "#89599b",
                                             "#c058c6",
                                             "#c058c6",
                                             "#c058c6",
                                             "#c058c6",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#5db8de",
                                             "#5db8de",
                                             "#5db8de",
                                             "#55b793",
                                             "#55b793",
                                             "#55b793",
                                             "#dadd48",
                                             "#dadd48",
                                             "#ba4c46",
                                             "#ba4c46",
                                             "#85c967",
                                             "#59388a",
                                             "#59388a",
                                             "#59388a",
                                             "#6d82d9",
                                             "#6d82d9",
                                             "#6d82d9",
                                             "#caad44",
                                             "#d95750",
                                             "#d95750",
                                             "#d95750",
                                             "#d95750","black"), 
    name="Functional Domains",
    breaks=c( "InterPro:IPR006703",
              "InterPro:IPR022103",
              "InterPro:IPR001370",
              "SSF57924",
              "G3DSA:1.10.1170.10",
              "InterPro:IPR016187",
              "InterPro:IPR001304",
              "InterPro:IPR016186",
              "G3DSA:1.10.533.10",
              "InterPro:IPR011010",
              "InterPro:IPR008593",
              "cd00397",
              "InterPro:IPR043502",
              "InterPro:IPR036691",
              "InterPro:IPR000305",
              "InterPro:IPR035901",
              "InterPro:IPR011335",
              "InterPro:IPR011604",
              "InterPro:IPR009057",
              "InterPro:IPR013098",
              "InterPro:IPR007110",
              "InterPro:IPR036179",
              "InterPro:IPR013783",
              "InterPro:IPR041588",
              "InterPro:IPR013762",
              "InterPro:IPR001584",
              "InterPro:IPR002104",
              "InterPro:IPR010998",
              "InterPro:IPR043128",
              "InterPro:IPR000477",
              "InterPro:IPR026960",
              "InterPro:IPR002156",
              "InterPro:IPR036397",
              "InterPro:IPR012337",
              "cd16713",
              "SSF57850",
              "cd09275",
              "cd09274",
              "cd09276",
              "cd03714",
              "cd01647",
              "cd01650",
              "InterPro:IPR038717",
              "InterPro:IPR027805",
              "InterPro:IPR002492",
              "InterPro:IPR019080",
              "InterPro:IPR001841",
              "InterPro:IPR013083",
              "G3DSA:4.10.60.10",
              "PF13920","gene_shortened"),
    labels=c("AIG1-type guanine nucleotide-binding (G) domain",
             "Baculoviral IAP repeat-containing protein 6",
             "BIR repeat",
             "Inhibitor of apoptosis (IAP) repeat",
             "Inhibitor Of Apoptosis Protein (2mihbC-IAP-1); Chain A",
             "C-type lectin fold",
             "C-type lectin-like",
             "C-type lectin-like/link domain superfamily",
             "Death Domain, Fas",
             "DNA breaking-rejoining enzyme, catalytic core",
             "DNA N-6-adenine-methyltransferase",
             "DNA_BRE_C ",
             "DNA/RNA polymerase superfamily",
             "Endonuclease/exonuclease/phosphatase superfamily",
             "GIY-YIG endonuclease",
             "GIY-YIG endonuclease superfamily",
             "Restriction endonuclease type II-like",
             "Exonuclease, phage-type/RecB, C-terminal",
             "Homeobox-like domain superfamily",
             "Immunoglobulin I-set",
             "Immunoglobulin-like domain",
             "Immunoglobulin-like domain superfamily",
             "Immunoglobulin-like fold",
             "Integrase zinc-binding domain",
             "Integrase-like, catalytic domain superfamily",
             "Integrase, catalytic core",
             "Integrase, catalytic domain",
             "Integrase/recombinase, N-terminal",
             "Reverse transcriptase/Diguanylate cyclase domain",
             "Reverse transcriptase domain",
             "Reverse transcriptase zinc-binding domain",
             "Ribonuclease H domain",
             "Ribonuclease H superfamily",
             "Ribonuclease H-like superfamily",
             "RING-HC_BIRC2_3_7",
             "RING/U-box",
             "RNase_HI_RT_DIRS1 ",
             "Ty3/Gypsy family of RNase HI in long-term repeat retroelements",
             "non-LTR RNase HI domain of reverse transcriptases",
             "RT_DIRS1 ",
             "RT_LTR: Reverse transcriptases (RTs) from retrotransposons and retroviruses ",
             "RT_nLTR_like ",
             "Tc1-like transposase, DDE domain",
             "Transposase, Helix-turn-helix domain",
             "Transposase, Tc1-like",
             "YqaJ viral recombinase",
             "Zinc finger, RING-type",
             "Zinc finger, RING/FYVE/PHD-type",
             "Zinc finger, CCHC-type",
             "Zinc finger, C3HC4 type (RING finger)",
             "Full Gene Length Shortened")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top")) 


## Create gene domain tree with BIR domains removed 
## March 1st, 2021, removing BIR domains from the domain annotation data frame tree in order to reduce tree clutter
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm <- IAP_GENE_Interpro_Domains_only_reorder_subset %>% 
  filter(Dbxref != "InterPro:IPR001370") %>% filter(Dbxref != "InterPro:IPR022103") %>% filter(Dbxref !="SSF57924") %>% filter(Dbxref != "G3DSA:1.10.1170.10") %>%
  filter(Dbxref !="G3DSA:1.10.533.10") %>% filter(Dbxref != "InterPro:IPR001841")  %>% filter(Dbxref != "InterPro:IPR013083") %>% filter(Dbxref != "G3DSA:4.10.60.10") %>%
        filter(Dbxref != "PF13920") %>% filter(Dbxref !="cd16713") %>% filter(Dbxref != "SSF57850")

#2 InterPro:IPR022103 "\"InterPro:IPR022103\"," Baculoviral IAP repeat-containing protein 6            "\"#5ba6a6\"," NA       
#3 InterPro:IPR001370 "\"InterPro:IPR001370\"," BIR repeat                                             "\"#524ed4\"," NA       
#4 SSF57924           "\"SSF57924\","           Inhibitor of apoptosis (IAP) repeat                    "\"#524ed4\"," NA       
#5 G3DSA:1.10.1170.10 "\"G3DSA:1.10.1170.10\"," Inhibitor Of Apoptosis Protein (2mihbC-IAP-1); Chain A "\"#524ed4\"," NA   


### Create main plot ###
IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data = IAP_GENE_Interpro_Domains_all_full_prot_shortened,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add protein domain boxes with geom_rect 
  geom_rect(data=IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm,inherit.aes = FALSE,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Dbxref)) +
  #add labels
  labs(y = NULL, x = "Nucleotide position") +
  # add theme
  theme_bw() + 
  # plot theme
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=14, family="sans"),
        legend.title = element_text(size=16, family="sans"),
        axis.text.x=element_text(size=18, family="sans"),
        axis.title.x.bottom = element_text(size=18, family="sans")) +
  # Change y axis ticks
  scale_x_continuous(breaks=c(0,500,1000,1500,1800), expand = c(0,0)) + 
  # Change domain labels 
  scale_fill_manual(values=c(                "#b2dbfe",
                                             "#b6eee2",
                                             "#b6eee2",
                                             "#b6eee2",
                                            # "brown", 
                                             "#53e1a1",
                                             "#d393d4",
                                             "#b77853",
                                             "#57e5d6",
                                             "#e1c2aa",
                                             "#d04d8f",
                                             "#d04d8f",
                                             "#d04d8f",
                                             "#495f3a",
                                             "#89599b",
                                             "#c058c6",
                                             "#c058c6",
                                             "#c058c6",
                                             "#c058c6",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#79e379",
                                             "#5db8de",
                                             "#5db8de",
                                             "#5db8de",
                                             "#55b793",
                                             "#55b793",
                                             "#55b793",
                                             #"#dadd48",
                                            # "#dadd48",
                                             "#ba4c46",
                                             "#ba4c46",
                                             "#85c967",
                                             "#59388a",
                                             "#59388a",
                                             "#59388a",
                                             "#6d82d9",
                                             "#6d82d9",
                                             "#6d82d9",
                                             "#caad44",
                                             #"#d95750",
                                             #"#d95750",
                                             #"#d95750",
                                             #"#d95750",
                                            "black"), 
                    name="Functional Domains",
                    breaks=c( "InterPro:IPR006703",
                              "InterPro:IPR016187",
                              "InterPro:IPR001304",
                              "InterPro:IPR016186",
                            #  "G3DSA:1.10.533.10",
                              "InterPro:IPR011010",
                              "InterPro:IPR008593",
                              "cd00397",
                              "InterPro:IPR043502",
                              "InterPro:IPR036691",
                              "InterPro:IPR000305",
                              "InterPro:IPR035901",
                              "InterPro:IPR011335",
                              "InterPro:IPR011604",
                              "InterPro:IPR009057",
                              "InterPro:IPR013098",
                              "InterPro:IPR007110",
                              "InterPro:IPR036179",
                              "InterPro:IPR013783",
                              "InterPro:IPR041588",
                              "InterPro:IPR013762",
                              "InterPro:IPR001584",
                              "InterPro:IPR002104",
                              "InterPro:IPR010998",
                              "InterPro:IPR043128",
                              "InterPro:IPR000477",
                              "InterPro:IPR026960",
                              "InterPro:IPR002156",
                              "InterPro:IPR036397",
                              "InterPro:IPR012337",
                              #"cd16713",
                              #"SSF57850",
                              "cd09275",
                              "cd09274",
                              "cd09276",
                              "cd03714",
                              "cd01647",
                              "cd01650",
                              "InterPro:IPR038717",
                              "InterPro:IPR027805",
                              "InterPro:IPR002492",
                              "InterPro:IPR019080",
                              #"InterPro:IPR001841",
                              #"InterPro:IPR013083",
                              #"G3DSA:4.10.60.10",
                              #"PF13920",
                            "gene_shortened"),
                    labels=c("AIG1-type guanine nucleotide-binding (G) domain",
                             "C-type lectin fold",
                             "C-type lectin-like",
                             "C-type lectin-like/link domain superfamily",
                            # "Death Domain, Fas",
                             "DNA breaking-rejoining enzyme, catalytic core",
                             "DNA N-6-adenine-methyltransferase",
                             "DNA_BRE_C ",
                             "DNA/RNA polymerase superfamily",
                             "Endonuclease/exonuclease/phosphatase superfamily",
                             "GIY-YIG endonuclease",
                             "GIY-YIG endonuclease superfamily",
                             "Restriction endonuclease type II-like",
                             "Exonuclease, phage-type/RecB, C-terminal",
                             "Homeobox-like domain superfamily",
                             "Immunoglobulin I-set",
                             "Immunoglobulin-like domain",
                             "Immunoglobulin-like domain superfamily",
                             "Immunoglobulin-like fold",
                             "Integrase zinc-binding domain",
                             "Integrase-like, catalytic domain superfamily",
                             "Integrase, catalytic core",
                             "Integrase, catalytic domain",
                             "Integrase/recombinase, N-terminal",
                             "Reverse transcriptase/Diguanylate cyclase domain",
                             "Reverse transcriptase domain",
                             "Reverse transcriptase zinc-binding domain",
                             "Ribonuclease H domain",
                             "Ribonuclease H superfamily",
                             "Ribonuclease H-like superfamily",
                             #"RING-HC_BIRC2_3_7",
                             #"RING/U-box",
                             "RNase_HI_RT_DIRS1 ",
                             "Ty3/Gypsy family of RNase HI in long-term repeat retroelements",
                             "non-LTR RNase HI domain of reverse transcriptases",
                             "RT_DIRS1 ",
                             "RT_LTR: Reverse transcriptases (RTs) from retrotransposons and retroviruses ",
                             "RT_nLTR_like ",
                             "Tc1-like transposase, DDE domain",
                             "Transposase, Helix-turn-helix domain",
                             "Transposase, Tc1-like",
                             "YqaJ viral recombinase",
                             #"Zinc finger, RING-type",
                             #"Zinc finger, RING/FYVE/PHD-type",
                             #"Zinc finger, CCHC-type",
                             #"Zinc finger, C3HC4 type (RING finger)",
                             "Full Gene Length Shortened")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top"))


## Create separate plot only with the genes of interest, separated by species
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro <- IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm %>% filter(!is.na(node)) %>%
  # order by species and gene
  arrange(desc(Species,gene)) %>%
  # remove existing height
  dplyr::select(!c(height_start,height_end)) %>%
  # 5/10/21 - remove genes LOC111100407 and LOC111103155 which only have Ig and C-type lectin domains respectively
  filter(gene != "LOC111100407") %>% filter(gene != "LOC111103155" )

  
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro_height <- IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro %>% distinct(gene) %>%
  # add in the height start and end 
mutate(height_start = rev(as.numeric(row.names(.)) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(.)) + .5))

# join back in height 
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro <- left_join(IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro, IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro_height)
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro$gene <- factor(IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro$gene, levels =unique(IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro$gene))

# make dataframe for the protein lines
IAP_GENE_Interpro_Domains_all_full_prot_shortened_only_retro <- IAP_GENE_Interpro_Domains_all_full_prot_shortened[IAP_GENE_Interpro_Domains_all_full_prot_shortened$gene %in%IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro$gene,]
IAP_GENE_Interpro_Domains_all_full_prot_shortened_only_retro <- IAP_GENE_Interpro_Domains_all_full_prot_shortened_only_retro %>% mutate(end = 1500) %>% 
  arrange(desc(Species,gene)) %>% mutate(node = rev(as.numeric(row.names(.))))

# find gene names to print
IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro_ID <- IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro %>% distinct(gene, Species) %>% mutate(species_gene = paste(Species, gene, sep = " "))
View(rev(IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro_ID$species_gene))

IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_only_retro <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data = IAP_GENE_Interpro_Domains_all_full_prot_shortened_only_retro,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add protein domain boxes with geom_rect 
  geom_rect(data=IAP_GENE_Interpro_Domains_only_reorder_subset_BIR_rm_only_retro,inherit.aes = FALSE,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Dbxref)) +
  #add labels
  labs(y = NULL, x = "Nucleotide position") +
  # add theme
  theme_bw() + 
  # plot theme
  theme(
    #axis.ticks.y = element_blank(), 
        axis.text.y = ggtext::element_markdown(size=18, family="sans"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=18, family="sans"),
        legend.title = element_text(size=18, family="sans"),
        axis.text.x=element_text(size=18, family="sans"),
        axis.title.x.bottom = element_text(size=18, family="sans")) +
  # Change y axis ticks
  scale_x_continuous(breaks=c(0,500,1000,1500)) + 
  # change y axis labels
 scale_y_discrete(name = NULL, limits = c("*C. gigas* LOC105339790",
                                               "*C. virginica* LOC111100017",
                                               #"*C. virginica* LOC111100407",
                                               "*C. virginica* LOC111100394",
                                               "*C. virginica* LOC111112532",
                                               "*C. virginica* LOC111103270",
                                               #"*C. virginica* LOC111103155",
                                               "*M. yessoensis* LOC110457936",
                                               "*M. yessoensis* LOC110458047",
                                               "*M. yessoensis* LOC110462612",
                                               "*M. yessoensis* LOC110465395",
                                               "*M. yessoensis* LOC110456394",
                                               "*M. yessoensis* LOC110456920",
                                               "*M. yessoensis* LOC110457934",
                                               "*M. yessoensis* LOC110456891",
                                               "*M. yessoensis* LOC110456458",
                                               "*M. yessoensis* LOC110451590",
                                               "*M. yessoensis* LOC110452306",
                                               "*M. yessoensis* LOC110460644"))  +
  # Change domain labels 
  scale_fill_manual(values=c("#a2843c",
                             
                             #"#5ec778",
                             #"#5ec778",
                             #"#5ec778",
                             
                             "#e18562",
                             
                             "#5ec778",
                             
                             "#d7ba7d",
                             
                             "#cfa73d",
                             
                             "#677fd8",
                             
                             "#5a388b",
                             "#5a388b",
                             "#5a388b",
                             
                             "#873015",
                             
                             "#cdd452",
                             
                             #"#d471b2",
                             #"#d471b2",
                             #"#d471b2",
                             #"#d471b2",
                             
                             "#7d9fbf",
                             "#7d9fbf",
                             "#7d9fbf",
                             "#7d9fbf",
                             "#7d9fbf",
                             
                             "#872762",
                             "#872762",
                             "#872762",
                             
                             "#e25180",
                             "#e25180",
                             "#e25180",
                             
                             "#b9475f",
                             "#b9475f",
                             "#b9475f",
                             
                             "#bf73cb",
                             
                             "#bf73cb",
                             "#bf73cb",
                             
                             "#5f8e3e",
                             "#5f8e3e",
                             "#5f8e3e",
                             
                             "#98b342",
    "black"),
                    name="Functional Domains",
                    breaks=c( "InterPro:IPR006703",
                              #"InterPro:IPR016187",
                              #"InterPro:IPR001304",
                              #"InterPro:IPR016186",
                              #  "G3DSA:1.10.533.10",
                              "InterPro:IPR011010",
                              "InterPro:IPR008593",
                              "cd00397",
                              "InterPro:IPR043502",
                              "InterPro:IPR036691",
                              "InterPro:IPR000305",
                              "InterPro:IPR035901",
                              "InterPro:IPR011335",
                              "InterPro:IPR011604",
                              "InterPro:IPR009057",
                              #"InterPro:IPR013098",
                              #"InterPro:IPR007110",
                              #"InterPro:IPR036179",
                              #"InterPro:IPR013783",
                              "InterPro:IPR041588",
                              "InterPro:IPR013762",
                              "InterPro:IPR001584",
                              "InterPro:IPR002104",
                              "InterPro:IPR010998",
                              "InterPro:IPR043128",
                              "InterPro:IPR000477",
                              "InterPro:IPR026960",
                              "InterPro:IPR002156",
                              "InterPro:IPR036397",
                              "InterPro:IPR012337",
                              #"cd16713",
                              #"SSF57850",
                              "cd09275",
                              "cd09274",
                              "cd09276",
                              "cd03714",
                              "cd01647",
                              "cd01650",
                              "InterPro:IPR038717",
                              "InterPro:IPR027805",
                              "InterPro:IPR002492",
                              "InterPro:IPR019080",
                              #"InterPro:IPR001841",
                              #"InterPro:IPR013083",
                              #"G3DSA:4.10.60.10",
                              #"PF13920",
                              "gene_shortened"),
                    labels=c("AIG1-type guanine nucleotide-binding (G) domain",
                             #"C-type lectin fold",
                             #"C-type lectin-like",
                             #"C-type lectin-like/link domain superfamily",
                             # "Death Domain, Fas",
                             "DNA breaking-rejoining enzyme, catalytic core",
                             "DNA N-6-adenine-methyltransferase",
                             "DNA_BRE_C ",
                             "DNA/RNA polymerase superfamily",
                             "Endonuclease/exonuclease/phosphatase superfamily",
                             "GIY-YIG endonuclease",
                             "GIY-YIG endonuclease superfamily",
                             "Restriction endonuclease type II-like",
                             "Exonuclease, phage-type/RecB, C-terminal",
                             "Homeobox-like domain superfamily",
                             #"Immunoglobulin I-set",
                             #"Immunoglobulin-like domain",
                             #"Immunoglobulin-like domain superfamily",
                             #"Immunoglobulin-like fold",
                             "Integrase zinc-binding domain",
                             "Integrase-like, catalytic domain superfamily",
                             "Integrase, catalytic core",
                             "Integrase, catalytic domain",
                             "Integrase/recombinase, N-terminal",
                             "Reverse transcriptase/Diguanylate cyclase domain",
                             "Reverse transcriptase domain",
                             "Reverse transcriptase zinc-binding domain",
                             "Ribonuclease H domain",
                             "Ribonuclease H superfamily",
                             "Ribonuclease H-like superfamily",
                             #"RING-HC_BIRC2_3_7",
                             #"RING/U-box",
                             "RNase_HI_RT_DIRS1 ",
                             "Ty3/Gypsy family of RNase HI in\n long-term repeat retroelements",
                             "non-LTR RNase HI domain of reverse transcriptases",
                             "RT_DIRS1 ",
                             "RT_LTR: Reverse transcriptases (RTs)\n from retrotransposons and retroviruses ",
                             "RT_nLTR_like ",
                             "Tc1-like transposase, DDE domain",
                             "Transposase, Helix-turn-helix domain",
                             "Transposase, Tc1-like",
                             "YqaJ viral recombinase",
                             #"Zinc finger, RING-type",
                             #"Zinc finger, RING/FYVE/PHD-type",
                             #"Zinc finger, CCHC-type",
                             #"Zinc finger, C3HC4 type (RING finger)",
                             "Full Gene Length Shortened")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top"))

###  Export and arrange domain plot with tree
IAP_GENE_raxml_treedata_vertical_legend <- cowplot::get_legend(IAP_GENE_raxml_treedata_vertical)
IAP_GENE_raxml_treedata_vertical_no_legend <- IAP_GENE_raxml_treedata_vertical + 
  theme(legend.position='none')

IAP_GENE_Interproscan_domain_plot_domain_subset_plot_no_legend <- IAP_GENE_Interproscan_domain_plot_domain_subset + theme(legend.position='none')
IAP_GENE_Interproscan_domain_plot_domain_subset_plot_legend <- cowplot::get_legend(IAP_GENE_Interproscan_domain_plot_domain_subset)
IAP_GENE_MY_CV_CG_tree <- IAP_GENE_raxml_treedata_vertical_no_legend + aplot::ylim2(IAP_GENE_Interproscan_domain_plot_domain_subset_plot_no_legend)

IAP_GENE_tr_dom_collapsed <- plot_grid(NULL,IAP_GENE_MY_CV_CG_tree, IAP_GENE_Interproscan_domain_plot_domain_subset_plot_no_legend, ncol=3, align='h', rel_widths = c(0.1, 0.75,0.9)) +
  # Add some space at top for labels
  theme(plot.margin = unit(c(1,0.5,1,0.0), "cm")) 
IAP_GENE_tr_dom_collapsed_legend <- plot_grid(NULL, IAP_GENE_raxml_treedata_vertical_legend, IAP_GENE_Interproscan_domain_plot_domain_subset_plot_legend,
                                         nrow = 1, align="hv", rel_widths  =c(0.5, 0.7,1)) 

## Create combined figure for publication
IAP_GENE_tr_dom_plus_legend <- plot_grid(IAP_GENE_tr_dom_collapsed, IAP_GENE_tr_dom_collapsed_legend,  ncol=1, rel_heights  = c(0.62, 0.1)) +
  # add labels for plot components
  draw_plot_label(c("A","B"), x= c(0.35, 0.48), y = c(1,1), size = 30, family = "sans") +
  # add margin at bottom
  theme(plot.margin = unit(c(0,0.0,1,0.0), "cm")) 

## Export plot with tree and domains aligned : USE THIS FOR PUBLICATION
ggsave(filename = "IAP_GENE_tr_dom_plus_legend_plot_10142020.tiff", plot=IAP_GENE_tr_dom_plus_legend, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
       width = 34,
       height = 27,
       units = "in",
       dpi=300)

## Mar 1st repeat export procedure for revised plot with BIR domains removed 
###  Export and arrange domain plot with tree

IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_plot_no_legend <- IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset + theme(legend.position='none')
IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_plot_legend <- cowplot::get_legend(IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset)

IAP_GENE_tr_dom_collapsed_BIR_rm <- plot_grid(NULL,IAP_GENE_MY_CV_CG_tree, IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_plot_no_legend, ncol=3, align='h',
                                              rel_widths = c(0.05, 0.7,0.9)) +
  # Add some space at top for labels
  theme(plot.margin = unit(c(1,0.5,1,0.0), "cm")) 
IAP_GENE_tr_dom_collapsed_BIR_rm_legend <- plot_grid(NULL, IAP_GENE_raxml_treedata_vertical_legend, IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_plot_legend,
                                              nrow = 1, align="hv", rel_widths  =c(0.5, 0.7,1)) 

## Create combined figure for publication
IAP_GENE_tr_dom_plus_legend_BIR_rm <- plot_grid(IAP_GENE_tr_dom_collapsed_BIR_rm, IAP_GENE_tr_dom_collapsed_BIR_rm_legend,  ncol=1, rel_heights  = c(0.62, 0.1)) +
  # add labels for plot components
  draw_plot_label(c("A","B"), x= c(0.35, 0.48), y = c(1,1), size = 30, family = "sans") +
  # add margin at bottom
  theme(plot.margin = unit(c(1,3,1,0.0), "cm")) 

## Export plot with tree and domains aligned : USE THIS FOR PUBLICATION
ggsave(filename = "IAP_GENE_tr_dom_plus_legend_BIR_rm_03092021.tiff", plot=IAP_GENE_tr_dom_plus_legend_BIR_rm, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
       width = 35.5,
       height = 27,
       units = "in",
       dpi=300)

### Export plot with just the genes that have retrotransposon machinery 
ggsave(filename = "IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_only_retro_05102021.tiff", plot=IAP_GENE_Interproscan_domain_plot_domain_BIR_rm_subset_only_retro , device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
       width = 23,
       height = 15,
       units = "in",
       dpi=300)


#### PLOT FULL IAP PROTEIN TREE ####
# Helpful online tutorial regarding tool: https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
# Tree data vignette https://yulab-smu.github.io/treedata-book/faq.html#different-x-labels-for-different-facet-panels
# Load and parse RAxML bipartitions bootstrapping file with treeio. File input is the bootstrapping analysis output
IAP_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_IAP_HMMER_Interpro_XP_list_all_MSA_RAxML")
IAP_raxml

# Get information in the features/attributes of the tree (the XP labels) with get_data
get.fields(IAP_raxml)
get.data(IAP_raxml)

# Convert to tibble tree dataframe object with tidytree to add external data
IAP_raxml_tibble <- as_tibble(IAP_raxml)

# Join protein product name,gene or locus, and species
colnames(IAP_raxml_tibble)[4] <- "protein_id"
IAP_raxml_tibble <- left_join(IAP_raxml_tibble, BIR_XP_gff_species_join, by = "protein_id")
colnames(IAP_raxml_tibble)[4] <- "label"

# Add combined gene and locus name column 
IAP_raxml_tibble$gene_locus_tag <- coalesce(IAP_raxml_tibble$gene, IAP_raxml_tibble$locus_tag)

# Remove text after isoform so I can collapse protein names into shorter list
IAP_raxml_tibble$product <- gsub("isoform.*", "", IAP_raxml_tibble$product)
IAP_raxml_tibble$product <- trimws(IAP_raxml_tibble$product , which = "both")

View(IAP_raxml_tibble %>% filter(!grepl("uncharacterized",product)) %>% distinct(product))

# Join with alias info
IAP_alias <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_shortened_product.csv")
IAP_raxml_tibble <- left_join(IAP_raxml_tibble, IAP_alias)
# export IAP_raxml_tibble
save(IAP_raxml_tibble, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_raxml_tibble.Rdata")

# Fill in blanks with uncharacterized locus name
IAP_raxml_tibble$alias[is.na(IAP_raxml_tibble$alias)] <- IAP_raxml_tibble$product[is.na(IAP_raxml_tibble$alias)]

# Remove uncharacterized protein and just keep gene name for those uncharacterized
IAP_raxml_tibble$alias <- gsub("uncharacterized protein", "",IAP_raxml_tibble$alias)

# Convert to treedata object to store tree plus outside data
IAP_raxml_treedata <- as.treedata(IAP_raxml_tibble)

# save treedata
save(IAP_raxml_treedata, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_raxml_treedata.Rdata")

## Plot circular protein tree to use for paper ##
IAP_raxml_treedata_circular_product <- ggtree(IAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=alias,angle=angle), size =1.8, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
theme(legend.position = "bottom", 
      legend.text = element_text(face = "italic", size=8, family="sans"),
      legend.title = element_text(size=12, family="sans")) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32",  "#a68340",
                                                 "#a3c763", "#c257b0", "#c083d0","#59a1cf","#c2134a","#ead76b"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis", 
"Elysia_chlorotica","Lottia_gigantea", "Octopus_bimaculoides", "Octopus_vulgaris", "Pomacea_canaliculata", "Biomphalaria_glabrata","Aplysia_californica"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis", 
                                 "Elysia chlorotica","Lottia gigantea", "Octopus bimaculoides", "Octopus vulgaris", "Pomacea canaliculata", "Biomphalaria glabrata", "Aplysia californica")) 
 
# Export plot to file to put together with gene tree for paper
ggsave(filename = "IAP_full_protein_circular_tree.tiff", plot=IAP_raxml_treedata_circular_product, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_protein_tree/",
       width = 10 ,
       height = 10,
       units = "in",
       dpi=300)

## Plot circular tree with gene labels ##
IAP_raxml_treedata_circular_gene <- ggtree(IAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=gene_locus_tag,angle=angle), size =1.8, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  #xlim(-100,100)  +
  # fix legend appearance
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32",  "#a68340",
                                                 "#a3c763", "#c257b0", "#c083d0","#59a1cf","#c2134a","#ead76b"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis", 
                                                                                                                                             "Elysia_chlorotica","Lottia_gigantea", "Octopus_bimaculoides", "Octopus_vulgaris", "Pomacea_canaliculata", "Biomphalaria_glabrata","Aplysia_californica"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis", 
                                 "Elysia chlorotica","Lottia gigantea", "Octopus bimaculoides", "Octopus vulgaris", "Pomacea canaliculata", "Biomphalaria glabrata", "Aplysia californica")) 

# Export plot to file for summary powerpoint - July 12th, 2020 
ggsave(filename = "IAP_full_gene_circular_tree.tiff", plot=IAP_raxml_treedata_circular_gene, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/",
       width = 10 ,
       height = 10,
       units = "in",
       dpi=300)

# PLOT AS PROTEIN TREE WITH GENE LABEL TO SEARCH FOR POTENTIAL ARTIFACTS
IAP_raxml_treedata_circular_gene <- ggtree(IAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=gene_locus_tag,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-80,80)  

IAP_raxml_treedata_circular_gene + scale_color_discrete(name = "Species", labels = c("Aplysia californica", 
"Biomphalaria glabrata", "Crassostrea gigas", "Crassostrea virginica","Elysia chlorotica","Lottia gigantea","Mizuhopecten yessoensis","Octopus bimaculoides",
"Octopus vulgaris", "Pomacea canaliculata","NA"))

# Where are the potential gene artifacts located?
# list of IAP potential gene artifacts 
#"LOC111111659", 
#"LOC111114013" ,
#"LOC111103682" ,
#"LOC111132589" , # node 454
#"LOC111102106" ,
#"LOC111114070"

IAP_raxml_treedata_artifact <- tree_subset(IAP_raxml_treedata, "XP_022336127.1", levels_back = 8)
IAP_raxml_treedata_artifact_tibble <- as.tibble(IAP_raxml_treedata_artifact)
IAP_raxml_treedata_artifact_tree <- ggtree(IAP_raxml_treedata_artifact, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=gene_locus_tag), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

#### PLOT GENE TREE AND PROTEIN TREE SIDE BY SIDE FOR FIGURE ####

combined_trees <- plot_grid(IAP_GENE_all_species_raxml_treedata_circular_gene, IAP_raxml_treedata_circular_product) +
                            draw_plot_label(c("A","B"), x= c(0, 0.5), y = c(0.8,0.8), size = 30, family = "sans") +
                  theme(plot.margin = unit(c(0,0,0,0),"cm"))

## Export plot with tree and domains aligned : USE THIS FOR PUBLICATION
ggsave(filename = "IAP_gene_prot_combined_trees_10232020.tiff", plot=combined_trees, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/",
       units = "in",
       height = 15, width = 18,
       dpi=300)

## Febraury 12th changing plot orientation and stacking the two trees 
# remove tree plot margins 
IAP_GENE_all_species_raxml_treedata_circular_gene_no_margin <- ggtree(IAP_GENE_all_species_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=label,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans"),
        plot.margin = margin(0,0,0,0)) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32",  "#a68340",
                                                 "#a3c763", "#c257b0", "#c083d0","#59a1cf","#c2134a","#ead76b"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis", 
                                                                                                                                             "Elysia_chlorotica","Lottia_gigantea", "Octopus_bimaculoides", "Octopus_vulgaris", "Pomacea_canaliculata", "Biomphalaria_glabrata","Aplysia_californica"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis", 
                                 "Elysia chlorotica","Lottia gigantea", "Octopus bimaculoides", "Octopus vulgaris", "Pomacea canaliculata", "Biomphalaria glabrata", "Aplysia californica")) 


IAP_raxml_treedata_circular_product_no_margin <- ggtree(IAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=alias,angle=angle), size =1.8, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans"),
        plot.margin = margin(0,0,0,0)) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=0.8) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=0.8) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=0.8) +
  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")) ) + # need to override aes to get rid of "a"
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32",  "#a68340",
                                                 "#a3c763", "#c257b0", "#c083d0","#59a1cf","#c2134a","#ead76b"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis", 
                                                                                                                                             "Elysia_chlorotica","Lottia_gigantea", "Octopus_bimaculoides", "Octopus_vulgaris", "Pomacea_canaliculata", "Biomphalaria_glabrata","Aplysia_californica"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis", 
                                 "Elysia chlorotica","Lottia gigantea", "Octopus bimaculoides", "Octopus vulgaris", "Pomacea canaliculata", "Biomphalaria glabrata", "Aplysia californica")) 

combined_stacked_trees <- ggpubr::ggarrange(IAP_GENE_all_species_raxml_treedata_circular_gene_no_margin,IAP_raxml_treedata_circular_product_no_margin,
                                    ncol = 1, labels = c("A","B"), font.label = list(size = 30, family = "sans"), common.legend = TRUE,
                                    legend = "bottom")

## add species tree so I can label the gene numbers 

ggsave(filename = "IAP_gene_prot_combined_trees_2_12_21.tiff", combined_stacked_trees, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/",
       units = "in",
       height = 18, width = 10,
       dpi=300)


#### PLOT IAP MY, CV, CG PROTEIN TREE WITH DOMAIN INFO ####

# Load and parse RAxML bipartitions bootstrapping file with treeio. File input is the bootstrapping analysis output
IAP_MY_CV_CG_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG_MSA_RAxML")
IAP_MY_CV_CG_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
IAP_MY_CV_CG_raxml_tibble <- as_tibble(IAP_MY_CV_CG_raxml)

# Join protein product name,gene or locus, and species
colnames(IAP_MY_CV_CG_raxml_tibble)[4] <- "protein_id"
IAP_MY_CV_CG_raxml_tibble <- left_join(IAP_MY_CV_CG_raxml_tibble, BIR_XP_gff_species_join, by = "protein_id")
colnames(IAP_MY_CV_CG_raxml_tibble)[4] <- "label"

# Add combined gene and locus name column 
IAP_MY_CV_CG_raxml_tibble$gene_locus_tag <- coalesce(IAP_MY_CV_CG_raxml_tibble$gene, IAP_MY_CV_CG_raxml_tibble$locus_tag)

# Remove text after isoform so I can collapse protein names into shorter list
IAP_MY_CV_CG_raxml_tibble$product <- gsub("isoform.*", "", IAP_MY_CV_CG_raxml_tibble$product)
IAP_MY_CV_CG_raxml_tibble$product <- trimws(IAP_MY_CV_CG_raxml_tibble$product , which = "both")

# Join with alias info
IAP_alias <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_shortened_product.csv")
IAP_MY_CV_CG_raxml_tibble <- left_join(IAP_MY_CV_CG_raxml_tibble, IAP_alias)

# Fill in blanks with uncharacterized locus name
IAP_MY_CV_CG_raxml_tibble$alias[is.na(IAP_MY_CV_CG_raxml_tibble$alias)] <- IAP_MY_CV_CG_raxml_tibble$product[is.na(IAP_MY_CV_CG_raxml_tibble$alias)]

# Remove uncharacterized protein and just keep gene name for those uncharacterized
IAP_MY_CV_CG_raxml_tibble$alias <- gsub("uncharacterized protein", "",IAP_MY_CV_CG_raxml_tibble$alias)

# Convert to treedata object to store tree plus outside data
IAP_MY_CV_CG_raxml_treedata <- as.treedata(IAP_MY_CV_CG_raxml_tibble)

# Drop tips to collapse Mizuhopecten tips (done by manually looking at tree and keeping the outgroup for each group)
IAP_to_drop <- c("236","237","238","239","240","205","206","115","116","117","120","121","122","200","201",
                 "196","185","186","187","194","190","192","183","179","189","161","162","141","142","146","143",
                 "144","145","149","154","128","129","130","13","17","15","16","20","18","19","21","26","22","23","24",
                 "25","50","52","53","54","59","60","61","62","63","70","71")

# get labels for nodes to drop
IAP_to_drop_label <- IAP_MY_CV_CG_raxml_tibble[IAP_MY_CV_CG_raxml_tibble$node %in% IAP_to_drop,]
IAP_to_drop_label <- IAP_to_drop_label$label

# labels for when to add shapes for dropped tips 
IAP_to_label <- c("241","207","118","119","202","184","195","188","193","189","191","182","181","163","147","148","153","127","14","12","49","51","63","72")
IAP_to_label_XP <- IAP_MY_CV_CG_raxml_tibble[IAP_MY_CV_CG_raxml_tibble$node %in% IAP_to_label,]
IAP_to_label_XP <- IAP_to_label_XP$label

#Collapse
IAP_MY_CV_CG_raxml_treedata_collapsed <- drop.tip(IAP_MY_CV_CG_raxml_treedata , IAP_to_drop_label)
View(as.tibble(IAP_MY_CV_CG_raxml_treedata_collapsed))
IAP_MY_CV_CG_raxml_treedata_collapsed # tips were dropped

# get new node numbers for where to add shapes 
IAP_collapsed_tibble <- as.tibble(IAP_MY_CV_CG_raxml_treedata_collapsed)
IAP_shape_label <- IAP_collapsed_tibble[IAP_collapsed_tibble$label %in% IAP_to_label_XP, ]
IAP_shape_node <- IAP_shape_label$node
length(IAP_shape_node)

# check nrows
IAP_collapsed_tibble %>% filter(!is.na(label)) %>% count() # 184 lines to plot 

# Add Domain type groupings for later statistics, remember the below data frame has been collapsed so that proteins with identical sequence have been removed 
IAP_domain_structure1 <- read_csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_Domain_Structure_groups.csv")
colnames(IAP_domain_structure1)[1] <- "label"
IAP_collapsed_tibble_domain_type <- left_join(IAP_collapsed_tibble, IAP_domain_structure1[,c("label","Domain_Name","Number")])

IAP_collapsed_tibble_domain_type %>% group_by(Species) %>% count()
#Species                     n
#<chr>                   <int>
#  1 Crassostrea_gigas          54
#2 Crassostrea_virginica      95
#3 Mizuhopecten_yessoensis    35

# Plot collapsed tree
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed <- 
  ggtree(IAP_MY_CV_CG_raxml_treedata_collapsed, aes(color=Species, fill=Species),  branch.length = "none") + 
  geom_tiplab(aes(label=gene), fontface="bold", size =3.5, offset=0) + # geom_tiplab2 flips the labels correctly
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=2.0) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=2.0) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=2.0) +
  # Add shape for tips removed (see IAP_shape_node above)
  geom_point2(aes(subset=(node==12)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==13)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==36)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==37)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==48)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==91)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==92)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==97)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==108)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==109)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==113)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==120)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==137)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==138)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==139)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==140)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==141)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==142)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==143)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==147)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==150)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==179)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  ## Add clade labels for the 21 domain groups  domain groups using the internal node number
  geom_cladelabel(261, label="1",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') + # get node order from below 
  geom_cladelabel(254, label="2",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(233, label="3",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(228, label="4",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(217, label="5",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(211, label="6",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(207, label="7",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(198, label="8",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(201, label="9",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(311, label="10", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(308, label="11", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(291, label="12", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(304, label="13", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(329, label="14", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(340, label="15", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(343, label="16", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(347, label="17", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(353, label="18", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(359, label="19", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(188, label="20", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(366, label="21", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black', extend = 0.5) +
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=14, family="sans"),
        legend.title = element_text(size=16, family="sans")) +
  #geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 2.0, fontface="bold") + # allows for subset
  xlim(-70,31.8) + #change scaling so branch lengths are smaller and all alias labels are showing
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis")) +
  guides(col = guide_legend(ncol =1, title.position = "top", override.aes = aes(label = "")) ) # need to override aes to get rid of "a"

# find internal node number to use for the clade labels above 
#IAP_MY_CV_CG_raxml_treedata_vertical_collapsed + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

#### PLOT ALL IAP DOMAINS WITHOUT BIR TYPE INFORMATION ####
# Use combination of geom_segment and geom_rect and combine plot with vertical tree using ggarrange from ggpubr
# Get only the Interproscan domains for my proteins of interest
IAP_MY_CV_CG_raxml_tibble_join <- IAP_collapsed_tibble %>% filter(!is.na(label)) # remove rows with just bootstrap information
colnames(IAP_MY_CV_CG_raxml_tibble_join)[4] <- "protein_id"
BIR_XP_gff_Interpro_Domains <-  left_join(IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")], BIR_XP_gff)
BIR_XP_gff_Interpro_Domains_only <- BIR_XP_gff_Interpro_Domains %>% 
  filter(source =="CDD" | grepl("InterPro:IPR", Dbxref) | grepl("G3DSA:1.10.533.10", Name)) # keep Interproscan domain lines, CDD NCBI lines, and death domain structure

BIR_XP_gff_Interpro_Domains_all <- BIR_XP_gff_Interpro_Domains

# Which domains were removed
BIR_XP_gff_Interpro_Domains_Name <- BIR_XP_gff_Interpro_Domains %>% distinct(Name) 
BIR_XP_gff_Interpro_Domains_only_Name <- BIR_XP_gff_Interpro_Domains_only%>% distinct(Name) 

BIR_XP_gff_Interpro_Domains_Name[!(BIR_XP_gff_Interpro_Domains_Name$Name %in% BIR_XP_gff_Interpro_Domains_only_Name$Name),] # 13 were removed 
#  2 SSF57924 : IAP superfamily          
#  3 G3DSA:1.10.1170.10: CATH Superfamily 1.10.1170.10 , Inhibitor Of Apoptosis Protein (2mihbC-IAP-1); Chain A
#  4 PF13920: Zinc finger, C3HC4 type (RING finger)           
#  5 mobidb-lite       
#  6 G3DSA:1.10.533.10: this includes the Death domains, Fas
#  7 Coil              
#  8 G3DSA:1.10.8.10: Ubiquitin-associated (UBA) domain
#  9 PIRSF036836 : E3 ubiquitin-protein ligase BOI-like      
#  10 G3DSA:3.30.710.10: Potassium Channel Kv1.1; Chain A 
#  11 SM00212: Ubiquitin-conjugating enzyme E2, catalytic domain homologues           
#  12 G3DSA:3.40.50.300 : P-loop containing nucleotide triphosphate hydrolases
#  13 G3DSA:1.10.10.2190: unnamed domain

# Which domains do I keep since some of these (other than the death domain) are also identified by Interproscan 

# filter out NA source lines to get full length of the protein
BIR_XP_gff_Interpro_Domains_fullprot <- BIR_XP_gff_Interpro_Domains %>% 
  filter(is.na(source))
BIR_XP_gff_Interpro_Domains_all_fullprot <- BIR_XP_gff_Interpro_Domains_all %>% 
  filter(is.na(source))

nrow(BIR_XP_gff_Interpro_Domains_fullprot %>% filter(is.na(source))) # 184
nrow(IAP_MY_CV_CG_raxml_tibble_join %>% filter(!is.na(protein_id))) # 184 - they agree, all proteins were found - number is less because mizuhopecten proteins were collapsed 
nrow(BIR_XP_gff_Interpro_Domains_all_full_prot) 

# Fill in the CDD rows that have NULL for DBxref with the Name column
BIR_XP_gff_Interpro_Domains_only$Dbxref[BIR_XP_gff_Interpro_Domains_only$Dbxref == "character(0)" ] <- "CDD"
BIR_XP_gff_Interpro_Domains_all$Dbxref[BIR_XP_gff_Interpro_Domains_all$Dbxref == "character(0)" ] <- "CDD"

# unlist
BIR_XP_gff_Interpro_Domains_only <- BIR_XP_gff_Interpro_Domains_only %>% unnest(Dbxref)
BIR_XP_gff_Interpro_Domains_all <- BIR_XP_gff_Interpro_Domains_all  %>% unnest(Dbxref)

# Change CDD rows to be the Name column
BIR_XP_gff_Interpro_Domains_only <- BIR_XP_gff_Interpro_Domains_only %>% mutate(Dbxref = ifelse(Dbxref == "CDD", Name, Dbxref))
BIR_XP_gff_Interpro_Domains_all<- BIR_XP_gff_Interpro_Domains_all  %>% mutate(Dbxref = ifelse(Dbxref == "CDD", Name, Dbxref))

# Get the node order from collapsed IAP tree
IAP_MY_CV_CG_raxml_treedata_tip  <- fortify(IAP_MY_CV_CG_raxml_treedata_collapsed) # not changing code from here down
IAP_MY_CV_CG_raxml_treedata_tip <- subset(IAP_MY_CV_CG_raxml_treedata_tip, isTip)
IAP_MY_CV_CG_raxml_treedata_tip_order <- IAP_MY_CV_CG_raxml_treedata_tip$label[order(IAP_MY_CV_CG_raxml_treedata_tip$y, decreasing=TRUE)]

# Reorder the protein and polygon
BIR_XP_gff_Interpro_Domains_fullprot <- BIR_XP_gff_Interpro_Domains_fullprot[match(IAP_MY_CV_CG_raxml_treedata_tip_order, BIR_XP_gff_Interpro_Domains_fullprot$protein_id),]
IAP_MY_CV_CG_raxml_treedata_tip_order <- as.data.frame(IAP_MY_CV_CG_raxml_treedata_tip_order)
colnames(IAP_MY_CV_CG_raxml_treedata_tip_order)[1] <- "protein_id"
BIR_XP_gff_Interpro_Domains_only <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, BIR_XP_gff_Interpro_Domains_only)

BIR_XP_gff_Interpro_Domains_all_fullprot <- BIR_XP_gff_Interpro_Domains_all_fullprot[match(IAP_MY_CV_CG_raxml_treedata_tip_order, BIR_XP_gff_Interpro_Domains_all_fullprot$protein_id),]
BIR_XP_gff_Interpro_Domains_all <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, BIR_XP_gff_Interpro_Domains_all)

# Add polygon height
BIR_XP_gff_Interpro_Domains_only_ID  <- BIR_XP_gff_Interpro_Domains_only  %>% distinct(protein_id) 
BIR_XP_gff_Interpro_Domains_only_ID <- BIR_XP_gff_Interpro_Domains_only_ID %>% 
  mutate(height_start = rev(as.numeric(row.names(BIR_XP_gff_Interpro_Domains_only_ID )) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(BIR_XP_gff_Interpro_Domains_only_ID)) + .5))

BIR_XP_gff_Interpro_Domains_all_ID  <- BIR_XP_gff_Interpro_Domains_all  %>% distinct(protein_id) 
BIR_XP_gff_Interpro_Domains_all_ID <- BIR_XP_gff_Interpro_Domains_all_ID %>% 
  mutate(height_start = rev(as.numeric(row.names(BIR_XP_gff_Interpro_Domains_all_ID )) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(BIR_XP_gff_Interpro_Domains_all_ID)) + .5))

# Join back in height
BIR_XP_gff_Interpro_Domains_only <- left_join(BIR_XP_gff_Interpro_Domains_only , BIR_XP_gff_Interpro_Domains_only_ID )
BIR_XP_gff_Interpro_Domains_all <- left_join(BIR_XP_gff_Interpro_Domains_all , BIR_XP_gff_Interpro_Domains_all_ID )

# Set factor level order of the nodes set levels in reverse order
BIR_XP_gff_Interpro_Domains_only$node <- factor(BIR_XP_gff_Interpro_Domains_only$node, levels = unique(BIR_XP_gff_Interpro_Domains_only$node))
BIR_XP_gff_Interpro_Domains_only$Dbxref <- factor(BIR_XP_gff_Interpro_Domains_only$Dbxref, levels = unique(BIR_XP_gff_Interpro_Domains_only$Dbxref))
BIR_XP_gff_Interpro_Domains_fullprot$node <- factor(BIR_XP_gff_Interpro_Domains_fullprot$node, levels = rev(BIR_XP_gff_Interpro_Domains_fullprot$node))

BIR_XP_gff_Interpro_Domains_all$node <-   factor(BIR_XP_gff_Interpro_Domains_all$node, levels = unique(BIR_XP_gff_Interpro_Domains_all$node))
BIR_XP_gff_Interpro_Domains_all$Dbxref <- factor(BIR_XP_gff_Interpro_Domains_all$Dbxref, levels = unique(BIR_XP_gff_Interpro_Domains_all$Dbxref))
BIR_XP_gff_Interpro_Domains_all_fullprot$node <- factor(BIR_XP_gff_Interpro_Domains_all_fullprot$node, levels = rev(BIR_XP_gff_Interpro_Domains_all_fullprot$node))


# Plotting all domains 
IAP_Interproscan_all_domain_plot <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data =BIR_XP_gff_Interpro_Domains_all_fullprot,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add boxes with geom_rect 
  geom_rect(data=BIR_XP_gff_Interpro_Domains_all,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Dbxref)) +
  #add labels
  labs(y = NULL, x = "Protein Domain Position (aa)") +
  # add text labels
  #geom_text(data=BIR_XP_gff_Interpro_Domains_fullprot,aes(x= end, y = node, label=alias),
  #          size=2.0, hjust=-.15, check_overlap = TRUE) + 
  # text theme
  theme_bw() + 
  # plot theme
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=8, family="sans"),
        legend.title = element_text(size=12, family="sans"))+
  # Change y axis ticks
  scale_x_continuous(breaks=c(0,500,1000,1500,2000,3000), expand = c(0,0)) + 
  # Change domain labels 
  scale_fill_manual(values=c("#d44172","#6d8dd7","#c972c4","#ca8853","#cd9c2e","#92b540",
                             "#da83b4","#45c097","#cba950","#65c874","#5b3788","#8a371d","#b1457b",
                             "#be4a5b","#6971d7","#50893b","#d55448","#c46a2f","#8a8a39","#d1766b"), 
                    name="Functional Domains",
                    breaks=c("cd16713",
                             "\"InterPro:IPR001370\"",
                             "\"InterPro:IPR013083\"",
                             "\"InterPro:IPR001841\"",
                             "\"InterPro:IPR015940\"",
                             "\"InterPro:IPR003131\"",
                             "cd18316",
                             "\"InterPro:IPR011333\"",
                             "\"InterPro:IPR000210\"",
                             "\"InterPro:IPR000608\"",
                             "\"InterPro:IPR016135\"",
                             "\"InterPro:IPR022103\"",
                             "\"InterPro:IPR036322\"",
                             "\"InterPro:IPR011047\"",
                             "\"InterPro:IPR019775\"",
                             "\"InterPro:IPR017907\"",
                             "\"InterPro:IPR038765\"",
                             "\"InterPro:IPR032171\"",
                             "\"InterPro:IPR027417\"",
                             "cd14321"),
                    labels=c("RING-HC_BIRC2_3_7",
                             "BIR repeat",
                             "Zinc finger, RING/FYVE/PHD-type",
                             "Zinc finger, RING-type",
                             "Ubiquitin-associated domain",
                             "Potassium channel tetramerisation-type BTB domain",
                             "BTB/POZ domain",
                             "SKP1/BTB/POZ domain superfamily",
                             "BTB/POZ domain",
                             "Ubiquitin-conjugating enzyme E2",
                             "Ubiquitin-conjugating enzyme/RWD-like",
                             "Baculoviral IAP repeat-containing protein 6",
                             "WD40-repeat-containing domain superfamily",
                             "Quinoprotein alcohol dehydrogenase-like superfamily",
                             "WD40 repeat, conserved site",
                             "Zinc finger, RING-type, conserved site",
                             "Papain-like cysteine peptidase superfamily",
                             "C-terminal of Roc (COR) domain",
                             "P-loop containing nucleoside triphosphate hydrolase",
                             "Ubiquitin-associated domain")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top"))


#### CHARACTERIZE TYPE 1 AND TYPE II BIR IAP REPEATS ####

# Fill with fasta protein sequences for each IAP can be found in
# ("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all_rm_dup.fa")

# Export the BIR cd00022, also used for coordinates for the model organism domain types 
BIR_XP_gff_Interpro_Domains_only_cd00022 <- BIR_XP_gff_Interpro_Domains_only %>% filter(Name =="cd00022")

# Add position column with the position of the BIR domain named, add Name column
BIR_XP_gff_Interpro_Domains_only_cd00022$alias <- str_trim(BIR_XP_gff_Interpro_Domains_only_cd00022$alias, side="both")

BIR_XP_gff_Interpro_Domains_only_cd00022 <- BIR_XP_gff_Interpro_Domains_only_cd00022 %>% 
  arrange(protein_id, start) %>% #arrange values by start 
  group_by(protein_id) %>% # group_by protein id
  mutate(Domain_Number = paste(signature_desc, row_number(), sep="")) %>% 
  ungroup() %>%
  mutate(Name=paste(protein_id, Domain_Number, alias, sep="_"))
# Add name column with concatenated protein_id and the BIR domain position
# remember that this file has the Mizuhopecten yessoensis sequences collapsed 

# Export the file in BED format 
#write.table(BIR_XP_gff_Interpro_Domains_only_cd00022[,c("protein_id","start","end","Name")], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_XP_gff_Interpro_Domains_only_cd00022.bed",
#            quote = FALSE,col.names = FALSE, row.names=FALSE, sep="\t")

## IN bluewaves, get the BIR sequence for each domain, make multiple sequence alignment in MAFFTT
## Load Sequence file from bluewaves back into R and search for subsequences of interest 
BIR_domain_model_MY_CV_CG <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq.fa")

# Locate string pattern for each type
#Type I = K/76 or R/76, H77, V/80 or L/80, C/84; Type II = E/76 or Q/76, H77, W/80 or H/80, C/84 
class(BIR_domain_model_MY_CV_CG$seq.text) # character

# Locate Type I and Type II 
# Use case_when for conditional mutate
BIR_domain_model_MY_CV_CG_type  <- BIR_domain_model_MY_CV_CG %>%
mutate(Type = case_when(
grepl("KH..V...C", seq.text) ~ "T1",
grepl("KH..L...C", seq.text) ~ "T1",
grepl("RH..V...C", seq.text)  ~ "T1",
grepl("RH..L...C", seq.text)  ~ "T1",
grepl("EH..W...C", seq.text) ~ "T2",
grepl("EH..H...C", seq.text) ~ "T2",
grepl("QH..W...C", seq.text) ~ "T2",
grepl("QH..H...C", seq.text) ~ "T2",
TRUE ~ NA_character_)) # final line is a catch all for things that don't match 
# quite a few NA's 

# Add in Species information (change some column names first for joining)
BIR_XP_gff_Interpro_Domains_only_cd00022_species <- BIR_XP_gff_Interpro_Domains_only_cd00022
colnames(BIR_XP_gff_Interpro_Domains_only_cd00022_species)[16] <-"seq.name"
IAP_collapsed_tibble_species <- IAP_collapsed_tibble[,c("label","Species")]
colnames(IAP_collapsed_tibble_species)[1] <- "protein_id"

BIR_domain_model_MY_CV_CG_type <- left_join(BIR_domain_model_MY_CV_CG_type,BIR_XP_gff_Interpro_Domains_only_cd00022_species[,c("seq.name","protein_id")])
BIR_domain_model_MY_CV_CG_type <- left_join(BIR_domain_model_MY_CV_CG_type, IAP_collapsed_tibble_species)                                           

# Export all the non Type 1 or Type II proteins to rerun alignment and Tree in MAFFT in order to better view multiple alignment patterns
BIR_domain_model_MY_CV_CG_non_T1_T2 <- BIR_domain_model_MY_CV_CG_type %>% filter(is.na(Type)) %>% distinct(seq.name, .keep_all = TRUE)
#phylotools::dat2fasta(BIR_domain_model_MY_CV_CG_non_T1_T2[,c("seq.name","seq.text")],
#            outfile="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_non_T1_T2.fa")

## Export only the C. vir and C. gigas non type I or Type II BIR repeats and run multiple alignment. This will allow me to View consensus percentages of AA for just C. vir and C. gig repeats
BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig <- BIR_domain_model_MY_CV_CG_type %>% filter(is.na(Type) & Species != "Mizuhopecten_yessoensis") %>% distinct(seq.name, .keep_all = TRUE)
#phylotools::dat2fasta(BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig[,c("seq.name","seq.text")],
#                      outfile="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig.fa")

## Export only the C. vir and C. gigas ALL BIR repeats and run multiple alignment. This will allow me to View consensus percentages of AA for ALL C. vir and C. gig repeats
BIR_domain_model_MY_CV_CG_Cvir_Cgig <- BIR_domain_model_MY_CV_CG_type %>% filter(Species != "Mizuhopecten_yessoensis") %>% distinct(seq.name, .keep_all = TRUE)
#phylotools::dat2fasta(BIR_domain_model_MY_CV_CG_Cvir_Cgig[,c("seq.name","seq.text")],
#                      outfile="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_Cvir_Cgig.fa")


# How many Conserved Type II domains in C gigas and C virginica?
BIR_domain_model_MY_CV_CG_type_distinct <- BIR_domain_model_MY_CV_CG_type %>%   distinct(seq.name, .keep_all = TRUE)
BIR_domain_model_MY_CV_CG_type_distinct %>% 
  group_by(Type, Species) %>% 
  count() 

### Based on alignments, adding more sequence patterns to look for 
BIR_domain_model_MY_CV_CG_type_updated  <- BIR_domain_model_MY_CV_CG_type_distinct %>%
  mutate(Type = case_when(
    grepl("KH..V...C", seq.text) ~ "T1",   # Conserved model org Type 1
    grepl("KH..L...C", seq.text) ~ "T1",  # Conserved model org Type 1
    grepl("RH..V...C", seq.text)  ~ "T1",  # Conserved model org Type 1
    grepl("RH..L...C", seq.text)  ~ "T1",  # Conserved model org Type 1
    grepl("H..L...C", seq.text) ~ "T1", # this pattern should be called conserved Type I
    grepl("EH..W...C", seq.text) ~ "T2",  # Conserved model org Type 2 
    grepl("EH..H...C", seq.text) ~ "T2",  # Conserved model org Type 2 
    grepl("QH..W...C", seq.text) ~ "T2",  # Conserved model org Type 2 
    grepl("QH..H...C", seq.text) ~ "T2",  # Conserved model org Type 2 
    grepl("SFCC", seq.text) ~ "Non_Zinc_binding", # lacks C for zinc binding
    grepl("TFCC", seq.text) ~ "Non_Zinc_binding", # lacks C for zinc binding 
    grepl("EH..GSR.C",seq.text) ~ "TX", # new type with Glycine and no proline
    grepl("EH..YKP",seq.text) ~ "T2-like_1",
    grepl("EHKN.FP",seq.text) ~ "T2-like_2",
    grepl("H.NMSP",seq.text) ~ "T1-like_1",
    grepl("IHRQQSP", seq.text) ~ "T1-like_2",
    grepl("TS.I.AIH..ISP", seq.text) ~ "T1_like_3",
    grepl("VHKENSP",seq.text) ~ "T1_like_4",
    grepl("CYSCHVVHEGW", seq.text) ~ "TY",
    grepl("RL..FK",seq.text) ~ "T2_like_3", 
    grepl("EHLDK", seq.text) ~ "T2_like_4",
    grepl("EH.KY", seq.text) ~ "T2_like_5",
    TRUE ~ "Unique")) # final line is a catch all for things that don't match 

# Fill in NAs after checking which species they come from 
BIR_domain_model_MY_CV_CG_type_updated %>% filter(is.na(Species))
model_org_species <- data.frame(seq.name = c("NP_001261918.1_BIR2_T2_DIAP1" ,"NP_001261918.1_BIR1_T2_DIAP1" ,"NP_477127.1_BIR3_T2_DIAP2" ,"NP_477127.1_BIR2_T2_DIAP2" ,
                                             "NP_477127.1_BIR1_T1_DIAP2" ,"NP_001243092.1_BIR3_T3_cIAP1" ,"NP_001243092.1_BIR2_T2_cIAP1" ,"NP_001243092.1_BIR1_T1_cIAP1"
                                             ,"NP_919376.1_BIR3_T2_cIAP1" ,"NP_919376.1_BIR2_T2_cIAP1" ,"NP_919376.1_BIR1_T1_cIAP1" ,"NP_892007.1_BIR3_T2_cIAP2" ,
                                             "NP_892007.1_BIR2_T2_cIAP2" ,"NP_892007.1_BIR1_T1_cIAP2" ,"NP_031490.2_BIR3_T2_cIAP2" ,"NP_031490.2_BIR2_T2_cIAP2" ,
                                             "NP_031490.2_BIR1_T1_cIAP2" ,"NP_001159.2_BIR1_T2_BIRC5" ,"NP_919378.1_BIR1_T2_BIRC5a" ,"NP_057336.3_BIR1_T2_BIRC6" ,
                                             "XP_009291311.1_BIR1_T2_BIRC6" ,"NP_001158.2_BIR3_T2_XIAP" ,"NP_001158.2_BIR2_T2_XIAP" ,"NP_001158.2_BIR1_T1_XIAP" ,
                                             "NP_919377.2_BIR3_T2_XIAP" ,"NP_919377.2_BIR2_T2_XIAP" ,"NP_919377.2_BIR1_T1_XIAP" ,"NP_647478.1_BIR1_T2_BIRC7a" ,
                                             "NP_071444.1_BIR1_T2_BIRC7b" ,"AAH39318.1_BIR1_T2_BIRC8"),
                                Species = c("Drosophila_melanogaster","Drosophila_melanogaster","Drosophila_melanogaster","Drosophila_melanogaster",
                                            "Drosophila_melanogaster","Homo_sapiens","Homo_sapiens","Homo_sapiens","Danio_rerio ","Danio_rerio ","Danio_rerio ",
                                            "Homo_sapiens","Homo_sapiens","Homo_sapiens","Mus_musculus","Mus_musculus","Mus_musculus","Homo_sapiens","Danio_rerio",
                                            "Homo_sapiens","Danio_rerio","Homo_sapiens","Homo_sapiens","Homo_sapiens","Danio_rerio","Danio_rerio","Danio_rerio",
                                            "Homo_sapiens","Homo_sapiens","Homo_sapiens "))

#NP_001261918.1_BIR2_T2_DIAP1     "Drosophila_melanogaster",
#NP_001261918.1_BIR1_T2_DIAP1     "Drosophila_melanogaster",
#NP_477127.1_BIR3_T2_DIAP2        "Drosophila_melanogaster",
#NP_477127.1_BIR2_T2_DIAP2        "Drosophila_melanogaster",
#NP_477127.1_BIR1_T1_DIAP2        "Drosophila_melanogaster",
#NP_001243092.1_BIR3_T3_cIAP1     "Homo_sapiens",
#NP_001243092.1_BIR2_T2_cIAP1     "Homo_sapiens",
#NP_001243092.1_BIR1_T1_cIAP1     "Homo_sapiens",
#NP_919376.1_BIR3_T2_cIAP1        "Danio_rerio ",
#NP_919376.1_BIR2_T2_cIAP1        "Danio_rerio ",
#NP_919376.1_BIR1_T1_cIAP1        "Danio_rerio ",
#NP_892007.1_BIR3_T2_cIAP2        "Homo_sapiens",
#NP_892007.1_BIR2_T2_cIAP2        "Homo_sapiens",
#NP_892007.1_BIR1_T1_cIAP2        "Homo_sapiens",
#NP_031490.2_BIR3_T2_cIAP2        "Mus_musculus",
#NP_031490.2_BIR2_T2_cIAP2        "Mus_musculus",
#NP_031490.2_BIR1_T1_cIAP2        "Mus_musculus",
#NP_001159.2_BIR1_T2_BIRC5        "Homo_sapiens",
#NP_919378.1_BIR1_T2_BIRC5a       "Danio_rerio",
#NP_057336.3_BIR1_T2_BIRC6        "Homo_sapiens",
#XP_009291311.1_BIR1_T2_BIRC6     "Danio_rerio",
#NP_001158.2_BIR3_T2_XIAP         "Homo_sapiens",
#NP_001158.2_BIR2_T2_XIAP         "Homo_sapiens",
#NP_001158.2_BIR1_T1_XIAP         "Homo_sapiens",
#NP_919377.2_BIR3_T2_XIAP         "Danio_rerio",
#NP_919377.2_BIR2_T2_XIAP         "Danio_rerio",
#NP_919377.2_BIR1_T1_XIAP         "Danio_rerio",
#NP_647478.1_BIR1_T2_BIRC7a       "Homo_sapiens",
#NP_071444.1_BIR1_T2_BIRC7b       "Homo_sapiens",
#AAH39318.1_BIR1_T2_BIRC8         "Homo_sapiens ",

# join this updated species info
BIR_domain_model_MY_CV_CG_type_updated <-  
BIR_domain_model_MY_CV_CG_type_updated %>%
  left_join(model_org_species, by = c("seq.name")) %>%
  mutate(Species = ifelse(is.na(Species.x), Species.y, Species.x)) %>% 
  select(-c(Species.x, Species.y))

# Look at stastics
BIR_domain_model_MY_CV_CG_type_updated_stats <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  group_by(Type, Species) %>% 
  count() 

# fewer NAs 
# Which C. vir and C. gig domains have still not been typed?
#BIR_domain_model_MY_CV_CG_type_updated_not_ID <- BIR_domain_model_MY_CV_CG_type_updated %>% filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>% filter(is.na(Type))
# Remove from multiple alignment to improve visualization
#BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig_MSA <- phylotools::read.fasta(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig_MSA.fa")
#BIR_domain_model_MY_CV_CG_type_updated_not_ID_MSA <- left_join(BIR_domain_model_MY_CV_CG_type_updated_not_ID[c("seq.name","Type")], BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig_MSA)
#dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_not_ID_MSA[,c("seq.name","seq.text")], outfile="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_non_T1_T2_Cvir_Cgig_MSA_notID.fa")
# Viewing now in UGENE

## View RAxML All BIR tree 
BIR_IAP_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_model_prot_IAP_prot_BIR_seq_MSA_RAxML")
BIR_IAP_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
BIR_IAP_raxml_tibble <- as_tibble(BIR_IAP_raxml)

# Add type 1 and type II as found above 
colnames(BIR_domain_model_MY_CV_CG_type_updated)[1] <- "label"
class(BIR_domain_model_MY_CV_CG_type_updated$label)
class(BIR_IAP_raxml_tibble$label)
BIR_IAP_raxml_tibble <- left_join(BIR_IAP_raxml_tibble, BIR_domain_model_MY_CV_CG_type_updated[,c("Species","label","Type")])

# Make a shortened species column
BIR_IAP_raxml_tibble$Species_shortened <- recode(BIR_IAP_raxml_tibble$Species, "Crassostrea_gigas"="C. gigas",
                                                 "Crassostrea_virginica"="C. virginica", "Mizuhopecten_yessoensis"="M. yessoensis",
                                                 "Homo_sapiens"="H. sapiens", "Drosophila_melanogaster" = "D. melanogaster", "Danio_rerio"="D. rerio",
                                                 "Mus_musculus"="M. musculus")
# add color to tibble in order to plot geom_hilight below
BIR_type_color <- data.frame(color= c("#d14e3a", 
                                      "#65c95d", 
                                      "#cb958a", 
                                      "#c4d648", 
                                      "#9944cd", 
                                      "#c6933b", 
                                      "#7095b9", 
                                      "#77cebd", 
                                      "#bcca90", 
                                      "#753a32", 
                                      "#ca96cd", 
                                      "#4a6139", 
                                      "#6257b3", 
                                      "#422d4f", 
                                      "#c6458b"),
                             Type = c("Non_Zinc_binding", 
                                      "T1",
                                      "T1-like_1", 
                                      "T1-like_2",
                                      "T1_like_3",
                                      "T1_like_4",
                                      "T2",
                                      "T2-like_1",
                                      "T2-like_2" ,
                                      "T2_like_3",       
                                      "T2_like_4",
                                      "T2_like_5",
                                      "TX",
                                      "TY",
                                      "Unique"))
BIR_IAP_raxml_tibble_color <- left_join(BIR_IAP_raxml_tibble, BIR_type_color)

# Convert to treedata
BIR_IAP_raxml_treedata <- as.treedata(BIR_IAP_raxml_tibble)

### PLOT BIR TYPE TREE WITH MSA ####
BIR_IAP_raxml_tree <- 
  ggtree(BIR_IAP_raxml_treedata, aes(color=Type, fill=Type), branch.length = "none") + 
  # label with the label
  geom_tiplab(aes(label=label), size = 2.3, offset = 0.1) + 
  # add label with species
  geom_tiplab(aes(label=Species_shortened), size = 2.3, offset = 6.3, color = "black", fontface = "italic") +
     #Edit theme
  theme(legend.position = c(0.25, 0.96), 
        legend.text = element_text(size=10, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  xlim(NA,45) + 
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=2.0) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=2.0) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=2.0) +
  # change position of legend title and spread across columns 
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = ""))) + # need to override aes to get rid of "a" 
  # change colors for species to match other trees 
  scale_colour_manual(name = "Type", values=c("#d14e3a", # "Non_Zinc_binding", 
                                              "#65c95d", # "T1",
                                              "#cb958a", # "T1-like_1", 
                                              "#c4d648", # "T1-like_2",
                                              "#9944cd", # "T1_like_3",
                                              "#c6933b", # "T1_like_4",
                                              "#7095b9", # "T2",
                                              "#77cebd", # "T2-like_1",
                                              "#bcca90", # "T2-like_2" ,
                                              "#753a32", # "T2_like_3",       
                                              "#ca96cd", # "T2_like_4",
                                              "#4a6139", # "T2_like_5",
                                              "#6257b3", # "TX",
                                              "#422d4f", # "TY",
                                              "#c6458b"),# "Unique"
                                              na.value="grey46",
                      breaks=c("Non_Zinc_binding", "T1","T1-like_1", "T1-like_2","T1_like_3","T1_like_4","T2","T2-like_1","T2-like_2" ,"T2_like_3",       
                               "T2_like_4","T2_like_5","TX","TY","Unique"),
                      labels = c("Non Zinc-binding", "Type 1","Type 1-like 1", "Type 1-like 2","Type 1-like 3","Type 1-like 4","Type 2","Type 2-like 1","Type 2-like 2" ,"Type 2-like 3",       
                                 "Type 2-like 4","Type 2-like 5","Type X","Type Y","Uncharacterized/\nModel Organism BIR")) 

# Use loop to clade label highlight to each node 
for(j in 1:dim(BIR_IAP_raxml_tibble_color)[1]){
  #Then add each clade label
  BIR_IAP_raxml_tree <- BIR_IAP_raxml_tree + geom_hilight(node=BIR_IAP_raxml_tibble_color$node[j], offset = .25, fill = BIR_IAP_raxml_tibble_color$color[j], alpha=.2, extend = 6.2)
}

# ggsave BIR tree only
#ggsave(filename = "BIR_tree_type.tiff", plot=BIR_IAP_raxml_tree, device="tiff",
#       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/",
#       width = 15 ,
#       height = 25,
#       units = "in",
#       dpi=300)
#

# save figure that has the correct species names updated on Feb. 4th, 2021
ggsave(filename = "BIR_tree_type_2_4_21.tiff", plot=BIR_IAP_raxml_tree, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/BIR_Type_MSA_figure/",
       width = 15 ,
       height = 25,
       units = "in",
       dpi=300)

### View the BIR tree of IAP domains that didn't match consensus Type 1 and Type II above ###
#BIR_IAP_non_T2_T1_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_domain_model_MY_CV_CG_non_T1_T2_MSA")
#BIR_IAP_non_T2_T1_raxml
#
## Convert to tibble
#BIR_IAP_non_T2_T1_raxml_tibble <- as.tibble(BIR_IAP_non_T2_T1_raxml)
#
## Convert to treedata
#BIR_IAP_non_T2_T1_raxml_treedata <- as.treedata(BIR_IAP_non_T2_T1_raxml_tibble)
#
## Plot tree
#BIR_IAP_non_T2_T1_raxml_tree <- 
#  ggtree(BIR_IAP_non_T2_T1_raxml_treedata, branch.length = "none") + 
#  geom_tiplab(aes(label=label), size = 2.0) + 
#  #Edit theme
#  theme(legend.position = "bottom", 
#        legend.text = element_text(face = "italic", size=6, family="sans"),
#        legend.title = element_text(size=12, family="sans")) +
#  xlim(NA,70) + 
#  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 3.0, fontface="bold") + # allows for subset
#  guides(col = guide_legend(ncol =1, title.position = "top", override.aes = aes(label = "")) ) # need to override aes to get rid of "a"
## View tree with Multiple alignment side by side 
#msaplot(BIR_IAP_non_T2_T1_raxml_tree, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_non_T1_T2_MSA.fa",
# offset=9)

## GENERATE MSA FIGURE  ##
# View MSA with ggmsa because there are more options for viewing sequence color and its easier than trying to make a geom_rect object again
# https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html#visualizing-multiple-sequence-alignments
# ggtree author mentions usage of this package: https://github.com/GuangchuangYu/ggtree-current-protocols, however, my plot grid way makes it 
    # easier to use
# additional protocols for ggmsa https://yulab-smu.github.io/ggmsa/articles/Extensions/extensions.html
# tips on adding layer to specific facet https://guangchuangyu.github.io/2016/12/add-layer-to-specific-panel-of-facet_plot-output/

# Load AAmultiple alignment using biostrings
#BIR_IAP_all_MSA <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA.fa", format = "fasta")
BIR_IAP_all_MSA <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA.fa")

# Reorder the multiple alignment to be the order from the RAxML tree 
# get tip order of tree 
# Get the node order from original GIMAP tree
BIR_IAP_raxml_treedata_tip  <- fortify(BIR_IAP_raxml_treedata )
BIR_IAP_raxml_treedata_tip = subset(BIR_IAP_raxml_treedata_tip, isTip)
BIR_IAP_raxml_treedata_tip_order <- BIR_IAP_raxml_treedata_tip$label[order(BIR_IAP_raxml_treedata_tip$y, decreasing=TRUE)]

# order MSA vector object using tip order 
BIR_IAP_all_MSA <- BIR_IAP_all_MSA[match(rev(BIR_IAP_raxml_treedata_tip_order), BIR_IAP_all_MSA[,1]),] # reverse the order to make it match when plotting 

# export back to then reload in correct order
dat2fasta(BIR_IAP_all_MSA, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA_treeorder.fa")

# reload as AAMultipleAlignment object 
BIR_IAP_all_MSA <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA_treeorder.fa", format = "fasta")

# Plot MSA with ggmsa
#BIR_IAP_all_MSA_treeorder <- ggmsa(BIR_IAP_all_MSA, start = 53, end = 85, 
#      color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
#      none_bg = TRUE, # keeps only the letters and not the full color background
#      posHighligthed = c(57,60, 67, 76, 77,80,82,84), # specify specific positions to highlight in the alignment 
#      seq_name = FALSE) +  # checked that the order is correct so I fixed this, add the sequence name so I can check its plotting in the right order
#     # increase text size
#   theme(text = element_text(size=10),
#         plot.margin = unit(c(0, 0, 0, 0), "cm")) # remove the y axis which just shows the counts of sequences 

# Plot revised MSA with the full BIR seqeuence
BIR_IAP_all_MSA_treeorder_revised <- ggmsa(BIR_IAP_all_MSA, 
                                   color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                   none_bg = TRUE, # keeps only the letters and not the full color background
                                   posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84), # specify specific positions to highlight in the alignment 
                                   seq_name = FALSE) +  # checked that the order is correct so I fixed this, add the sequence name so I can check its plotting in the right order
  # increase text size
  theme(text = element_text(size=10),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) # remove the y axis which just shows the counts of sequences 

## Plotting the tree with the MSA side by side (though it doesn't align well)
# MSA with same axis spacing as the tree
#BIR_tree <- BIR_IAP_raxml_tree +
#aplot::ylim2(BIR_IAP_all_MSA_treeorder)

# plot tree and the alignment using plot_grid 
BIR_tree_type_MSA <- plot_grid(BIR_tree, BIR_IAP_all_MSA_treeorder, ncol=2, align="h", axis ="tb",
                              labels = c("A","B"), rel_heights = c(0.8,1))

plot_grid(BIR_tree, BIR_IAP_all_MSA_treeorder, ncol=2, align="h", axis ="tb",
          labels = c("A","B"), rel_heights = c(0.8,1))

# Remove in between white space in plot in Inkscape
ggsave(filename = "BIR_tree_MSA_plot_2_4_21.tiff", plot=BIR_IAP_all_MSA_treeorder_revised, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/BIR_Type_MSA_figure",
       width = 16 ,
       height = 25,
       units = "in",
       dpi=300)

# final solution was to export the tree and MSA seperately and put them together in inkscape because I kept getting so many errors and the MSA
  # takes and extremely long time to render.


#### PLOT IAP COLLAPSED TREE WITH IAP PROTEIN DOMAINS ####

# Join domain type with label and full domains
colnames(BIR_XP_gff_Interpro_Domains_only_cd00022)[16] <-"Domain_Name"
# Join by start and end to make specific for domain
BIR_XP_gff_Interpro_Domains_only_BIR_type <- left_join(BIR_XP_gff_Interpro_Domains_only, BIR_XP_gff_Interpro_Domains_only_cd00022[,c("protein_id","Domain_Name", "start","end")], by = c("protein_id","start","end"))

# Join in the Type after recoding to reduce the number of levels, keeping T1, T1-like, T2, T2-like T3, T4
BIR_domain_model_MY_CV_CG_type_updated_Type_reduced <- BIR_domain_model_MY_CV_CG_type_updated
BIR_domain_model_MY_CV_CG_type_updated_Type_reduced$Type <- recode(BIR_domain_model_MY_CV_CG_type_updated_Type_reduced$Type,
                                                              "T1-like_1"="T1", 
                                                              "T1-like_2"="T1",
                                                              "T1_like_3"="T1",
                                                              "T1_like_4"="T1",      
                                                              "T2-like_1"="T2",
                                                              "T2-like_2"="T2" ,
                                                              "T2_like_3"="T2",
                                                              "T2_like_4"="T2",
                                                              "T2_like_5"="T2")

colnames(BIR_domain_model_MY_CV_CG_type_updated_Type_reduced)[1] <- "Domain_Name"
BIR_XP_gff_Interpro_Domains_only_BIR_type <- left_join(BIR_XP_gff_Interpro_Domains_only_BIR_type, BIR_domain_model_MY_CV_CG_type_updated_Type_reduced[,c("Domain_Name","Type","Species")])

# Mutate type to fill in NAs with Dbxref, but still keep T2 and T1 for BIR domains
BIR_XP_gff_Interpro_Domains_only_BIR_type <- BIR_XP_gff_Interpro_Domains_only_BIR_type %>%
  mutate(Type = case_when(
    Dbxref != "\"InterPro:IPR001370\"" & is.na(.$Type) ~ as.character(Dbxref), # exclude the non classified BIR domains because these are covering up my other colors 
    TRUE ~ as.character(.$Type))) 
# this leaves some NA lines for type but thats okay they will be removed when plotting 

# Add in generic BIR for proteins where there is now no BIR repeat, meaning that the BIR repeat was no identified with CDD database but with just other interproscan tools
BIR_XP_gff_Interpro_Domains_only_BIR_type %>% filter()
# fixing this later!

# shorten the full protein line by 3300
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened <- BIR_XP_gff_Interpro_Domains_fullprot %>%
  mutate(end = case_when(
    grepl("BIRC6",alias) & end >= 3000 ~ as.numeric(.$end - 3300),
    TRUE ~ as.numeric(end)
  ))

#  Reduce 2 BIR6 introns region, the first large region by 2500 base pairs and the second by 800  
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened <-   BIR_XP_gff_Interpro_Domains_only_BIR_type  %>% 
  # shorten for those domains where the start position is greater than 3000 (these are the only domains affected)
  ## fix start and end in separate mutates 
  mutate(start = case_when(
    grepl("BIRC6",alias) & start >= 3000 ~ as.numeric(.$start - 2500),
    TRUE ~ as.numeric(start)
  )) %>% 
  mutate(end = case_when(
    grepl("BIRC6",alias) & end >= 3000 ~ as.numeric(.$end - 2500),
    TRUE ~ as.numeric(end)
  )) %>% 
  # mutate again to shorten the second region 
  mutate(start = case_when(
    grepl("BIRC6",alias) & start >= 1500 ~ as.numeric(.$start - 800),
    TRUE ~ as.numeric(start)
  )) %>% 
  mutate(end = case_when(
    grepl("BIRC6",alias) & end >= 1500 ~ as.numeric(.$end - 800),
    TRUE ~ as.numeric(end)
  ))

# Add in two sets of railroad tracks 
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened %>% filter(grepl("BIRC6",alias)) %>%
  distinct(protein_id, alias, height_start, height_end) # get the list of protein IDs to add lines for 
BIRC6_railroad <- data.frame(protein_id = c("XP_022331415.1","XP_022331413.1","XP_022331412.1","XP_022331414.1","XP_022331416.1",
                                            "XP_022332917.1","XP_022332915.1","XP_022332914.1","XP_022332916.1","XP_022332918.1",
                                            "XP_019925481.1","XP_019925480.1","XP_011436597.1","XP_019925483.1","XP_019925482.1",
                                            "XP_021363252.1","XP_022331415.1","XP_022331413.1","XP_022331412.1","XP_022331414.1",
                                            "XP_022331416.1","XP_022332917.1","XP_022332915.1","XP_022332914.1","XP_022332916.1",
                                            "XP_022332918.1","XP_019925481.1","XP_019925480.1","XP_011436597.1","XP_019925483.1","XP_019925482.1","XP_021363252.1"),
                             alias = c("BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like",
                                       "BIRC6","BIRC6","BIRC6","BIRC6","BIRC6","BIRC6-like","BIRC6-like","BIRC6-like",
                                       "BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like","BIRC6-like",
                                       "BIRC6","BIRC6","BIRC6","BIRC6","BIRC6","BIRC6-like"),
                            start = c(960 ,960 ,960 ,960 ,960 ,960,960 ,960 ,960 ,960 ,960 ,960 ,960 ,960 ,960 ,960, 1220,1220,1220,1220,1220,1220,
                                      1220,1220,1220,1220,1220,1220,1220,1220,1220,1220), 
                            end = c(962 ,962 ,962 ,962 ,962 ,962,962 ,962 ,962 ,962 ,962 ,962 ,962 ,962 ,962 ,962 ,1223,1223,1223,1223,1223,1223,
                                    1223,1223,1223,1223,1223,1223,1223,1223,1223,1223), 
                            height_start = c(101.75 ,100.75 ,99.75 ,98.75 ,97.75 ,96.75 ,95.75 ,94.75 ,93.75 ,92.75 ,91.75 ,90.75 ,89.75 ,88.75 ,
                                             87.75 ,86.75 ,101.75,100.75,99.75 ,98.75 ,97.75 ,96.75 ,95.75 ,94.75 ,93.75 ,92.75 ,91.75 ,90.75 ,89.75 ,
                                             88.75 ,87.75 ,86.75),
                            height_end = c(102.5,101.5,100.5,99.5,98.5,97.5,96.5,95.5,94.5,93.5,92.5,91.5,90.5,89.5,88.5,87.5,102.5,101.5,100.5,
                                           99.5,98.5,97.5,96.5,95.5,94.5,93.5,92.5,91.5,90.5,89.5,88.5,87.5),
                            Type = c("Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened",
                            "Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened",
                            "Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened",
                            "Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened",
                            "Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened","Intron_shortened",
                            "Intron_shortened","Intron_shortened","Intron_shortened"))

# Add in new rows 
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill <- plyr::rbind.fill(BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened, BIRC6_railroad)

# Set factor level order of the nodes set levels in reverse order (already set)
#BIR_XP_gff_Interpro_Domains_only_BIR_type$node <- factor(BIR_XP_gff_Interpro_Domains_only_BIR_type$node, levels = unique(BIR_XP_gff_Interpro_Domains_only_BIR_type$node))

BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill$Type <- factor(BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill$Type, 
                                                         levels = c( "\"InterPro:IPR022103\"",
                                                                    "Non_Zinc_binding", "T1","T2","TX","TY","Unique",  
                                                                    "Intron_shortened",
                                                                    "\"InterPro:IPR036322\"",
                                                                    "\"InterPro:IPR019775\"",
                                                                    "cd16713",
                                                                    "\"InterPro:IPR013083\"",
                                                                    "\"InterPro:IPR001841\"",
                                                                    "\"InterPro:IPR000608\"",
                                                                    "\"InterPro:IPR016135\"",
                                                                    "\"InterPro:IPR015940\"",
                                                                    "cd14321",
                                                                    "\"InterPro:IPR032171\"",
                                                                    "\"InterPro:IPR027417\"",
                                                                    "G3DSA:1.10.533.10",
                                                                    "\"InterPro:IPR003131\"",
                                                                    "cd18316",
                                                                    "\"InterPro:IPR011333\"",
                                                                    "\"InterPro:IPR000210\"",
                                                                    "\"InterPro:IPR011047\"",
                                                                    "\"InterPro:IPR038765\"",
                                                                    "\"InterPro:IPR001370\""))

### Plot select domains with Domain Information ###
# get rows we want to add star to - these do not have a type designation
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_with_Type <- 
  BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>% 
  filter(!is.na(Type)) %>% 
  distinct(protein_id)
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_no_Type <- BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill[!(BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill$protein_id %in% 
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_with_Type$protein_id), ] %>% distinct(protein_id)
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_no_Type
# Get coordinates for adding star
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_coord <- BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened[BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened$protein_id %in% BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_no_Type$protein_id,] %>%
  select(protein_id, end, node)

## Join domain structure designations and use to create dataframe for plotting shading
# Join with domain structure designations with proteins collapsed to add annotation
IAP_domain_structure <- read_csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_Domain_Structure_groups.csv")
class(BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened$protein_id)
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect <- left_join(BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened, IAP_domain_structure)

# Join with height for each protein domain
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect <- left_join(BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect, 
unique(BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill[,c("protein_id","height_start","height_end")]))

# find xmin and xmax for plot
min(BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened$start) # 1
max(BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened$end) # 1832

# get domain rect coordinates
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect <- BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect %>%
  # group by number
  group_by(Number) %>%
  # keep min height start and max height end within each group
  mutate(min_start_per_group = min(height_start),
            max_end_per_group = max(height_end)) %>% 
  filter(!is.na(Domain_Name)) %>% 
  # keep distinct rows
  distinct(Number, min_start_per_group, max_end_per_group,
           Domain_Name, Color_group, Bold_group) %>%
  ungroup() %>%
  mutate(xmin = 1,
         xmax = 2200) # extend a bit to allow room for plotting text 

# Create geom text dataframe
BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_text <- BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect %>%
  # plot text at midpoint in each group, which would be the mean value
  mutate(y_difference = (max_end_per_group - min_start_per_group)/2,
         y_midpoint =  min_start_per_group + y_difference,
         text_name = 1970,
         text_number = 2400) 

## Create main plot
IAP_Interproscan_domain_plot_BIR_type_domain_subset <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data =BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add stars for all lines where there is not a "Type" listed  - meaning that these have BIR domains found by interproscan but not CDD search
  # do this after the initial geom_segment so that the levels are preserved below for the domains
  geom_point(data = BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened[BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened$protein_id %in% BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_no_Type$protein_id,],
             aes(x=as.numeric(end), y = as.numeric(node)),
             color="black",
             size=6,
             shape= "*") +
  # add protein domain boxes with geom_rect 
  geom_rect(data=BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill,inherit.aes = FALSE,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Type)) +
  #add labels
  labs(y = NULL, x = "Protein Domain Position (aa)") +
  # add theme
  theme_bw() + 
  # plot theme
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=14, family="sans"),
        legend.title = element_text(size=16, family="sans"),
        axis.text.x=element_text(size=18, family="sans"),
        axis.title.x.bottom = element_text(size=18, family="sans")) +
  # Change y axis ticks
  scale_x_continuous(breaks=c(0,500,1000,1500,1800), expand = c(0,0)) + 
  # Change domain labels 
  scale_fill_manual(values=c( 
                             "#5db8de",
                             "#b2dbfe",
                             "#5ba6a6",
                             "#524ed4",
                             "#b6eee2",
                             "#53e1a1",
                             "#d393d4",
                             # "Non_Zinc_binding", "T1","T2","T3","T4","Unique",
                             "black",
                             "#cd6137",
                             "#d27b3d",
                             "#d04d8f",
                             "#89599b",
                             "#c058c6",
                             "#8aaa75",
                             "#50803d",
                             "#a6b348",
                             "#85c967",
                             "#b77853",
                             "#55b793",
                             "#dadd48",
                             "grey", # dummy values
                             "grey",
                             "grey",
                             "grey",
                             "grey",
                             "grey",
                             "grey",
                             "grey"), 
                    name="Functional Domains",
                    breaks=c("\"InterPro:IPR022103\"",
                             "Non_Zinc_binding", "T1","T2","TX","TY","Unique",
                             "Intron_shortened",
                             "\"InterPro:IPR036322\"",
                             "\"InterPro:IPR019775\"",
                             "cd16713",
                             "\"InterPro:IPR013083\"",
                             "\"InterPro:IPR001841\"",
                             "\"InterPro:IPR000608\"",
                             "\"InterPro:IPR016135\"",
                             "\"InterPro:IPR015940\"",
                             "cd14321",
                             "\"InterPro:IPR032171\"",
                             "\"InterPro:IPR027417\"",
                             "G3DSA:1.10.533.10"),
                    labels=c("Baculoviral IAP repeat-containing protein 6",
                             "Non-Zinc-binding BIR", "BIR Type I","BIR Type II","BIR Type X","BIR Type Y","Uncharacterized BIR", # switch to roman numeral
                             "Intron not to scale",
                             "WD40-repeat-containing domain superfamily",
                             "WD40 repeat, conserved site",
                             "RING-HC BIRC2/3/7",
                             "Zinc finger, RING/FYVE/PHD-type",
                             "Zinc finger, RING-type",
                             "Ubiquitin-conjugating enzyme E2",
                             "Ubiquitin-conjugating enzyme/RWD-like",
                             "Ubiquitin-associated domain",
                             "Ubiquitin-association domain of IAPs",
                             "C-terminal of Roc (COR) domain",
                             "P-loop containing nucleoside triphosphate hydrolase",
                             "Death domain, Fas")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top")) 

# add shading domain group boxes and text
IAP_Interproscan_domain_plot_BIR_type_domain_subset_shaded_text <- 
  IAP_Interproscan_domain_plot_BIR_type_domain_subset + 
  geom_rect(data=BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect, inherit.aes=FALSE,
          aes(ymin=min_start_per_group ,ymax=max_end_per_group,xmin=xmin,xmax=xmax,
              group=Color_group), 
          # fill works if you put outside of aes and inlude the colors in the data itseld
              fill = BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_rect$Color_group,
          color = "gray86", # add border color
          size=0.2, # set border line thickness 
          alpha=0.1)  + # make translucent 
  # add text for domain name
 geom_text(data = BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_text, inherit.aes = FALSE,
           aes(x=text_name, y = y_midpoint, label = Domain_Name, family= "sans", fontface = Bold_group),
           size = 7)  
  # add text for Number, DECIDED TO REMOVE THE NUMBERS FROM HERE AND JUST KEEP WITH THE TREE PLOT
  #geom_text(data = BIR_XP_gff_Interpro_Domains_fullprot_BIR6_shortened_domain_text, inherit.aes = FALSE,
  #          aes(x=text_number, y = y_midpoint, label = Number, family= "sans", fontface = Bold_group),
  #          size = 6)  

###  Export and arrange domain plot with tree
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_legend <- cowplot::get_legend(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed)
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_no_legend <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed + 
  theme(legend.position='none')

IAP_MY_CV_CG_tree <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_no_legend + aplot::ylim2(IAP_Interproscan_domain_plot_no_legend)
IAP_Interproscan_domain_plot_no_legend <- IAP_Interproscan_domain_plot_BIR_type_domain_subset_shaded_text + theme(legend.position='none')
IAP_Interproscan_domain_plot_legend <- cowplot::get_legend(IAP_Interproscan_domain_plot_BIR_type_domain_subset_shaded_text)

IAP_tr_dom_collapsed <- plot_grid(NULL,IAP_MY_CV_CG_tree, IAP_Interproscan_domain_plot_no_legend, ncol=3, align='h', rel_widths = c(0.2, 0.7,0.8)) +
  # Add some space at top for labels
  theme(plot.margin = unit(c(1,0.0,0.0,0.0), "cm")) 
IAP_tr_dom_collapsed_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_legend, IAP_Interproscan_domain_plot_legend,
  nrow = 1, align="hv", rel_widths  =c(0.7, 0.7,1)) 

## Create combined figure for publication
IAP_tr_dom_plus_legend <- plot_grid(IAP_tr_dom_collapsed, IAP_tr_dom_collapsed_legend,  ncol=1, rel_heights  = c(0.8, 0.1)) +
  # add labels for plot components
  draw_plot_label(c("A","B","C"), x= c(0.38, 0.53, 0.9), y = c(1,1,1), size = 30, family = "sans")

## Export plot with tree and domains aligned : USE THIS FOR PUBLICATION
#ggsave(filename = "IAP_tr_dom_plus_legend_plot_07302020.tiff", plot=IAP_tr_dom_plus_legend, device="tiff",
#       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
#       width = 34,
#       height = 27,
#       units = "in",
#       dpi=300)

ggsave(filename = "IAP_tr_dom_plus_legend_plot_09172020.tiff", plot=IAP_tr_dom_plus_legend, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
       width = 34,
       height = 27,
       units = "in",
       dpi=300)


#### BIR TYPE MULTIPANEL FIGURE ####

## GOAL: build a multipanel figure showing the sequence of model organism BIR types, the sequence of C. vir and C. gig types, the number of domains, and the number of total proteins those domains are in. 
# this will be the top panel of the figure (A). Panel A will also display above the sequences what their predicted secondary structure is. 
# Panel B will be a table with the number of BIR proteins that have one two or three BIR repeats.

### Statistics for BIR domains ###
# code not correct, edited and fixed on 2/5/21
## The code below is only applicable for C. gigas and C. virginica. In my previous analysis when I collapsed the proteins for M. yessonsis, I did not carry through with analyzing 
  # these proteins. In order to get accurate numbers for M. yessonensis, I will need to go back and add these back in and repeat this analysis. This will be relatively straightforward to
  # do to get the correct number of BIR domains each protein has, but not so straightforward for analyzing new BIR types.
## Also, the BIR domain analysis was only conducted for BIR domains that were confirmed by CDD, which is why the full gene lists don't match the whole gene list

## Make table of statistics for how many genes of each type and variant type for C. gigas and C. virginica
# join table with gene info 
BIR_domain_model_MY_CV_CG_type_updated_gene <- left_join(BIR_domain_model_MY_CV_CG_type_updated, unique(All_molluscs_CDS_gff[,c("gene","protein_id")]), by = "protein_id") 
BIR_domain_model_MY_CV_CG_type_updated_gene %>% distinct(gene, Type, Species) %>% count(Type, Species) %>% 
  # remember the M. yessoensis is not totally accurate because it was collapsed before analysis 
  filter(Species== "Crassostrea_virginica" | Species == "Crassostrea_gigas")

BIR_domain_model_MY_CV_CG_type_updated_gene %>% distinct(gene, Species) %>% count(Species) 
#Species  n
#1       Crassostrea_gigas 35  ..... 5 genes missing from C. gig
#2   Crassostrea_virginica 59  ..... 10 missing from C. vir...why?

## Testing to determine where my gene numbers changed and why
# were genes lost in my calling of types? (check the DF right before the BIRs were called into types)
left_join(BIR_domain_model_MY_CV_CG_type_distinct, BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG) %>% distinct(gene,Species) %>% count(Species) 
    # no seems like they were lost before this 
    # Species  n
    # 1       Crassostrea_gigas 35
    # 2   Crassostrea_virginica 59
    # 3 Mizuhopecten_yessoensis 27
    # 4                    <NA>  1

# are all the genes listed in this original Interproscan DF? YES
left_join(BIR_XP_gff_Interpro_Domains, BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG) %>% distinct(gene,Species) %>% count(Species) 
#Species                     n
#<chr>                   <int>
#1 Crassostrea_gigas          40 # YES
#2 Crassostrea_virginica      68 # yes all but one
#3 Mizuhopecten_yessoensis    34 # yes

# were they already lost here? YES -- figured out because I subset for only those BIRs that were confirmed by CDD
left_join(BIR_XP_gff_Interpro_Domains_only_cd00022,BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG) %>% distinct(gene,Species) %>% count(Species) 
    #Joining, by = "protein_id"
    ## A tibble: 3 x 2
    #Species                     n
    #<chr>                   <int>
    #  1 Crassostrea_gigas          35
    #2 Crassostrea_virginica      59
    #3 Mizuhopecten_yessoensis    27
# were thye lost here
left_join(BIR_XP_gff_Interpro_Domains_only, BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG) %>% distinct(gene,Species) %>% count(Species) 
# Joining, by = "protein_id"
# # A tibble: 3 x 2
# Species                     n
# <chr>                   <int>
#   1 Crassostrea_gigas          40
# 2 Crassostrea_virginica      68
# 3 Mizuhopecten_yessoensis    34

# Number of proteins with a certain number of BIR domains
## The use of this Interpro DF that was created for the purpose of plotting is okay, but remember because Mizuhopecten was collapsed for plotting, the code below is
  # not valid for getting Mizuhopecten statistics 
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>%
  # filter by BIR types that were plotted 
  filter(Type == "\"InterPro:IPR022103\"" | Type == "Non_Zinc_binding" | Type == "T1" |
           Type == "T2" | Type == "TX" | Type =="TY" | Type == "Unique") %>% filter(source == "CDD") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  group_by(protein_id) %>% mutate(total_CDD_BIR_per_protein = n()) %>% 
  # Keep only one row for each protein
  distinct(protein_id, .keep_all = TRUE) %>%
  ungroup() %>% 
  count(total_CDD_BIR_per_protein, Species) %>% filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  arrange(desc(total_CDD_BIR_per_protein, n, Species))

## A tibble: 6 x 3
#total_CDD_BIR_per_protein Species                   n
#<int> <chr>                 <int>
#  1                         3 Crassostrea_gigas         1
#2                         3 Crassostrea_virginica     2
#3                         2 Crassostrea_gigas        26
#4                         2 Crassostrea_virginica    35
#5                         1 Crassostrea_gigas        18
#6                         1 Crassostrea_virginica    47

BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR <- BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>%
  # filter by BIR types that were plotted 
  filter(Type == "\"InterPro:IPR022103\"" | Type == "Non_Zinc_binding" | Type == "T1" |
           Type == "T2" | Type == "TX" | Type =="TY" | Type == "Unique") %>% filter(source == "CDD") %>%
   filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  group_by(protein_id) %>% mutate(total_CDD_BIR_per_protein = n()) %>% 
  # Keep only one row for each protein
  distinct(protein_id, .keep_all = TRUE)

# number of Proteins with each type (summarized for all the Type I or Type I-like or Type II and Type II-like)
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>%
  # filter by BIR types that were plotted 
  filter(Type == "\"InterPro:IPR022103\"" | Type == "Non_Zinc_binding" | Type == "T1" |
           Type == "T2" | Type == "TX" | Type =="TY" | Type == "Unique") %>% filter(source == "CDD") %>%
  filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  distinct(protein_id, Type, .keep_all = TRUE) %>% 
  count(Type, Species)

# Number of genes containing each type (summarized for all the Type I or Type I-like or Type II and Type II-like)
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>%
  # filter by BIR types that were plotted 
  filter(Type == "\"InterPro:IPR022103\"" | Type == "Non_Zinc_binding" | Type == "T1" |
           Type == "T2" | Type == "TX" | Type =="TY" | Type == "Unique") %>% filter(source == "CDD") %>%
  filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  distinct(protein_id, Type, .keep_all = TRUE) %>% 
  left_join(., unique(All_molluscs_CDS_gff[,c("gene","protein_id")]), by = "protein_id") %>%
  ungroup() %>%
  distinct(gene,Type, Species) %>%
  count(Type, Species)

BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_type <- BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill %>%
  # filter by BIR types that were plotted 
  filter(Type == "\"InterPro:IPR022103\"" | Type == "Non_Zinc_binding" | Type == "T1" |
           Type == "T2" | Type == "TX" | Type =="TY" | Type == "Unique") %>% filter(source == "CDD") %>%
  filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  distinct(protein_id, Type, .keep_all = TRUE) %>% 
  left_join(., unique(All_molluscs_CDS_gff[,c("gene","protein_id")]), by = "protein_id") %>%
  ungroup() %>%
  distinct(gene,Type, Species)

# How many genes have one two or three BIRs? Including this in my summary table 
# join with gene info
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR %>% 
  # join with gene info
  left_join(., unique(All_molluscs_CDS_gff[,c("gene","protein_id")]), by = "protein_id") %>%
  ungroup() %>%
  distinct(gene,total_CDD_BIR_per_protein, Species) %>%
   filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  count(Species) # C. virginica is missing 10 and C. gigas is missing 5..why is that? 

BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR_gene <- BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR %>% 
  # join with gene info
  left_join(., unique(BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG[,c("gene","protein_id")]), by = "protein_id") %>%
  ungroup() %>%
  distinct(gene,total_CDD_BIR_per_protein, Species)

BIR_XP_gff_Interpro_Domains_gene_count <- BIR_XP_gff_Interpro_Domains %>%
  # filter by BIR types that were plotted 
  filter(Name =="cd00022")  %>%
  left_join(., unique(BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG[,c("gene","protein_id","Species")]), by = "protein_id") %>% 
  filter( Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas" | Species == "Mizuhopecten_yessoensis") %>%
  # filter out lines that have the exact same Target (meaning had two Interproscan BIR hits)
  distinct(Target, .keep_all = TRUE) %>% 
  group_by(protein_id) %>% mutate(total_CDD_BIR_per_protein = n()) %>% 
  # Keep only one row for each protein
  distinct(protein_id, .keep_all = TRUE) %>%
  ungroup() %>%
  distinct(gene,total_CDD_BIR_per_protein, Species)

BIR_XP_gff_Interpro_Domains_gene_count %>% ungroup() %>% count(Species, total_CDD_BIR_per_protein)

# compare to C. gig and C. vir BIR domain counts from two different starting dataframes using anti_join 
  # anti_join() return all rows from x without a match in y.
anti_join(BIR_XP_gff_Interpro_Domains_gene_count,BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR_gene) %>% View()
anti_join(BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_number_BIR_gene,BIR_XP_gff_Interpro_Domains_gene_count) %>% View()

# Now I can add the full M. yessoensis counts for the BIR numbers with the full gene list! 

### Make tree showing the number of BIR domains in each  
IAP_GENE_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene_MSA_RAxML")
IAP_GENE_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
IAP_GENE_raxml_tibble_BIR <- as_tibble(IAP_GENE_raxml)
colnames(IAP_GENE_raxml_tibble_BIR)[4] <- "gene"
# add species data and change name to include the species data so I can add a highlight for BIR number and keep the text black
#  CV_CG_MY_gene_species
IAP_GENE_raxml_tibble_BIR_join <- left_join(IAP_GENE_raxml_tibble_BIR, BIR_XP_gff_Interpro_Domains_gene_count) %>% mutate(label = case_when(
  Species == "Crassostrea_gigas" ~ paste("C.gig",gene, sep = "-"),
  Species == "Crassostrea_virginica" ~ paste("C.vir",gene, sep = "-"))) 

# join in to fill in the rest of the genes that were not analyzed for BIR
IAP_GENE_raxml_tibble_BIR_join <- left_join(IAP_GENE_raxml_tibble_BIR_join, unique(BIR_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG[,c("gene","Species")]), by = "gene") %>%
  select(-Species.x) %>% rename(Species = Species.y) %>%  mutate(label = case_when(
    Species == "Crassostrea_gigas" ~ paste("C.gig",gene, sep = "-"),
    Species == "Crassostrea_virginica" ~ paste("C.vir",gene, sep = "-"),
    Species == "Mizuhopecten_yessoensis" ~ paste("M.yes",gene, sep = "-"))) %>% 
    # fill in NA total_CDD_BIR_per_protein_label  
  mutate(total_CDD_BIR_per_protein = case_when(
    !is.na(total_CDD_BIR_per_protein) ~ as.character(total_CDD_BIR_per_protein),
    is.na(total_CDD_BIR_per_protein) & !is.na(gene) ~ "non_CDD",
    !is.na(gene) & !is.na(bootstrap) ~ NA_character_,
    gene == "LOC111099688" ~ "2")) %>%
  rename(species_label = label) %>% rename(label = gene) %>% 
# Fix the three genes that are still blank, C. vir LOC111099688 recoded to  LOC111105597..the other genes are M.yessoensis 
 mutate(species_label = case_when(
  label == "LOC111099688" ~ "C_vir-LOC111099688",
    label == "LOC110443976" ~ "M.yes-LOC110443976",
    label == "LOC110440802" ~ "M.yes-LOC110440802",
  TRUE ~ as.character(species_label))) %>%
  mutate(Species = case_when(
    label == "LOC111099688" ~ "Crassostrea_virginica",
    label == "LOC110443976" ~ "Mizuhopecten_yessoensis",
    label == "LOC110440802" ~ "Mizuhopecten_yessoensis",
    TRUE ~ as.character(Species))) 

# Join the IAPs that are Type X,Type Y, NZBIR
BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_type_new <- BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_type %>%
  filter(Type == "TX" | Type == "Non_Zinc_binding" | Type == "TY") %>% rename(label = gene)
IAP_GENE_raxml_tibble_BIR_join <- left_join(IAP_GENE_raxml_tibble_BIR_join, BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_type_new)

#Convert to treedata object to store tree plus outside data
IAP_GENE_raxml_BIR_treedata <- as.treedata(IAP_GENE_raxml_tibble_BIR_join)

# make dataframe with color added 
# add color to tibble in order to plot geom_hilight below
BIR_type_color_gene <- data.frame(color= c("#d14e3a",
                                           "#65c95d",
                                           "#cb958a",
                                           "#c4d648",
                                           "#9944cd"),
                             total_CDD_BIR_per_protein = c("1","2","3","non_CDD", "4"))
IAP_GENE_raxml_tibble_BIR_join_color <- left_join(IAP_GENE_raxml_tibble_BIR_join, BIR_type_color_gene)

## Plot IAP gene sequence tree as circular 
# basic plot to get color legend
IAP_GENE_BIR_number_circular_gene_legend <- ggtree(IAP_GENE_raxml_BIR_treedata, layout="circular", 
                                                   aes(color=total_CDD_BIR_per_protein), 
                                                   branch.length = "none") + 
  scale_colour_manual(name = "Species", values=c("#d14e3a",
                                                 "#65c95d",
                                                 "#cb958a",
                                                 "#c4d648",
                                                 "#9944cd"), na.value="grey46", breaks=c("1", "2","3","non_CDD","4"),
                      labels = c("1 BIR Domain", "2 BIR Domains","3 BIR Domains","Non-CDD BIR","4 BIR Domains")) 
IAP_GENE_BIR_number_circular_gene_legend <- get_legend(IAP_GENE_BIR_number_circular_gene_legend)

# add geom_highlight for BIR number
IAP_GENE_BIR_number_circular_gene <- 
  ggtree(IAP_GENE_raxml_BIR_treedata, layout="circular", 
         #aes(color=Type), 
         branch.length = "none") + 
  geom_tiplab2(aes(label=species_label,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  #xlim(-100,100)  +
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=2.0) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=2.0) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=2.0) +
                  # fix legend appearance
  guides(col = guide_legend(ncol =3, title.position = "top", override.aes = aes(label = "")))  # need to override aes to get rid of "a"

# Use loop to clade label highlight to each node 
for(j in 1:dim(IAP_GENE_raxml_tibble_BIR_join_color)[1]){
  #Then add each clade label
  IAP_GENE_BIR_number_circular_gene <- IAP_GENE_BIR_number_circular_gene + geom_hilight(node=IAP_GENE_raxml_tibble_BIR_join_color$node[j], offset = .25, fill = IAP_GENE_raxml_tibble_BIR_join_color$color[j], extend = 12.0)
}

# Put together plot and other legend in cowplot
IAP_GENE_BIR_number_circular_gene_plus_legend <- plot_grid(IAP_GENE_BIR_number_circular_gene, IAP_GENE_BIR_number_circular_gene_legend,  ncol=1, rel_heights  = c(0.62, 0.1)) 

# Export plot to file 
ggsave(filename = "IAP_GENE_BIR_number_circular_gene.tiff", plot=IAP_GENE_BIR_number_circular_gene_plus_legend, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/BIR_Type_MSA_figure/",
       width = 10 ,
       height = 10,
       units = "in",
       dpi=300)

## MSA of model organism sequences ###
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens <- BIR_domain_model_MY_CV_CG_type_updated %>% filter(Species == "Homo_sapiens" | Species == "Drosophila_melanogaster") 

# only keep a few of each type 
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens

H_sapiens_keep <- data.frame(label = c("NP_001243092.1_BIR1_T1_cIAP1","NP_892007.1_BIR1_T1_cIAP2","NP_001158.2_BIR1_T1_XIAP","NP_001261918.1_BIR2_T2_DIAP1","NP_477127.1_BIR3_T2_DIAP2",
                                       "NP_001243092.1_BIR2_T2_cIAP1","NP_001159.2_BIR1_T2_BIRC5","NP_001158.2_BIR2_T2_XIAP","AAH39318.1_BIR1_T2_BIRC8"))
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens <- BIR_domain_model_MY_CV_CG_type_updated_H_sapiens[BIR_domain_model_MY_CV_CG_type_updated_H_sapiens$label %in% H_sapiens_keep$label,] %>% arrange(Type)

# Subset MSA
BIR_IAP_all_MSA_subset <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA.fa")
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_H_sapiens$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset$seq.name
#[1] "AAH39318.1_BIR1_T2_BIRC8"     H sapiens
#"NP_001158.2_BIR2_T2_XIAP"     H sapiens
#"NP_001159.2_BIR1_T2_BIRC5"   H. sapiens
#"NP_001243092.1_BIR2_T2_cIAP1" H. sapiens
#"NP_477127.1_BIR3_T2_DIAP2"   D. melanogaster
#[6] "NP_001261918.1_BIR2_T2_DIAP1" D. melanogaster
#"NP_001158.2_BIR1_T1_XIAP"   H. sapiens  
#"NP_892007.1_BIR1_T1_cIAP2"   H. sapiens
#"NP_001243092.1_BIR1_T1_cIAP1" H. sapiens

BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "NP_001243092.1_BIR1_T1_cIAP1" ~ "H. sapiens cIAP1 BIR1 Type I",
    seq.name ==   "NP_892007.1_BIR1_T1_cIAP2"    ~ "H. sapiens cIAP2 BIR1 Type I",
    seq.name ==   "NP_001158.2_BIR1_T1_XIAP"     ~  "H. sapiens XIAP BIR1 Type I",
    seq.name ==  "NP_001261918.1_BIR2_T2_DIAP1" ~ "D. melanogaster DIAP1 BIR2 Type II",
    seq.name ==  "NP_477127.1_BIR3_T2_DIAP2"    ~ "D. melanogaster DIAP2 BIR3 Type II",
    seq.name ==  "NP_001243092.1_BIR2_T2_cIAP1" ~ "H. sapiens cIAP1 BIR2 Type II",
    seq.name ==  "NP_001158.2_BIR2_T2_XIAP"  ~ "H. sapiens XIAP BIR2 Type II",
    seq.name ==  "NP_001159.2_BIR1_T2_BIRC5"  ~ "H. sapiens BIRC5 BIR1 Type II",
    seq.name ==  "AAH39318.1_BIR1_T2_BIRC8" ~ "H. sapiens BIRC8 BIR1 Type II"))

# put in reverse order
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename[c(6,7,8,1,2,3,5,4),]

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename, 
                                                                                color = "Shapely_AA", 
                                                                                none_bg = TRUE, # keeps only the letters and not the full color background
                                                                                posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84), # specify specific positions to highlight in the alignment 
                                                                                seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.5,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct 

### plot Type I C_vir C_gig
BIR_domain_model_MY_CV_CG_type_updated_Type_1 <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  # filter out just the type 1 and the type 1 variants
  filter(grepl("T1", Type)) %>%
  # keep C. vir and C. gigas and  M yessoensis
  filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # only keep one distinct example from each type and species, do this by just keeping the top row in each group randomly for an example
  group_by(Type, Species) %>%
  arrange(Type) %>%
  filter(row_number() ==1)

# Subset MSA
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_Type_1$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset$seq.name
#"XP_011423762.1_BIR1_LOC105325768" 
#"XP_021360744.1_BIR1_BIRC8-like"   
#"XP_022287996.1_BIR1_BIRC7-like"   
#"XP_011444161.1_BIR1_LOC105340029" 
#"XP_022288684.1_BIR1_LOC111100858"
#"XP_021363602.1_BIR1_BIRC3-like"
#"XP_019919899.1_BIR1_BIRC7"        
#"XP_022288097.1_BIR1_BIRC3-like" 

BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "XP_011423762.1_BIR1_LOC105325768" ~   "C. gigas Type I",
    #seq.name == "XP_021360744.1_BIR1_BIRC8-like"   ~  "M. yessoensis Type I",
    seq.name == "XP_022287996.1_BIR1_BIRC7-like"   ~  "C. virginica Type I",
    seq.name == "XP_011444161.1_BIR1_LOC105340029" ~  "C. gigas Type I-like 3",
    seq.name == "XP_022288684.1_BIR1_LOC111100858" ~  "C. virginica Type I-like 3",
    seq.name == "XP_021363602.1_BIR1_BIRC3-like" ~  "M. yessoensis Type I-like 4",
    seq.name == "XP_019919899.1_BIR1_BIRC7"        ~  "C. gigas Type I-like 1",
    seq.name == "XP_022288097.1_BIR1_BIRC3-like"  ~  "C. virginica Type I-like 2"))

# put in reverse and correct order
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename[c(3,4,6,5,1,2),] 

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename, 
                                                                         color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                                                         none_bg = TRUE, # keeps only the letters and not the full color background
                                                                         posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84), # specify specific positions to highlight in the alignment 
                                                                         seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct 

### plot Type II C_vir C_gig
BIR_domain_model_MY_CV_CG_type_updated_Type_2 <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  # filter out just the type 2 and the type 2 variants
  filter(grepl("T2", Type)) %>%
  # keep C. vir and C. gigas and  M yessoensis
  filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas" ) %>%
  # only keep one distinct example from each type and species, do this by just keeping the top row in each group randomly for an example
  group_by(Type, Species) %>%
  arrange(Type) %>%
  filter(row_number() ==1)

# Subset MSA
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_Type_2$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset$seq.name
#"XP_011412926.1_BIR2_BIRC7A"  
#"XP_021341279.1_BIR3_BIRC3-like"  
#"XP_022287929.1_BIR1_LOC111100400" 
#"XP_011416423.1_BIR1_BIRC7B"      
#"XP_022287912.1_BIR2_LOC111100394"
#"XP_022292343.1_BIR1_BIRC7A-like"  
#"XP_021350197.1_BIR2_BIRC3-like"   
#"XP_011427116.1_BIR1_BIRC7A"      
#"XP_022291629.1_BIR1_BIRC2-like"   
#"XP_011431980.1_BIR2_BIRC7"       
#"XP_022287996.1_BIR2_BIRC7-like"   

BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "XP_011412926.1_BIR2_BIRC7A"  ~ "C. gigas Type II",
    #seq.name == "XP_021341279.1_BIR3_BIRC3-like"  ~ "M. yessoensis Type II",
    seq.name == "XP_022287929.1_BIR1_LOC111100400"~ "C. virginica Type II",
    seq.name == "XP_011416423.1_BIR1_BIRC7B"      ~ "C. gigas Type II-like 3",
    seq.name == "XP_022287912.1_BIR2_LOC111100394"~ "C. virginica Type II-like 3",
    seq.name == "XP_022292343.1_BIR1_BIRC7A-like" ~ "C. virginica Type II-like 4",
    #seq.name == "XP_021350197.1_BIR2_BIRC3-like"  ~ "M. yessoensis Type II-like 5",
    seq.name == "XP_011427116.1_BIR1_BIRC7A"      ~ "C. gigas BIR1 Type II-like 1",
    seq.name == "XP_022291629.1_BIR1_BIRC2-like"  ~ "C. virginica Type II-like 1",
    seq.name == "XP_011431980.1_BIR2_BIRC7"       ~ "C. gigas BIR2 Type II-like 2",
    seq.name == "XP_022287996.1_BIR2_BIRC7-like"  ~ "C. virginica Type II-like 2"))

# put in reverse and correct order
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename[c(5,3,4,8,9,6,7,1,2),] 

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename, 
                                                                         color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                                                         none_bg = TRUE, # keeps only the letters and not the full color background
                                                                         posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84),  # specify specific positions to highlight in the alignment 
                                                                         seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct  

### plot Type X 
BIR_domain_model_MY_CV_CG_type_updated_Type_X <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  filter(grepl("TX", Type)) %>%
  # keep C. vir and C. gigas and  M yessoensis
  filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas" ) %>%
  # only keep one distinct example from each type and species, do this by just keeping the top row in each group randomly for an example
  group_by(Type, Species) %>%
  arrange(Type) %>%
  filter(row_number() ==1)

# Subset MSA
BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_Type_X$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset$seq.name
#"XP_022287912.1_BIR1_LOC111100394"  

BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "XP_022287912.1_BIR1_LOC111100394"  ~ "C. virginica Type X"))

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename, 
                                                                         color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                                                         none_bg = TRUE, # keeps only the letters and not the full color background
                                                                         posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84),  # specify specific positions to highlight in the alignment 
                                                                         seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct  

### plot Type Y 
BIR_domain_model_MY_CV_CG_type_updated_Type_Y <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  filter(grepl("TY", Type)) %>%
  # keep C. vir and C. gigas and  M yessoensis
  filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas" ) %>%
  # only keep one distinct example from each type and species, do this by just keeping the top row in each group randomly for an example
  group_by(Type, Species) %>%
  arrange(Type) %>%
  filter(row_number() ==1)

# Subset MSA
BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_Type_Y$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset$seq.name
#"XP_011437445.1_BIR1_BIRC2"       
#"XP_022287934.1_BIR1_LOC111100402" 

BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "XP_011437445.1_BIR1_BIRC2"  ~ "C. gigas Type Y",
    seq.name == "XP_022287934.1_BIR1_LOC111100402" ~ "C. virginica Type Y"))

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename, 
                                                                         color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                                                         none_bg = TRUE, # keeps only the letters and not the full color background
                                                                         posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84),  # specify specific positions to highlight in the alignment 
                                                                         seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct 

### plot NZBIR 
BIR_domain_model_MY_CV_CG_type_updated_Type_NZ <- BIR_domain_model_MY_CV_CG_type_updated %>% 
  filter(grepl("Non", Type)) %>%
  # keep C. vir and C. gigas and  M yessoensis
  filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas") %>%
  # only keep one distinct example from each type and species, do this by just keeping the top row in each group randomly for an example
  group_by(Type, Species) %>%
  arrange(Type) %>%
  filter(row_number() ==1)

# Subset MSA
BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset <- BIR_IAP_all_MSA_subset[match(BIR_domain_model_MY_CV_CG_type_updated_Type_NZ$label, BIR_IAP_all_MSA_subset[,1]),]

# edit seq name for figure
BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset$seq.name
#"XP_019922709.1_BIR1_DIAP2"      
#"XP_022291916.1_BIR1_BIRC3-like" 

BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename <- BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset %>%
  mutate(seq.name = case_when(
    seq.name == "XP_019922709.1_BIR1_DIAP2"  ~ "C. gigas Non Zinc Binding",
    seq.name == "XP_022291916.1_BIR1_BIRC3-like" ~ "C. virginica Non Zinc Binding"))

# export back to then reload in correct order
dat2fasta(BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename.fa")

# reload as AAMultipleAlignment object 
BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename_msa <- ggmsa(BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename,  
                                                                         color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
                                                                         none_bg = TRUE, # keeps only the letters and not the full color background
                                                                         posHighligthed = c(2,30,31,33,34,48,53,55,57,58, 60, 67, 77,80,82,84),  # specify specific positions to highlight in the alignment 
                                                                         seq_name = TRUE) +
  # increase text size
  theme(text = element_text(size=12),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct 

## Plot all MSAs on a grid
# arrange with ggarrange
MSA_arrange <- ggarrange(BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset_rename_msa,
                         BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset_rename_msa,
                         BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset_rename_msa,
                         BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset_rename_msa,
                         BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset_rename_msa,
                         BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset_rename_msa,
                         ncol = 1)

# Export the table to save
ggsave(filename = "BIR_MSA_by_type.tiff", plot=MSA_arrange, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_full_tree_MSA",
       width = 13 ,
       height = 6,
       units = "in",
       dpi=300)
# save version with the M. yessoensis sequences removed 
ggsave(filename = "BIR_MSA_by_type_2_8_21.tiff", plot=MSA_arrange, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_full_tree_MSA",
       width = 13 ,
       height = 6,
       units = "in",
       dpi=300)

## Get example sequences to run in secondary protein prediction software RaptorX
BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset$seq.name
BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset$seq.name
BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset$seq.name
BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset$seq.name
BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset$seq.name
BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset$seq.name

BIR_example_list <- rbind(BIR_domain_model_MY_CV_CG_type_updated_H_sapiens_MSA_subset,
                          BIR_domain_model_MY_CV_CG_type_updated_T1_MSA_subset,
                          BIR_domain_model_MY_CV_CG_type_updated_T2_MSA_subset,
                          BIR_domain_model_MY_CV_CG_type_updated_TY_MSA_subset,
                          BIR_domain_model_MY_CV_CG_type_updated_NZ_MSA_subset,
                          BIR_domain_model_MY_CV_CG_type_updated_TX_MSA_subset)

# get the non-MSA version of these sequences 
BIR_example_list_sequences <- BIR_domain_model_MY_CV_CG_type_distinct[BIR_domain_model_MY_CV_CG_type_distinct$seq.name %in% BIR_example_list$seq.name,]

# export as fasta 
dat2fasta(BIR_example_list_sequences[,c(1,2)], "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_example_list_sequences.fa")


#### PLOT IAP GENE TREE WITH INTRONLESS GENE AND GENE DENSITY INFO ####

# export gene features for use in plotting and bedtools 
C_vir_rtracklayer_gene_bed <- C_vir_rtracklayer %>% filter(type == "gene") %>% dplyr::select(seqid, start, end, gene) 

# Write table for use in any external software 
write.table(C_vir_rtracklayer_gene_bed, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_rtracklayer_gene.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE)

C_vir_rtracklayer %>% filter(type == "exon") %>% distinct(start, end, .keep_all = TRUE) %>% count(gene)

# export chromosome lengths for use in plotting and bedtools 
C_vir_rtracklayer_chromosome_bed <- C_vir_rtracklayer %>% filter(type == "region") %>% dplyr::select(seqid, start, end)

# Write table for use in any external software 
write.table(C_vir_rtracklayer_chromosome_bed, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_rtracklayer_chromosome.bed",
            quote = FALSE,col.names = FALSE, row.names=FALSE)

## Import gene density information that was calculated in bedops 
C_vir_rtracklayer_gene_100kb_density <- read.table("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_rtracklayer_gene_100kb_density.bed",
                                                   col.names = c("seqid","start","end"))
C_vir_rtracklayer_gene_1Mb_density <- read.table("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_rtracklayer_gene_1Mb_density.bed",
                                                   col.names = c("seqid","start","end"))

# separate the counts into its own column for histogram plotting 
C_vir_rtracklayer_gene_100kb_density_sep <- C_vir_rtracklayer_gene_100kb_density %>% separate(end , into = c("stop","count"))
C_vir_rtracklayer_gene_1Mb_density_sep <- C_vir_rtracklayer_gene_1Mb_density %>% separate(end , into = c("stop","count"))

# Identify intronless genes in  C. virginica
C_vir_rtracklayer_intronless <- C_vir_rtracklayer %>% filter(type == "exon") %>% distinct(start, end, .keep_all = TRUE) %>% count(gene) %>% filter(n == 1)

C_vir_rtracklayer_intronless_IAP <- C_vir_rtracklayer_intronless[C_vir_rtracklayer_intronless$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,]

# gene n
# 245  LOC111105137 1
# 463  LOC111109152 1
# 853  LOC111116378 1
# 898  LOC111116826 1
# 906  LOC111117137 1
# 941  LOC111117856 1
# 1569 LOC111132301 1
# 1591 LOC111132489 1

# see how  these overlap with BIR types
left_join(C_vir_rtracklayer_intronless_IAP, BIR_XP_gff_Interpro_Domains_only_BIR_type_BIR6_shortened_fill_type)

# gene n Type               Species
# 1 LOC111105137 1 <NA>                  <NA> # non-cdd
# 2 LOC111109152 1   T2 Crassostrea_virginica
# 3 LOC111116378 1   T2 Crassostrea_virginica
# 4 LOC111116826 1   T2 Crassostrea_virginica
# 5 LOC111117137 1   T2 Crassostrea_virginica
# 6 LOC111117856 1   T2 Crassostrea_virginica
# 7 LOC111132301 1 <NA>                  <NA>  # non-cdd
#   8 LOC111132489 1 <NA>                  <NA> # non-cdd

## get coordinates for intronless genes 
C_vir_rtracklayer_intronless_IAP_coord <- C_vir_rtracklayer_gene_bed[C_vir_rtracklayer_gene_bed$gene %in% C_vir_rtracklayer_intronless_IAP$gene,]

# identify in C. gig
C_gig_rtracklayer_intronless <- C_gig_rtracklayer %>% filter(type == "exon") %>% distinct(start, end, .keep_all = TRUE) %>% count(gene) %>% filter(n == 1)
C_gig_rtracklayer_intronless_IAP <- C_gig_rtracklayer_intronless[C_gig_rtracklayer_intronless$gene %in% BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY$gene,]

# only three in C. gigas 
  #gene n
  #771  LOC105338159 1
  #1012 LOC105345723 1
  #1187 LOC109617982 1

# change chromosome names to be all character
C_vir_chromsome_char <- data.frame(seqid = c(   "NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1","NC_035786.1","NC_035787.1","NC_035788.1","NC_035789.1","NC_007175.2"),
                                   Chromosome = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","MT"))

## identify in all mollusc Gff
All_mollusc_exon_gff_intronless <- All_mollusc_exon_gff %>% distinct(start, end, .keep_all = TRUE) %>% count(gene) %>% filter(n == 1)
All_mollusc_exon_gff_intronless_locus <- All_mollusc_exon_gff %>% distinct(start, end, .keep_all = TRUE) %>% count(locus_tag) %>% filter(n == 1)
All_mollusc_exon_gff_intronless_IAP <- All_mollusc_exon_gff_intronless[All_mollusc_exon_gff_intronless$gene %in% All_gene_species$gene,]
All_mollusc_exon_gff_intronless_IAP_locus <- All_mollusc_exon_gff_intronless_locus[All_mollusc_exon_gff_intronless_locus$locus_tag %in% All_gene_species$gene,]


# match with species
left_join(All_mollusc_exon_gff_intronless_IAP, All_gene_species[,c("gene","Species")])
All_mollusc_exon_gff_intronless_IAP_locus %>% rename(gene = locus_tag) %>% left_join(., All_gene_species[,c("gene","Species")])

#Joining, by = "gene"
#gene n           Species
#1       EGW08_017511 1 Elysia_chlorotica
#2  LOTGIDRAFT_111196 1   Lottia_gigantea
#3  LOTGIDRAFT_113794 1   Lottia_gigantea
#4  LOTGIDRAFT_119876 1   Lottia_gigantea
#5  LOTGIDRAFT_125617 1   Lottia_gigantea
#6  LOTGIDRAFT_139361 1   Lottia_gigantea
#7   LOTGIDRAFT_59237 1   Lottia_gigantea
#8   LOTGIDRAFT_59297 1   Lottia_gigantea
#9   LOTGIDRAFT_59316 1   Lottia_gigantea
#10  LOTGIDRAFT_59317 1   Lottia_gigantea
#11  LOTGIDRAFT_59566 1   Lottia_gigantea
#12  LOTGIDRAFT_69569 1   Lottia_gigantea
#13  LOTGIDRAFT_69570 1   Lottia_gigantea


# Joining, by = "gene"
# gene n                 Species
# 1  LOC105338159 1       Crassostrea_gigas
# 2  LOC105345723 1       Crassostrea_gigas
# 3  LOC106075513 1   Biomphalaria_glabrata
# 4  LOC106077001 1   Biomphalaria_glabrata
# 5  LOC109617982 1       Crassostrea_gigas
# 6  LOC110462844 1 Mizuhopecten_yessoensis
# 7  LOC111105137 1   Crassostrea_virginica
# 8  LOC111109152 1   Crassostrea_virginica
# 9  LOC111116378 1   Crassostrea_virginica
# 10 LOC111116826 1   Crassostrea_virginica
# 11 LOC111117137 1   Crassostrea_virginica
# 12 LOC111117856 1   Crassostrea_virginica
# 13 LOC111132301 1   Crassostrea_virginica
# 14 LOC111132489 1   Crassostrea_virginica
# 15 LOC115223103 1        Octopus_vulgaris

## get coodinates for all Cvir IAPs
  # already calculated this on line 1166
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label <- BIR_XP_gff_species_join_haplotig_collapsed_CV_gene %>% 
  left_join(.,C_vir_chromsome_char) %>%
  dplyr::select(Chromosome, start, end,gene) %>%
  rename(ChromStart = start, ChromEnd = end, Gene = gene)

# turn the gene and chr info into factors and start and end into numeric
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$Gene <- as.factor(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$Gene)
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$Chromosome <- as.factor(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$Chromosome)
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$ChromStart <- as.numeric(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$ChromStart)
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$ChromEnd <- as.numeric(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label$ChromEnd)

# How many genes on chromosomes 6 and 7
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label_count <-  BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label
BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label_count$Chromosome <- as.character(BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label_count$Chromosome)

BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label_count %>% filter(Chromosome == "chr6" | Chromosome == "chr7") %>% dplyr::group_by(Chromosome) %>%
 dplyr::count()
    # Chromosome     n
    # <chr>      <int>
    #   1 chr6          27
    # 2 chr7          27

### Generate plot using RCircos 

# change ideogram labels
C_vir_rtracklayer_chromosome_bed_label <- C_vir_rtracklayer_chromosome_bed %>% 
  left_join(.,C_vir_chromsome_char) %>%
  dplyr::select(Chromosome, start, end) %>% 
  rename(ChromStart = start, ChromStart = start, ChromEnd = end) %>%
  # add band and stain columns 
  mutate(Band = "all", Stain = "gvar")

C_vir_rtracklayer_chromosome_bed_label$Chromosome <- as.factor(C_vir_rtracklayer_chromosome_bed_label$Chromosome)
class(C_vir_rtracklayer_chromosome_bed_label$ChromStart)
C_vir_rtracklayer_chromosome_bed_label$ChromStart <- as.numeric(C_vir_rtracklayer_chromosome_bed_label$ChromStart)
C_vir_rtracklayer_chromosome_bed_label$ChromEnd <- as.numeric(C_vir_rtracklayer_chromosome_bed_label$ChromEnd)

# turn on graphics device
out.file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/Cvir_rtracklayer_RCircos11.pdf";
pdf(file=out.file, height = 7, width = 7);
RCircos.Set.Plot.Area();

# Initialize RCircos
RCircos.Set.Core.Components(cyto.info = C_vir_rtracklayer_chromosome_bed_label, 
                            # add three tracks inside: gene density, IAP location, intronless IAPs
                            tracks.inside = 3, tracks.outside = 0) ;

# modifying plot parameters
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 3000;
rcircos.params$text.size <- 0.2
#rcircos.params$chr.name.pos <- 1
rcircos.params$char.width <- 50
rcircos.params$line.color <- "black"
rcircos.params$track.background <- "white"
rcircos.params$hist.color <- "blue"
rcircos.params$grid.line.color <- "black"
rcircos.params$chrom.width <- 0.05
#rcircos.params$highlight.width <- 0.1
RCircos.Reset.Plot.Parameters(rcircos.params);

# plot ideogram 
RCircos.Chromosome.Ideogram.Plot();

# Add histogram of gene density 
C_vir_rtracklayer_gene_100kb_density_sep_label <- C_vir_rtracklayer_gene_100kb_density_sep %>%
  left_join(.,C_vir_chromsome_char) %>%
  dplyr::select(Chromosome, start, stop,count) %>%
  rename(ChromStart = start, ChromEnd = stop, Data = count)
C_vir_rtracklayer_gene_100kb_density_sep_label$Chromosome <- as.factor(C_vir_rtracklayer_gene_100kb_density_sep_label$Chromosome)
C_vir_rtracklayer_gene_100kb_density_sep_label$ChromStart <- as.numeric(C_vir_rtracklayer_gene_100kb_density_sep_label$ChromStart)
C_vir_rtracklayer_gene_100kb_density_sep_label$ChromEnd  <- as.numeric(C_vir_rtracklayer_gene_100kb_density_sep_label$ChromEnd )
C_vir_rtracklayer_gene_100kb_density_sep_label$Data <- as.numeric(C_vir_rtracklayer_gene_100kb_density_sep_label$Data)

C_vir_rtracklayer_gene_1Mb_density_sep_label <- C_vir_rtracklayer_gene_1Mb_density_sep %>%
  left_join(.,C_vir_chromsome_char) %>%
  dplyr::select(Chromosome, start, stop,count) %>%
  rename(ChromStart = start, ChromEnd = stop, Data = count)
C_vir_rtracklayer_gene_1Mb_density_sep_label$Chromosome <- as.factor(C_vir_rtracklayer_gene_1Mb_density_sep_label$Chromosome)
C_vir_rtracklayer_gene_1Mb_density_sep_label$ChromStart <- as.numeric(C_vir_rtracklayer_gene_1Mb_density_sep_label$ChromStart)
C_vir_rtracklayer_gene_1Mb_density_sep_label$ChromEnd  <- as.numeric(C_vir_rtracklayer_gene_1Mb_density_sep_label$ChromEnd )
C_vir_rtracklayer_gene_1Mb_density_sep_label$Data <- as.numeric(C_vir_rtracklayer_gene_1Mb_density_sep_label$Data)

# add IAP gene labels 
gene.data <- BIR_XP_gff_species_join_haplotig_collapsed_CV_gene_label
RCircos.Gene.Connector.Plot(gene.data, track.num = 2, side = "in", is.sorted = FALSE, genomic.columns = 3);
#RCircos.Gene.Name.Plot(gene.data,track.num = 2, side = "out", name.col = 4, is.sorted = FALSE);
RCircos.Histogram.Plot(C_vir_rtracklayer_gene_1Mb_density_sep_label, data.col = 4, track.num = 1, side = "in");

# add label for intronless genes
C_vir_rtracklayer_intronless_IAP_coord_label <- C_vir_rtracklayer_intronless_IAP_coord %>% left_join(.,C_vir_chromsome_char) %>%
  dplyr::select(Chromosome, start, end,gene) %>%
  rename(ChromStart = start, ChromEnd = end) %>% 
  # change text to be star for intronless genes
  mutate(gene_label = "I")
RCircos.Gene.Name.Plot(C_vir_rtracklayer_intronless_IAP_coord_label,track.num = 3, side = "in", name.col = 4, is.sorted = FALSE);

# remember  to turn off at the end
dev.off()

### Plot IAP gene positions using new C_gigas assembly that is at the chromosome level - March 12th, 2021

# Load the new assembly  GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff as rtracklayer file, and make bed file
C_gigas_chr_assembly <- rtracklayer::readGFF("/Users/erinroberts/Downloads/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff")
C_gigas_chr_assembly <- as.data.frame(C_gigas_chr_assembly)

# export chromosome lengths for use in plotting and bedtools 
C_gigas_chr_assembly_chromosome_bed <- C_gigas_chr_assembly %>% filter(type == "region" & genome == "chromosome") %>% dplyr::select(seqid, start, end)  

## Generate plot using RCircos 

# change ideogram labels
C_gigas_chr_assembly_chromosome_bed_label <- C_gigas_chr_assembly_chromosome_bed %>% 
  rename(Chromosome = seqid, ChromStart = start, ChromEnd = end) %>%
  # add band and stain columns 
  mutate(Band = "all", Stain = "gvar")

C_gigas_chr_assembly_chromosome_bed_label$Chromosome <- as.factor(C_gigas_chr_assembly_chromosome_bed_label$Chromosome)
class(C_gigas_chr_assembly_chromosome_bed_label$ChromStart)
C_gigas_chr_assembly_chromosome_bed_label$ChromStart <- as.numeric(C_gigas_chr_assembly_chromosome_bed_label$ChromStart)
C_gigas_chr_assembly_chromosome_bed_label$ChromEnd <-   as.numeric(C_gigas_chr_assembly_chromosome_bed_label$ChromEnd)

## get coodinates for all Cvir IAPs
# already calculated this on line 1166
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label <- BIR_XP_gff_species_join_haplotig_collapsed_CG_gene %>% 
  rename(Chromosome = seqid, ChromStart = start, ChromEnd = end, Gene = gene)

# turn the gene and chr info into factors and start and end into numeric
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$Gene <- as.factor(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$Gene)
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$Chromosome <- as.factor(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$Chromosome)
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$ChromStart <- as.numeric(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$ChromStart)
BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$ChromEnd <- as.numeric(BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label$ChromEnd)

# turn on graphics device
out.file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/Cgig_rtracklayer_RCircos.pdf";
pdf(file=out.file, height = 7, width = 7);
RCircos.Set.Plot.Area();

# Initialize RCircos
RCircos.Set.Core.Components(cyto.info = C_gigas_chr_assembly_chromosome_bed_label, 
                            # add three tracks inside: gene density, IAP location, intronless IAPs
                            tracks.inside = 1, tracks.outside = 0) ;

# modifying plot parameters
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 3000;
rcircos.params$text.size <- 0.2
#rcircos.params$chr.name.pos <- 1
rcircos.params$char.width <- 50
rcircos.params$line.color <- "black"
rcircos.params$track.background <- "white"
rcircos.params$hist.color <- "blue"
rcircos.params$grid.line.color <- "black"
rcircos.params$chrom.width <- 0.05
#rcircos.params$highlight.width <- 0.1
RCircos.Reset.Plot.Parameters(rcircos.params);

# plot ideogram 
RCircos.Chromosome.Ideogram.Plot();

# add IAP gene labels 
gene.data <- BIR_XP_gff_species_join_haplotig_collapsed_CG_gene_label
RCircos.Gene.Connector.Plot(gene.data, track.num = 2, side = "in", is.sorted = FALSE, genomic.columns = 3);
#RCircos.Gene.Name.Plot(gene.data,track.num = 2, side = "out", name.col = 4, is.sorted = FALSE);

dev.off()

## Can't make the plot because these gene locations are not on the actual chromosomes..

#### PLOT IAP TREE WITH DESEQ2 INFORMATION ####
# for plotting here, proteins with identical sequence are collaped for the purpose of plotting 
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_vir_apop_LFC_IAP.Rdata")
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_gig_apop_LFC_IAP.Rdata")
# remember this is the IAPs that were significantly different with challenge, not expression of all IAPs!

C_vir_apop_LFC_IAP_OG <- C_vir_apop_LFC_IAP
C_gig_apop_LFC_IAP_OG <- C_gig_apop_LFC_IAP

# Keep tables separate for plotting 
C_vir_apop_LFC_IAP$Species <- "Crassostrea_virginica"
C_gig_apop_LFC_IAP$Species <- "Crassostrea_gigas"

# Full Join the list of transcripts so the missing labels are in there and can be ordered for plotting 
C_gig_apop_LFC_IAP_full_XP <- full_join(C_gig_apop_LFC_IAP, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])
# Many NA's meaning that many of those proteins were ones that were collapsed as duplicates in CD-Hit
# For NAs meaning that those were exactly identical to another protein and were collapsed. Need to replace these NA's with their uncollapsed parent protein
C_gig_apop_LFC_IAP_full_XP_collapsed <- C_gig_apop_LFC_IAP_full_XP %>% filter(is.na(node) & !is.na(log2FoldChange))  %>% 
  distinct(protein_id, .keep_all = TRUE)
C_gig_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6 <- BIR_seq_rm_dup_clstr6 %>% filter(protein_id %in% (C_gig_apop_LFC_IAP_full_XP_collapsed$protein_id))
# Find parent proteins in these clusters
C_gig_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6_cluster <- BIR_seq_rm_dup_clstr6[BIR_seq_rm_dup_clstr6$cluster %in% C_gig_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6$cluster,]

# Recode these proteins which were collapsed for the purpose of plotting
C_gig_apop_LFC_IAP$protein_id <- recode(C_gig_apop_LFC_IAP$protein_id, 
                                        "XP_011445380.1" = "XP_011445382.1" ,
                                        "XP_011436808.1" = "XP_011436809.1" ,
                                        "XP_011445383.1" = "XP_011445382.1" ,
                                        "XP_019925515.1" = "XP_019925513.1" ,
                                        "XP_019925516.1" = "XP_019925513.1" ,
                                        "XP_019925512.1" = "XP_019925513.1" ,
                                        "XP_011418792.1" = "XP_011418791.1" ,
                                        "XP_019925514.1" = "XP_019925513.1" ,
                                        "XP_011445381.1" = "XP_011445382.1",
                                        "XP_011437419.1" = "XP_011437418.1" ,
                                        "XP_011428386.1" = "XP_011428384.1")
# Rejoin the Full list of transcript and check for fixed NA
C_gig_apop_LFC_IAP_full_XP <- full_join(C_gig_apop_LFC_IAP, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

C_vir_apop_LFC_IAP_full_XP <- full_join(C_vir_apop_LFC_IAP, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])
# For NAs for GIMAP proteins, meaning that those were exactly identical to another protein and were collapsed. Need to replace these NA's with their uncollapsed parent protein
C_vir_apop_LFC_IAP_full_XP_collapsed <- C_vir_apop_LFC_IAP_full_XP %>% filter(is.na(node) & !is.na(log2FoldChange)) %>% 
  distinct(protein_id, .keep_all = TRUE)
C_vir_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6 <- BIR_seq_rm_dup_clstr6 %>% filter(protein_id %in% (C_vir_apop_LFC_IAP_full_XP_collapsed$protein_id))
# Find parent proteins in these clusters
C_vir_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6_cluster <- BIR_seq_rm_dup_clstr6[BIR_seq_rm_dup_clstr6$cluster %in% C_vir_apop_LFC_IAP_full_XP_collapsed_BIR_seq_rm_dup_clstr6$cluster,]


# Recode these proteins for the purpose of plotting
C_vir_apop_LFC_IAP$protein_id <- recode(C_vir_apop_LFC_IAP$protein_id, 
                                       "XP_022287921.1" = "XP_022287919.1",
                                       "XP_022288976.1" = "XP_022288977.1",
                                       "XP_022292112.1" = "XP_022292108.1",
                                       "XP_022295524.1" = "XP_022295527.1",
                                       "XP_022288687.1" = "XP_022288684.1",
                                       "XP_022287932.1" = "XP_022287934.1",
                                       "XP_022291031.1" = "XP_022291030.1",
                                       "XP_022292110.1" = "XP_022292108.1",
                                       "XP_022292109.1"= "XP_022292108.1",
                                       "XP_022292969.1" = "XP_022292970.1",
                                       "XP_022288681.1" = "XP_022288684.1",
                                       "XP_022293781.1" = "XP_022293782.1",
                                       "XP_022288683.1" = "XP_022288684.1",
                                       "XP_022286792.1" = "XP_022295668.1",
                                       "XP_022292414.1" = "XP_022292412.1",
                                       "XP_022291628.1" = "XP_022291629.1",
                                       "XP_022289978.1" = "XP_022289977.1")
                                        
# Rejoin the Full list of transcript and check for fixed NA
C_vir_apop_LFC_IAP_full_XP <- full_join(C_vir_apop_LFC_IAP, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

# Reorder both to be the order of the GIMAP tree XPs
# Get the node order from original GIMAP tree (done in code chunk regarding domain information above)
# Reorder the proteins
C_vir_apop_LFC_IAP_full_XP <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, C_vir_apop_LFC_IAP_full_XP)
C_gig_apop_LFC_IAP_full_XP <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, C_gig_apop_LFC_IAP_full_XP)

# Add in grouping designations for domain structure groups 

# Set factor level order of the nodes set levels in reverse order
C_vir_apop_LFC_IAP_full_XP$protein_id <- factor(C_vir_apop_LFC_IAP_full_XP$protein_id, levels = rev(unique(C_vir_apop_LFC_IAP_full_XP$protein_id)))
C_gig_apop_LFC_IAP_full_XP$protein_id <- factor(C_gig_apop_LFC_IAP_full_XP$protein_id, levels = rev(unique(C_gig_apop_LFC_IAP_full_XP$protein_id)))
C_vir_apop_LFC_IAP_full_XP$node <- factor(C_vir_apop_LFC_IAP_full_XP$node, levels = rev(unique(C_vir_apop_LFC_IAP_full_XP$node)))
C_gig_apop_LFC_IAP_full_XP$node <- factor(C_gig_apop_LFC_IAP_full_XP$node, levels = rev(unique(C_gig_apop_LFC_IAP_full_XP$node)))

# Create geom_rect object for LFC and const. plots - need to have the nodes in the correct order first
# Join in Domain annotation document (with proteins collapsed) in order to get the node number for geom_strip
IAP_domain_structure <- read_csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_Domain_Structure_groups.csv")
IAP_domain_structure_label <- IAP_domain_structure
colnames(IAP_domain_structure_label)[1] <- "label" # change the name
IAP_domain_structure_node <-  left_join(IAP_MY_CV_CG_raxml_tibble[,c("label","node")], IAP_domain_structure_label)
IAP_domain_structure_node <- IAP_domain_structure_node[!is.na(IAP_domain_structure_node$label),]
# reorder the nodes based on tree to correctly call nodes for each group
IAP_domain_structure_node_order <- IAP_domain_structure_node[match(IAP_MY_CV_CG_raxml_treedata_tip_order$protein_id, IAP_domain_structure_node$label),] %>%
  mutate(order = as.numeric(row.names(IAP_domain_structure_node_order)) )
# rename as protien_id
colnames(IAP_domain_structure_node_order)[1] <- "protein_id"

# get coordinates
  IAP_domain_structure_node_order_coordinates <-   IAP_domain_structure_node_order %>%
    # put order in reverse and this will be y coordinates 
    mutate(ycord =  as.numeric(rev(order))) %>%
    group_by(as.character(Number), Color_group)%>%
    # get coordinate labels for each 
    summarize(ymax = max(ycord),
      ymin = min(ycord)) %>%
  # collapse the rows with summarize each
    # subtract half from bottom and add half to top 
    mutate(ymax = as.numeric(ymax) + 0.5,
           ymin = as.numeric(ymin) - 0.5)

## Plot C. virginica
# set factor levels 
C_vir_apop_LFC_IAP_full_XP$experiment <- factor(C_vir_apop_LFC_IAP_full_XP$experiment, levels=c("Hatchery_Probiotic_RI",  "Lab_Pro_RE22","ROD", "Dermo", "NA"), 
                                                labels= c("CVBAC-B",  "CVBAC-A","CVBAC-C", "CVPMA", "NA"))
C_vir_apop_LFC_IAP_full_XP$group_by_sim <- factor(C_vir_apop_LFC_IAP_full_XP$group_by_sim, levels=c("Hatchery_Probiotic_RI" ,"Lab_RI_6hr" , "Lab_RI_RI_24hr", "Lab_S4_6hr","Lab_S4_24hr", "Lab_RE22" ,
                          "ROD_susceptible_seed","ROD_resistant_seed", "Dermo_Susceptible_36hr", "Dermo_Susceptible_7d", "Dermo_Susceptible_28d","Dermo_Tolerant_36hr",   
                          "Dermo_Tolerant_7d","Dermo_Tolerant_28d" ),
                          labels= c("RI" ,"RI 6h", "RI 24h", "S4 6h","S4 24h", "RE22" ,"S. ROD","R. ROD", "S. 36h", "S. 7d", "S. 28d","T. 36h",   
                      "T. 7d","T. 28d"))

C_vir_apop_LFC_IAP_tile_plot_COLLAPSED <- ggplot(C_vir_apop_LFC_IAP_full_XP[!(is.na(C_vir_apop_LFC_IAP_full_XP$experiment)),]) + 
  geom_tile(aes(x=group_by_sim, y = protein_id, fill=log2FoldChange))  + 
  # make translucent 
  #scale_fill_viridis_c(breaks = seq(min(C_vir_apop_LFC_GIMAP_full_XP$log2FoldChange, na.rm = TRUE),max(C_vir_apop_LFC_GIMAP_full_XP$log2FoldChange, na.rm=TRUE),length.out = 15), 
  #                  option="magma", guide=guide_legend()) +
  scale_fill_viridis_c(name = "Log2 Fold Change", 
                       limits = c(-11,10),
                       breaks = c(-10,-5,-1,0.5,1,3,5,7,10), 
                       option="plasma",
                       guide=guide_legend(), na.value = "transparent") +
  facet_grid(.~experiment, scales="free",space="free", drop=TRUE) + 
  labs(
    #x="Treatment", 
    x= NULL,
    y = NULL,
       title = "*C. virginica* Experiment"
    ) +
  theme_minimal() + 
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size=20, family="sans"),
        plot.title = ggtext::element_markdown(size=20, family="sans", hjust= 0.5),
        axis.text.x.bottom = element_text(size=16, family="sans", face = "bold", angle = 90, vjust=0.5, hjust =1),
        legend.position = "bottom",
        legend.title = element_text(size=16, family="sans"), 
        legend.text = element_text(size=14, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(size=0.2, color="gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(
          size = 16, face = "bold", family="sans"),
        strip.background = element_rect(colour = "white", fill="white")) +
  guides(fill=guide_legend(ncol=3, title.position="top")) +
  # remove NA row from the list (get list from rev(unique(C_vir_apop_LFC_IAP_full_XP$protein_id)))
  scale_y_discrete(limits=c(  "XP_011431980.1", "XP_022287996.1", "XP_022287997.1", "XP_021341279.1", "XP_021361235.1", "XP_011428387.1", "XP_011428384.1", "XP_022288977.1", "XP_019920103.1", "XP_021350197.1",
                              "XP_011416423.1", "XP_022295668.1", "XP_022291030.1", "XP_022287919.1", "XP_022287917.1", "XP_022287913.1", "XP_022287914.1", "XP_022287912.1", "XP_019919109.1", "XP_019919110.1", "XP_022292970.1",
                              "XP_022294409.1", "XP_022294410.1", "XP_021369350.1", "XP_011452303.1", "XP_022312791.1", "XP_022311552.1", "XP_011441452.1", "XP_022300945.1", "XP_022311926.1", "XP_022311074.1", "XP_021356418.1",
                              "XP_022293362.1", "XP_011455007.1", "XP_021355384.1", "XP_022287287.1", "XP_022287292.1", "XP_019927576.1", "XP_011449366.1", "XP_022320574.1", "XP_022320313.1", "XP_021362910.1", "XP_021372726.1",
                              "XP_021365179.1", "XP_022287983.1", "XP_022287991.1", "XP_021355347.1", "XP_021345058.1", "XP_021360744.1", "XP_022293034.1", "XP_011418793.1", "XP_011418791.1", "XP_021365162.1", "XP_021365133.1",
                              "XP_021345059.1", "XP_021365104.1", "XP_021365115.1", "XP_021366455.1", "XP_022305768.1", "XP_021365097.1", "XP_021353464.1", "XP_021347771.1", "XP_011427116.1", "XP_022291629.1", "XP_022293782.1",
                              "XP_022295018.1", "XP_022335805.1", "XP_022336007.1", "XP_021355270.1", "XP_022293068.1", "XP_022292343.1", "XP_022292342.1", "XP_022293067.1", "XP_019922713.1", "XP_019922711.1", "XP_019922715.1",
                              "XP_019920663.1", "XP_019922714.1", "XP_019922712.1", "XP_021344072.1", "XP_021369394.1", "XP_011435625.1", "XP_022315256.1", "XP_021363603.1", "XP_021363602.1", "XP_011435383.1", "XP_021363252.1",
                              "XP_019925482.1", "XP_019925483.1", "XP_011436597.1", "XP_019925480.1", "XP_019925481.1", "XP_022332918.1", "XP_022332916.1", "XP_022332914.1", "XP_022332915.1", "XP_022332917.1", "XP_022331416.1",
                              "XP_022331414.1", "XP_022331412.1", "XP_022331413.1", "XP_022331415.1", "XP_021345008.1", "XP_021348428.1", "XP_011455592.1", "XP_022291483.1", "XP_011452445.1", "XP_011412926.1", "XP_022322280.1",
                              "XP_022322279.1", "XP_022322278.1", "XP_022337137.1", "XP_021366760.1", "XP_011436809.1", "XP_011445382.1", "XP_022292415.1", "XP_022292412.1", "XP_022291345.1", "XP_019929019.1", "XP_011433586.1",
                              "XP_019921695.1", "XP_019921696.1", "XP_011423762.1", "XP_011423764.1", "XP_011427115.1", "XP_022291924.1", "XP_022293780.1", "XP_022295543.1", "XP_022290906.1", "XP_021360358.1", "XP_011433457.1",
                              "XP_022292108.1", "XP_022288420.1", "XP_021349549.1", "XP_021344945.1", "XP_022292341.1", "XP_022291916.1", "XP_019922709.1", "XP_019922710.1", "XP_021362783.1", "XP_011444161.1", "XP_022289969.1",
                              "XP_022288685.1", "XP_022288684.1", "XP_022288686.1", "XP_021372320.1", "XP_022288032.1", "XP_011437420.1", "XP_022288056.1", "XP_022287930.1", "XP_022287929.1", "XP_011437418.1", "XP_019925758.1",
                              "XP_019926829.1", "XP_011437445.1", "XP_022287977.1", "XP_022287934.1", "XP_022287975.1", "XP_022287971.1", "XP_022287973.1", "XP_022287974.1", "XP_011413590.1", "XP_022295036.1", "XP_022293846.1",
                              "XP_022295029.1", "XP_022293847.1", "XP_022295527.1", "XP_022295528.1", "XP_022295529.1", "XP_019925513.1", "XP_022288097.1", "XP_022289977.1", "XP_022288105.1", "XP_022288100.1", "XP_022288104.1",
                              "XP_019924100.1", "XP_019919899.1", "XP_011421261.1", "XP_022290206.1", "XP_022288934.1", "XP_022288750.1", "XP_022288931.1", "XP_022288746.1", "XP_022288748.1"),
                   labels=c( "XP_011431980.1", "XP_022287996.1", "XP_022287997.1", "XP_021341279.1", "XP_021361235.1", "XP_011428387.1", "XP_011428384.1", "XP_022288977.1", "XP_019920103.1", "XP_021350197.1",
                             "XP_011416423.1", "XP_022295668.1", "XP_022291030.1", "XP_022287919.1", "XP_022287917.1", "XP_022287913.1", "XP_022287914.1", "XP_022287912.1", "XP_019919109.1", "XP_019919110.1", "XP_022292970.1",
                             "XP_022294409.1", "XP_022294410.1", "XP_021369350.1", "XP_011452303.1", "XP_022312791.1", "XP_022311552.1", "XP_011441452.1", "XP_022300945.1", "XP_022311926.1", "XP_022311074.1", "XP_021356418.1",
                             "XP_022293362.1", "XP_011455007.1", "XP_021355384.1", "XP_022287287.1", "XP_022287292.1", "XP_019927576.1", "XP_011449366.1", "XP_022320574.1", "XP_022320313.1", "XP_021362910.1", "XP_021372726.1",
                             "XP_021365179.1", "XP_022287983.1", "XP_022287991.1", "XP_021355347.1", "XP_021345058.1", "XP_021360744.1", "XP_022293034.1", "XP_011418793.1", "XP_011418791.1", "XP_021365162.1", "XP_021365133.1",
                             "XP_021345059.1", "XP_021365104.1", "XP_021365115.1", "XP_021366455.1", "XP_022305768.1", "XP_021365097.1", "XP_021353464.1", "XP_021347771.1", "XP_011427116.1", "XP_022291629.1", "XP_022293782.1",
                             "XP_022295018.1", "XP_022335805.1", "XP_022336007.1", "XP_021355270.1", "XP_022293068.1", "XP_022292343.1", "XP_022292342.1", "XP_022293067.1", "XP_019922713.1", "XP_019922711.1", "XP_019922715.1",
                             "XP_019920663.1", "XP_019922714.1", "XP_019922712.1", "XP_021344072.1", "XP_021369394.1", "XP_011435625.1", "XP_022315256.1", "XP_021363603.1", "XP_021363602.1", "XP_011435383.1", "XP_021363252.1",
                             "XP_019925482.1", "XP_019925483.1", "XP_011436597.1", "XP_019925480.1", "XP_019925481.1", "XP_022332918.1", "XP_022332916.1", "XP_022332914.1", "XP_022332915.1", "XP_022332917.1", "XP_022331416.1",
                             "XP_022331414.1", "XP_022331412.1", "XP_022331413.1", "XP_022331415.1", "XP_021345008.1", "XP_021348428.1", "XP_011455592.1", "XP_022291483.1", "XP_011452445.1", "XP_011412926.1", "XP_022322280.1",
                             "XP_022322279.1", "XP_022322278.1", "XP_022337137.1", "XP_021366760.1", "XP_011436809.1", "XP_011445382.1", "XP_022292415.1", "XP_022292412.1", "XP_022291345.1", "XP_019929019.1", "XP_011433586.1",
                             "XP_019921695.1", "XP_019921696.1", "XP_011423762.1", "XP_011423764.1", "XP_011427115.1", "XP_022291924.1", "XP_022293780.1", "XP_022295543.1", "XP_022290906.1", "XP_021360358.1", "XP_011433457.1",
                             "XP_022292108.1", "XP_022288420.1", "XP_021349549.1", "XP_021344945.1", "XP_022292341.1", "XP_022291916.1", "XP_019922709.1", "XP_019922710.1", "XP_021362783.1", "XP_011444161.1", "XP_022289969.1",
                             "XP_022288685.1", "XP_022288684.1", "XP_022288686.1", "XP_021372320.1", "XP_022288032.1", "XP_011437420.1", "XP_022288056.1", "XP_022287930.1", "XP_022287929.1", "XP_011437418.1", "XP_019925758.1",
                             "XP_019926829.1", "XP_011437445.1", "XP_022287977.1", "XP_022287934.1", "XP_022287975.1", "XP_022287971.1", "XP_022287973.1", "XP_022287974.1", "XP_011413590.1", "XP_022295036.1", "XP_022293846.1",
                             "XP_022295029.1", "XP_022293847.1", "XP_022295527.1", "XP_022295528.1", "XP_022295529.1", "XP_019925513.1", "XP_022288097.1", "XP_022289977.1", "XP_022288105.1", "XP_022288100.1", "XP_022288104.1",
                             "XP_019924100.1", "XP_019919899.1", "XP_011421261.1", "XP_022290206.1", "XP_022288934.1", "XP_022288750.1", "XP_022288931.1", "XP_022288746.1", "XP_022288748.1")) +
  guides(fill=guide_legend(ncol=4, title.position="top"))  

  # add geom rect boxes 
C_vir_apop_LFC_IAP_tile_plot_COLLAPSED <- C_vir_apop_LFC_IAP_tile_plot_COLLAPSED + 
  # color odd numbered groups as blue and even groups 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "gray86",], inherit.aes = FALSE, 
            aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
            # fill works if you put outside of aes and inlude the colors in the data itself # add border color
            fill = "gray86",
            color = "gray86", # add border color
            size=0.2, # set border line thickness 
            alpha=0.1) +  # make translucent 
# color odd numbered groups as blue and even groups 
geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "slategray1",], inherit.aes = FALSE, 
          aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
          # fill works if you put outside of aes and inlude the colors in the data itself # add border color
          color = "gray86", # add border color
          fill = "slategray1",
          size=0.2, # set border line thickness 
          alpha=0.1)   # make translucent 

# check number of rows matches with tree 
Cvir_LFC_collapsed <- C_vir_apop_LFC_IAP_full_XP[!(is.na(C_vir_apop_LFC_IAP_full_XP$experiment)),]
Cvir_LFC_collapsed %>% distinct(protein_id) %>% count() # 29 to plot 

rev(unique(C_vir_apop_LFC_IAP_full_XP$protein_id)) # 184 levels correct

## Repeat procedure for C. gigas
# Change factor level order of experiment for ggplot 
C_gig_apop_LFC_IAP_full_XP$experiment <-  factor(C_gig_apop_LFC_IAP_full_XP$experiment, levels = c("Zhang", "Rubio","He","deLorgeril_sus", "deLorgeril_res","NA"), 
       labels = c("CVBAC-A" ,"CVBAC-B" , "CGOSHV1-B" ,"CGOSHV1-A Sus.", "CGOSHV1-A Res.","CGOSHV1-A Res."))
        # make sure to use Rmarkdown newline notation since using element_markdown for strip text formatting
C_gig_apop_LFC_IAP_full_XP$group_by_sim <- factor(C_gig_apop_LFC_IAP_full_XP$group_by_sim, levels= c("Zhang_Valg"          ,"Zhang_Vtub"          ,"Zhang_LPS"          , "Rubio_J2_8"          ,"Rubio_J2_9"          ,"Rubio_LGP32"         ,"Rubio_LMG20012T"     ,"He_6hr"             ,
           "He_12hr"             ,"He_24hr"             ,"He_48hr"            , "He_120hr"            ,"deLorgeril_res_6hr"  ,"deLorgeril_res_12hr" ,"deLorgeril_res_24hr" ,"deLorgeril_res_48hr",
           "deLorgeril_res_60hr" ,"deLorgeril_res_72hr" ,"deLorgeril_sus_6hr" , "deLorgeril_sus_12hr" ,"deLorgeril_sus_24hr" ,"deLorgeril_sus_48hr" ,"deLorgeril_sus_60hr" ,"deLorgeril_sus_72hr"), 
labels= c("*V. alg*","*V.tub, V. ang*","LPS *M. Lut*", "*V. cras* NVir","*V. cras* Vir" ,"*V. tas* Vir","*V. tas* NVir","6hr",
          "12hr", "24hr", "48hr", "120hr","6hr","12hr","24hr" ,"48hr",
          "60hr", "72hr" ,"6hr",  "12hr","24hr","48hr" ,
          "60hr", "72hr"))

C_gig_apop_LFC_IAP_tile_plot_COLLAPSED <- ggplot(C_gig_apop_LFC_IAP_full_XP[!(is.na(C_gig_apop_LFC_IAP_full_XP$experiment)),], aes(x=group_by_sim, y=protein_id, fill=log2FoldChange)) + 
  geom_tile() + 
  #scale_fill_viridis_c(breaks = seq(min(C_gig_apop_LFC_GIMAP$log2FoldChange, na.rm = TRUE),max(C_gig_apop_LFC_GIMAP$log2FoldChange, na.rm=TRUE),length.out = 15), 
  #                     option="plasma", guide=guide_legend()) +
  #scale_fill_gradientn(colors = viridis_pal, limits=c(-10, 10),breaks = c(-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10),
  #                     guide=guide_legend(), na.value = "transparent") + 
  scale_fill_viridis_c(name = "Log2 Fold Change", 
                       limits = c(-11,10),
                       breaks = c(-10,-8,-6,-4,-2,-1,0,0.5,1,2,3,4,6,8,10), 
                       option="plasma",
                       guide=guide_legend(), na.value = "transparent") +
  # add facet to plot 
  facet_grid(.~experiment, scales="free",space="free", drop=TRUE) + 
  labs(
    #x="Treatment", 
    x = NULL,
    y = NULL, 
      title = "*C. gigas* Experiment"
  ) +
  theme_minimal() + 
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size=20, family="sans"),
        plot.title = ggtext::element_markdown(size=20, family="sans", hjust= 0.5),
        axis.text.x.bottom = ggtext::element_markdown(size=16, family="sans", face = "bold", angle = 90, vjust=0.5, hjust =1),
        legend.position = "bottom",
        legend.title = element_text(size=16, family="sans"), 
        legend.text = element_text(size=14, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(size=0.2, color="gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = ggtext::element_markdown(
          size = 16, face = "bold", family="sans"),
        strip.background = element_rect(colour = "white", fill="white")) +
  guides(fill=guide_legend(ncol=3, title.position="top")) +
  # remove NA row from the list 
  scale_y_discrete(limits=c("XP_011431980.1", "XP_022287996.1", "XP_022287997.1", "XP_021341279.1", "XP_021361235.1", "XP_011428387.1", "XP_011428384.1", "XP_022288977.1", "XP_019920103.1", "XP_021350197.1",
                            "XP_011416423.1", "XP_022295668.1", "XP_022291030.1", "XP_022287919.1", "XP_022287917.1", "XP_022287913.1", "XP_022287914.1", "XP_022287912.1", "XP_019919109.1", "XP_019919110.1", "XP_022292970.1",
                            "XP_022294409.1", "XP_022294410.1", "XP_021369350.1", "XP_011452303.1", "XP_022312791.1", "XP_022311552.1", "XP_011441452.1", "XP_022300945.1", "XP_022311926.1", "XP_022311074.1", "XP_021356418.1",
                            "XP_022293362.1", "XP_011455007.1", "XP_021355384.1", "XP_022287287.1", "XP_022287292.1", "XP_019927576.1", "XP_011449366.1", "XP_022320574.1", "XP_022320313.1", "XP_021362910.1", "XP_021372726.1",
                            "XP_021365179.1", "XP_022287983.1", "XP_022287991.1", "XP_021355347.1", "XP_021345058.1", "XP_021360744.1", "XP_022293034.1", "XP_011418793.1", "XP_011418791.1", "XP_021365162.1", "XP_021365133.1",
                            "XP_021345059.1", "XP_021365104.1", "XP_021365115.1", "XP_021366455.1", "XP_022305768.1", "XP_021365097.1", "XP_021353464.1", "XP_021347771.1", "XP_011427116.1", "XP_022291629.1", "XP_022293782.1",
                            "XP_022295018.1", "XP_022335805.1", "XP_022336007.1", "XP_021355270.1", "XP_022293068.1", "XP_022292343.1", "XP_022292342.1", "XP_022293067.1", "XP_019922713.1", "XP_019922711.1", "XP_019922715.1",
                            "XP_019920663.1", "XP_019922714.1", "XP_019922712.1", "XP_021344072.1", "XP_021369394.1", "XP_011435625.1", "XP_022315256.1", "XP_021363603.1", "XP_021363602.1", "XP_011435383.1", "XP_021363252.1",
                            "XP_019925482.1", "XP_019925483.1", "XP_011436597.1", "XP_019925480.1", "XP_019925481.1", "XP_022332918.1", "XP_022332916.1", "XP_022332914.1", "XP_022332915.1", "XP_022332917.1", "XP_022331416.1",
                            "XP_022331414.1", "XP_022331412.1", "XP_022331413.1", "XP_022331415.1", "XP_021345008.1", "XP_021348428.1", "XP_011455592.1", "XP_022291483.1", "XP_011452445.1", "XP_011412926.1", "XP_022322280.1",
                            "XP_022322279.1", "XP_022322278.1", "XP_022337137.1", "XP_021366760.1", "XP_011436809.1", "XP_011445382.1", "XP_022292415.1", "XP_022292412.1", "XP_022291345.1", "XP_019929019.1", "XP_011433586.1",
                            "XP_019921695.1", "XP_019921696.1", "XP_011423762.1", "XP_011423764.1", "XP_011427115.1", "XP_022291924.1", "XP_022293780.1", "XP_022295543.1", "XP_022290906.1", "XP_021360358.1", "XP_011433457.1",
                            "XP_022292108.1", "XP_022288420.1", "XP_021349549.1", "XP_021344945.1", "XP_022292341.1", "XP_022291916.1", "XP_019922709.1", "XP_019922710.1", "XP_021362783.1", "XP_011444161.1", "XP_022289969.1",
                            "XP_022288685.1", "XP_022288684.1", "XP_022288686.1", "XP_021372320.1", "XP_022288032.1", "XP_011437420.1", "XP_022288056.1", "XP_022287930.1", "XP_022287929.1", "XP_011437418.1", "XP_019925758.1",
                            "XP_019926829.1", "XP_011437445.1", "XP_022287977.1", "XP_022287934.1", "XP_022287975.1", "XP_022287971.1", "XP_022287973.1", "XP_022287974.1", "XP_011413590.1", "XP_022295036.1", "XP_022293846.1",
                            "XP_022295029.1", "XP_022293847.1", "XP_022295527.1", "XP_022295528.1", "XP_022295529.1", "XP_019925513.1", "XP_022288097.1", "XP_022289977.1", "XP_022288105.1", "XP_022288100.1", "XP_022288104.1",
                            "XP_019924100.1", "XP_019919899.1", "XP_011421261.1", "XP_022290206.1", "XP_022288934.1", "XP_022288750.1", "XP_022288931.1", "XP_022288746.1", "XP_022288748.1"),
                   labels=c("XP_011431980.1", "XP_022287996.1", "XP_022287997.1", "XP_021341279.1", "XP_021361235.1", "XP_011428387.1", "XP_011428384.1", "XP_022288977.1", "XP_019920103.1", "XP_021350197.1",
                            "XP_011416423.1", "XP_022295668.1", "XP_022291030.1", "XP_022287919.1", "XP_022287917.1", "XP_022287913.1", "XP_022287914.1", "XP_022287912.1", "XP_019919109.1", "XP_019919110.1", "XP_022292970.1",
                            "XP_022294409.1", "XP_022294410.1", "XP_021369350.1", "XP_011452303.1", "XP_022312791.1", "XP_022311552.1", "XP_011441452.1", "XP_022300945.1", "XP_022311926.1", "XP_022311074.1", "XP_021356418.1",
                            "XP_022293362.1", "XP_011455007.1", "XP_021355384.1", "XP_022287287.1", "XP_022287292.1", "XP_019927576.1", "XP_011449366.1", "XP_022320574.1", "XP_022320313.1", "XP_021362910.1", "XP_021372726.1",
                            "XP_021365179.1", "XP_022287983.1", "XP_022287991.1", "XP_021355347.1", "XP_021345058.1", "XP_021360744.1", "XP_022293034.1", "XP_011418793.1", "XP_011418791.1", "XP_021365162.1", "XP_021365133.1",
                            "XP_021345059.1", "XP_021365104.1", "XP_021365115.1", "XP_021366455.1", "XP_022305768.1", "XP_021365097.1", "XP_021353464.1", "XP_021347771.1", "XP_011427116.1", "XP_022291629.1", "XP_022293782.1",
                            "XP_022295018.1", "XP_022335805.1", "XP_022336007.1", "XP_021355270.1", "XP_022293068.1", "XP_022292343.1", "XP_022292342.1", "XP_022293067.1", "XP_019922713.1", "XP_019922711.1", "XP_019922715.1",
                            "XP_019920663.1", "XP_019922714.1", "XP_019922712.1", "XP_021344072.1", "XP_021369394.1", "XP_011435625.1", "XP_022315256.1", "XP_021363603.1", "XP_021363602.1", "XP_011435383.1", "XP_021363252.1",
                            "XP_019925482.1", "XP_019925483.1", "XP_011436597.1", "XP_019925480.1", "XP_019925481.1", "XP_022332918.1", "XP_022332916.1", "XP_022332914.1", "XP_022332915.1", "XP_022332917.1", "XP_022331416.1",
                            "XP_022331414.1", "XP_022331412.1", "XP_022331413.1", "XP_022331415.1", "XP_021345008.1", "XP_021348428.1", "XP_011455592.1", "XP_022291483.1", "XP_011452445.1", "XP_011412926.1", "XP_022322280.1",
                            "XP_022322279.1", "XP_022322278.1", "XP_022337137.1", "XP_021366760.1", "XP_011436809.1", "XP_011445382.1", "XP_022292415.1", "XP_022292412.1", "XP_022291345.1", "XP_019929019.1", "XP_011433586.1",
                            "XP_019921695.1", "XP_019921696.1", "XP_011423762.1", "XP_011423764.1", "XP_011427115.1", "XP_022291924.1", "XP_022293780.1", "XP_022295543.1", "XP_022290906.1", "XP_021360358.1", "XP_011433457.1",
                            "XP_022292108.1", "XP_022288420.1", "XP_021349549.1", "XP_021344945.1", "XP_022292341.1", "XP_022291916.1", "XP_019922709.1", "XP_019922710.1", "XP_021362783.1", "XP_011444161.1", "XP_022289969.1",
                            "XP_022288685.1", "XP_022288684.1", "XP_022288686.1", "XP_021372320.1", "XP_022288032.1", "XP_011437420.1", "XP_022288056.1", "XP_022287930.1", "XP_022287929.1", "XP_011437418.1", "XP_019925758.1",
                            "XP_019926829.1", "XP_011437445.1", "XP_022287977.1", "XP_022287934.1", "XP_022287975.1", "XP_022287971.1", "XP_022287973.1", "XP_022287974.1", "XP_011413590.1", "XP_022295036.1", "XP_022293846.1",
                            "XP_022295029.1", "XP_022293847.1", "XP_022295527.1", "XP_022295528.1", "XP_022295529.1", "XP_019925513.1", "XP_022288097.1", "XP_022289977.1", "XP_022288105.1", "XP_022288100.1", "XP_022288104.1",
                            "XP_019924100.1", "XP_019919899.1", "XP_011421261.1", "XP_022290206.1", "XP_022288934.1", "XP_022288750.1", "XP_022288931.1", "XP_022288746.1", "XP_022288748.1")) +
  guides(fill=guide_legend(ncol=3, title.position="top"))

C_gig_apop_LFC_IAP_tile_plot_COLLAPSED <- C_gig_apop_LFC_IAP_tile_plot_COLLAPSED + 
# color odd numbered groups as blue and even groups 
geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "gray86",], inherit.aes = FALSE, 
          aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
          # fill works if you put outside of aes and inlude the colors in the data itself # add border color
          fill = "gray86",
          color = "gray86", # add border color
          size=0.2, # set border line thickness 
          alpha=0.1) +  # make translucent 
  # color odd numbered groups as blue and even groups 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "slategray1",], inherit.aes = FALSE, 
            aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
            # fill works if you put outside of aes and inlude the colors in the data itself # add border color
            color = "gray86", # add border color
            fill = "slategray1",
            size=0.2, # set border line thickness 
            alpha=0.1)   # make translucent 

#### PLOT CONSTITUTIVELY EXPRESSED IAPS ####
# As with the LFC plot, the proteins with identical sequence have been collapsed for the purpose of plotting with the output tree from RAxML

# Load data from Transcriptome dataframes
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_gig_vst_common_df_all_mat_limma_IAP_gather_avg.RData")
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_vst_common_df_all_mat_limma_IAP_gather_avg.RData")

# Keep tables separate for plotting 
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg$Species <- "Crassostrea_virginica"
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg$Species <- "Crassostrea_gigas"

# Full Join the list of transcripts so the missing labels are in there and can be ordered for plotting 
C_gig_vst_common_df_all_mat_limma_IAP_XP <- full_join(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

# Several NA's meaning that those proteins were ones that were collapsed as duplicates in CD-Hit
C_gig_vst_common_df_all_mat_limma_IAP_XP_collapsed <- C_gig_vst_common_df_all_mat_limma_IAP_XP %>% filter(is.na(node) & !is.na(avg_vst_counts_per_treatment))  %>% 
  dplyr::distinct(protein_id)
C_gig_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6 <- BIR_seq_rm_dup_clstr6 %>% filter(protein_id %in% (C_gig_vst_common_df_all_mat_limma_IAP_XP_collapsed$protein_id))
# Find parent proteins in these clusters
C_gig_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6_cluster <- BIR_seq_rm_dup_clstr6[BIR_seq_rm_dup_clstr6$cluster %in% C_gig_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6$cluster,]

# Recode these proteins for the purpose of plotting
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg$protein_id <- recode(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg$protein_id, 
                                                                      # "XP_019925512.1"= "XP_019925513.1",
                                                                      # "XP_011437419.1"= "XP_011437418.1",
                                                                      "XP_011428385.1"= "XP_011428384.1",  #these proteins only ones to recode after LFC genes removed from list 
                                                                      "XP_011414430.1"= "XP_019919109.1" #these proteins only ones to recode after LFC genes removed from list 
                                                                      # "XP_011445380.1"= "XP_011445382.1",
                                                                      # "XP_011445381.1"= "XP_011445382.1",
                                                                      # "XP_011445383.1"= "XP_011445382.1",
                                                                      # "XP_011436808.1"= "XP_011436809.1"
                                                                      )
# Rejoin the Full list of transcript and check for fixed NA
C_gig_vst_common_df_all_mat_limma_IAP_XP <- full_join(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

# Now check Cvir proteins
# Full Join the list of transcripts so the missing labels are in there and can be ordered for plotting 
C_vir_vst_common_df_all_mat_limma_IAP_XP <- full_join(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])
#several NAs

C_vir_vst_common_df_all_mat_limma_IAP_XP_collapsed <- C_vir_vst_common_df_all_mat_limma_IAP_XP %>% filter(is.na(node) & !is.na(avg_vst_counts_per_treatment))  %>% 
  dplyr::distinct(protein_id)
C_vir_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6 <- BIR_seq_rm_dup_clstr6 %>% filter(protein_id %in% (C_vir_vst_common_df_all_mat_limma_IAP_XP_collapsed$protein_id))
# Find parent proteins in these clusters
C_vir_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6_cluster <- BIR_seq_rm_dup_clstr6[BIR_seq_rm_dup_clstr6$cluster %in% C_vir_vst_common_df_all_mat_limma_IAP_XP_collapsed_BIR_seq_rm_dup_clstr6$cluster,]

# Recode these proteins for the purpose of plotting
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg$protein_id <- recode(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg$protein_id, 
                                                                      "XP_022288682.1"="XP_022288684.1",
                                                                      "XP_022288101.1"="XP_022288100.1", 
                                                                      "XP_022288102.1"="XP_022288100.1", 
                                                                     # "XP_022289978.1"="XP_022289977.1", # ones commented out were removed when the LFC genes were removed 
                                                                      "XP_022287965.1"="XP_022287971.1",
                                                                      "XP_022287969.1"="XP_022287971.1",
                                                                      "XP_022290205.1"="XP_022290206.1",
                                                                      "XP_022288933.1"="XP_022288934.1",
                                                                      # "XP_022295524.1"="XP_022295527.1",
                                                                      "XP_022288031.1"="XP_022288032.1",
                                                                      #"XP_022292109.1"="XP_022292108.1",
                                                                      "XP_022288975.1"="XP_022288977.1",
                                                                      "XP_022293361.1"="XP_022293362.1",
                                                                      #"XP_022292414.1"="XP_022292412.1",
                                                                      #"XP_022292969.1"="XP_022292970.1",
                                                                      #"XP_022293781.1"="XP_022293782.1",
                                                                      #"XP_022291628.1"="XP_022291629.1",
                                                                      "XP_022286791.1"="XP_022295668.1"
                                                                      #"XP_022287921.1"="XP_022287919.1"
                                                                     )

# Rejoin the Full list of transcript and check for fixed NA
C_vir_vst_common_df_all_mat_limma_IAP_XP <- full_join(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg, IAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

# Reorder both to be the order of the tree XPs
# Get the node order from original tree (done in code chunk regarding domain information above)
# Reorder the proteins
C_vir_vst_common_df_all_mat_limma_IAP_XP <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, C_vir_vst_common_df_all_mat_limma_IAP_XP)
C_gig_vst_common_df_all_mat_limma_IAP_XP <- full_join(IAP_MY_CV_CG_raxml_treedata_tip_order, C_gig_vst_common_df_all_mat_limma_IAP_XP)

# Set factor level order of the nodes set levels in reverse order
C_vir_vst_common_df_all_mat_limma_IAP_XP$protein_id <- factor(C_vir_vst_common_df_all_mat_limma_IAP_XP$protein_id, levels = rev(unique(C_vir_vst_common_df_all_mat_limma_IAP_XP$protein_id)))
C_gig_vst_common_df_all_mat_limma_IAP_XP$protein_id <- factor(C_gig_vst_common_df_all_mat_limma_IAP_XP$protein_id, levels = rev(unique(C_gig_vst_common_df_all_mat_limma_IAP_XP$protein_id)))
C_vir_vst_common_df_all_mat_limma_IAP_XP$node <- factor(C_vir_vst_common_df_all_mat_limma_IAP_XP$node, levels = rev(unique(C_vir_vst_common_df_all_mat_limma_IAP_XP$node)))
C_gig_vst_common_df_all_mat_limma_IAP_XP$node <- factor(C_gig_vst_common_df_all_mat_limma_IAP_XP$node, levels = rev(unique(C_gig_vst_common_df_all_mat_limma_IAP_XP$node)))

# Edit factor labels# need to get rid of ROD susceptible vs resistant level and the Dermo tolerant and Dermo susceptible levels
C_vir_vst_common_df_all_mat_limma_IAP_XP$Experiment <- recode_factor(C_vir_vst_common_df_all_mat_limma_IAP_XP$Experiment,
                                                              "ROD_Susceptible"="ROD", 
                                                              "ROD_Resistant"="ROD", 
                                                              "Dermo_Susceptible"="Dermo",
                                                              "Dermo_Tolerant"="Dermo")
C_vir_vst_common_df_all_mat_limma_IAP_XP$Experiment <- factor(C_vir_vst_common_df_all_mat_limma_IAP_XP$Experiment, levels=c("Hatchery_Probiotic", "Pro_RE22","ROD","Dermo", "NA"),
                                                              labels=c("CVBAC-B",  "CVBAC-A","CVBAC-C", "CVPMA", "NA"))

C_vir_vst_common_df_all_mat_limma_IAP_XP$Condition <- factor(C_vir_vst_common_df_all_mat_limma_IAP_XP$Condition, levels=c("Untreated_control","Bacillus_pumilus_RI0695", "Pro_RE22_Control_no_treatment", "Bacillus_pumilus_RI06_95_exposure_6h","Bacillus_pumilus_RI06_95_exposure_24h",
                                                                                                "Phaeobacter_inhibens_S4_exposure_6h", "Phaeobacter_inhibens_S4_exposure_24h", "Vibrio_coralliilyticus_RE22_exposure_6h",
                                                                                                "ROD_Res_Control","ROD_Res_Challenge","ROD_Sus_Control","ROD_Sus_Challenge","Dermo_Sus_36h_Control","Dermo_Sus_28d_Control",
                                                                                                "Dermo_Sus_7d_Control","Dermo_Sus_36h_Injected","Dermo_Sus_7d_Injected","Dermo_Sus_28d_Injected","Dermo_Tol_36h_Control",
                                                                                                "Dermo_Tol_7d_Control","Dermo_Tol_28d_Control","Dermo_Tol_36h_Injected","Dermo_Tol_7d_Injected","Dermo_Tol_28d_Injected"),
                                                  labels=c("Con" , "RI", "Con","RI 6h", "RI 24h", "S4 6h","S4 24h", "RE22" ,
                                                           "S. Con", "S. ROD", "R. Con","R. ROD", "S. Con. 36h", "S. Con. 7d", "S. Con 28d", 
                                                           "S. Pm 36h", "S. Pm 7d", "S. Pm 28d", "T. Con 36h", "T. Con 7d","T. Con 28d",
                                                           "T. Pm 36h", "T. Pm 7d", "T. Pm 28d"))
# Plot Cvir IAP collapsed
# Fix scale (check yellow)
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED <- ggplot(C_vir_vst_common_df_all_mat_limma_IAP_XP, aes(x=Condition, y=protein_id, fill=avg_vst_counts_per_treatment)) + 
  geom_tile() + 
  scale_fill_viridis_c(name = "Average Read Count", 
                       # limits = c(-11,10),
                       breaks = c(0,0.5,1,2,3,4,6,8,11,12), 
                       option="plasma",
                       guide=guide_legend(), na.value = "transparent") +
  facet_grid(.~Experiment, scales="free",space="free", drop= TRUE) + 
  labs(
    #x="Treatment", 
    x= NULL,
    y = NULL,
    title = "*C. virginica* Experiment"
  ) +
  theme_minimal() + 
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size=20, family="sans"),
        plot.title = ggtext::element_markdown(size=20, family="sans", hjust= 0.5),
        axis.text.x.bottom = element_text(size=16, family="sans", face = "bold", angle = 90, vjust=0.5, hjust =1),
        legend.position = "bottom",
        legend.title = element_text(size=16, family="sans"), 
        legend.text = element_text(size=14, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(size=0.2, color="gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(
          size = 16, face = "bold", family="sans"),
        strip.background = element_rect(colour = "white", fill="white")) +
  guides(fill=guide_legend(ncol=3, title.position="top")) 

# add geom_rect boxes
# color odd numbered groups as blue and even groups 
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED <- C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED+ 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "gray86",], inherit.aes = FALSE, 
          aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
          # fill works if you put outside of aes and inlude the colors in the data itself # add border color
          fill = "gray86",
          color = "gray86", # add border color
          size=0.2, # set border line thickness 
          alpha=0.1) +  # make translucent 
  # color odd numbered groups as blue and even groups 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "slategray1",], inherit.aes = FALSE, 
            aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
            # fill works if you put outside of aes and inlude the colors in the data itself # add border color
            color = "gray86", # add border color
            fill = "slategray1",
            size=0.2, # set border line thickness 
            alpha=0.1)   # make translucent 

# drop NA column with gtable (see https://stackoverflow.com/questions/40141684/suppress-na-column-when-faceting)
Cvir_const_IAP_gt <- ggplot_gtable(ggplot_build(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED + theme(legend.position = "none")))
#find column to drop
gtable_show_layout(Cvir_const_IAP_gt) # drop 15
lemon::gtable_show_names(Cvir_const_IAP_gt)
rm_grobs_cvir <- Cvir_const_IAP_gt$layout$name %in% c("strip-t-5","panel-1-5","axis-b-5", "axis-t-5", "axis-r-1", "ylab-r")
# remove grobs
Cvir_const_IAP_gt$grobs[rm_grobs_cvir] <- NULL
Cvir_const_IAP_gt$layout <- Cvir_const_IAP_gt$layout[!rm_grobs_cvir, ]
lemon::gtable_show_names(Cvir_const_IAP_gt)
# remove extra width
Cvir_const_IAP_gt$widths
Cvir_const_IAP_gt$widths[13] = unit(0, "cm")
Cvir_const_IAP_gt$widths[17] = unit(0, "cm")
# check result again
lemon::gtable_show_names(Cvir_const_IAP_gt)
#save plot as ggplot object
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED_NArm <- ggplotify::as.ggplot(Cvir_const_IAP_gt)

# Edit factor levels for C. gig
C_gig_vst_common_df_all_mat_limma_IAP_XP$Experiment <- factor(C_gig_vst_common_df_all_mat_limma_IAP_XP$Experiment, levels = c("Zhang", "Rubio","He","deLorgeril_Susceptible", "deLorgeril_Resistant"), 
                                                              labels = c("CVBAC-A" ,"CVBAC-B" , "CGOSHV1-B" ,"CGOSHV1-A Sus.", "CGOSHV1-A Res."))

C_gig_vst_common_df_all_mat_limma_IAP_XP$Condition <- factor(C_gig_vst_common_df_all_mat_limma_IAP_XP$Condition, levels = c( "Zhang_Control","V_aes_V_alg1_V_alg2","V_tub_V_ang","LPS_M_lut","Rubio_Control","Vcrass_J2_8","Vcrass_J2_9","Vtasma_LGP32",
                                                                              "Vtasma_LMG20012T","Time0_control","6h_control","6h_OsHV1","12h_control","12h_OsHV1","24h_control","24h_OsHV1","48h_control",
                                                                              "48h_OsHV1","120hr_control","120hr_OsHV1","AF21_Resistant_control_0h","AF21_Resistant_6h","AF21_Resistant_12h","AF21_Resistant_24h",
                                                                              "AF21_Resistant_48h","AF21_Resistant_60h","AF21_Resistant_72h","AF11_Susceptible_control_0h","AF11_Susceptible_6h","AF11_Susceptible_12h","AF11_Susceptible_24h","AF11_Susceptible_48h","AF11_Susceptible_60h","AF11_Susceptible_72h"), 
                                                                  labels= c("Con","*V. alg*","*V.tub, V. ang*","LPS, *M. Lut*", "Con","*V. cras* NVir","*V. cras* Vir" ,"*V. tas* Vir","*V. tas* NVir",
                                                                            "0h Con","6h Con", "6h OsHv-1","12h Con","12h OsHv-1", "24h Con","24h OsHv-1",
                                                                            "48h Con","48h OsHv-1","120h Con", "120h OsHv-1","0h Con","6h OsHv-1","12h OsHv-1","24h OsHv-1" ,"48h OsHv-1",
                                                                            "60h OsHv-1","72h OsHv-1" ,"0h Con","6h OsHv-1", "12h OsHv-1","24h OsHv-1","48h OsHv-1" ,
                                                                            "60h OsHv-1","72h OsHv-1"))
# check nrow
C_gig_vst_common_df_all_mat_limma_IAP_XP %>% distinct(protein_id) %>% count() # 184 correct

# Plot C_gig IAP collapsed tree
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED <- ggplot(C_gig_vst_common_df_all_mat_limma_IAP_XP, aes(x=Condition, y=node, fill=avg_vst_counts_per_treatment)) + 
  geom_tile() + 
  scale_fill_viridis_c(name = "Average Read Count", 
                       # limits = c(-11,10),
                       breaks = c(0,0.5,1,2,3,4,6,8,11,12,13,14), 
                       option="plasma",
                       guide=guide_legend(), na.value = "transparent") +
  facet_grid(.~Experiment, scales="free",space="free", drop= TRUE) + 
  labs(
    #x="Treatment", 
    x= NULL,
    y = NULL,
    title = "*C. gigas* Experiment"
  ) +
  theme_minimal() + 
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size=20, family="sans"),
        plot.title = ggtext::element_markdown(size=20, family="sans", hjust= 0.5),
        axis.text.x.bottom = ggtext::element_markdown(size=16, family="sans", face = "bold", angle = 90, vjust=0.5, hjust =1),
        legend.position = "bottom",
        legend.title = element_text(size=16, family="sans"), 
        legend.text = element_text(size=14, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(size=0.2, color="gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = ggtext::element_markdown(
          size = 16, face = "bold", family="sans"),
        strip.background = element_rect(colour = "white", fill="white")) +
  guides(fill=guide_legend(ncol=3, title.position="top")) 
# add geom_rect boxes
# color odd numbered groups as blue and even groups 
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED <- C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED  + 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "gray86",], inherit.aes = FALSE, 
            aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
            # fill works if you put outside of aes and inlude the colors in the data itself # add border color
            fill = "gray86",
            color = "gray86", # add border color
            size=0.2, # set border line thickness 
            alpha=0.1) +  # make translucent 
  # color odd numbered groups as blue and even groups 
  geom_rect(data= IAP_domain_structure_node_order_coordinates[IAP_domain_structure_node_order_coordinates$Color_group == "slategray1",], inherit.aes = FALSE, 
            aes(ymin=as.numeric(ymin),ymax=as.numeric(ymax),xmin=-Inf,xmax=Inf), 
            # fill works if you put outside of aes and inlude the colors in the data itself # add border color
            color = "gray86", # add border color
            fill = "slategray1",
            size=0.2, # set border line thickness 
            alpha=0.1)   # make translucent 

# drop NA column with gtable (see https://stackoverflow.com/questions/40141684/suppress-na-column-when-faceting)
Cgig_const_IAP_gt <- ggplot_gtable(ggplot_build(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED + theme(legend.position = "none")))
#find column to drop
gtable_show_layout(Cgig_const_IAP_gt) # drop 15
lemon::gtable_show_names(Cgig_const_IAP_gt)
rm_grobs <- Cgig_const_IAP_gt$layout$name %in% c("strip-t-6","panel-1-6","axis-b-6", "axis-t-6", "axis-r-1", "ylab-r")
# remove grobs
Cgig_const_IAP_gt$grobs[rm_grobs] <- NULL
Cgig_const_IAP_gt$layout <- Cgig_const_IAP_gt$layout[!rm_grobs, ]
lemon::gtable_show_names(Cgig_const_IAP_gt)
# remove extra width
Cgig_const_IAP_gt$widths
Cgig_const_IAP_gt$widths[15] = unit(0, "cm")
Cgig_const_IAP_gt$widths[19] = unit(0, "cm")
# check result again
lemon::gtable_show_names(Cgig_const_IAP_gt)
#save plot as ggplot object
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED_NArm <- ggplotify::as.ggplot(Cgig_const_IAP_gt)

## Plot vst side by side 
# get legend
C_vir_vst_legend <- get_legend(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED)
C_gig_vst_legend <- get_legend(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_tile_plot_COLLAPSED)

# get legend plot using plot_grid and add on tree legend
vst_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_legend, C_vir_vst_legend, C_gig_vst_legend,
                        nrow = 1, ncol= 4, align="hv", rel_widths  =c(0.5,0.57,1,1))

vst_plots <- plot_grid(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_no_legend, Cvir_const_IAP_gt,
                       Cgig_const_IAP_gt ,
          ncol =3, align="h", axis="tb")

vst_combined <- plot_grid(vst_plots, vst_legend, ncol=1, rel_heights  = c(0.8, 0.1)) + 
  # Add some space at top for labels
  theme(plot.margin = unit(c(1.2,0.0,0.0,0.0), "cm")) +
  draw_plot_label(c("A","B","C"), x= c(0.21, 0.33, 0.66), y = c(1,1,1), size = 30, family = "sans", vjust = 0.2) +
  draw_label("Treatment", x=0.67, y=  0.09, vjust=-0.5, size = 20, fontfamily = "sans", angle= 0) 

# Export Const. expression plot to use in publication 
#ggsave(filename = "IAP_Const_C_vir_C_gig_07302020.tiff", plot=vst_combined, device="tiff",
#       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_Const",
#       width = 40,
#       height = 25,
#       units = "in",
#       dpi=300, limitsize = FALSE)

ggsave(filename = "IAP_Const_C_vir_C_gig_09172020.tiff", plot=vst_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_Const",
       width = 40,
       height = 25,
       units = "in",
       dpi=300, limitsize = FALSE)

ggsave(filename = "IAP_Const_C_vir_C_gig_03102021.tiff", plot=vst_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_Const",
       width = 40,
       height = 25,
       units = "in",
       dpi=300, limitsize = FALSE)

## Plot LFC expression side by side
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_LFC_axis <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_no_legend + aplot::ylim2(C_vir_apop_LFC_IAP_tile_plot_COLLAPSED)
C_vir_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend <- C_vir_apop_LFC_IAP_tile_plot_COLLAPSED + theme(legend.position = "none")
C_gig_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend <- C_gig_apop_LFC_IAP_tile_plot_COLLAPSED + theme(legend.position = "none")

LFC_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_legend,  get_legend(C_vir_apop_LFC_IAP_tile_plot_COLLAPSED), get_legend(C_gig_apop_LFC_IAP_tile_plot_COLLAPSED),
                        nrow = 1, ncol= 4, align="hv", rel_widths  =c(0.5,0.57,1,1))

LFC_plots <- plot_grid(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_LFC_axis, C_vir_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend, C_gig_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend, 
                          ncol=3, align="h", axis="tb")

LFC_combined <- plot_grid(LFC_plots, LFC_legend, ncol=1, rel_heights  = c(0.8, 0.1)) + 
  # Add some space at top for labels
  theme(plot.margin = unit(c(1.2,0.0,0.0,0.0), "cm")) +
  draw_plot_label(c("A","B","C"), x= c(0.21, 0.33, 0.66), y = c(1,1,1), size = 30, family = "sans", vjust = 0.2) +
  draw_label("Treatment", x=0.67, y=  0.09, vjust=-0.5, size = 20, fontfamily = "sans", angle= 0)

# Export LFC plot ## USE FOR PUBLICATION##
#ggsave(filename = "IAP_LFC_C_vir_C_gig_07302020.tiff", plot=LFC_combined, device="tiff",
#       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_LFC",
#       width = 40,
#       height = 30,
#       units = "in",
#       dpi=300)

ggsave(filename = "IAP_LFC_C_vir_C_gig_09172020.tiff", plot=LFC_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_LFC",
       width = 40,
       height = 30,
       units = "in",
       dpi=300)

ggsave(filename = "IAP_LFC_C_vir_C_gig_03102021.tiff", plot=LFC_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_LFC",
       width = 40,
       height = 30,
       units = "in",
       dpi=300)


#### PLOT IAP PROTEIN LFC AND CONST TREE WITH GENE NAMES COLLAPSED ####

class(IAP_MY_CV_CG_raxml_treedata_collapsed)
IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name <- as_tibble(IAP_MY_CV_CG_raxml_treedata_collapsed)

# For repeated gene names, change the duplicated nodes to be a dash
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name$gene_collapse <-  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name$gene
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name <-   IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name %>%  arrange(gene, node)
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name$gene_collapse <- ave(
    # use base R ave function to replace the duplicated values with NA
    IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name$gene_collapse, 
    IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name$gene, 
    FUN = function(a) replace(a, duplicated(a), NA_integer_)
  )

  # change all NAs to be dashes
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name <-   IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name %>%
    mutate( gene_collapse = replace_na(gene_collapse, "---------------------"))
  
  # check the duplicated values are now dashes
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name %>% arrange(gene, node) %>% group_by(gene) %>% 
    mutate(number_proteins = n()) %>% filter(number_proteins > 1 )

  # convert back to tree object
  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name <- as.treedata(  IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name)
  
# Plot collapsed tree
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME <- 
  ggtree(IAP_MY_CV_CG_raxml_treedata_collapsed_gene_name, aes(color=Species, fill=Species),  branch.length = "none") + 
  geom_tiplab(aes(label=gene_collapse), fontface="bold", size =3.5, offset=0) + # geom_tiplab2 flips the labels correctly
  # add circle for 90-100 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 90), color = "black", fill="black", shape=21, size=2.0) +
  # add triangle for 70-89 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 70 & as.numeric(bootstrap) < 90),color = "black", fill="black", shape=24, size=2.0) +
  # add upside down traingle for 50-69 instead of bootstrap values
  geom_nodepoint(aes(subset = as.numeric(bootstrap) >= 50  &  as.numeric(bootstrap) < 70 ), color = "black",fill="black", shape=25, size=2.0) +
  # Add shape for tips removed (see IAP_shape_node above)
  geom_point2(aes(subset=(node==12)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==13)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==36)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==37)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==48)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==91)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==92)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==97)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==108)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==109)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==113)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==120)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==137)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==138)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==139)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==140)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==141)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==142)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==143)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==147)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==150)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  geom_point2(aes(subset=(node==179)), shape=22, size=2.0, color = '#c55d32', fill='#c55d32') +
  ## Add clade labels for the 21 domain groups  domain groups using the internal node number
  geom_cladelabel(261, label="1",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') + # get node order from below 
  geom_cladelabel(254, label="2",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(233, label="3",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(228, label="4",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(217, label="5",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(211, label="6",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(207, label="7",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(198, label="8",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(201, label="9",  offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(311, label="10", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(308, label="11", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(291, label="12", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(304, label="13", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(329, label="14", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(340, label="15", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(343, label="16", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(347, label="17", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(353, label="18", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(359, label="19", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(188, label="20", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black') +
  geom_cladelabel(366, label="21", offset = 9.5, offset.text=0.5, family="sans", fontsize = 7, barsize=2, color='black', extend = 0.5) +
  #Edit theme
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=14, family="sans"),
        legend.title = element_text(size=16, family="sans")) +
  #geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 2.0, fontface="bold") + # allows for subset
  xlim(-70,31.8) + #change scaling so branch lengths are smaller and all alias labels are showing
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis"),
                      labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis")) +
  guides(col = guide_legend(ncol =1, title.position = "top", override.aes = aes(label = "")) ) # need to override aes to get rid of "a"

# RE-EXPORT LFC AND CONST FIGURES WITH THIS NEW COLLAPSED TREE BEING USED 
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_legend <- cowplot::get_legend(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME)
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_no_legend <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME + 
  theme(legend.position='none')

## Plot const with dashed gene names
# get legend plot using plot_grid and add on tree legend
vst_gene_name_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_legend, C_vir_vst_legend, C_gig_vst_legend,
                        nrow = 1, ncol= 4, align="hv", rel_widths  =c(0.5,0.57,1,1))

vst_gene_name_plots <- plot_grid(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_no_legend, Cvir_const_IAP_gt,
                       Cgig_const_IAP_gt ,
                       ncol =3, align="h", axis="tb")

vst_gene_name_combined <- plot_grid(vst_gene_name_plots, vst_gene_name_legend, ncol=1, rel_heights  = c(0.8, 0.1)) + 
  # Add some space at top for labels
  theme(plot.margin = unit(c(1.2,0.0,0.0,0.0), "cm")) +
  draw_plot_label(c("A","B","C"), x= c(0.21, 0.33, 0.66), y = c(1,1,1), size = 30, family = "sans", vjust = 0.2) +
  draw_label("Treatment", x=0.67, y=  0.09, vjust=-0.5, size = 20, fontfamily = "sans", angle= 0) 

# Export Const. expression plot
ggsave(filename = "IAP_Const_C_vir_C_gig_dashed_gene_01122021.tiff", plot=vst_gene_name_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_Const",
       width = 40,
       height = 25,
       units = "in",
       dpi=300, limitsize = FALSE)

ggsave(filename = "IAP_Const_C_vir_C_gig_dashed_gene_03102021.tiff", plot=vst_gene_name_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_Const",
       width = 40,
       height = 25,
       units = "in",
       dpi=300, limitsize = FALSE)

## Plot LFC expression side by side
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_LFC_axis <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_no_legend + aplot::ylim2(C_vir_apop_LFC_IAP_tile_plot_COLLAPSED)

LFC_gene_name_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_legend,  get_legend(C_vir_apop_LFC_IAP_tile_plot_COLLAPSED), get_legend(C_gig_apop_LFC_IAP_tile_plot_COLLAPSED),
                        nrow = 1, ncol= 4, align="hv", rel_widths  =c(0.5,0.57,1,1))

LFC_gene_name_plots <- plot_grid(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_gene_name_LFC_axis, C_vir_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend, C_gig_apop_LFC_IAP_tile_plot_COLLAPSED_no_legend, 
                       ncol=3, align="h", axis="tb")

LFC_gene_name_combined <- plot_grid(LFC_gene_name_plots, LFC_gene_name_legend, ncol=1, rel_heights  = c(0.8, 0.1)) + 
  # Add some space at top for labels
  theme(plot.margin = unit(c(1.2,0.0,0.0,0.0), "cm")) +
  draw_plot_label(c("A","B","C"), x= c(0.21, 0.33, 0.66), y = c(1,1,1), size = 30, family = "sans", vjust = 0.2) +
  draw_label("Treatment", x=0.67, y=  0.09, vjust=-0.5, size = 20, fontfamily = "sans", angle= 0)

# Export LFC plot
ggsave(filename = "IAP_LFC_C_vir_C_gig_gene_name_01122021.tiff", plot=LFC_gene_name_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_LFC",
       width = 40,
       height = 30,
       units = "in",
       dpi=300)

ggsave(filename = "IAP_LFC_C_vir_C_gig_gene_name_03102021.tiff", plot=LFC_gene_name_combined, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_LFC",
       width = 40,
       height = 30,
       units = "in",
       dpi=300)
#### EXPORT PROTEIN DOMAIN TREE PLOT WITH COLLAPSED GENE NAMES ####

# Use data frame with collapsed --- for transcripts from the same gene 
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME

###  Export and arrange domain plot with tree
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME_legend <- cowplot::get_legend(IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME)
IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME_no_legend <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME + 
  theme(legend.position='none')

IAP_MY_CV_CG_treecollapsed_GENE_NAME <- IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME_no_legend + aplot::ylim2(IAP_Interproscan_domain_plot_no_legend)

IAP_tr_dom_collapsed_GENE_NAME <- plot_grid(NULL,IAP_MY_CV_CG_treecollapsed_GENE_NAME, IAP_Interproscan_domain_plot_no_legend, ncol=3, align='h', rel_widths = c(0.2, 0.7,0.8)) +
  # Add some space at top for labels
  theme(plot.margin = unit(c(1,0.0,0.0,0.0), "cm")) 
IAP_tr_dom_collapsed__GENE_NAME_legend <- plot_grid(NULL, IAP_MY_CV_CG_raxml_treedata_vertical_collapsed_GENE_NAME_legend, IAP_Interproscan_domain_plot_legend,
                                         nrow = 1, align="hv", rel_widths  =c(0.7, 0.7,1)) 

## Create combined figure for publication
IAP_tr_dom_plus_legend_collapsed_GENE_NAME <- plot_grid(IAP_tr_dom_collapsed_GENE_NAME, IAP_tr_dom_collapsed__GENE_NAME_legend,  ncol=1, rel_heights  = c(0.8, 0.1)) +
  # add labels for plot components
  draw_plot_label(c("A","B","C"), x= c(0.38, 0.53, 0.9), y = c(1,1,1), size = 30, family = "sans")

## Export plot with tree and domains aligned 

ggsave(filename = "IAP_tr_dom_plus_legend_collapsed_GENE_NAME_03102021.tiff", plot=IAP_tr_dom_plus_legend_collapsed_GENE_NAME, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_tree_domain",
       width = 34,
       height = 27,
       units = "in",
       dpi=300)


#### CREATE TABLE WITH IAP DOMAIN STRUCTURE FOR ALL IAP PROTEINS ####

# This data frames is currently collapsed where duplicate proteins have been removed, need to add back in those removed into the correct groups
IAP_domain_structure <- read_csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_Domain_Structure_groups.csv")

# Examine those that have been removed 
BIR_XP_gff_species_join_haplotig_collapsed
BIR_XP_gff_species_join_haplotig_collapsed_CV_CG <- BIR_XP_gff_species_join_haplotig_collapsed %>% filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas")

setdiff(IAP_domain_structure$protein_id, BIR_XP_gff_species_join_haplotig_collapsed_CV_CG$protein_id) 
    # 35 proteins missing in second data frame because Mizuhopecten yessoensis proteins are also included in the original domain group table
length(setdiff(BIR_XP_gff_species_join_haplotig_collapsed_CV_CG$protein_id, IAP_domain_structure$protein_id)) # 83 proteins that were removed from being duplicated sequences for RAxML

# Examine the cluster file where it shows the groupings of clusters and those that were removed
BIR_seq_rm_dup_clstr6_CV_CG <- BIR_seq_rm_dup_clstr6 %>% filter(Species == "Crassostrea_virginica" | Species == "Crassostrea_gigas")

# Join collapsed IAP domain structure data frame to this to put them into their clusters 
BIR_seq_rm_dup_clstr6_CV_CG_domain <- left_join(BIR_seq_rm_dup_clstr6_CV_CG, IAP_domain_structure[,c("protein_id","Domain_Name","Number")])

# Do all have the same domain desigations within clusters
BIR_seq_rm_dup_clstr6_CV_CG_domain %>% group_by(cluster, Domain_Name) %>% distinct(cluster, Domain_Name) %>%  filter(n() >1) # 0, all have only 1 domain name added

# Group by cluster and paste the values for those that are empty
BIR_seq_rm_dup_clstr6_CV_CG_domain_group_fill <- BIR_seq_rm_dup_clstr6_CV_CG_domain %>% group_by(cluster) %>% fill(Domain_Name, Number) %>% 
  # change NAs to be called not classified
  mutate(Domain_Name = case_when(
    is.na(Domain_Name) ~ "not_classified",
    TRUE ~ Domain_Name
  ))

# subset and rename data frame 
IAP_domain_structure_no_dup_rm <- BIR_seq_rm_dup_clstr6_CV_CG_domain_group_fill[,c("protein_id", "gene","product", "Species","Domain_Name", "Number")]

# check that all original IAP proteins identified are there
setdiff(IAP_domain_structure_no_dup_rm$protein_id, BIR_XP_gff_species_join_haplotig_collapsed_CV_CG$protein_id) # 6 proteins now extra in the IAP domain structure file, these are the protein haplotigs that were removed 
# remove those proteins that are haplotigs

haplotid_protein_id_rm <- c("XP_022308010.1", "XP_022292821.1", "XP_022336127.1", "XP_022290466.1", "XP_022308067.1", "XP_022304464.1")
IAP_domain_structure_no_dup_rm <- IAP_domain_structure_no_dup_rm %>% filter(!(protein_id %in% haplotid_protein_id_rm))

# check again 
setdiff(IAP_domain_structure_no_dup_rm$protein_id, BIR_XP_gff_species_join_haplotig_collapsed_CV_CG$protein_id) # files match 
setdiff(BIR_XP_gff_species_join_haplotig_collapsed_CV_CG$protein_id, IAP_domain_structure_no_dup_rm $protein_id) # no proteins missing from origincal haplotig collapsed file

# check number of total proteins 
IAP_domain_structure_no_dup_rm %>% distinct(protein_id, Species) %>% count(Species) # this is correct! 
 # Species                   n
 # <chr>                 <int>
 #   1 Crassostrea_gigas        74
 # 2 Crassostrea_virginica   158

# export IAP_domain_structure_no_dup_rm for use in DEG workspace 
save(IAP_domain_structure_no_dup_rm, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/IAP_domain_structure_no_dup_rm.RData")

#### TOTAL ANNOTATED IAP DOMAIN STATS ####

IAP_domain_structure_no_dup_rm

# What percent of the total annotated in genome transcripts in each is in each domain type
# perform separately so I can create a wide table for data visualization
C_vir_IAP_domain_structure_no_dup_rm_count <- IAP_domain_structure_no_dup_rm %>%
  filter(Species == "Crassostrea_virginica") %>% 
  # add up number of unique proteins in each experiment
  distinct(protein_id , .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(CV_total_per_domain = n()) %>%
  mutate(CV_percent_of_total = CV_total_per_domain/sum(CV_total_per_domain)) %>% arrange(desc(CV_total_per_domain))

C_gig_IAP_domain_structure_no_dup_rm_count <- IAP_domain_structure_no_dup_rm %>%
  filter(Species == "Crassostrea_gigas") %>% 
  # add up number of unique proteins in each experiment
  distinct(protein_id , .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(CG_total_per_domain = n()) %>%
  mutate(CG_percent_of_total = CG_total_per_domain/sum(CG_total_per_domain)) %>% arrange(desc(CG_total_per_domain))

# join 
C_vir_C_gig_IAP_domain_structure_no_dup_rm_count <- left_join(C_vir_IAP_domain_structure_no_dup_rm_count, C_gig_IAP_domain_structure_no_dup_rm_count)
C_vir_C_gig_IAP_domain_structure_no_dup_rm_count$CG_total_per_domain <- replace_na(C_vir_C_gig_IAP_domain_structure_no_dup_rm_count$CG_total_per_domain,0)
C_vir_C_gig_IAP_domain_structure_no_dup_rm_count$CG_percent_of_total <- replace_na(C_vir_C_gig_IAP_domain_structure_no_dup_rm_count$CG_percent_of_total,0)

# create table for Paper Figure 
C_vir_C_gig_IAP_domain_structure_no_dup_rm_count_table <- C_vir_C_gig_IAP_domain_structure_no_dup_rm_count %>% 
  mutate(Domain_Name = case_when(Domain_Name == "not_classified" ~ "Non-Domain Grouped", TRUE ~ Domain_Name)) %>%
  gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Total Identified Transcripts with Each Domain Structure**")) %>%
  tab_spanner(label = md("*Crassostrea virginica*"), columns = c("CV_total_per_domain", "CV_percent_of_total")) %>%
  tab_spanner(label = md("*Crassostrea gigas*"), columns = c("CG_total_per_domain", "CG_percent_of_total")) %>%
  cols_label( CV_total_per_domain  =md("**Transcripts<br>Per Type**"),
              CV_percent_of_total = md("**Percent<br>of Total**"),
              CG_total_per_domain  =md("**Transcripts<br>Per Type**"),
              CG_percent_of_total = md("**Percent of Total**")) %>%
  fmt_percent(columns = vars(CV_percent_of_total), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(CG_percent_of_total), decimals = 2) %>% # format as percent
  summary_rows(fns=list(Total = "sum"), 
               columns = c("CV_total_per_domain","CG_total_per_domain"), 
               formatter = fmt_number) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_options(table.font.color = "black", table.font.size = 20) 
  
# save as png
gtsave(C_vir_C_gig_IAP_domain_structure_no_dup_rm_count_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/C_vir_C_gig_IAP_domain_structure_no_dup_rm_count_table.png")

# get same statistics but how many genes in each group 
# perform separately so I can create a wide table for data visualization
C_vir_IAP_domain_structure_no_dup_rm_count_GENE <- IAP_domain_structure_no_dup_rm %>%
  filter(Species == "Crassostrea_virginica") %>% 
  # add up number of unique proteins in each experiment
  distinct(gene, .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(CV_total_per_domain = n()) %>%
  mutate(CV_percent_of_total = CV_total_per_domain/sum(CV_total_per_domain)) %>% arrange(desc(CV_total_per_domain))

C_gig_IAP_domain_structure_no_dup_rm_count_GENE <- IAP_domain_structure_no_dup_rm %>%
  filter(Species == "Crassostrea_gigas") %>% 
  # add up number of unique proteins in each experiment
  distinct(gene , .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(CG_total_per_domain = n()) %>%
  mutate(CG_percent_of_total = CG_total_per_domain/sum(CG_total_per_domain)) %>% arrange(desc(CG_total_per_domain))


### Export Data frame of DEG IAPs and domain groupsing to in WGCNA

# export data frame
save(C_vir_C_gig_apop_LFC_IAP_OG_domain_structure, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/C_vir_C_gig_apop_LFC_IAP_OG_domain_structure.RData")

#### IAP SUMMARY TABLE STATS LFC, CONST, DOMAINS ####
### THIS CODE MAKES INITIAL TABLES AND STATISTICS, MORE DETAILED PLOTS AND PATHWAY ANALYSIS ARE IN THE DESEQ2 ANALYSIS SCRIPT 
# Join original IAP stats from C_gig and C_vir experiments with the domain info
C_vir_apop_LFC_IAP_OG_domain_structure <- left_join(C_vir_apop_LFC_IAP_OG, IAP_domain_structure_no_dup_rm) %>% filter(!is.na(protein_id))
C_vir_apop_LFC_IAP_OG_domain_structure$Species <- "Crassostrea_virginica"
# change ID column to transcript_id for joining purposes
colnames(C_vir_apop_LFC_IAP_OG_domain_structure)[6] <- "transcript_id"
C_gig_apop_LFC_IAP_OG_domain_structure <- left_join(C_gig_apop_LFC_IAP_OG, IAP_domain_structure_no_dup_rm)
C_gig_apop_LFC_IAP_OG_domain_structure$Species <- "Crassostrea_gigas"

# combine into one
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure <- rbind(C_vir_apop_LFC_IAP_OG_domain_structure, C_gig_apop_LFC_IAP_OG_domain_structure)

# save for used in other script
save(C_vir_C_gig_apop_LFC_IAP_OG_domain_structure, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/C_vir_C_gig_apop_LFC_IAP_OG_domain_structure")

# Were alternativly spliced forms used across different experiments
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_alternative <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure %>% distinct(Species, gene,transcript_id, experiment) %>% ungroup() %>%
  group_by(gene) %>% mutate(gene_count = n()) %>% filter( gene_count >1) %>% distinct(gene, transcript_id, .keep_all = TRUE) 

C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_alternative <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_alternative %>%
  ungroup() %>% group_by(gene) %>% 
  mutate(transcript_type = n()) %>% filter(transcript_type >1)


# How many transcripts not used at all from any experiment
length(unique(C_vir_apop_LFC_IAP_OG_domain_structure$protein_id)) # 37 total used 
length(unique(C_gig_apop_LFC_IAP_OG_domain_structure$protein_id)) # 38 total used 

BIR_XP_gff_species_join_haplotig_collapsed_CV <- BIR_XP_gff_species_join_haplotig_collapsed_CV_CG  %>% filter(Species =="Crassostrea_virginica")
BIR_XP_gff_species_join_haplotig_collapsed_CG <- BIR_XP_gff_species_join_haplotig_collapsed_CV_CG %>% filter(Species =="Crassostrea_gigas")

length(setdiff(BIR_XP_gff_species_join_haplotig_collapsed_CV$protein_id, unique(C_vir_apop_LFC_IAP_OG_domain_structure$protein_id))) # 121 proteins not used at all in any C. vir experiment
length(setdiff(BIR_XP_gff_species_join_haplotig_collapsed_CG$protein_id, unique(C_gig_apop_LFC_IAP_OG_domain_structure$protein_id))) # 36 proteins not used at all in any C. gig experiment

# What percent of the total LFC transcripts in each is in each domain type
# perform separately so I can create a wide table for data visualization
C_vir_apop_LFC_IAP_OG_domain_structure_count <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure %>%
  filter(Species == "Crassostrea_virginica") %>% 
  # add up number of TOTAL proteins in each experiment
  group_by(Domain_Name) %>%
  dplyr::summarise(CV_total_per_domain = n()) %>%
  mutate(CV_percent_of_total = CV_total_per_domain/sum(CV_total_per_domain)) %>% arrange(desc(CV_total_per_domain)) 

    # add in missing row for TII-DD-RING
Cvir_TII_DD_RING <- data.frame(Domain_Name = c("TII-DD-RING"), CV_total_per_domain = 0, CV_percent_of_total = 0)
C_vir_apop_LFC_IAP_OG_domain_structure_count <- rbind(C_vir_apop_LFC_IAP_OG_domain_structure_count, Cvir_TII_DD_RING)

C_gig_apop_LFC_IAP_OG_domain_structure_count <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure %>%
  filter(Species == "Crassostrea_gigas") %>% 
  # add up number of TOTAL proteins in each experiment
  group_by(Domain_Name) %>%
  dplyr::summarise(CG_total_per_domain = n()) %>%
  mutate(CG_percent_of_total = CG_total_per_domain/sum(CG_total_per_domain)) %>% arrange(desc(CG_total_per_domain))
    # add in missing row for NZBIR-TII-UBA-DD-RING
Cgig_NZBIR <- data.frame(Domain_Name = c("NZBIR-TII-UBA-DD-RING"), CG_total_per_domain = 0, CG_percent_of_total = 0)
C_gig_apop_LFC_IAP_OG_domain_structure_count <- rbind(C_gig_apop_LFC_IAP_OG_domain_structure_count, Cgig_NZBIR)

# join 
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count <- left_join(C_vir_apop_LFC_IAP_OG_domain_structure_count, C_gig_apop_LFC_IAP_OG_domain_structure_count)
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count$CG_total_per_domain <- replace_na(C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count$CG_total_per_domain,0)
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count$CG_percent_of_total <- replace_na(C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count$CG_percent_of_total,0)

# create table
C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count_table <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count %>% 
  mutate(Domain_Name = case_when(Domain_Name == "not_classified" ~ "Non-Domain Grouped", TRUE ~ Domain_Name)) %>%
  gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Total Significantly Differentially Expressed Transcripts with Each Domain Structure**")) %>%
  tab_spanner(label = md("*Crassostrea virginica*"), columns = c("CV_total_per_domain", "CV_percent_of_total")) %>%
  tab_spanner(label = md("*Crassostrea gigas*"), columns = c("CG_total_per_domain", "CG_percent_of_total")) %>%
  cols_label( CV_total_per_domain  =md("**Transcripts<br>Per Type**"),
              CV_percent_of_total = md("**Percent<br>of Total**"),
              CG_total_per_domain  =md("**Transcripts<br>Per Type**"),
              CG_percent_of_total = md("**Percent of Total**")) %>%
  fmt_percent(columns = vars(CV_percent_of_total), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(CG_percent_of_total), decimals = 2) %>% # format as percent
  summary_rows(fns=list(Total = "sum"), 
               columns = c("CV_total_per_domain","CG_total_per_domain"), 
               formatter = fmt_number) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_options(table.font.color = "black", table.font.size = 20) 

# save as png
gtsave(C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/C_vir_C_gig_apop_LFC_IAP_OG_domain_structure_count_table.png")

# Join IAP domain information with constitutively expressed IAP information
C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_domain <- left_join(C_vir_vst_common_df_all_mat_limma_IAP_gather_avg,  IAP_domain_structure_no_dup_rm)
C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_domain <- left_join(C_gig_vst_common_df_all_mat_limma_IAP_gather_avg,  IAP_domain_structure_no_dup_rm)

# Combine tables to get stats about which domain types are present 
# First combine LFC and const
C_vir_C_gig_apop_LFC_IAP_full_XP_stats <- C_vir_C_gig_apop_LFC_IAP_OG_domain_structure %>% 
  mutate(Data_Type = "LFC") %>% select(protein_id, experiment, Species,product, Data_Type, group_by_sim)

C_vir_vst_common_df_all_mat_limma_IAP_XP_stats <- C_vir_vst_common_df_all_mat_limma_IAP_gather_avg_domain  %>% ungroup() %>%
  mutate(Data_Type = "Const") %>% select(protein_id, Experiment, Species,product, Data_Type, Condition)
colnames(C_vir_vst_common_df_all_mat_limma_IAP_XP_stats)[2] <- "experiment"
colnames(C_vir_vst_common_df_all_mat_limma_IAP_XP_stats)[6] <- "group_by_sim"

C_gig_vst_common_df_all_mat_limma_IAP_XP_stats <- C_gig_vst_common_df_all_mat_limma_IAP_gather_avg_domain %>% ungroup() %>%
  mutate(Data_Type = "Const") %>% select(protein_id, Experiment, Species,product, Data_Type, Condition)
colnames(C_gig_vst_common_df_all_mat_limma_IAP_XP_stats)[2] <- "experiment"
colnames(C_gig_vst_common_df_all_mat_limma_IAP_XP_stats)[6] <- "group_by_sim"

# Make all experiment and group_by_sim levels the same - will save headaches later 
levels(factor(C_vir_C_gig_apop_LFC_IAP_full_XP_stats$experiment))
    #"deLorgeril_res", "deLorgeril_sus", "Dermo","Hatchery_Probiotic_RI", "He","Lab_Pro_RE22","ROD","Rubio","Zhang"
levels(factor(C_vir_C_gig_apop_LFC_IAP_full_XP_stats$group_by_sim))
    #[1] "deLorgeril_res_12hr"    "deLorgeril_res_24hr"    "deLorgeril_res_48hr"    "deLorgeril_res_60hr"    "deLorgeril_res_6hr"     "deLorgeril_res_72hr"    "deLorgeril_sus_12hr"   
    #[8] "deLorgeril_sus_24hr"    "deLorgeril_sus_48hr"    "deLorgeril_sus_60hr"    "deLorgeril_sus_6hr"     "deLorgeril_sus_72hr"    "Dermo_Susceptible_28d"  "Dermo_Susceptible_36hr"
    #[15] "Dermo_Susceptible_7d"   "Dermo_Tolerant_28d"     "Dermo_Tolerant_36hr"    "Dermo_Tolerant_7d"      "Hatchery_Probiotic_RI"  "He_120hr"               "He_12hr"               
    #[22] "He_24hr"                "He_48hr"                "He_6hr"                 "Lab_RE22"               "Lab_RI_6hr"             "Lab_RI_RI_24hr"         "Lab_S4_24hr"           
    #[29] "Lab_S4_6hr"             "ROD_susceptible_seed"   "Rubio_J2_8"             "Rubio_J2_9"             "Rubio_LGP32"            "Rubio_LMG20012T"        "Zhang_LPS"             
    #[36] "Zhang_Valg"             "Zhang_Vtub" 

levels(factor(C_vir_vst_common_df_all_mat_limma_IAP_XP_stats$experiment))
 #"Dermo_Susceptible"  "Dermo_Tolerant"     "Hatchery_Probiotic" "Pro_RE22"           "ROD_Resistant"      "ROD_Susceptible" 
levels(factor(C_gig_vst_common_df_all_mat_limma_IAP_XP_stats$experiment))
 # "deLorgeril_Resistant"   "deLorgeril_Susceptible" "He"                     "Rubio"                  "Zhang" 

levels(factor(C_vir_vst_common_df_all_mat_limma_IAP_XP_stats$group_by_sim))
    # [1] "Dermo_Tol_7d_Control"                    "Dermo_Tol_28d_Control"                   "Dermo_Tol_36h_Control"                   "Dermo_Tol_28d_Injected"                 
    # [5] "Dermo_Tol_36h_Injected"                  "Dermo_Tol_7d_Injected"                   "Dermo_Sus_36h_Control"                   "Dermo_Sus_28d_Control"                  
    # [9] "Dermo_Sus_7d_Control"                    "Dermo_Sus_28d_Injected"                  "Dermo_Sus_36h_Injected"                  "Dermo_Sus_7d_Injected"                  
    # [13] "Bacillus_pumilus_RI0695"                 "Untreated_control"                       "ROD_Res_Control"                         "ROD_Res_Challenge"                      
    # [17] "ROD_Sus_Control"                         "ROD_Sus_Challenge"                       "Pro_RE22_Control_no_treatment"           "Bacillus_pumilus_RI06_95_exposure_24h"  
    # [21] "Bacillus_pumilus_RI06_95_exposure_6h"    "Phaeobacter_inhibens_S4_exposure_24h"    "Phaeobacter_inhibens_S4_exposure_6h"     "Vibrio_coralliilyticus_RE22_exposure_6h"

levels(factor(C_gig_vst_common_df_all_mat_limma_IAP_XP_stats$group_by_sim))
    # [1] "Zhang_Control"               "LPS_M_lut"                   "V_aes_V_alg1_V_alg2"         "V_tub_V_ang"                 "Rubio_Control"               "Vcrass_J2_8"                
    # [7] "Vcrass_J2_9"                 "Vtasma_LGP32"                "Vtasma_LMG20012T"            "AF21_Resistant_60h"          "AF21_Resistant_48h"          "AF21_Resistant_24h"         
    # [13] "AF21_Resistant_72h"          "AF21_Resistant_6h"           "AF21_Resistant_12h"          "AF21_Resistant_control_0h"   "AF11_Susceptible_72h"        "AF11_Susceptible_60h"       
    # [19] "AF11_Susceptible_48h"        "AF11_Susceptible_24h"        "AF11_Susceptible_12h"        "AF11_Susceptible_6h"         "AF11_Susceptible_control_0h" "Time0_control"              
    # [25] "6h_control"                  "6h_OsHV1"                    "12h_control"                 "12h_OsHV1"                   "24h_control"                 "24h_OsHV1"                  
    # [31] "48h_control"                 "48h_OsHV1"                   "120hr_control"               "120hr_OsHV1"                

# Keeping levels from the LFC experiment and changing levels for the vst const
C_vir_vst_common_df_all_mat_limma_IAP_XP_stats$experiment <- plyr::revalue(C_vir_vst_common_df_all_mat_limma_IAP_XP_stats$experiment,
              c("Dermo_Susceptible" = "Dermo", "Dermo_Tolerant" = "Dermo", "Hatchery_Probiotic"="Hatchery_Probiotic_RI", "Pro_RE22"="Lab_Pro_RE22", "ROD_Resistant"="ROD","ROD_Susceptible"="ROD"))
C_gig_vst_common_df_all_mat_limma_IAP_XP_stats$experiment <- plyr::revalue(C_gig_vst_common_df_all_mat_limma_IAP_XP_stats$experiment,
               c("deLorgeril_Resistant"="deLorgeril_res", "deLorgeril_Susceptible"= "deLorgeril_sus","He"="He","Rubio"="Rubio","Zhang"="Zhang"))

# Rbind all together
LFC_cont_comb <- rbind(C_vir_C_gig_apop_LFC_IAP_full_XP_stats, C_vir_vst_common_df_all_mat_limma_IAP_XP_stats,
                       C_gig_vst_common_df_all_mat_limma_IAP_XP_stats)

# check experiment levels
levels(factor(LFC_cont_comb$experiment))

# Join with Domain type information
LFC_cont_comb_domain_type <- left_join(LFC_cont_comb, IAP_domain_structure_no_dup_rm[, c("protein_id","Domain_Name","Number")])

# Export recoded LFC data frames with transcript IDs to search for modules with matchin transcripts in the WGCNA data
save(LFC_cont_comb_domain_type, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/LFC_cont_comb_domain_type.RData")

# Number of LFC and Const IAP in each experiment, group, and Domain_Name
LFC_cont_comb_summary_count <- LFC_cont_comb_domain_type  %>% 
  distinct(protein_id, experiment,Species, Data_Type, group_by_sim, Domain_Name) %>%
  group_by(experiment, Species, Data_Type, group_by_sim, Domain_Name) %>%
  count() 

# Number of LFC total per species
LFC_comb_summary_count_species <- LFC_cont_comb_summary_count %>% filter(Data_Type == "LFC") %>% 
  dplyr::group_by(Species, Domain_Name, Data_Type)  %>% dplyr::summarize(total_IAP = sum(n)) %>% ungroup() %>%
  group_by(Species) %>% dplyr::mutate(percent = total_IAP/sum(total_IAP)*100) %>% arrange(Data_Type, -total_IAP)

# Are any Domain types unique to 1 experiment?
LFC_cont_comb_domain_type_protein_number <- LFC_cont_comb_domain_type %>% 
  group_by(Domain_Name, experiment, Data_Type, Species) %>% 
  distinct(protein_id) %>% 
  count()

# create ggplot plot to better visualize patterns in domain usage
LFC_cont_comb_domain_type_protein_number_plot <- ggplot(LFC_cont_comb_domain_type_protein_number, aes(x=Domain_Name, y=n, fill=experiment, color=experiment)) + geom_col() + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust =1)) + facet_grid(Species~Data_Type) 

### Make aggregate table regarding how many DEG IAP transcripts used out of total in genome, percent shared, percent unique
# Calculate number of times transcripts are shared across experiments
LFC_cont_comb_summary_unique_shared <- LFC_cont_comb %>% 
  filter(Data_Type == "LFC") %>% 
  rowwise() %>%
    mutate(num_groups = sum(group_by(., experiment, Species) %>%
                              distinct(protein_id, experiment) %>%
                              ungroup() %>%
                              select(protein_id) %>%
                              unlist() %in% protein_id)) %>% 
    arrange(num_groups) 

# check numbers
LFC_cont_comb_summary_unique_shared %>% distinct(protein_id, group_by_sim, experiment, .keep_all = TRUE) %>% count(experiment) 
  
# Calculate the total number of proteins used in each experiment (only counting shared ones once) and the proportion of total identified IAPs
LFC_cont_comb_proportion_IAP_used <-  LFC_cont_comb %>%
  filter(Data_Type == "LFC") %>% 
  ungroup() %>% 
  distinct(experiment, group_by_sim, protein_id, .keep_all = TRUE) %>% # include group by sim when getting the full count
  group_by(experiment) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>%
  distinct(experiment, total_count, experiment_count, Species) %>%
  mutate(proportion_used = as.numeric(case_when(
           Species == "Crassostrea_virginica"~ as.character(experiment_count/158*100),
           Species == "Crassostrea_gigas"~ as.character(experiment_count/74*100),
           TRUE~NA_character_))) %>%
  arrange(desc(proportion_used))
  
# Calculate how many unique transcripts in each experiment and proportion unique
LFC_cont_comb_summary_unique <- LFC_cont_comb_summary_unique_shared %>%
  # first make count for the individual proteins across all levels of each experiment rather than the total count
  distinct(protein_id, experiment, .keep_all = TRUE) %>% 
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>% 
  ungroup() %>%
  filter(num_groups == 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  dplyr::mutate(total_unique = n()) %>%
  distinct(experiment, total_unique, experiment_count) %>%
  mutate(percent_not_shared_each_exp = total_unique/experiment_count*100)

# Calculate percent of shared transcripts between experiments in each experiment 
LFC_cont_comb_summary_shared <- LFC_cont_comb_summary_unique_shared %>%
  # first make count for the individual proteins across all levels of each experiment rather than the total count
  distinct(protein_id, experiment, .keep_all = TRUE) %>% 
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>% 
  ungroup() %>%
  filter(num_groups > 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  dplyr::mutate(total_shared = n()) %>%
  distinct(experiment, total_shared, experiment_count) %>%
  mutate(percent_shared_each_exp = total_shared/experiment_count*100)

## final df with count, %unique % shared % percent used out of all IAPS
LFC_cont_comb_summary_shared_used <- left_join(LFC_cont_comb_proportion_IAP_used, LFC_cont_comb_summary_shared)
LFC_cont_comb_summary_shared_unique <-   left_join(LFC_cont_comb_summary_shared_used, LFC_cont_comb_summary_unique)

# Make table of summary IAP shared and unique transcripts
# change C. gigas experiment labels
levels(factor(LFC_cont_comb_summary_shared_unique$experiment))
LFC_cont_comb_summary_shared_unique$experiment <- factor(LFC_cont_comb_summary_shared_unique$experiment, levels = c("Hatchery_Probiotic_RI","Lab_Pro_RE22" ,  "ROD","Dermo" ,
                                                                                                                    "Zhang","Rubio","He","deLorgeril_sus", "deLorgeril_res"),
                                                         labels= c("Hatchery\nPro. RI","Lab Pro. S4, RI\n or RE22",  "ROD","Dermo","Zhang Vibrio spp." ,   
                                                                                                                   "Rubio Vibrio spp.","He OsHV-1","de Lorgeril Sus. OsHV-1", "de Lorgeril Res. OsHV-1"))
LFC_cont_comb_summary_shared_unique_table <- LFC_cont_comb_summary_shared_unique %>%
  select(Species, experiment,total_count, proportion_used, experiment_count,total_shared,percent_shared_each_exp, total_unique, percent_not_shared_each_exp) %>%
  # replace all NA's
  mutate(
    across(everything(), ~replace_na(.x, 0)),
    proportion_used = proportion_used*0.01,
    percent_shared_each_exp= percent_shared_each_exp*0.01,
    percent_not_shared_each_exp = percent_not_shared_each_exp*0.01) %>%
  gt::gt(rowname_col = "experiment", groupname = "Species") %>%
  tab_header(title = gt::md("**Unique and Shared Differentially Expressed IAP Transcripts Across Experiments**")) %>%
  cols_label( total_count = md("**Total Differentially<br>Expressed IAP<br>Transcripts**"),
              experiment_count = md("**Number of Different Differentially<br>Expressed IAP<br>Transcripts in Each Experiment**"),
              proportion_used  =md("**Percent Differentially<br>Expressed of Total<br>Genome IAP Transcripts**"),
              total_shared = md("**Number of Different IAP Transcripts<br>Significant Across<br>Multiple Experiments**"),
              total_unique = md("**Number of Different <br>IAP Transcripts Only Expressed<br>in One Experiment**"),
              percent_shared_each_exp = md("**Percent Different IAP Transcripts<br>Significant Across<br>Multiple Experiments**"),
              percent_not_shared_each_exp = md("**Percent Different <br>IAP Transcripts Only Expressed<br>in One Experiment**")) %>%
  fmt_percent(columns = vars(proportion_used), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(percent_shared_each_exp), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(percent_not_shared_each_exp), decimals = 2) %>% # format as percent
tab_row_group(group = "Crassostrea virginica",rows = c(5,6,8,9)) %>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups(groups = "Crassostrea virginica")) %>%
  tab_row_group(group = "Crassostrea gigas",rows = c(1:4,7)) %>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups(groups = "Crassostrea gigas")) %>%
  summary_rows(groups = TRUE, fns = list(Average = "mean", SD = "sd"), formatter = fmt_percent, columns = c("proportion_used","percent_shared_each_exp","percent_not_shared_each_exp")) %>%
  tab_options(table.font.color = "black")

# save as png
gtsave(LFC_cont_comb_summary_shared_unique_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/LFC_cont_comb_summary_shared_unique_table.png")

### Analyze domain structure of unique and shared DEG transcripts ###
LFC_cont_comb_summary_unique_shared_domain_type <- left_join(LFC_cont_comb_summary_unique_shared, IAP_domain_structure_no_dup_rm)

# Count domain structure of unique transcripts in each experiments
LFC_cont_comb_summary_unique_shared_domain_type %>%
filter(num_groups == 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(experiment, Domain_Name) %>%
  dplyr::summarise(total_unique_per_domain = n()) %>% View()

# count domain structure of unique transcripts across all experiments 
LFC_cont_comb_summary_unique_shared_domain_type_unique <- LFC_cont_comb_summary_unique_shared_domain_type %>%
  filter(num_groups == 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(total_unique_per_domain = n()) %>%
  arrange(desc(total_unique_per_domain)) 

# Count domain structure of shared transcripts in each experiments
LFC_cont_comb_summary_unique_shared_domain_type %>%
  filter(num_groups > 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(experiment, Domain_Name) %>%
  filter(!is.na(Domain_Name)) %>%
  dplyr::summarise(total_shared_per_domain = n()) %>% 
  arrange(desc(total_shared_per_domain)) %>% View()

# count domain structure of shared transcripts across all experiments 
LFC_cont_comb_summary_unique_shared_domain_type_shared <- LFC_cont_comb_summary_unique_shared_domain_type %>%
  filter(num_groups > 1) %>%
  # add up number of unique proteins in each experiment
  distinct(experiment, protein_id, .keep_all = TRUE) %>%
  group_by(Domain_Name) %>%
  dplyr::summarise(total_shared_per_domain = n()) %>%
  filter(!is.na(Domain_Name)) %>%
  arrange(desc(total_shared_per_domain))

# Join two together and create table
unique_shared_by_domain <- left_join(LFC_cont_comb_summary_unique_shared_domain_type_shared, LFC_cont_comb_summary_unique_shared_domain_type_unique) %>%
  mutate(total_unique_per_domain = replace_na(total_unique_per_domain, 0))

# create table
unique_shared_by_domain_table <- unique_shared_by_domain %>% 
gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Number Unique and Shared Transcripts with Each Domain Structure Across Experiments**")) %>%
  cols_label( Domain_Name = md("**Domain Structure Type**"),
              total_unique_per_domain  =md("**Total Unique**"),
              total_shared_per_domain = md("**Total Shared**")) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*"))

# save as png
gtsave(unique_shared_by_domain_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/unique_shared_by_domain_table.png")

### Analyze Unique and Shared GENES across experiments
# Join gene info to LFC and cont table 
LFC_cont_comb_gene <- left_join(LFC_cont_comb, unique(IAP_domain_structure_no_dup_rm[,c("protein_id","gene","Domain_Name")]))

# Calculate number of times genes are shared across experiments
LFC_cont_comb_summary_unique_shared_GENE <- LFC_cont_comb_gene %>% 
  filter(Data_Type == "LFC") %>% 
  rowwise() %>%
  mutate(num_groups = sum(group_by(., experiment, Species) %>%
                            distinct(gene, experiment) %>%
                            ungroup() %>%
                            select(gene) %>%
                            unlist() %in% gene)) %>% 
  arrange(num_groups) 

# check numbers
LFC_cont_comb_summary_unique_shared_GENE %>% distinct(gene, group_by_sim, experiment, .keep_all = TRUE) %>% count(experiment) 

# Calculate the total number of genes used in each experiment (only counting shared ones once) and the proportion of total identified IAP genes
LFC_cont_comb_proportion_IAP_used_GENE <-  LFC_cont_comb_gene %>%
  filter(Data_Type == "LFC") %>% 
  ungroup() %>% 
  distinct(experiment, group_by_sim, gene, .keep_all = TRUE) %>% # include group by sim when getting the full count
  group_by(experiment) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  distinct(experiment, gene, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>%
  distinct(experiment, total_count, experiment_count, Species) %>%
  mutate(proportion_used = as.numeric(case_when(
    Species == "Crassostrea_virginica"~ as.character(experiment_count/69*100),
    Species == "Crassostrea_gigas"~ as.character(experiment_count/40*100),
    TRUE~NA_character_))) %>%
  arrange(desc(proportion_used))

# Calculate how many unique genes in each experiment and proportion unique
LFC_cont_comb_summary_unique_GENE <- LFC_cont_comb_summary_unique_shared_GENE %>%
  # first make count for the individual genes across all levels of each experiment rather than the total count
  distinct(gene, experiment, .keep_all = TRUE) %>% 
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>% 
  ungroup() %>%
  filter(num_groups == 1) %>%
  # add up number of unique genes in each experiment
  distinct(experiment, gene, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  dplyr::mutate(total_unique = n()) %>%
  distinct(experiment, total_unique, experiment_count) %>%
  mutate(percent_not_shared_each_exp = total_unique/experiment_count*100)

# Calculate percent of shared genes between experiments in each experiment 
LFC_cont_comb_summary_shared_GENE <- LFC_cont_comb_summary_unique_shared_GENE %>%
  # first make count for the individual proteins across all levels of each experiment rather than the total count
  distinct(gene, experiment, .keep_all = TRUE) %>% 
  group_by(experiment) %>%
  mutate(experiment_count = n()) %>% 
  ungroup() %>%
  filter(num_groups > 1) %>%
  # add up number of shared genes in each experiment
  distinct(experiment, gene, .keep_all = TRUE) %>%
  group_by(experiment) %>%
  dplyr::mutate(total_shared = n()) %>%
  distinct(experiment, total_shared, experiment_count) %>%
  mutate(percent_shared_each_exp = total_shared/experiment_count*100)

## final df with count, %unique % shared % percent used out of all IAPS
LFC_cont_comb_summary_shared_used_GENE <- left_join(LFC_cont_comb_proportion_IAP_used_GENE, LFC_cont_comb_summary_shared_GENE)
LFC_cont_comb_summary_shared_unique_GENE <-   left_join(LFC_cont_comb_summary_shared_used_GENE, LFC_cont_comb_summary_unique_GENE)

# Make table of summary IAP shared and unique transcripts
# change C. gigas experiment labels
levels(factor(LFC_cont_comb_summary_shared_unique_GENE$experiment))
LFC_cont_comb_summary_shared_unique_GENE$experiment <- factor(LFC_cont_comb_summary_shared_unique_GENE$experiment, levels = c("Hatchery_Probiotic_RI","Lab_Pro_RE22" ,  "ROD","Dermo" ,
                                                                                                                    "Zhang","Rubio","He","deLorgeril_sus", "deLorgeril_res"),
                                                         labels= c("Hatchery\nPro. RI","Lab Pro. S4, RI\n or RE22",  "ROD","Dermo","Zhang Vibrio spp." ,   
                                                                   "Rubio Vibrio spp.","He OsHV-1","de Lorgeril Sus. OsHV-1", "de Lorgeril Res. OsHV-1"))
LFC_cont_comb_summary_shared_unique_table_GENE <- LFC_cont_comb_summary_shared_unique_GENE %>%
  select(Species, experiment,total_count, proportion_used, experiment_count,total_shared,percent_shared_each_exp, total_unique, percent_not_shared_each_exp) %>%
  # replace all NA's
  mutate(
    across(everything(), ~replace_na(.x, 0)),
    proportion_used = proportion_used*0.01,
    percent_shared_each_exp= percent_shared_each_exp*0.01,
    percent_not_shared_each_exp = percent_not_shared_each_exp*0.01) %>%
  gt::gt(rowname_col = "experiment", groupname = "Species") %>%
  tab_header(title = gt::md("**Unique and Shared Differentially Expressed IAP Genes Across Experiments**")) %>%
  cols_label( total_count = md("**Total Differentially<br>Expressed IAP<br>Genes**"),
              experiment_count = md("**Number of Different Differentially<br>Expressed IAP<br>Genes in Each Experiment**"),
              proportion_used  =md("**Percent Differentially<br>Expressed of Total<br>Genome IAP Genes**"),
              total_shared = md("**Number of Different IAP Genes<br>Significant Across<br>Multiple Experiments**"),
              total_unique = md("**Number of Different <br>IAP Genes Only Expressed<br>in One Experiment**"),
              percent_shared_each_exp = md("**Percent Different IAP Genes<br>Significant Across<br>Multiple Experiments**"),
              percent_not_shared_each_exp = md("**Percent Different <br>IAP Genes Only Expressed<br>in One Experiment**")) %>%
  fmt_percent(columns = vars(proportion_used), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(percent_shared_each_exp), decimals = 2) %>% # format as percent
  fmt_percent(columns = vars(percent_not_shared_each_exp), decimals = 2) %>% # format as percent
  tab_row_group(group = "Crassostrea virginica",rows = c(5,6,8,9)) %>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups(groups = "Crassostrea virginica")) %>%
  tab_row_group(group = "Crassostrea gigas",rows = c(1:4,7)) %>%
  tab_style(style = cell_text(style = "italic"),locations = cells_row_groups(groups = "Crassostrea gigas")) %>%
  summary_rows(groups = TRUE, fns = list(Average = "mean", SD = "sd"), formatter = fmt_percent, columns = c("proportion_used","percent_shared_each_exp","percent_not_shared_each_exp")) %>%
  tab_options(table.font.color = "black")

# save as png
gtsave(LFC_cont_comb_summary_shared_unique_table_GENE, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/LFC_cont_comb_summary_shared_unique_table_GENE.png")

## Generate const IAP TABLE
LFC_cont_comb_summary_count_CONST_IAP_TABLE <- LFC_cont_comb_summary_count %>%
  ungroup() %>%
  distinct(Domain_Name, Species, .keep_all = TRUE) %>% 
  filter(Data_Type == "Const") %>% 
  # change NA to be ungrouped
  mutate(Domain_Name = case_when(
    Domain_Name == "not_classified" ~ "Non-Domain Grouped",
    TRUE ~ Domain_Name
  )) %>%
  select(Domain_Name, n, Species) %>% 
  spread(Species, n, fill = 0) %>% 
  gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Domain Structure of Constitutively Expressed IAP Transcripts**")) %>%
  cols_label( Domain_Name = md("**Domain Structure Type**"),
              Crassostrea_gigas  =md("*Crassostrea gigas*"),
              Crassostrea_virginica = md("*Crassostrea virginica*")) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_footnote(footnote = md("Note: 3 *C. gigas* and 9 *C. virginica* transcripts were not placed into domain structure groups based on bootstrap support"),
                             location = cells_stub("Non-Domain Grouped")) %>%
  #add total column
  grand_summary_rows(fns = list(Total="sum"), formatter = fmt_number, decimals = 0)%>%
  tab_options(table.font.color = "black")
# save as png
gtsave(LFC_cont_comb_summary_count_CONST_IAP_TABLE, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/LFC_cont_comb_summary_count_CONST_IAP_TABLE.png")

# Generate LFC IAP TABLE
levels(factor(LFC_cont_comb_summary_count$experiment))
#[1] "deLorgeril_res"        "deLorgeril_sus"        "Dermo"                 "Hatchery_Probiotic_RI" "He"                    "Lab_Pro_RE22"          "ROD"                  
# [8] "Rubio"                 "Zhang"
LFC_cont_comb_summary_count_LFC_IAP_TABLE <- LFC_cont_comb_summary_count %>%
  ungroup() %>%
  filter(Data_Type == "LFC") %>% 
  dplyr::count(experiment, Domain_Name) %>%
  mutate(Domain_Name = case_when(Domain_Name == "not_classified" ~ "Non-Domain Grouped", TRUE ~ Domain_Name)) %>%
  select(experiment, Domain_Name, n) %>% 
  # spread table into wide format
  spread(experiment, n, fill = 0) %>%
  gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Domain Structure of Significantly Differentially Expressed IAP Transcripts**")) %>%
  cols_label( Domain_Name = md("**Domain Structure Type**"),
              Hatchery_Probiotic_RI ="Hatchery\nProbiotic RI" ,
              Lab_Pro_RE22 ="Lab Pro. S4,\nRI or RE22",
              Zhang =md("Zhang *Vibrio* spp."),
              Rubio =md("Rubio *Vibrio* spp."),
              deLorgeril_sus  = "de Lorgeril\nSus. OsHV-1",
              deLorgeril_res = "de Lorgeril\nRes. OsHV-1" ,
              He = "He OsHV-1") %>%
  #add total column
  grand_summary_rows(fns = list(Total="sum"), formatter = fmt_number, decimals = 0) %>%
  # add spanner over experiments for each species
  tab_spanner(
    label = md("*Crassostrea gigas*"),
    columns = vars(He, Zhang,Rubio,deLorgeril_sus, deLorgeril_res)) %>%
  tab_spanner(
    label = md("*Crassostrea virginica*"),
      columns = vars(Hatchery_Probiotic_RI, Lab_Pro_RE22, ROD, Dermo)) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_options(table.font.color = "black")
# save as png
gtsave(LFC_cont_comb_summary_count_LFC_IAP_TABLE, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/LFC_cont_comb_summary_count_LFC_IAP_TABLE.png")

## Create table with statistics on DEGS, apop DEGs, percent apop DEGs
# get statistics on percent of IAPs out of all apoptosis transcript by joining with summary LFC 
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_gig_sig_table.RData")
nrow(C_vir_gig_sig_table) # 38
# levels don't match 
levels(factor(LFC_comb$group_by_sim))
levels(factor(C_vir_gig_sig_table$group_by_sim))
levels(factor(C_vir_gig_sig_table$experiment)) # "deLorgeril" "Dermo"      "He"         "Pro_RE22"   "Probiotic"  "ROD"        "Rubio"      "Zhang" 
C_vir_gig_sig_table$group_by_sim <- as.factor(C_vir_gig_sig_table$group_by_sim)
C_vir_gig_sig_table$experiment <- as.factor(C_vir_gig_sig_table$experiment)

# export table of apoptosis data statistics with gtable for C. virginica
C_vir_sig_table_apop_table <- C_vir_gig_sig_table %>%
  # change labels and add species
  mutate(Species = case_when(
    experiment == "deLorgeril" |experiment ==  "He"| experiment == "Rubio"| experiment == "Zhang" ~"Crassostrea_gigas",
    experiment == "Dermo" |experiment ==  "Pro_RE22"| experiment == "Probiotic"| experiment == "ROD"   ~"Crassostrea_virginica",
    TRUE ~ NA_character_))  %>% 
  mutate(apop_percent = apop_percent*0.01) %>% 
  filter(Species == "Crassostrea_virginica") %>%
  # change order
  select(group_by_sim, experiment, sig_total, num_sig_apop, apop_percent) %>%
  gt::gt(rowname_col = "group_by_sim", groupname_col = "experiment") %>%
  # use md function to invoke markdown format and add title
  tab_header(title = gt::md("**Significant Differentially Expressed *C. virginica* Transcripts**")) %>% 
  fmt_percent(columns = vars(apop_percent), decimals = 2) %>% # format as percent
  # fix column labels
  cols_label( group_by_sim = md("**Experiment**"),
              sig_total  =md("**Total Significant DEGs**"),
              num_sig_apop = md("**Significant Apoptosis DEGs**"),
              apop_percent = md("**Percent Apoptosis DEGs**")) %>%
  # change group labels and put in bold
  tab_row_group(group = "ROD",rows = 1:2) %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_row_groups(groups = "ROD")) %>%
  tab_row_group(group = "Hatchery Probiotic RI",rows = 3) %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_row_groups(groups = "Hatchery Probiotic RI")) %>% 
  tab_row_group(group = "Dermo",rows = 4:9) %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_row_groups(groups = "Dermo")) %>%
  tab_row_group(group = "Lab Probiotic S4,RI, RE22",rows =10:14) %>%
  tab_style(style = cell_text(weight = "bold"),locations =cells_row_groups(groups = "Lab Probiotic S4,RI, RE22")) %>%
  summary_rows(groups = c("ROD","Dermo","Lab Probiotic S4,RI, RE22"), fns = list(Average = "mean", SD = "sd"), formatter = fmt_percent, columns = "apop_percent") %>%
  summary_rows(groups = c("ROD","Dermo","Lab Probiotic S4,RI, RE22"), fns = list(Total = "sum"), columns = "num_sig_apop") %>%
  tab_options(table.font.color = "black")
# save as png
gtsave(C_vir_sig_table_apop_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/C_vir_sig_table_apop_table.png")

# export table of apoptosis data statistics with gtable for C. gigas
C_gig_sig_table_apop_table <- C_vir_gig_sig_table %>%
  # change labels and add species
  mutate(Species = case_when(
    experiment == "deLorgeril" |experiment ==  "He"| experiment == "Rubio"| experiment == "Zhang" ~"Crassostrea_gigas",
    experiment == "Dermo" |experiment ==  "Pro_RE22"| experiment == "Probiotic"| experiment == "ROD"   ~"Crassostrea_virginica",
    TRUE ~ NA_character_))  %>%
  mutate(apop_percent = apop_percent*0.01) %>% 
  filter(Species == "Crassostrea_gigas") %>%
  # change order
  select(group_by_sim, experiment, sig_total, num_sig_apop, apop_percent) %>%
  gt::gt(rowname_col = "group_by_sim", groupname_col = "experiment") %>%
  # use md function to invoke markdown format and add title
  tab_header(title = gt::md("**Significant Differentially Expressed *C. gigas* Transcripts**")) %>% 
  fmt_percent(columns = vars(apop_percent), decimals = 2) %>% # format as percent
  # fix column labels
  cols_label( group_by_sim = md("**Experiment**"),
              sig_total  =md("**Total Significant DEGs**"),
              num_sig_apop = md("**Significant Apoptosis DEGs**"),
              apop_percent = md("**Percent Apoptosis DEGs**")) %>%
  # change group labels and put in bold
  tab_row_group(group = "Zhang Vibrio spp.",rows = 1:3) %>%
  tab_style( style = cell_text(weight = "bold"), locations = cells_row_groups(groups = "Zhang Vibrio spp.")) %>% 
  tab_row_group(group = "Rubio Vibrio spp.",rows = 4:7) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups(groups = "Rubio Vibrio spp.")) %>%
  tab_row_group(group = "He OsHV-1",rows = 8:12) %>%
  tab_style(style = cell_text(weight = "bold"),locations = cells_row_groups(groups = "He OsHV-1"))   %>%
  tab_row_group(group = "de Lorgeril OsHV-1",rows = 13:24) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups(groups = "de Lorgeril OsHV-1")) %>%
  summary_rows(groups = TRUE, fns = list(Avgerage = "mean", SD = "sd"), formatter = fmt_percent, columns = "apop_percent") %>%
  summary_rows(groups = TRUE, fns = list(Total = "sum"), columns = "num_sig_apop")%>%
  tab_options(table.font.color = "black")
# save as png
gtsave(C_gig_sig_table_apop_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/C_gig_sig_table_apop_table.png")

#### IAP SUMMARY TABLES INTO HEATMAPS ####
# Use the package complex heatmap to create better looking graphic than including a lot of separate tables 
# Structure:
    # ROWS = each DEG and const experimental challenge
    # COLUMNS:
      # Heatmap one
        # 1. Total DEGs: heatmap
        # 2. Apoptosis DEGs: heatmap
        # 3. % apoptosis DEGs: bar graph 
        # - also include text annotation with the sample size for each group
      # Heatmap 2
        # 4. Number IAP DEGs
        # 5-19: column of numbers for each domain structure
      # Heatmap 3
        # 20. # of genes in each challenge that are shared in other experiments
        # 21. # of genes in each challenge that are unique to that experiment
        # 22. % shared genes
        # 23. % unique genes
      # Heatmap 4
        # Same as heatmap 3 but with transcripts instead of genes 

## Heatmap 1 (axes as groups in each experiment)
# Generate the matrix for each data set
total_DEG_mat <- C_vir_gig_sig_table[,c(1,4)]
total_DEG_mat <- total_DEG_mat %>% column_to_rownames(., var = "group_by_sim")
total_DEG_mat <- as.matrix(total_DEG_mat)

apop_DEG_mat <- C_vir_gig_sig_table[,c(1,3)]
apop_DEG_mat <- apop_DEG_mat %>% column_to_rownames(., var = "group_by_sim")
apop_DEG_mat <- as.matrix(apop_DEG_mat)

apop_perc_DEG_mat <- C_vir_gig_sig_table[,c(1,5)]
apop_perc_DEG_mat <- apop_perc_DEG_mat %>% column_to_rownames(., var = "group_by_sim")
apop_perc_DEG_mat <- as.matrix(apop_DEG_mat)

# generate domain structure matrix
LFC_group_domain_count <- LFC_cont_comb_summary_count %>% ungroup() %>% filter(Data_Type == "LFC") %>% distinct(Domain_Name, group_by_sim, n) %>%
  spread(Domain_Name, n, fill = 0)   %>% column_to_rownames(., var = "group_by_sim")
LFC_group_domain_count_mat <- as.matrix(LFC_group_domain_count)

# generate color matrices for each 

# define plotting options
ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)

# plot heatmap list
ht_list = 
  
## Heatmap 2 (axes as individual experiments)
# generate matrices of shared and unique TRANSCRIPTs in each experiment as compared with other experiments
  LFC_cont_comb_summary_shared_unique

Total_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, count) %>% column_to_rownames(., var = "experiment") %>% as.matrix()
shared_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, total_shared) %>% column_to_rownames(., var = "experiment") %>% as.matrix()
unique_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, total_unique) %>% column_to_rownames(., var = "experiment") %>% as.matrix()

propotion_used_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, proportion_used) %>% column_to_rownames(., var = "experiment") %>% as.matrix()
propotion_shared_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, percent_shared_each_exp) %>% column_to_rownames(., var = "experiment") %>% as.matrix()
propotion_unique_IAP_trans_mat <- LFC_cont_comb_summary_shared_unique %>% select(experiment, percent_not_shared_each_exp) %>% column_to_rownames(., var = "experiment") %>% as.matrix()

# generate matrices of shared and unique GENES in each experiment as compared with other experiments
# Calculate number of times GENES are shared across experiments



# create matrix for each row

## define plotting options
ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)

# plot heatmap list

ht_list = 

experiment = as.data.frame(C_vir_gig_sig_table[,c(1,2)])
class(experiment)
ha = ComplexHeatmap::HeatmapAnnotation(type = experiment, annotation_name_side = "left")

ComplexHeatmap::Heatmap(total_DEG_mat, name = "total significant transcripts", row_names_side = "left", row_order = rownames(total_DEG_mat), top_annotation = ha)
                        
                      


# Will come back and finish this code later if its necessary 

#### COMPARE BIR TYPE VARIANTS BETWEEN IAP DOMAIN GROUPS ####

# Comparing file with the type 1 BIR variants with the full IAP non collapsed list 
BIR_domain_model_MY_CV_CG_type_updated
IAP_domain_structure_no_dup_rm

# Join BIR variant info by protein_id
IAP_domain_structure_no_dup_rm_BIR_variant <- left_join(IAP_domain_structure_no_dup_rm, BIR_domain_model_MY_CV_CG_type_updated[,c("protein_id","Type","Species")])
# some NAs are present due to the BIR protein sequences that were collapsed for being identical 

# Do different clusters with the same domain structure have different BIR variant types?
IAP_domain_structure_no_dup_rm_BIR_variant_distinct_type <- IAP_domain_structure_no_dup_rm_BIR_variant %>% distinct(Domain_Name, Number, Type) %>% arrange(Domain_Name, Type)
  # results show that where there are two separate clusters with the same overall domain type they mostly have the same type variant. More within cluster diversity than between cluster 

# how many BIR of each type
BIR_domain_model_MY_CV_CG_type_updated %>% distinct(protein_id, Type, Species) %>% count(Type, Species) %>% arrange(Type, desc(n)) %>% ungroup() %>%
  group_by(Species) %>% mutate(total = sum(n)) %>% mutate(percent = n/total*100) %>% View()

#### IAP AND CONST GENE AND TRANSCRIPT USAGE ####

# Remember I found 69 unique genes in C. virginica and 40 unique genes in C. gigas
IAP_domain_structure_no_dup_rm %>% distinct(gene, Species) %>% count(Species)
#Species                   n
#<chr>                 <int>
#  1 Crassostrea_gigas        40
#2 Crassostrea_virginica    69

IAP_domain_structure_no_dup_rm %>% distinct(protein_id, Species) %>% count(Species)
#Species                   n
#<chr>                 <int>
#  1 Crassostrea_gigas        74
#2 Crassostrea_virginica   158

# Join gene info to LFC and cont table 
LFC_cont_comb_gene <- left_join(LFC_cont_comb, unique(IAP_domain_structure_no_dup_rm[,c("protein_id","gene","Domain_Name")]))

# How many genes of the total IAP genes are being used in each species and data type across experiments - this includes ones that overlap between!
LFC_cont_comb_gene %>% distinct(gene, Data_Type, Species) %>% count(Data_Type, Species)
  #Data_Type               Species  n
  #1     Const     Crassostrea_gigas 13
  #2     Const Crassostrea_virginica 38
  #3       LFC     Crassostrea_gigas 25 (62.5%)
  #4       LFC Crassostrea_virginica 25
  # LOC111099688 recoded to  LOC111105597 for graphing purposes on the tree b/c of cd identical sequences
# How many transcripts of the total IAP transcripts are being used in each species and data type across experiments
LFC_cont_comb_gene %>% distinct(protein_id, Data_Type, Species) %>% count(Data_Type, Species)
#Data_Type               Species  n
#1     Const     Crassostrea_gigas 16
#2     Const Crassostrea_virginica 40
#3       LFC     Crassostrea_gigas 38
#4       LFC Crassostrea_virginica 37

## How many unique genes and transcripts being used across both LFC and cont (how much of total gene diversity is being used)
LFC_cont_comb_gene %>% distinct(gene, Species) %>% count(Species)
  #Species  n
  #1     Crassostrea_gigas 33
  #2 Crassostrea_virginica 53

LFC_cont_comb_gene %>% distinct(protein_id, Species) %>% count(Species)
#Species  n
#1     Crassostrea_gigas 53 transcripts
#2 Crassostrea_virginica 76

## how many LFC transcripts for each gene being used?
LFC_cont_comb_gene_transcript_usage <- LFC_cont_comb_gene %>% filter(Data_Type == "LFC") %>% distinct(protein_id, experiment, Species, group_by_sim, gene) %>%
  group_by(experiment, Species, group_by_sim, gene) %>% count() %>% arrange(desc(n)) 

#### Are different LFC transcripts within genes being used between experiments?
LFC_cont_comb_gene_prot_comb <- LFC_cont_comb_gene %>% filter(Data_Type == "LFC") %>% distinct(protein_id, experiment, Species, group_by_sim, gene)  %>%
  group_by(experiment, group_by_sim, gene) %>% 
  # paste together combo of prot id's
  mutate(comb_prot_id = paste(protein_id, collapse = "_")) %>% distinct(gene, experiment, comb_prot_id, Species)

# find which genes have several different transcript combos across more than one experiment
LFC_cont_comb_gene_prot_comb_diff <- LFC_cont_comb_gene_prot_comb %>% ungroup() %>% distinct(comb_prot_id, gene, Species, experiment) %>% group_by(gene) %>% filter(n()>1) %>% 
  ungroup() %>% distinct(gene, Species,experiment) %>% group_by(gene) %>% filter(n()>1)

# How many genes have differential trancript usage between experiments and how many do not?
LFC_cont_comb_gene_prot_comb_diff %>% ungroup() %>% distinct(gene, Species) %>% count(Species)
# 6 genes in each have differential transcript usage - a small percentage of the total genes express a different set of transcripts between experiments
#Species                   n
#<chr>                 <int>
#  1 Crassostrea_gigas         17
#2 Crassostrea_virginica     6

# Which genes are those that have differential transcript usage?
LFC_cont_comb_gene_prot_comb_differential_usage <- LFC_cont_comb_gene_prot_comb[LFC_cont_comb_gene_prot_comb$gene %in% LFC_cont_comb_gene_prot_comb_diff$gene,]

# view these genes but with the product and domain type info
LFC_cont_comb_gene_differential_usage <- LFC_cont_comb_gene %>% filter(Data_Type == "LFC") %>% filter(gene %in% LFC_cont_comb_gene_prot_comb_diff$gene)

# Should any of these be collapsed based on having duplicate sequence
LFC_cont_comb_gene_differential_usage_collapse <- left_join(LFC_cont_comb_gene_differential_usage, BIR_seq_rm_dup_clstr6_CV_CG[,c("protein_id","stat","cluster")])
     
### Which annotated IAP genes and transcripts not utilized in any experiment
Not_used_IAP_prot_gene <- IAP_domain_structure_no_dup_rm[!(IAP_domain_structure_no_dup_rm$protein_id %in% LFC_cont_comb_gene$protein_id),]
Not_used_IAP_gene <- IAP_domain_structure_no_dup_rm[!(IAP_domain_structure_no_dup_rm$gene %in% LFC_cont_comb_gene$gene),]

# How many genes not used in any experiment?
Not_used_IAP_gene_uniq <- Not_used_IAP_gene %>% distinct(Species, gene)
Not_used_IAP_gene %>% distinct(Species, gene) %>% count(Species)
#Species                   n
#<chr>                 <int>
#  1 Crassostrea_gigas         7
#2 Crassostrea_virginica    16

# How many transcripts not used in any experiment?
Not_used_IAP_prot_gene %>% distinct(Species, protein_id) %>% count(Species)
  #Species                   n
  #<chr>                 <int>
  #  1 Crassostrea_gigas        21
  #2 Crassostrea_virginica    82

# For those genes that shared between experiments- are different transcripts being used in the different experiments?
LFC_cont_comb_summary_unique_shared_GENE %>%
  # first make count for the individual proteins across all levels of each experiment rather than the total count
  distinct(gene, experiment, .keep_all = TRUE) %>% 
  ungroup() %>%
  filter(num_groups > 1) %>% 
  # concatenate gene and protein ID to compare 
  group_by(gene,protein_id,experiment) %>%
  mutate(gene_prot_comb = paste(protein_id, gene, sep = "_")) %>%
  # get distinct combos
  ungroup() %>%
  distinct(gene,gene_prot_comb, Species) %>%
  # any time where there is more than one gene row means there are two different transcript combos
  group_by(gene, Species) %>%
  summarize(count = n()) %>% arrange(Species, desc(count)) %>% View()

## Are any genes both cont and DEG?
LFC_cont_comb_gene_shared_LFC_cont_gene <- LFC_cont_comb_gene %>% distinct(Species, gene, Data_Type) %>% group_by(gene, Species) %>% filter(n() > 1) %>% ungroup() %>% distinct(gene, Species)
LFC_cont_comb_gene %>% distinct(Species, gene, Data_Type) %>% group_by(gene, Species) %>% filter(n() > 1) %>% ungroup() %>% distinct(gene, Species) %>% count(Species)
#Species                   n
#<chr>                 <int>
#  1 Crassostrea_gigas         5
#2 Crassostrea_virginica    10

# Which transcripts do these represent
LFC_cont_comb_gene_shared_LFC_cont_gene_prot <- left_join(LFC_cont_comb_gene_shared_LFC_cont_gene, IAP_domain_structure_no_dup_rm) %>% 
  # also join which proteins were cont vs. LFC for each gene
  left_join(., unique(LFC_cont_comb_gene[,c("protein_id","Data_Type", "experiment")]))

# how many const transcripts vs. LFC transcripts used 

#### PLOT IAP FULL PROTEIN TREE WITH BIR MSA AND IAP TABLE ####
# tutorial for review #https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html
# object with protein tree
IAP_raxml_treedata_circular_product 

# Subset tree to get an example of each domain type in C. virginica and C. gigas
BIR_IAP_all_MSA_subset <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_model_prot_IAP_prot_BIR_seq_MSA.fa")
BIR_IAP_raxml_tibble_figure_subset <- BIR_IAP_raxml_tibble %>% filter(Species == "Crassostrea_gigas" | Species == "Crassostrea_virginica") %>% distinct(Type, Species, .keep_all = TRUE) %>% 
  # remove the "unique" types
  filter(Type != "Unique") %>% arrange(Type)
BIR_IAP_raxml_tibble_figure_subset_label <- BIR_IAP_raxml_tibble_figure_subset$label

# Make table to go alongside MSA
BIR_IAP_raxml_tibble_figure_subset_table <- BIR_IAP_raxml_tibble_figure_subset %>%
  dplyr::select(Type) %>%
  mutate(Type = case_when(
    Type == "Non_Zinc_binding"~ "Non Zinc-binding",
    Type == "T1"              ~ "T1"              ,
    Type == "T1_like_3"       ~ "T1-like-3"       ,
    Type == "T1-like_1"       ~ "T1-like-1"       ,
    Type == "T1-like_2"       ~ "T1-like-2"       ,
    Type == "T2"              ~ "T2"              ,
    Type == "T2_like_3"       ~ "T2-like-3"       ,
    Type == "T2_like_4"       ~ "T2-like-4"       ,
    Type == "T2-like_1"       ~ "T2-like-1"       ,
    Type == "T2-like_2"       ~ "T2-like-2"       ,
    Type == "TX"~ "TX",
    Type == "TY"  ~ "TY")) %>%
    #Species = case_when(
    #  Species == "Crassostrea_gigas" ~ "C. gigas",
    #  Species == "Crassostrea_virginica" ~ "C. virginica")) %>% 
  gt::gt() %>%
  cols_label(Type = md("**BIR Type**")) %>%
  tab_options(table.font.size = 20,
              table.font.color = "black")
# save as png
gtsave(BIR_IAP_raxml_tibble_figure_subset_table, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/BIR_IAP_raxml_tibble_figure_subset_table.png")
# this table does not look good though next to the MSA, not going to use for figure 
# need to change the labels to be the Type 1 and 2 designations

# subset MSA
class(BIR_IAP_all_MSA)
BIR_IAP_all_MSA_figure_subset <- BIR_IAP_all_MSA_subset[match(rev(BIR_IAP_raxml_tibble_figure_subset_label), BIR_IAP_all_MSA_subset[,1]),]
colnames(BIR_IAP_raxml_tibble_figure_subset)[4] <- "seq.name"
# make label the type info
BIR_IAP_all_MSA_figure_subset_join <- left_join(BIR_IAP_all_MSA_figure_subset, BIR_IAP_raxml_tibble_figure_subset[,c("seq.name","Type","Species")], by = "seq.name")
BIR_IAP_all_MSA_figure_subset_join <- BIR_IAP_all_MSA_figure_subset_join[-1]
BIR_IAP_all_MSA_figure_subset_join <- BIR_IAP_all_MSA_figure_subset_join %>% 
  dplyr::mutate(Species = case_when(
      Species == "Crassostrea_gigas" ~ "C. gig.",
      Species == "Crassostrea_virginica" ~ "C. vir."),
      Type = case_when(
        Type == "Non_Zinc_binding"~ "Non Zinc-binding",
        Type == "T1"              ~ "T1"              ,
        Type == "T1_like_3"       ~ "T1-like-3"       ,
        Type == "T1-like_1"       ~ "T1-like-1"       ,
        Type == "T1-like_2"       ~ "T1-like-2"       ,
        Type == "T2"              ~ "T2"              ,
        Type == "T2_like_3"       ~ "T2-like-3"       ,
        Type == "T2_like_4"       ~ "T2-like-4"       ,
        Type == "T2-like_1"       ~ "T2-like-1"       ,
        Type == "T2-like_2"       ~ "T2-like-2"       ,
        Type == "TX"~ "TX",
        Type == "TY"  ~ "TY"), 
      seq.names = paste(Type, Species, sep = ", "))
BIR_IAP_all_MSA_figure_subset_join <- BIR_IAP_all_MSA_figure_subset_join[,c(4,1)]

# export back to then reload in correct order
dat2fasta(BIR_IAP_all_MSA_figure_subset_join, outfile ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_all_MSA_figure_subset.fa")

# reload as AAMultipleAlignment object 
BIR_IAP_all_MSA_figure_subset <- Biostrings::readAAMultipleAlignment("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_all_MSA_figure_subset.fa", format = "fasta")

# Plot MSA with ggmsa
BIR_IAP_all_MSA_figure_subset_msa <- ggmsa(BIR_IAP_all_MSA_figure_subset, start = 53, end = 85, 
      color = "Zappo_AA",  # Zappo colors by amino acid chemical characteristics 
      none_bg = TRUE, # keeps only the letters and not the full color background
      posHighligthed = c(57,60, 67, 76, 77,80,82,84), # specify specific positions to highlight in the alignment 
      seq_name = TRUE) +
     # increase text size
   theme(text = element_text(size=10),
         plot.margin = unit(c(0, 0, 0, 0), "cm")) # remove the y axis which just shows the counts of sequences 
# checked and the order is correct 

## Plot MSA and tree
# align with plot_grid 
BIR_IAP_all_MSA_figure_subset_grid <- plot_grid(IAP_raxml_treedata_circular_product , BIR_IAP_all_MSA_figure_subset_msa, ncol=1,
                               rel_heights = c(1,0.2))

# Because gt tables can't be converted into grobs, add those in using Inkscape
ggsave(filename = "BIR_IAP_all_MSA_figure_subset_grid.tiff", plot=BIR_IAP_all_MSA_figure_subset_grid, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_full_tree_MSA",
       width = 13 ,
       height = 13,
       units = "in",
       dpi=300)

#### PLOT FULL GIMAP PROTEIN TREE ####
# Load and parse RAxML bipartitions bootstrapping file with treeio. File input is the bootstrapping analysis output
GIMAP_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.AIG_GIMAP_HMMER_Interpro_XP_list_all_MSA_RaxML")
GIMAP_raxml

# Get information in the features/attributes of the tree (the XP labels) with get_data
get.fields(GIMAP_raxml)
get.data(GIMAP_raxml)

# Convert to tibble tree dataframe object with tidytree to add external data
GIMAP_raxml_tibble <- as_tibble(GIMAP_raxml)

# Join protein product name,gene or locus, and species
colnames(GIMAP_raxml_tibble)[4] <- "protein_id"
GIMAP_raxml_tibble <- left_join(GIMAP_raxml_tibble, AIG1_XP_ALL_gff_GIMAP_species_join, by = "protein_id")
colnames(GIMAP_raxml_tibble)[4] <- "label"

# Add combined gene and locus name column 
GIMAP_raxml_tibble$gene_locus_tag <- coalesce(GIMAP_raxml_tibble$gene, GIMAP_raxml_tibble$locus_tag)

# Convert to treedata object to store tree plus outside data
GIMAP_raxml_treedata <- as.treedata(GIMAP_raxml_tibble)

# save treedata
save(GIMAP_raxml_treedata, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_raxml_treedata.Rdata")

# Plot circular tree
GIMAP_raxml_treedata_circular_product <- ggtree(GIMAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=product,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-130,130)  

GIMAP_raxml_treedata_circular_product + scale_color_discrete(name = "Species", labels = c("Aplysia californica", 
                             "Biomphalaria glabrata", "Crassostrea gigas", "Crassostrea virginica","Elysia chlorotica","Lottia gigantea","Mizuhopecten yessoensis",
                             "Pomacea canaliculata","NA"))


# specific color for each species 
                          
# Plot normal tree to use with heatmap
GIMAP_raxml_treedata_tree_product <- ggtree(GIMAP_raxml_treedata, aes(color=Species)) + 
  geom_tiplab(aes(label=product)) + 
  theme(legend.position = "right", legend.text = element_text(face = "italic"))  

#Subset tree after biomphalaria to zoom in on C. vir and C.gig
GIMAP_raxml_treedata_pomacea_down_subset <- tree_subset(GIMAP_raxml_treedata, 60, levels_back = 7)
ggtree(GIMAP_raxml_treedata_pomacea_down_subset, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=product,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-70,70)  

# PLOT AS GENE TREES TO SEARCH FOR POTENTIAL ARTIFACTS
GIMAP_raxml_treedata_circular_gene <- ggtree(GIMAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=gene_locus_tag,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-80,80)  

GIMAP_raxml_treedata_circular_gene + scale_color_discrete(name = "Species", labels = c("Aplysia californica", 
                                                                                          "Biomphalaria glabrata", "Crassostrea gigas", "Crassostrea virginica","Elysia chlorotica","Lottia gigantea","Mizuhopecten yessoensis",
                                                                                          "Pomacea canaliculata","NA"))

#Subset tree after biomphalaria to zoom in on C. vir and C.gig
GIMAP_raxml_treedata_pomacea_down_subset_genes <- ggtree(GIMAP_raxml_treedata_pomacea_down_subset, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=gene_locus_tag,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-70,70)  

## PLOT GENE TREE OF ONLY THE THREE INTERESTING CVIR CGIG BRANCHES
GIMAP_raxml_treedata_CV_CG_mixed_branches <- tree_subset(GIMAP_raxml_treedata, "XP_021367963.1", levels_back = 3)
GIMAP_raxml_treedata_CV_CG_mixed_branches_tibble <- as.tibble(GIMAP_raxml_treedata_CV_CG_mixed_branches)
GIMAP_raxml_treedata_CV_CG_mixed_branches_tree <- ggtree(GIMAP_raxml_treedata_CV_CG_mixed_branches, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=gene_locus_tag), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

# Susbet for individual clusters
GIMAP_raxml_treedata_CV_CG_mixed_branch_B <- tree_subset(GIMAP_raxml_treedata, "XP_021367963.1", levels_back = 2)
GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble <- as.tibble(GIMAP_raxml_treedata_CV_CG_mixed_branch_B)
GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tree <- ggtree(GIMAP_raxml_treedata_CV_CG_mixed_branch_B, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=gene_locus_tag), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

ggtree(GIMAP_raxml_treedata_CV_CG_mixed_branch_B, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=label), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list <- GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble %>% filter(Species=="Crassostrea_virginica") %>% select(label)

GIMAP_raxml_treedata_CV_CG_mixed_branch_C <- tree_subset(GIMAP_raxml_treedata, "XP_009059109.1", levels_back = 2)
GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble <- as.tibble(GIMAP_raxml_treedata_CV_CG_mixed_branch_C)
GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tree <- ggtree(GIMAP_raxml_treedata_CV_CG_mixed_branch_C, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=gene_locus_tag), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list <- GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble %>% filter(Species=="Crassostrea_virginica") %>% select(label)

GIMAP_raxml_treedata_CV_CG_mixed_branch_D <- tree_subset(GIMAP_raxml_treedata, "XP_021353565.1", levels_back = 2)
GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble <- as.tibble(GIMAP_raxml_treedata_CV_CG_mixed_branch_D)
GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tree <- ggtree(GIMAP_raxml_treedata_CV_CG_mixed_branch_D, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=gene_locus_tag), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-40,40)  

GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list <- GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble %>% filter(Species=="Crassostrea_virginica") %>% select(label)

# Pull out sequences for each branch 
# split seq name and product so I can look up 
AIG_seq_rm_dup_phylo_split <- separate(AIG_seq_rm_dup_phylo, seq.name, into = c("protein_id", "product"), sep = "\\s",
                                       extra = "merge")

colnames(GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list)[1] <- "protein_id" 
colnames(GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list)[1] <- "protein_id"
colnames(GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list)[1] <- "protein_id"

GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list_seq <- left_join(GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list, AIG_seq_rm_dup_phylo_split)
GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list_seq <- left_join(GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list, AIG_seq_rm_dup_phylo_split)
GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list_seq <- left_join(GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list, AIG_seq_rm_dup_phylo_split)

# column names must be seq.name and seq.text, combining back columns
GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list_seq <- unite(GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list_seq, seq.name, sep =" ", 1,2)  
GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list_seq <- unite(GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list_seq, seq.name, sep =" ", 1,2) 
GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list_seq <- unite(GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list_seq, seq.name, sep =" ", 1,2) 

# Export fasta to bluewaves to align with MAFFT (haven't done that yet)
dat2fasta(GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list_seq, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_Artifact_Investigation/GIMAP_raxml_treedata_CV_CG_mixed_branch_B_tibble_XP_list_seq.fa")
dat2fasta(GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list_seq, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_Artifact_Investigation/GIMAP_raxml_treedata_CV_CG_mixed_branch_C_tibble_XP_list_seq.fa")
dat2fasta(GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list_seq, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_Artifact_Investigation/GIMAP_raxml_treedata_CV_CG_mixed_branch_D_tibble_XP_list_seq.fa")

### PLOT FULL MOLLUSC GIMAP TREE WITH THE 1 GENE ARTIFACT REMOVED ###

# Drop tip to remove specific haplotig tip
GIMAP_drop <- c("XP_022296317.1") 
GIMAP_raxml_haplotig_rm_treedata  <- drop.tip(GIMAP_raxml_treedata ,GIMAP_drop)

#Convert to tibble 
GIMAP_raxml_haplotig_rm_tibble <- as_tibble(GIMAP_raxml_haplotig_rm_treedata )

# Remove text after isoform so I can collapse protein names into shorter list
GIMAP_raxml_haplotig_rm_tibble$product <- gsub("isoform.*", "", GIMAP_raxml_haplotig_rm_tibble$product)
GIMAP_raxml_haplotig_rm_tibble$product <- trimws(GIMAP_raxml_haplotig_rm_tibble$product , which = "both")

# Shorten product names by joinined edited excel spreadsheet
#View(unique(GIMAP_raxml_haplotig_rm_tibble$product))
GIMAP_alias <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_shortened_product.csv")

GIMAP_raxml_haplotig_rm_tibble <- left_join(GIMAP_raxml_haplotig_rm_tibble, GIMAP_alias)

# Fill in blanks with uncharacterized locus name
GIMAP_raxml_haplotig_rm_tibble$alias[is.na(GIMAP_raxml_haplotig_rm_tibble$alias)] <- GIMAP_raxml_haplotig_rm_tibble$product[is.na(GIMAP_raxml_haplotig_rm_tibble$alias)]

# Remove uncharacterized protein and just keep gene name for those uncharacterized
GIMAP_raxml_haplotig_rm_tibble$alias <- gsub("uncharacterized protein", "",GIMAP_raxml_haplotig_rm_tibble$alias)

# Convert back to treedata object 
GIMAP_raxml_haplotig_rm_treedata <- as.treedata(GIMAP_raxml_haplotig_rm_tibble)

# plot with bootstrap values
ggtree(GIMAP_raxml_haplotig_rm_treedata, layout="fan", aes(color=Species),  branch.length = "none") + 
  geom_tiplab2(aes(label=alias), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 3, fontface="bold") # allows for subset

## NOTE: Collapse expanded clades is not valid to do on this protein tree, would be valid on the gene tree (which I can run too)

## Calculate total gene counts across all species for each protein name
# Are there genes with two different protein names? find non-matching alias rows in each gene group
GIMAP_raxml_haplotig_rm_tibble_gene_dup <- GIMAP_raxml_haplotig_rm_tibble %>% group_by(gene) %>% filter(n()>1) %>% filter(!is.na(product)) %>% 
  distinct(gene,alias, .keep_all = TRUE) %>% ungroup() %>% group_by(gene) %>% filter(n()>1) # YES there are seven genes across all species (1 in C. gigas, 2 in C. virginica)
# where the different protein isoforms have different product annotations (Weird!)

# remove "-like" to facilitate collapsing protein names
GIMAP_raxml_haplotig_rm_tibble$alias_likerm <- str_remove(GIMAP_raxml_haplotig_rm_tibble$alias, "-like") 

# Keep in the ones that have two different named proteins for a single gene by adding distinct with gene and alias
GIMAP_raxml_haplotig_rm_tibble_gene_type_count <- GIMAP_raxml_haplotig_rm_tibble %>% 
  distinct(gene,alias_likerm, .keep_all = TRUE) %>% # because some dupcliated gene names still 
  group_by(Species, alias_likerm) %>%  
  summarize(gene_count_alias = n())%>% 
  ungroup() %>% #ungroup
  filter(!grepl("LOC",alias_likerm)) 

# Keep in the ones that have two different named proteins for a single gene by adding distinct with gene and alias
GIMAP_raxml_haplotig_rm_tibble_gene_type_count_spread <- GIMAP_raxml_haplotig_rm_tibble %>% 
  distinct(gene,alias_likerm, .keep_all = TRUE) %>% # because some dupcliated gene names still 
  group_by(Species, alias_likerm) %>%  
  summarize(gene_count_alias = n()) %>% 
  ungroup() %>% #ungroup
  filter(!grepl("LOC",alias_likerm)) %>% # remove rows with LOC
  spread(Species, gene_count_alias) # spread by species 

# replace_nas 
GIMAP_raxml_haplotig_rm_tibble_gene_type_count_spread_na <- GIMAP_raxml_haplotig_rm_tibble_gene_type_count_spread %>%
  replace(is.na(.), 0)

# plot as heatmap
ggplot(GIMAP_raxml_haplotig_rm_tibble_gene_type_count, aes(x=Species,y=alias_likerm, fill=gene_count_alias)) + 
  geom_tile() + scale_fill_viridis(discrete=FALSE) +
  labs(title="Gene Counts by Product in Each species",
       y = "Percent of Genes", x="Product")

# plot as columns
ggplot(GIMAP_raxml_haplotig_rm_tibble_gene_type_count, aes(x=Species,y=gene_count_alias, fill=alias_likerm)) + 
  geom_col(position="dodge")

# plot only GIMAP columns
GIMAP_raxml_haplotig_rm_tibble_gene_type_count %>% filter(grepl("GIMAP",alias_likerm)) %>% 
  ggplot(aes(x=Species,y=gene_count_alias, fill=alias)) + 
  geom_col(position="dodge")

GIMAP_raxml_haplotig_rm_prot_collapse_tibble_gene_type_count %>% filter(grepl("GIMAP",alias)) %>% 
  ggplot(aes(x=Species,y=gene_count_alias, fill=alias)) + 
  geom_col(position="fill")

# Calculate frequency of genes with particular protein annotations
GIMAP_raxml_haplotig_rm_tibble_gene_type_count_perspecies_freq <- GIMAP_raxml_haplotig_rm_tibble_gene_type_count %>%
  group_by(Species) %>% mutate(gene_in_species_total = sum(gene_count_alias)) %>% mutate(gene_prot_percent = (gene_count_alias / gene_in_species_total)*100) %>%
  filter(!is.na(Species)) 

# Plot percent (remember this still includes 7 genes with double hits in the list)
GIMAP_raxml_haplotig_rm_tibble_gene_type_count_perspecies_freq %>% 
  # keep only the named ones so I can look at diversity
  filter(grepl("IAN",alias_likerm) | grepl("GIMAP", alias_likerm)) %>%
  ggplot(aes(x=alias_likerm,y=gene_prot_percent, fill=Species)) + 
  geom_col(position="dodge") + coord_flip() +
  labs(title="Gene Frequency by Product in Each species",
       y = "Percent of Genes", x="Product")

GIMAP_raxml_haplotig_rm_tibble_gene_type_count_perspecies_freq %>% 
  # keep only the named ones so I can look at diversity
  filter(grepl("IAN",alias_likerm) | grepl("GIMAP", alias_likerm)) %>%
  ggplot(aes(x=Species,y=gene_prot_percent, fill=alias_likerm)) + 
  geom_col(position="dodge")  +
  labs(title="Gene Frequency by Product in Each species",
       y = "Percent of Genes", x="Product")

# Find protein with highest percent in each organism 
GIMAP_raxml_haplotig_rm_tibble_gene_type_count_perspecies_freq %>% 
  group_by(Species) %>% top_n(n=1, wt = gene_prot_percent) %>%
  ggplot(aes(x=alias_likerm,y=gene_prot_percent, fill=Species)) + 
  geom_col(position="dodge") + 
  labs(title="Most Abundant Product in Each species",
       y = "Percent of Genes", x="Product")

### PLOT FULL MOLLUSC GIMAP TREE WITH THE 1 GENE ARTIFACT REMOVED AND GENE COUNTS COLLAPSED ###
# Need to check the code below.. doing this is problematic because some of the genes have two differently annotated proteins 
# Identify protein isoforms from the same gene within a parent branch to be dropped with droptip 
GIMAP_raxml_tibble_haplotig_rm_prot_to_keep <- distinct_at(GIMAP_raxml_tibble, vars("parent","gene"), .keep_all = TRUE) 
GIMAP_raxml_tibble_haplotig_rm_prot_to_keep %>% group_by(gene_locus_tag) %>% filter(n()>1) %>% filter(!is.na(gene_locus_tag)) %>% View()
GIMAP_raxml_tibble_haplotig_rm_prot_to_keep_label <- GIMAP_raxml_tibble_haplotig_rm_prot_to_keep %>% filter(!is.na(label)) %>% select(label)
GIMAP_raxml_tibble_label <- GIMAP_raxml_tibble %>% filter(!is.na(label)) %>% select(label) 
GIMAP_raxml_tibble_haplotig_rm_prot_to_drop <- setdiff(GIMAP_raxml_tibble_label$label, GIMAP_raxml_tibble_haplotig_rm_prot_to_keep_label$label)

# Drop extra protein isoforms in parent branches and haplotig from the original GIMAP tree_data
GIMAP_raxml_haplotig_rm_prot_collapse_treedata <- drop.tip(GIMAP_raxml_treedata, GIMAP_raxml_tibble_haplotig_rm_prot_to_drop)

# Drop tip to remove specific haplotig tip
GIMAP_drop <- c("XP_022296317.1") 
GIMAP_raxml_haplotig_rm_prot_collapse_treedata  <- drop.tip(GIMAP_raxml_haplotig_rm_prot_collapse_treedata ,GIMAP_drop)

#Convert to tibble 
GIMAP_raxml_haplotig_rm_prot_collapse_tibble <- as_tibble(GIMAP_raxml_haplotig_rm_prot_collapse_treedata )

# Remove text after isoform so I can collapse protein names into shorter list
GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product <- gsub("isoform.*", "", GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product)
GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product <- trimws(GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product , which = "both")

# Shorten product names by joinined edited excel spreadsheet
#View(unique(GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product))
GIMAP_alias <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_shortened_product.csv")

GIMAP_raxml_haplotig_rm_prot_collapse_tibble <- left_join(GIMAP_raxml_haplotig_rm_prot_collapse_tibble, GIMAP_alias)

# Fill in blanks with uncharacterized locus name
GIMAP_raxml_haplotig_rm_prot_collapse_tibble$alias[is.na(GIMAP_raxml_haplotig_rm_prot_collapse_tibble$alias)] <- GIMAP_raxml_haplotig_rm_prot_collapse_tibble$product[is.na(GIMAP_raxml_haplotig_rm_prot_collapse_tibble$alias)]

# Remove uncharacterized protein and just keep gene name for those uncharacterized
GIMAP_raxml_haplotig_rm_prot_collapse_tibble$alias <- gsub("uncharacterized protein", "",GIMAP_raxml_haplotig_rm_prot_collapse_tibble$alias)

# Convert back to treedata object 
GIMAP_raxml_haplotig_rm_prot_collapse_treedata <- as.treedata(GIMAP_raxml_haplotig_rm_prot_collapse_tibble)

# Plot circular tree with protein list collapsed by shared genes
ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=label,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-70,70)  

# Plot circular tree with protein name aliases
ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=alias,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-70,70) +
  geom_text(aes(label = bootstrap), hjust = 2, vjust= -3, size = 2)

# Plot as vertical tree with protein name aliases and bootstrap values 
ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, aes(color=Species), branch.length = "none") + 
  geom_tiplab(aes(label=alias), size =2.0, offset=.5) + # geom_tiplab1 for vertical trees
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-70,70) +
  geom_text(aes(label = bootstrap), hjust = 3, vjust = -0.2, size = 2)

# Plot fan tree with protein name aliases
GIMAP_protein_haplotig_rm_fantree <- ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, layout="fan", aes(color=Species),  branch.length = "none") + 
  geom_tiplab2(aes(label=alias), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
 geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.2, size = 3, fontface="bold")  

GIMAP_protein_haplotig_rm_fantree_subset_boot_50 <- ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, layout="fan", aes(color=Species),  branch.length = "none") + 
  geom_tiplab2(aes(label=alias), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 3, fontface="bold") # allows for subset

ggtree(GIMAP_raxml_haplotig_rm_prot_collapse_treedata, layout="fan", aes(color=Species),  branch.length = "none") + 
  geom_tiplab2(aes(label=gene_locus_tag), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 3, fontface="bold")

#### PLOT GIMAP MY, CV, CG PROTEIN TREE ####

# Load and parse RAxML bipartitions bootstrapping file with treeio. File input is the bootstrapping analysis output
GIMAP_MY_CV_CG_raxml <- read.raxml(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/RAxML/RAxML_bipartitionsBranchLabels.AIG1_dup_seq_rm_kept_haplotig_collapsed_MY_CV_CG_MSA_RaxML")
GIMAP_MY_CV_CG_raxml

# Convert to tibble tree dataframe object with tidytree to add external data
GIMAP_MY_CV_CG_raxml_tibble <- as_tibble(GIMAP_MY_CV_CG_raxml)

# Join protein product name,gene or locus, and species
colnames(GIMAP_MY_CV_CG_raxml_tibble)[4] <- "protein_id"
GIMAP_MY_CV_CG_raxml_tibble <- left_join(GIMAP_MY_CV_CG_raxml_tibble, AIG1_XP_ALL_gff_GIMAP_species_join, by = "protein_id")
colnames(GIMAP_MY_CV_CG_raxml_tibble)[4] <- "label"

# Add combined gene and locus name column 
GIMAP_MY_CV_CG_raxml_tibble$gene_locus_tag <- coalesce(GIMAP_MY_CV_CG_raxml_tibble$gene, GIMAP_MY_CV_CG_raxml_tibble$locus_tag)

# Remove text after isoform so I can collapse protein names into shorter list
GIMAP_MY_CV_CG_raxml_tibble$product <- gsub("isoform.*", "", GIMAP_MY_CV_CG_raxml_tibble$product)
GIMAP_MY_CV_CG_raxml_tibble$product <- trimws(GIMAP_MY_CV_CG_raxml_tibble$product , which = "both")

# Join with alias info
GIMAP_MY_CV_CG_raxml_tibble <- left_join(GIMAP_MY_CV_CG_raxml_tibble, GIMAP_alias)

# Fill in blanks with uncharacterized locus name
GIMAP_MY_CV_CG_raxml_tibble$alias[is.na(GIMAP_MY_CV_CG_raxml_tibble$alias)] <- GIMAP_MY_CV_CG_raxml_tibble$product[is.na(GIMAP_MY_CV_CG_raxml_tibble$alias)]

# Remove uncharacterized protein and just keep gene name for those uncharacterized
GIMAP_MY_CV_CG_raxml_tibble$alias <- gsub("uncharacterized protein ", "",GIMAP_MY_CV_CG_raxml_tibble$alias)

# fill species NA with a value
GIMAP_MY_CV_CG_raxml_tibble[is.na(GIMAP_MY_CV_CG_raxml_tibble$Species), ] <- as.character("none")

# Convert to treedata object to store tree plus outside data
GIMAP_MY_CV_CG_raxml_treedata <- as.treedata(GIMAP_MY_CV_CG_raxml_tibble)

# Plot fan tree
ggtree(GIMAP_MY_CV_CG_raxml_treedata, aes(color=Species), layout="fan",  branch.length = "none") + 
  geom_tiplab2(aes(label=alias), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 3, fontface="bold") # allows for subset

# Plot vertical tree and edit colors
GIMAP_MY_CV_CG_raxml_treedata_vertical <- 
  ggtree(GIMAP_MY_CV_CG_raxml_treedata, aes(color=Species, fill=Species), branch.length = "none") + 
  geom_tiplab(aes(label=alias), size =2.2, offset=0) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  geom_text2(aes(label=bootstrap, subset = as.numeric(bootstrap) > 50), hjust = 1, vjust = -0.2, size = 1.8, fontface="bold") + # allows for subset
  xlim(-70,31.8) + #change scaling so branch lengths are smaller and all alias labels are showing
  scale_colour_manual(name = "Species", values=c("#0a8707","#6a70d8", "#c55d32"), na.value="grey46", breaks=c("Crassostrea_gigas", "Crassostrea_virginica","Mizuhopecten_yessoensis"),
                    labels = c("Crassostrea gigas", "Crassostrea virginica","Mizuhopecten yessoensis")) +
  guides(col = guide_legend(ncol =1, title.position = "top",  override.aes = aes(label = "")))

#### PLOT GIMAP DOMAIN STRUCTURE AND COMBINE WITH TREE ####
# Use combination of geom_segment and geom_rect and combine plot with vertical tree using ggarrange from ggpubr
# Get only the Interproscan domains for my proteins of interest
GIMAP_MY_CV_CG_raxml_tibble_join <- GIMAP_MY_CV_CG_raxml_tibble %>% filter(!is.na(label)) # remove rows with just bootstrap information
colnames(GIMAP_MY_CV_CG_raxml_tibble_join)[4] <- "protein_id"
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains <-  left_join(GIMAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")], AIG1_XP_ALL_gff_GIMAP)
#AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains %>% 
  filter(grepl("InterPro:IPR", Dbxref) | source == "Coils") # keep Interproscan domain lines and coiled coil lines 

# REMOVE THIS SUBSETTING 


AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains %>% 
  filter(is.na(source))

nrow(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot %>% filter(is.na(source))) # 144
nrow(GIMAP_MY_CV_CG_raxml_tibble %>% filter(!is.na(label))) # 144 - they agree, all proteins were found 

AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$Dbxref[AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$Dbxref =="character(0)"] <- "Coil"
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only  <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only %>% unnest(Dbxref)

# count most common
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only %>% group_by(Dbxref) %>% count() %>% View()

# Get the node order from original GIMAP tree
GIMAP_MY_CV_CG_raxml_treedata_tip  <- fortify(GIMAP_MY_CV_CG_raxml_treedata)
GIMAP_MY_CV_CG_raxml_treedata_tip = subset(GIMAP_MY_CV_CG_raxml_treedata_tip, isTip)
GIMAP_MY_CV_CG_raxml_treedata_tip_order <- GIMAP_MY_CV_CG_raxml_treedata_tip$label[order(GIMAP_MY_CV_CG_raxml_treedata_tip$y, decreasing=TRUE)]

# Reorder the protein and polygon
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot[match(GIMAP_MY_CV_CG_raxml_treedata_tip_order, AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot$protein_id),]
GIMAP_MY_CV_CG_raxml_treedata_tip_order <- as.data.frame(GIMAP_MY_CV_CG_raxml_treedata_tip_order)
colnames(GIMAP_MY_CV_CG_raxml_treedata_tip_order)[1] <- "protein_id"
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only <- full_join(GIMAP_MY_CV_CG_raxml_treedata_tip_order, AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only)

# Add polygon height
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID  <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only  %>% distinct(protein_id) 
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID <- AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID %>% 
  mutate(height_start = rev(as.numeric(row.names(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID )) - 0.25)) %>%
  mutate(height_end = rev(as.numeric(row.names(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID)) + .5))

# Join back in height
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only <- left_join(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only , AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only_ID )

# Set factor level order of the nodes set levels in reverse order
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$node <- factor(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$node, levels = unique(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$node))
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$Dbxref <- factor(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$Dbxref, levels = unique(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only$Dbxref))
AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot$node <- factor(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot$node, levels = rev(AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot$node))

# Plot the line segments of the NA source lines (which have the full protein start and end)
GIMAP_Interproscan_domain_plot <- ggplot() + 
  # plot length of each protein as line
  geom_segment(data =AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot,
               aes(x=as.numeric(start), xend=as.numeric(end), y=node, yend=node), color = "grey") +
  # add boxes with geom_rect 
  geom_rect(data=AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_only,
            aes(xmin=start, xmax=end, ymin=height_start, ymax=height_end, fill= Dbxref)) +
  #add labels
  labs(y = NULL, x = "Protein Domain Position (aa)") +
  # add text labels
  #geom_text(data=AIG1_XP_ALL_gff_GIMAP_Interpro_Domains_fullprot,aes(x= end, y = node, label=alias),
  #          size=2.2, hjust=-.15, check_overlap = TRUE) + 
  # text theme
  theme_bw() + 
  # plot theme
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size=8, family="sans"),
        legend.title = element_text(size=12, family="sans")) +
  # Change y axis ticks, turn off expand so that extra white space is removed
  scale_x_continuous(breaks=c(0,500,1000), expand = c(0,0)) + 
  # Change domain labels 
  scale_fill_manual(values=c("#6d83da","#49b9d3","#c24c6e","#a68742","#bab237","#9bad47","#c2464c","#b35535","#be6ec6","#568538","#45c097","#892863","#d56cad","#cb8130","#57398c","#62c36f"), 
                    name="Functional Domains",
                    breaks=c("\"InterPro:IPR006703\"", "Coil", "\"InterPro:IPR027417\"", "\"InterPro:IPR029071\"", "\"InterPro:IPR013783\"",
                             "\"InterPro:IPR036179\"", "\"InterPro:IPR007110\"", "\"InterPro:IPR003598\"", "\"InterPro:IPR013151\"", "\"InterPro:IPR003599\"",
                             "\"InterPro:IPR001876\"", "\"InterPro:IPR036443\"", "\"InterPro:IPR013761\"", "\"InterPro:IPR001660\"", "\"InterPro:IPR011029\"",
                             "\"InterPro:IPR001315\""),
                    labels=c("AIG1-type guanine nucleotide-binding (G) domain","Coil","P-loop containing nucleoside triphosphate hydrolase",
                             "Ubiquitin-like domain superfamily","Immunoglobulin-like fold","Immunoglobulin-like domain superfamily",
                             "Immunoglobulin-like domain","Immunoglobulin subtype 2","Immunoglobulin","Immunoglobulin subtype","Zinc finger, RanBP2-type",
                             "Zinc finger, RanBP2-type superfamily","Sterile alpha motif/pointed domain superfamily","Sterile alpha motif domain",
                             "Death-like domain superfamily","CARD domain")) +
  # change number of legend columns and put the legend title on top
  guides(fill=guide_legend(ncol=3, title.position="top"))

# Plot Tree and domains together
# rescale the ylim of the tree to the ylim of the domain plot using ggtree ylim2 
p2 <- GIMAP_MY_CV_CG_raxml_treedata_vertical +ylim2(GIMAP_Interproscan_domain_plot)

# Use cowplot to align the plots
GIMAP_tr_dom_plot <- cowplot::plot_grid(p2, GIMAP_Interproscan_domain_plot, ncol = 2,
                   align = "h", axis="b")

#### PLOT GIMAP TREES WITH DESEQ2 AND DOMAIN INFORMATION ####
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_vir_apop_LFC_GIMAP.Rdata")
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_gig_apop_LFC_GIMAP.Rdata")

# Keep tables separate for plotting 
C_vir_apop_LFC_GIMAP$Species <- "Crassostrea_virginica"
C_gig_apop_LFC_GIMAP$Species <- "Crassostrea_gigas"

# Full Join the list of transcripts so the missing labels are in there and can be ordered for plotting 
C_gig_apop_LFC_GIMAP_full_XP <- full_join(C_gig_apop_LFC_GIMAP, GIMAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])
# No sig GIMAPs have an NA in the tibble tree meaning that none of those proteins were ones that were collapsed as duplicates in CD-Hit

C_vir_apop_LFC_GIMAP_full_XP <- full_join(C_vir_apop_LFC_GIMAP, GIMAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])
# For NAs for GIMAP proteins, meaning that those were exactly identical to another protein and were collapsed. Need to replace these NA's with their uncollapsed parent protein
C_vir_apop_LFC_GIMAP_full_XP_collapsed <- C_vir_apop_LFC_GIMAP_full_XP %>% filter(is.na(node) & !is.na(log2FoldChange)) 
C_vir_apop_LFC_GIMAP_full_XP_collapsed_AIG_seq_rm_dup_clstr6 <- AIG_seq_rm_dup_clstr6 %>% filter(protein_id %in% (C_vir_apop_LFC_GIMAP_full_XP_collapsed$protein_id))
# Find parent proteins in these clusters
C_vir_apop_LFC_GIMAP_full_XP_collapsed_AIG_seq_rm_dup_clstr6_cluster <- AIG_seq_rm_dup_clstr6[AIG_seq_rm_dup_clstr6$cluster %in% C_vir_apop_LFC_GIMAP_full_XP_collapsed_AIG_seq_rm_dup_clstr6$cluster,]

# Recode original DESeq2 with correct protein name (they are all exactly the same length, same sequence and from the same gene. 
    # They are likely just erroneously called isoforms )
C_vir_apop_LFC_GIMAP$protein_id <- recode(C_vir_apop_LFC_GIMAP$protein_id, "XP_022301588.1" = "XP_022301589.1", "XP_022295278.1"="XP_022295277.1")

# Rejoin the Full list of transcript and check for fixed NA
C_vir_apop_LFC_GIMAP_full_XP <- full_join(C_vir_apop_LFC_GIMAP, GIMAP_MY_CV_CG_raxml_tibble_join[,c("protein_id","node","alias")])

# Reorder both to be the order of the GIMAP tree XPs
# Get the node order from original GIMAP tree (done in code chunk regarding domain information above)
GIMAP_MY_CV_CG_raxml_treedata_tip_order 

# Reorder the proteins
C_vir_apop_LFC_GIMAP_full_XP <- full_join(GIMAP_MY_CV_CG_raxml_treedata_tip_order, C_vir_apop_LFC_GIMAP_full_XP)
C_gig_apop_LFC_GIMAP_full_XP <- full_join(GIMAP_MY_CV_CG_raxml_treedata_tip_order, C_gig_apop_LFC_GIMAP_full_XP)

# Set factor level order of the nodes set levels in reverse order
C_vir_apop_LFC_GIMAP_full_XP$protein_id <- factor(C_vir_apop_LFC_GIMAP_full_XP$protein_id, levels = rev(unique(C_vir_apop_LFC_GIMAP_full_XP$protein_id)))
C_gig_apop_LFC_GIMAP_full_XP$protein_id <- factor(C_gig_apop_LFC_GIMAP_full_XP$protein_id, levels = rev(unique(C_gig_apop_LFC_GIMAP_full_XP$protein_id)))
C_vir_apop_LFC_GIMAP_full_XP$node <- factor(C_vir_apop_LFC_GIMAP_full_XP$node, levels = rev(unique(C_vir_apop_LFC_GIMAP_full_XP$node)))
C_gig_apop_LFC_GIMAP_full_XP$node <- factor(C_gig_apop_LFC_GIMAP_full_XP$node, levels = rev(unique(C_gig_apop_LFC_GIMAP_full_XP$node)))

# Plot LFC data
C_vir_apop_LFC_GIMAP_tile_plot <- ggplot(C_vir_apop_LFC_GIMAP_full_XP, aes(x=group_by_sim, y = protein_id, fill=log2FoldChange, na.rm= TRUE)) + 
  geom_tile()  + 
  #scale_fill_viridis_c(breaks = seq(min(C_vir_apop_LFC_GIMAP_full_XP$log2FoldChange, na.rm = TRUE),max(C_vir_apop_LFC_GIMAP_full_XP$log2FoldChange, na.rm=TRUE),length.out = 15), 
  #                   option="plasma", guide=guide_legend()) +
  scale_fill_viridis_c(name = "Log2 Fold Change", breaks = c(-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10), 
                       option="plasma", guide=guide_legend(), na.value = "transparent") +
  labs(x="Treatment", y =NULL) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x.top = element_text(size=8, family="sans"),
        #axis.text.y.left = element_text(family ="sans"),
        axis.title = element_text(size=12, family="sans"),
        legend.position = "bottom",
        legend.title = element_text(size=12, family="sans"), 
        legend.text = element_text(size=8, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x = element_line(size=0.2, color="gray"),
        panel.grid.major.y = element_line(size=0.2, color="gray")) +
  # remove NA row from the list 
  scale_y_discrete(limits=c("XP_022295827.1 ","XP_022296700.1 ","XP_022296665.1", "XP_022298533.1", "XP_022301840.1", "XP_011422393.2", "XP_011428816.1", "XP_011455059.1", "XP_011455056.1", "XP_011455058.1", "XP_011441429.1",
                            "XP_011456413.1 ","XP_019923700.1 ","XP_019917937.1", "XP_022302920.1", "XP_022297620.1", "XP_022297618.1", "XP_022297619.1", "XP_019921664.1", "XP_011453355.1", "XP_022300251.1", "XP_021343321.1", "XP_022303234.1",
                            "XP_011422350.2 ","XP_011421513.2 ","XP_011437597.1", "XP_011434056.2", "XP_011414633.1", "XP_011432416.1", "XP_011432417.1", "XP_019919146.1", "XP_022304619.1", "XP_022304621.1", "XP_022304617.1", "XP_022304615.1",
                            "XP_022304620.1 ","XP_022304616.1 ","XP_022304618.1", "XP_022304624.1", "XP_022307700.1", "XP_022304628.1", "XP_022304626.1", "XP_022303805.1", "XP_022303804.1", "XP_022303800.1", "XP_022297006.1", "XP_022298015.1",
                            "XP_022298016.1 ","XP_022299409.1 ","XP_022299405.1", "XP_022299418.1", "XP_022299408.1", "XP_022299416.1", "XP_022299407.1", "XP_022299410.1", "XP_022296312.1", "XP_011431833.1", "XP_019925442.1", "XP_022301589.1",
                            "XP_022301590.1 ","XP_022296314.1 ","XP_022301668.1", "XP_022301132.1", "XP_022301131.1", "XP_022301812.1", "XP_022301145.1", "XP_022300526.1", "XP_022297635.1", "XP_022335699.1", "XP_022302183.1", "XP_022310198.1",
                            "XP_022304950.1 ","XP_022304158.1 ","XP_022291053.1", "XP_022291807.1", "XP_022290458.1", "XP_022291746.1", "XP_021376515.1", "XP_021357452.1", "XP_011439007.1", "XP_022299725.1", "XP_011422444.1", "XP_022302477.1",
                            "XP_022302475.1 ","XP_022301435.1 ","XP_022301436.1", "XP_021364561.1", "XP_021356527.1", "XP_021352274.1", "XP_021356530.1", "XP_021359550.1", "XP_021359549.1", "XP_021352273.1", "XP_022302348.1", "XP_022302161.1",
                            "XP_022339582.1 ","XP_011420627.1 ","XP_011414635.1", "XP_022332608.1", "XP_022332192.1", "XP_011456539.1", "XP_011443783.1", "XP_021365591.1", "XP_021344399.1", "XP_021365588.1", "XP_021365603.1", "XP_021365598.1",
                            "XP_021365602.1", "XP_021370037.1 ","XP_021370077.1", "XP_021370080.1", "XP_021370085.1", "XP_021370086.1", "XP_021370081.1", "XP_021370087.1", "XP_021370075.1", "XP_021370082.1", "XP_021370084.1", "XP_021370078.1",
                            "XP_021370079.1", "XP_021370083.1 ","XP_021370088.1", "XP_021353565.1", "XP_011449219.1", "XP_022301090.1", "XP_022296099.1", "XP_022301089.1", "XP_022296137.1", "XP_011449218.1", "XP_022299663.1", "XP_022301104.1",
                            "XP_022299662.1", "XP_022296100.1 ","XP_021367963.1", "XP_021367951.1", "XP_022292453.1", "XP_022291927.1", "XP_011427681.1", "XP_011424928.1", "XP_022295277.1", "XP_022295280.1", "XP_022295281.1", "XP_022295274.1",
                            "XP_022295286.1"),
                   labels=c("XP_022295827.1 ","XP_022296700.1 ","XP_022296665.1", "XP_022298533.1", "XP_022301840.1", "XP_011422393.2", "XP_011428816.1", "XP_011455059.1", "XP_011455056.1", "XP_011455058.1", "XP_011441429.1",
                            "XP_011456413.1 ","XP_019923700.1 ","XP_019917937.1", "XP_022302920.1", "XP_022297620.1", "XP_022297618.1", "XP_022297619.1", "XP_019921664.1", "XP_011453355.1", "XP_022300251.1", "XP_021343321.1", "XP_022303234.1",
                            "XP_011422350.2 ","XP_011421513.2 ","XP_011437597.1", "XP_011434056.2", "XP_011414633.1", "XP_011432416.1", "XP_011432417.1", "XP_019919146.1", "XP_022304619.1", "XP_022304621.1", "XP_022304617.1", "XP_022304615.1",
                            "XP_022304620.1 ","XP_022304616.1 ","XP_022304618.1", "XP_022304624.1", "XP_022307700.1", "XP_022304628.1", "XP_022304626.1", "XP_022303805.1", "XP_022303804.1", "XP_022303800.1", "XP_022297006.1", "XP_022298015.1",
                            "XP_022298016.1 ","XP_022299409.1 ","XP_022299405.1", "XP_022299418.1", "XP_022299408.1", "XP_022299416.1", "XP_022299407.1", "XP_022299410.1", "XP_022296312.1", "XP_011431833.1", "XP_019925442.1", "XP_022301589.1",
                            "XP_022301590.1 ","XP_022296314.1 ","XP_022301668.1", "XP_022301132.1", "XP_022301131.1", "XP_022301812.1", "XP_022301145.1", "XP_022300526.1", "XP_022297635.1", "XP_022335699.1", "XP_022302183.1", "XP_022310198.1",
                            "XP_022304950.1 ","XP_022304158.1 ","XP_022291053.1", "XP_022291807.1", "XP_022290458.1", "XP_022291746.1", "XP_021376515.1", "XP_021357452.1", "XP_011439007.1", "XP_022299725.1", "XP_011422444.1", "XP_022302477.1",
                            "XP_022302475.1 ","XP_022301435.1 ","XP_022301436.1", "XP_021364561.1", "XP_021356527.1", "XP_021352274.1", "XP_021356530.1", "XP_021359550.1", "XP_021359549.1", "XP_021352273.1", "XP_022302348.1", "XP_022302161.1",
                            "XP_022339582.1 ","XP_011420627.1 ","XP_011414635.1", "XP_022332608.1", "XP_022332192.1", "XP_011456539.1", "XP_011443783.1", "XP_021365591.1", "XP_021344399.1", "XP_021365588.1", "XP_021365603.1", "XP_021365598.1",
                            "XP_021365602.1", "XP_021370037.1 ","XP_021370077.1", "XP_021370080.1", "XP_021370085.1", "XP_021370086.1", "XP_021370081.1", "XP_021370087.1", "XP_021370075.1", "XP_021370082.1", "XP_021370084.1", "XP_021370078.1",
                            "XP_021370079.1", "XP_021370083.1 ","XP_021370088.1", "XP_021353565.1", "XP_011449219.1", "XP_022301090.1", "XP_022296099.1", "XP_022301089.1", "XP_022296137.1", "XP_011449218.1", "XP_022299663.1", "XP_022301104.1",
                            "XP_022299662.1", "XP_022296100.1 ","XP_021367963.1", "XP_021367951.1", "XP_022292453.1", "XP_022291927.1", "XP_011427681.1", "XP_011424928.1", "XP_022295277.1", "XP_022295280.1", "XP_022295281.1", "XP_022295274.1",
                            "XP_022295286.1")) +
  scale_x_discrete(limits = c("Hatchery_Probiotic_RI" ,"Lab_RI_6hr" , "Lab_RI_RI_24hr", "Lab_S4_6hr","Lab_S4_24hr", "Lab_RE22" ,
                              "ROD_susceptible_seed","ROD_resistant_seed", "Dermo_Susceptible_36hr", "Dermo_Susceptible_7d", "Dermo_Susceptible_28d","Dermo_Tolerant_36hr",   
                              "Dermo_Tolerant_7d","Dermo_Tolerant_28d" ), 
                   labels= c("Hatchery\n RI" ,"Lab RI 6hr", "Lab RI 24hr", "Lab S4 6hr","Lab S4 24hr", "Lab RE22" ,
                             "ROD Sus.\n seed","ROD Res.\n seed", "Dermo\n Sus. 36hr", "Dermo\n Sus. 7d", "Dermo\n Sus. 28d","Dermo\n Tol. 36hr",   
                             "Dermo\n Tol. 7d","Dermo\n Tol. 28d"), position="top") +
  guides(fill=guide_legend(ncol=2, title.position="top"))

C_gig_apop_LFC_GIMAP_tile_plot <- ggplot(C_gig_apop_LFC_GIMAP_full_XP, aes(x=group_by_sim, y=protein_id, fill=log2FoldChange)) + 
  geom_tile() + 
  #scale_fill_viridis_c(breaks = seq(min(C_gig_apop_LFC_GIMAP$log2FoldChange, na.rm = TRUE),max(C_gig_apop_LFC_GIMAP$log2FoldChange, na.rm=TRUE),length.out = 15), 
  #                     option="plasma", guide=guide_legend()) +
  scale_fill_viridis_c(name = "Log2 Fold Change", breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), 
                       option="plasma", guide=guide_legend(), na.value = "transparent") +
  labs(x="Treatment", y =NULL) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x.top = element_text(size=8, family="sans"),
        axis.title.x.top = element_text(size=12, family="sans"),
        legend.position = "bottom",
        legend.title = element_text(size=12, family="sans"), 
        legend.text = element_text(size=8, family="sans"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x = element_line(size=0.2, color="gray"),
        panel.grid.major.y = element_line(size=0.2, color="gray")) +
  # remove NA row from the list 
  scale_y_discrete(limits=c("XP_022295827.1 ","XP_022296700.1 ","XP_022296665.1", "XP_022298533.1", "XP_022301840.1", "XP_011422393.2", "XP_011428816.1", "XP_011455059.1", "XP_011455056.1", "XP_011455058.1", "XP_011441429.1",
                            "XP_011456413.1 ","XP_019923700.1 ","XP_019917937.1", "XP_022302920.1", "XP_022297620.1", "XP_022297618.1", "XP_022297619.1", "XP_019921664.1", "XP_011453355.1", "XP_022300251.1", "XP_021343321.1", "XP_022303234.1",
                            "XP_011422350.2 ","XP_011421513.2 ","XP_011437597.1", "XP_011434056.2", "XP_011414633.1", "XP_011432416.1", "XP_011432417.1", "XP_019919146.1", "XP_022304619.1", "XP_022304621.1", "XP_022304617.1", "XP_022304615.1",
                            "XP_022304620.1 ","XP_022304616.1 ","XP_022304618.1", "XP_022304624.1", "XP_022307700.1", "XP_022304628.1", "XP_022304626.1", "XP_022303805.1", "XP_022303804.1", "XP_022303800.1", "XP_022297006.1", "XP_022298015.1",
                            "XP_022298016.1 ","XP_022299409.1 ","XP_022299405.1", "XP_022299418.1", "XP_022299408.1", "XP_022299416.1", "XP_022299407.1", "XP_022299410.1", "XP_022296312.1", "XP_011431833.1", "XP_019925442.1", "XP_022301589.1",
                            "XP_022301590.1 ","XP_022296314.1 ","XP_022301668.1", "XP_022301132.1", "XP_022301131.1", "XP_022301812.1", "XP_022301145.1", "XP_022300526.1", "XP_022297635.1", "XP_022335699.1", "XP_022302183.1", "XP_022310198.1",
                            "XP_022304950.1 ","XP_022304158.1 ","XP_022291053.1", "XP_022291807.1", "XP_022290458.1", "XP_022291746.1", "XP_021376515.1", "XP_021357452.1", "XP_011439007.1", "XP_022299725.1", "XP_011422444.1", "XP_022302477.1",
                            "XP_022302475.1 ","XP_022301435.1 ","XP_022301436.1", "XP_021364561.1", "XP_021356527.1", "XP_021352274.1", "XP_021356530.1", "XP_021359550.1", "XP_021359549.1", "XP_021352273.1", "XP_022302348.1", "XP_022302161.1",
                            "XP_022339582.1 ","XP_011420627.1 ","XP_011414635.1", "XP_022332608.1", "XP_022332192.1", "XP_011456539.1", "XP_011443783.1", "XP_021365591.1", "XP_021344399.1", "XP_021365588.1", "XP_021365603.1", "XP_021365598.1",
                            "XP_021365602.1", "XP_021370037.1 ","XP_021370077.1", "XP_021370080.1", "XP_021370085.1", "XP_021370086.1", "XP_021370081.1", "XP_021370087.1", "XP_021370075.1", "XP_021370082.1", "XP_021370084.1", "XP_021370078.1",
                            "XP_021370079.1", "XP_021370083.1 ","XP_021370088.1", "XP_021353565.1", "XP_011449219.1", "XP_022301090.1", "XP_022296099.1", "XP_022301089.1", "XP_022296137.1", "XP_011449218.1", "XP_022299663.1", "XP_022301104.1",
                            "XP_022299662.1", "XP_022296100.1 ","XP_021367963.1", "XP_021367951.1", "XP_022292453.1", "XP_022291927.1", "XP_011427681.1", "XP_011424928.1", "XP_022295277.1", "XP_022295280.1", "XP_022295281.1", "XP_022295274.1",
                            "XP_022295286.1"),
                   labels=c("XP_022295827.1 ","XP_022296700.1 ","XP_022296665.1", "XP_022298533.1", "XP_022301840.1", "XP_011422393.2", "XP_011428816.1", "XP_011455059.1", "XP_011455056.1", "XP_011455058.1", "XP_011441429.1",
                            "XP_011456413.1 ","XP_019923700.1 ","XP_019917937.1", "XP_022302920.1", "XP_022297620.1", "XP_022297618.1", "XP_022297619.1", "XP_019921664.1", "XP_011453355.1", "XP_022300251.1", "XP_021343321.1", "XP_022303234.1",
                            "XP_011422350.2 ","XP_011421513.2 ","XP_011437597.1", "XP_011434056.2", "XP_011414633.1", "XP_011432416.1", "XP_011432417.1", "XP_019919146.1", "XP_022304619.1", "XP_022304621.1", "XP_022304617.1", "XP_022304615.1",
                            "XP_022304620.1 ","XP_022304616.1 ","XP_022304618.1", "XP_022304624.1", "XP_022307700.1", "XP_022304628.1", "XP_022304626.1", "XP_022303805.1", "XP_022303804.1", "XP_022303800.1", "XP_022297006.1", "XP_022298015.1",
                            "XP_022298016.1 ","XP_022299409.1 ","XP_022299405.1", "XP_022299418.1", "XP_022299408.1", "XP_022299416.1", "XP_022299407.1", "XP_022299410.1", "XP_022296312.1", "XP_011431833.1", "XP_019925442.1", "XP_022301589.1",
                            "XP_022301590.1 ","XP_022296314.1 ","XP_022301668.1", "XP_022301132.1", "XP_022301131.1", "XP_022301812.1", "XP_022301145.1", "XP_022300526.1", "XP_022297635.1", "XP_022335699.1", "XP_022302183.1", "XP_022310198.1",
                            "XP_022304950.1 ","XP_022304158.1 ","XP_022291053.1", "XP_022291807.1", "XP_022290458.1", "XP_022291746.1", "XP_021376515.1", "XP_021357452.1", "XP_011439007.1", "XP_022299725.1", "XP_011422444.1", "XP_022302477.1",
                            "XP_022302475.1 ","XP_022301435.1 ","XP_022301436.1", "XP_021364561.1", "XP_021356527.1", "XP_021352274.1", "XP_021356530.1", "XP_021359550.1", "XP_021359549.1", "XP_021352273.1", "XP_022302348.1", "XP_022302161.1",
                            "XP_022339582.1 ","XP_011420627.1 ","XP_011414635.1", "XP_022332608.1", "XP_022332192.1", "XP_011456539.1", "XP_011443783.1", "XP_021365591.1", "XP_021344399.1", "XP_021365588.1", "XP_021365603.1", "XP_021365598.1",
                            "XP_021365602.1", "XP_021370037.1 ","XP_021370077.1", "XP_021370080.1", "XP_021370085.1", "XP_021370086.1", "XP_021370081.1", "XP_021370087.1", "XP_021370075.1", "XP_021370082.1", "XP_021370084.1", "XP_021370078.1",
                            "XP_021370079.1", "XP_021370083.1 ","XP_021370088.1", "XP_021353565.1", "XP_011449219.1", "XP_022301090.1", "XP_022296099.1", "XP_022301089.1", "XP_022296137.1", "XP_011449218.1", "XP_022299663.1", "XP_022301104.1",
                            "XP_022299662.1", "XP_022296100.1 ","XP_021367963.1", "XP_021367951.1", "XP_022292453.1", "XP_022291927.1", "XP_011427681.1", "XP_011424928.1", "XP_022295277.1", "XP_022295280.1", "XP_022295281.1", "XP_022295274.1",
                            "XP_022295286.1")) +
  scale_x_discrete(limits = c("Zhang_Valg"          ,"Zhang_Vtub"          ,"Zhang_LPS"          , "Rubio_J2_8"          ,"Rubio_J2_9"          ,"Rubio_LGP32"         ,"Rubio_LMG20012T"     ,"He_6hr"             ,
                              "He_12hr"             ,"He_24hr"             ,"He_48hr"            , "He_120hr"            ,"deLorgeril_res_6hr"  ,"deLorgeril_res_12hr" ,"deLorgeril_res_24hr" ,"deLorgeril_res_48hr",
                              "deLorgeril_res_60hr" ,"deLorgeril_res_72hr" ,"deLorgeril_sus_6hr" , "deLorgeril_sus_12hr" ,"deLorgeril_sus_24hr" ,"deLorgeril_sus_48hr" ,"deLorgeril_sus_60hr" ,"deLorgeril_sus_72hr"), 
                   labels= c("Zhang\n V. alg","Zhang\n V.tub\n V. ang","Zhang\n LPS\nM. Lut", "Rubio\nV. crass\n J2_8\n NVir","Rubio\nV. crass\n J2_9\n Vir" ,"Rubio\nV. tasma\n LGP32\n Vir","Rubio\nV. tasma\n LMG20012T\n NVir","He OsHv-1\n 6hr",
                             "He OsHv-1\n 12hr", "He OsHv-1\n24hr", "He OsHv-1\n48hr", "He OsHv-1\n 120hr","deLorgeril\nOsHV-1\n Res. 6hr","deLorgeril\nOsHV-1\n Res. 12hr","deLorgeril\nOsHV-1\n Res. 24hr" ,"deLorgeril\nOsHV-1\n Res. 48hr",
                             "deLorgeril\nOsHV-1\n Res. 60hr","deLorgeril\nOsHV-1\n Res. 72hr" ,"deLorgeril\nOsHV-1\n Sus. 6hr", "deLorgeril\nOsHV-1\n Sus. 12hr","deLorgeril\nOsHV-1\n Sus. 24hr","deLorgeril\nOsHV-1\n Sus. 48hr" ,
                             "deLorgeril\nOsHV-1\n Sus. 60hr","deLorgeril\nOsHV-1\n Sus. 72hr"), position="top") +
  guides(fill=guide_legend(ncol=3, title.position="top"))

## Plot DESeq2 alongside tree and domain structure

# Use cowplot to extract legends and then add separately
GIMAP_Interproscan_domain_plot_legend <- cowplot::get_legend(GIMAP_Interproscan_domain_plot)
GIMAP_Interproscan_domain_plot_no_legend <- GIMAP_Interproscan_domain_plot + theme(legend.position='none')

GIMAP_MY_CV_CG_raxml_treedata_vertical_legend <- cowplot::get_legend(GIMAP_MY_CV_CG_raxml_treedata_vertical)
GIMAP_MY_CV_CG_raxml_treedata_vertical_no_legend <- GIMAP_MY_CV_CG_raxml_treedata_vertical + theme(legend.position='none')
p2_no_legend <- GIMAP_MY_CV_CG_raxml_treedata_vertical_no_legend + aplot::ylim2(GIMAP_Interproscan_domain_plot_no_legend)

C_vir_apop_LFC_GIMAP_tile_plot_legend <- cowplot::get_legend(C_vir_apop_LFC_GIMAP_tile_plot)
C_vir_apop_LFC_GIMAP_tile_plot_no_legend <- C_vir_apop_LFC_GIMAP_tile_plot + theme(legend.position='none')

C_gig_apop_LFC_GIMAP_tile_plot_legend <- cowplot::get_legend(C_gig_apop_LFC_GIMAP_tile_plot)
C_gig_apop_LFC_GIMAP_tile_plot_no_legend <- C_gig_apop_LFC_GIMAP_tile_plot + theme(legend.position='none')

# Now plots are aligned vertically with the legend in one row underneath
#C vir plot
Cvir_GIMAP_tr_dom_LFC <- plot_grid(p2_no_legend, GIMAP_Interproscan_domain_plot_no_legend, C_vir_apop_LFC_GIMAP_tile_plot_no_legend, ncol=3, align='h',
                                   labels = c('A', 'B', 'C'), label_size = 12)
Cvir_GIMAP_tr_dom_LFC_legend <- plot_grid(GIMAP_MY_CV_CG_raxml_treedata_vertical_legend, GIMAP_Interproscan_domain_plot_legend,C_vir_apop_LFC_GIMAP_tile_plot_legend, 
                                          ncol = 3, align="hv")
Cvir_GIMAP_tr_dom_LFC_plus_legend <- plot_grid(Cvir_GIMAP_tr_dom_LFC, Cvir_GIMAP_tr_dom_LFC_legend, ncol=1, rel_heights =c(1, 0.2))

# C gig plots
Cgig_GIMAP_tr_dom_LFC <- plot_grid(p2_no_legend, GIMAP_Interproscan_domain_plot_no_legend, C_gig_apop_LFC_GIMAP_tile_plot_no_legend, ncol=3, align='h',
                                   labels = c('A', 'B', 'C'), label_size = 12)
Cgig_GIMAP_tr_dom_LFC_legend <- plot_grid(GIMAP_MY_CV_CG_raxml_treedata_vertical_legend, GIMAP_Interproscan_domain_plot_legend,C_gig_apop_LFC_GIMAP_tile_plot_legend, 
                                          ncol = 3, align="hv")
Cgig_GIMAP_tr_dom_LFC_plus_legend <- plot_grid(Cgig_GIMAP_tr_dom_LFC, Cgig_GIMAP_tr_dom_LFC_legend, ncol=1, rel_heights =c(1, 0.2))


#### PLOT CONSTITUTIVELY EXPRESSED GIMAPS ####


# Load data from Transcriptome dataframes
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_gig_vst_common_df_all_mat_limma_GIMAP.RData")
load(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_vst_common_df_all_mat_limma_GIMAP.RData")

C_gig_vst_common_df_all_mat_limma_GIMAP
C_vir_vst_common_df_all_mat_limma_GIMAP


#### PLOT ORTHOFINDER SPECIES TREE  ####
Mollusc_Species_Tree_text <-"((Octopus_bimaculoides:0.0710909,Octopus_sinensis:0.056727)N1:0.21781,((Mizuhopecten_yessoensis:0.315015,(Crassostrea_gigas:0.0955031,Crassostrea_virginica:0.0982277)N5:0.236348)N3:0.0835452,(Lottia_gigantea:0.31253,(Pomacea_canaliculata:0.34807,(Elysia_chlorotica:0.303751,(Biomphalaria_glabrata:0.296022,Aplysia_californica:0.248891)N8:0.0608488)N7:0.129889)N6:0.0520687)N4:0.0492055)N2:0.21781)N0;"
Mollusc_Species_Tree <- read.newick(text=Mollusc_Species_Tree_text)
Mollusc_Species_tibble <- as.tibble(Mollusc_Species_Tree)

# add gene numbers to this table 
All_mollusc_IAP_gene_list_after_haplotig_collapsed_sinensis <- All_mollusc_IAP_gene_list_after_haplotig_collapsed %>% 
  mutate(Species = case_when(Species == "Octopus_vulgaris" ~ "Octopus_sinensis",
         TRUE ~ Species))

Mollusc_Species_tibble_join <- Mollusc_Species_tibble %>% rename(Species = label) %>% 
  left_join(., All_mollusc_IAP_gene_list_after_haplotig_collapsed_sinensis) %>% rename(label = Species, Genes = n)
Mollusc_Species_tibble_join$Genes <- as.numeric(Mollusc_Species_tibble_join$Genes)

Mollusc_Species_Treedata <- as.treedata(Mollusc_Species_tibble_join)

# write out species and genus
genus <- c('Octopus', 'Octopus','Mizuhopecten','Crassostrea','Crassostrea','Lottia','Pomacea','Elysia','Biomphalaria','Aplysia')
species <- c('bimaculoides','sinensis','yessoensis','gigas','virginica','gigantea','canaliculata','chlorotica',
             'glabrata','californica')
d <- data.frame(label = Mollusc_Species_Tree$tip.label,
                genus = genus,
                species = species)

# Plot species tree only
Mollusc_Species_Tree <- ggtree(Mollusc_Species_Treedata, branch.length = "none") %<+% d +
  geom_tiplab(align=TRUE, aes(label=paste0('italic(', genus,')~italic(', species, ')')), parse=T) + # italicize species labels 
   xlim(0,17) + theme(plot.margin = unit(c(0,0,0,0), "cm"), text = element_text(size = 14))

# Create simple heatmap to plot the gene number to add next to the tree 
no_y_axis <- function () 
  theme(axis.line.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

no_x_axis <- function () 
  theme(axis.line.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# put in order of the species tree 
All_mollusc_IAP_gene_list_after_haplotig_collapsed_sinensis$Species <- factor(All_mollusc_IAP_gene_list_after_haplotig_collapsed_sinensis$Species, 
          levels = c("Octopus_bimaculoides","Octopus_sinensis","Mizuhopecten_yessoensis","Crassostrea_gigas","Crassostrea_virginica","Lottia_gigantea","Pomacea_canaliculata","Elysia_chlorotica","Biomphalaria_glabrata","Aplysia_californica"),
          labels = c("Octopus bimaculoides","Octopus sinensis","Mizuhopecten yessoensis","Crassostrea gigas","Crassostrea virginica","Lottia gigantea","Pomacea canaliculata","Elysia chlorotica","Biomphalaria glabrata","Aplysia californica"))

species_gene_heatmap <- All_mollusc_IAP_gene_list_after_haplotig_collapsed_sinensis %>% 
  # add dummy category so I can only plot one row
  mutate(category = "all") %>%
ggplot(., aes(y=Species, x =category)) + geom_tile(aes(fill=n, width = 0.2)) + no_y_axis() + no_x_axis() + 
  theme(panel.background = element_rect(fill = "transparent"),
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size=10, family="sans"),
        legend.title = element_text(size=12, family="sans"),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin = margin(0,0,0,0),
        axis.ticks.margin = unit(0,"null"),
        axis.ticks.length = unit(0, "pt"), 
        panel.spacing = unit(0,"null"),
    #    axis.text.y = element_text(size = 14, family = "sans", face = "italic")) +
    axis.text.y = element_blank(),
    plot.margin = unit(c(-10,-10,-10,-10), "cm")) +
  # fill tiles with the number
  geom_text(size = 6, aes(label = round(n, 1))) +
  scale_fill_viridis_c(name = "IAP Gene Number", 
                       limits = c(0,90),
                       breaks = c(10,20,30,40,50,60,70,80,90), 
                       option="plasma",
                       guide=guide_legend(), na.value = "transparent") + 
  guides(fill=guide_legend(ncol=1, title.position="right")) 

# Plot the tree and the heatmap next to each other 
Mollusc_Species_Tree_heatmap <- Mollusc_Species_Tree + species_gene_heatmap

ggsave(filename = "Mollusc_Species_Tree_heatmap.tiff", plot= Mollusc_Species_Tree_heatmap, device="tiff",
       path="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/IAP_gene_tree/",
       width = 8.7 ,
       height = 7,
       units = "in",
       dpi=300)

#### SCRAP CODE ####

### Collapse and label duplicated sequences ###
# Load full sequence file so labels can be collapsed 
AIG_full_seq <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all.fa")
BIR_full_seq <- phylotools::read.fasta("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_HMMER_Interpro_XP_list_all.fa")
nrow(AIG_full_seq) # 403
nrow(BIR_full_seq)# 1082 - why were sequences added? when looking up sequences
# get rid of rows with duplicate headers
BIR_full_seq <- BIR_full_seq[!duplicated(BIR_full_seq$seq.name),]
colnames(AIG_full_seq)

# Test manually removing exact duplicates with base R
AIG_full_seq_man_rm_dup <- AIG_full_seq[!duplicated(AIG_full_seq$seq.text),]
nrow(AIG_full_seq_man_rm_dup) # 326
BIR_full_seq_man_rm_dup <- BIR_full_seq[!duplicated(BIR_full_seq$seq.text),]
nrow(BIR_full_seq_man_rm_dup) # 528

# Collapse duplicates using alakazam package
AIG_full_seq_collapse <- collapseDuplicates(AIG_full_seq, id = "seq.name", seq= "seq.text", text_fields = "seq.name", 
                                            add_count = TRUE, # adds column showing how many collapsed
                                            sep="_",
                                            verbose=TRUE)
#FUNCTION> collapseDuplicates
#FIRST_ID> XP_013073227.1 PREDICTED: protein AIG1-like [Biomphalaria glabrata]
#TOTAL> 403
#UNIQUE> 325
#COLLAPSED> 78
#DISCARDED> 0

BIR_full_seq_collapse <- collapseDuplicates(BIR_full_seq, id = "seq.name", seq= "seq.text", text_fields = "seq.name", 
                                            add_count = TRUE, # adds column showing how many collapsed
                                            sep="_",
                                            verbose=TRUE)
#FUNCTION> collapseDuplicates
#FIRST_ID> XP_022287996.1 baculoviral IAP repeat-containing protein 7-like [Crassostrea virginica]
#TOTAL> 791
#UNIQUE> 527
#COLLAPSED> 264
#DISCARDED> 0

# There is a discrepancy in duplicates collapsed by CD-HIT vs this software. CD-HIT produced and AIG list of 313 sq vs 325, and BIR a list of 499 vs. 527 from collapseDuplicates. This is still different from the 
# 

# Which are additional in collapseDuplicates that were removed by CD-HIT?
AIG_extra_collapse <- AIG_full_seq_collapse[!(AIG_full_seq_collapse$seq.name %in% names(AIG_seq_rm_dup)),]
BIR_extra_collapse <- BIR_full_seq_collapse[!(BIR_full_seq_collapse$seq.name %in% names(BIR_seq_rm_dup)),]

# Test if a few should have been removed by collapse duplicates
AIG_extra_collapse_test <- AIG_extra_collapse[1,]
AIG_full_seq[AIG_full_seq$seq.text %in%  AIG_extra_collapse_test$seq.text,] # only 1 hit in the original list
AIG_extra_collapse_test <- AIG_extra_collapse[2,]
AIG_full_seq[AIG_full_seq$seq.text %in%  AIG_extra_collapse_test$seq.text,]


# Potentially erroneous proteins in GIMAP:
# calponin homology domain-containing protein DDB_G0272472-like
# myb-like protein X
# PREDICTED: dentin sialophosphoprotein-like
# PREDICTED: GTPase Era, mitochondrial-like
# PREDICTED: mitochondrial ribosome-associated GTPase 1-like
# PREDICTED: reticulocyte-binding protein 2 homolog a isoform X2
# PREDICTED: putative protein PHLOEM PROTEIN 2-LIKE A3
# PREDICTED: rho-associated protein kinase let-502-like
# reticulocyte-binding protein 2 homolog a-like
# transmembrane GTPase fzo-like
# trichohyalin-like
# vicilin-like seed storage protein At2g18540 isoform X1
# Reviewing Domains and comparing list to Lu et al., 2020
# Lu et al. 2020 counted rho-associated protein kinase let-502-like , putative protein PHLOEM PROTEIN 2-LIKE A3 as a GIMAP, 
# Lu et al. 2020 called PREDICTED: dentin sialophosphoprotein-like
# the P-loop_NTPase Superfamily seems to be what is bringing up these extra hits
# Are there other domains that are specific to the erroneous proteins 

# Examine domains from these proteins
AIG1_CDD_GIMAP_questioned_hits <- AIG1_XP_ALL_gff_GIMAP_species[grepl("calponin",AIG1_XP_ALL_gff_GIMAP_species$product, ignore.case = TRUE) | grepl("myb", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) |
                                                                  grepl("dentin sialophosphoprotein", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | grepl("GTPase Era", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | 
                                                                  grepl("mitochondrial", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | grepl("trichohyalin", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE)|
                                                                  grepl("reticulocyte", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | grepl("PHLOEM", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) |
                                                                  grepl("rho", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | grepl("fzo", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE) | grepl("vicilin", AIG1_XP_ALL_gff_GIMAP_species$product,ignore.case = TRUE),]

# What domains do the GIMAPs have?
AIG1_CDD_GIMAP_only <- AIG1_XP_ALL_gff_GIMAP_species[grepl("IMAP",AIG1_XP_ALL_gff_GIMAP_species$product) | grepl("immune", AIG1_XP_ALL_gff_GIMAP_species$product),]
AIG1_CDD_GIMAP_only_domains <- AIG1_CDD_GIMAP_only %>% group_by(signature_desc) %>% summarise(count = n())
length(unique(AIG1_CDD_GIMAP_only$protein_id)) # 242 proteins 
AIG1_CDD_GIMAP_only_domains$type <-"GIMAP"

# What domains are in non GIMAPs?
AIG1_CDD_GIMAP_non_gimap <- AIG1_XP_ALL_gff_GIMAP_species[!grepl("IMAP", AIG1_XP_ALL_gff_GIMAP_species$product) | !grepl("uncharacterized",AIG1_XP_ALL_gff_GIMAP_species$product) | !grepl("immune", AIG1_XP_ALL_gff_GIMAP_species$product),]
AIG1_CDD_GIMAP_non_gimap_domains <- AIG1_CDD_GIMAP_non_gimap  %>% group_by(Short.name) %>% summarise(count = n())
AIG1_CDD_GIMAP_non_gimap_domains$type <-"non-GIMAP"

# What domains are unique to the non-GIMAPS?
AIG1_CDD_GIMAP_non_gimap_domains_NOT_SHARED <-AIG1_CDD_GIMAP_non_gimap_domains[!(AIG1_CDD_GIMAP_non_gimap_domains$Short.name %in%AIG1_CDD_GIMAP_only_domains$Short.name),]
# Pkc domains, STKc, FN3

# try removing the Pkc and STK containing domain proteins to see if this removes erroneous non GIMAP hits
AIG1_CDD_GIMAP_no_Pkc_STK <- AIG1_CDD_GIMAP %>% group_by(protein_id) %>% filter(all(!grepl("PKc",ignore.case = TRUE, Short.name) | !grepl("PKc",ignore.case = TRUE, Short.name))) 
View(unique(AIG1_CDD_GIMAP_no_Pkc_STK$product))

# compare to both gimap and non gimap lists
setdiff(AIG1_CDD_GIMAP_only$protein_id, AIG1_CDD_GIMAP_no_Pkc_STK$protein_id) #0 none we wanted were remove
setdiff(AIG1_CDD_GIMAP_no_Pkc_STK$protein_id, AIG1_CDD_GIMAP_only$protein_id) # 184 remain that weren't removed



#### SESSION INFO ####
sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.0.0        ggpubr_0.3.0         viridis_0.5.1        viridisLite_0.3.0    alakazam_1.0.1       chopper_0.1.8        data.table_1.12.8    rtracklayer_1.44.4  
# [9] GenomicRanges_1.36.1 GenomeInfoDb_1.20.0  forcats_0.5.0        stringr_1.4.0        dplyr_1.0.0          purrr_0.3.4          readr_1.3.1          tidyr_1.1.0         
# [17] tibble_3.0.1         tidyverse_1.3.0      plyr_1.8.6           ggimage_0.2.8        tidytree_0.3.3       treeio_1.8.2         phylotools_0.2.2     ggrepel_0.8.2       
# [25] ggplot2_3.3.2        ggtree_1.16.6        Biostrings_2.52.0    XVector_0.24.0       IRanges_2.18.3       S4Vectors_0.22.1     BiocGenerics_0.30.0  ape_5.4             
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1            seqinr_3.6-1                ggsignif_0.6.0              ellipsis_0.3.1              rio_0.5.16                  fs_1.4.1                   
# [7] rstudioapi_0.11             farver_2.0.3                fansi_0.4.1                 lubridate_1.7.9             xml2_1.3.2                  ade4_1.7-15                
# [13] jsonlite_1.6.1              Rsamtools_2.0.3             broom_0.5.6                 dbplyr_1.4.4                BiocManager_1.30.10         compiler_3.6.1             
# [19] httr_1.4.1                  rvcheck_0.1.8               backports_1.1.8             assertthat_0.2.1            Matrix_1.2-18               lazyeval_0.2.2             
# [25] cli_2.0.2                   prettyunits_1.1.1           tools_3.6.1                 igraph_1.2.5                gtable_0.3.0                glue_1.4.1                 
# [31] GenomeInfoDbData_1.2.1      Rcpp_1.0.3                  carData_3.0-4               Biobase_2.44.0              cellranger_1.1.0            vctrs_0.3.1                
# [37] nlme_3.1-148                openxlsx_4.1.5              rvest_0.3.5                 lifecycle_0.2.0             rstatix_0.6.0               XML_3.99-0.3               
# [43] zlibbioc_1.30.0             MASS_7.3-51.6               scales_1.1.1                hms_0.5.3                   SummarizedExperiment_1.14.1 yaml_2.2.1                 
# [49] curl_4.3                    gridExtra_2.3               stringi_1.4.6               zip_2.0.4                   BiocParallel_1.18.1         rlang_0.4.6                
# [55] pkgconfig_2.0.3             matrixStats_0.56.0          bitops_1.0-6                lattice_0.20-41             labeling_0.3                GenomicAlignments_1.20.1   
# [61] tidyselect_1.1.0            magrittr_1.5                R6_2.4.1                    magick_2.4.0                generics_0.0.2              DelayedArray_0.10.0        
# [67] DBI_1.1.0                   pillar_1.4.4                haven_2.3.1                 foreign_0.8-72              withr_2.2.0                 abind_1.4-5                
# [73] RCurl_1.98-1.2              janeaustenr_0.1.5           modelr_0.1.8                crayon_1.3.4                car_3.0-8                   utf8_1.1.4                 
# [79] progress_1.2.2              grid_3.6.1                  readxl_1.3.1                blob_1.2.1                  digest_0.6.25               reprex_0.3.0               
# [85] gridGraphics_0.5-0          munsell_0.5.0               ggplotify_0.0.5   
