#### R script to identify and test gene family expansion
# Erin Roberts, 2020
# PhD Candidate University of Rhode Island 

#### Load packages ####
library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree) # install the dev version to get the get.tree function
library(ggrepel)
library(phylotools)
library(treeio)
library(tidytree)
library(ggimage)
library(plyr)
library(tidyverse)
library(tidytext)
library(rtracklayer)
library(data.table)
library(chopper)
library(alakazam)
library(phylotools)

#### Import Genomes and Annotations for each species in order to facilitate loookup #####
#load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
#load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
# Load the gff file for all mollusc genomes
All_molluscs_CDS_gff <- readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/All_molluscs_CDS.gff")
All_molluscs_CDS_gff <- as.data.frame(All_molluscs_CDS_gff)
All_mollusc_gene_gff <- readGFF(file="/Volumes/My Passport for Mac/OrthoFinder_Genomes_Mar_2020_Paper1/GFF3/All_mollusc_gene.gff")
All_mollusc_gene_gff <- as.data.frame(All_mollusc_gene_gff)

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
length(unique(Cgig_gff_IAP_family_XP$protein_id)) #65
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

#### Load May15th Orthogroup Analysis of 10 Mollusc species from Orthogroup.tsv file ####
# Load tsv
Orthogroups <- read_tsv("/Volumes/My Passport for Mac/OrthoFinder_3_25_2020_Bluewaves_Backup/Results_May15/Orthogroups/Orthogroups.tsv",
                        col_names = c("Orthogroup","Elysia_chlorotica", "Aplysia_californica", "Crassostrea_gigas", "Lottia_gigantea", 
                                      "Biomphalaria_glabrata", "Octopus_bimaculoides",
                                      "C_virginica", "Mizuhopecten_yessoensis",	"Pomacea_canaliculata",	"Octopus_sinensis"))

#### ASSESS ORTHOGROUPS USING ONLY THE CV AND CG REFERENCE ANNOTATIONS ####
CV_CG_IAP_list <- as.list(CV_CG_IAP)
CV_CG_GIMAP_list <- as.list(CV_CG_GIMAP)

#CV_CG_IAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_IAP_list, collapse="|"), i))),]
length(CV_CG_IAP_list_lookup$Orthogroup)
# 27 orthogroups

#CV_CG_GIMAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_GIMAP_list, collapse="|"), i))),]
length(CV_CG_GIMAP_list_lookup$Orthogroup)
#9

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

#Count IAP genes across species
BIR_XP_gff_species_gene_count <- BIR_XP_gff_species %>% group_by(Species) %>% filter(is.na(locus_tag)) %>% distinct(gene)  %>% summarise(gene_count = n())
BIR_XP_gff_species_locus_tag_count <- BIR_XP_gff_species %>% group_by(Species) %>% filter(is.na(gene)) %>%  distinct(locus_tag) %>% summarise(locus_tag_count = n())
colnames(BIR_XP_gff_species_locus_tag_count)[2] <- "gene_count"
BIR_XP_gff_species_gene_locus_tag_count <- rbind(BIR_XP_gff_species_gene_count, BIR_XP_gff_species_locus_tag_count)

#Count GIMAP and IAN genes across species to compare with Lu et al. 2020 paper
AIG1_XP_ALL_gff_GIMAP_species_gene_count <- AIG1_XP_ALL_gff_GIMAP_species %>% group_by(Species) %>% filter(is.na(locus_tag)) %>% distinct(gene)  %>% summarise(gene_count = n())
AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count <- AIG1_XP_ALL_gff_GIMAP_species %>% group_by(Species) %>% filter(is.na(gene)) %>%  distinct(locus_tag) %>% summarise(locus_tag_count = n())
colnames(AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count)[2] <- "gene_count"
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_count <- rbind(AIG1_XP_ALL_gff_GIMAP_species_gene_count, AIG1_XP_ALL_gff_GIMAP_species_locus_tag_count)
# All are between five and 1 over

## EXPORT GENE LISTS PER SPECIES TO EXAMINE POTENTIAL ARTIFACTS
BIR_XP_gff_species_genes <- BIR_XP_gff_species %>% filter(is.na(locus_tag)) %>% distinct(gene)
BIR_XP_gff_species_locus_tag <- BIR_XP_gff_species %>% filter(is.na(gene)) %>%  distinct(locus_tag) 
colnames(BIR_XP_gff_species_locus_tag)[1] <- "gene"
BIR_XP_gff_species_gene_locus_tag <- rbind(BIR_XP_gff_species_genes, BIR_XP_gff_species_locus_tag)

length(BIR_XP_gff_species_gene_locus_tag$gene) # 380
write.table(BIR_XP_gff_species_gene_locus_tag$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_genes_HMMER_Interpro_BIR.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

AIG1_XP_ALL_gff_GIMAP_species_gene <- AIG1_XP_ALL_gff_GIMAP_species %>% filter(is.na(locus_tag)) %>% distinct(gene)
AIG1_XP_ALL_gff_GIMAP_species_locus_tag <- AIG1_XP_ALL_gff_GIMAP_species  %>% filter(is.na(gene)) %>%  distinct(locus_tag)
colnames(AIG1_XP_ALL_gff_GIMAP_species_locus_tag )[1] <- "gene"
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag <- rbind(AIG1_XP_ALL_gff_GIMAP_species_gene, AIG1_XP_ALL_gff_GIMAP_species_locus_tag)
length(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$gene) #252

write.table(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$gene, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/GIMAP_genes_HMMER_Interpro_AIG.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Convert locus tag and gene LOC to Extrez ID to get sequences with Batch join by Name column 
colnames(BIR_XP_gff_species_gene_locus_tag)[1] <- "Name"
BIR_XP_gff_species_gene_locus_tag_convert <- left_join(BIR_XP_gff_species_gene_locus_tag, All_mollusc_gene_gff[,c("Name","Dbxref","start","end")])
BIR_XP_gff_species_gene_locus_tag_convert$Dbxref <- str_remove(BIR_XP_gff_species_gene_locus_tag_convert$Dbxref,"GeneID:")
BIR_XP_gff_species_gene_locus_tag_convert <- BIR_XP_gff_species_gene_locus_tag_convert %>% filter(Dbxref != "character(0)")

colnames(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag)[1] <- "Name"
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
by(BIR_XP_gff_species_gene_locus_tag, BIR_XP_gff_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/BIR_Gene_list_", i$Species[1], ".txt"), 
quote = FALSE,col.names = FALSE, row.names=FALSE))

by(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag, AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/AIG_GIMAP_Gene_list_", i$Species[1], ".txt"), 
quote = FALSE,col.names = FALSE, row.names=FALSE))

## EXPORT ONLY C. VIRGINICA GENE LIST AS BED FILE WITH THE START AND END COORDINATES TO LOOK AT MAPPING COVERAGE AND COMPARE IDENTITY
BIR_XP_gff_species_gene_locus_tag_Cvir <- BIR_XP_gff_species_gene_locus_tag %>% filter(Species =="Crassostrea_virginica")
AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir <- AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag %>% filter(Species =="Crassostrea_virginica")
colnames(BIR_XP_gff_species_gene_locus_tag_Cvir)[1] <- "gene"
colnames(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_C_vir)[1] <- "gene"

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


## Review Matches
BIR_XP_gff_species 
AIG1_XP_ALL_gff_GIMAP_species
View(unique(BIR_XP_gff_species$product)) # has one phosphatase and actin regulator 4-B-like but does have the IAP repeats 
View(unique(AIG1_XP_ALL_gff_GIMAP_species$product))
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
# the P-loop_NTPase Superfamily seems to be what is brining up these extra hits
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

##################### ORTHOGROUP ANALYSIS ###############################

#### USE FULL IAP AND GIMAP LISTS TO PULL OUT ALL MOLLUSC ORTHOGROUPS ####
BIR_XP_gff_species_list <- as.list(unique(BIR_XP_gff_species$protein_id))
#BIR_XP_gff_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_species_list, collapse="|"), i))),]
length(BIR_XP_gff_species_list_lookup$Orthogroup)
# 35 orthogroups (got one extra orthogroup when setting hmmer to eval -3)

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_species_list_lookup$Orthogroup) # 4 missed "OG0003807" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
setdiff(BIR_XP_gff_species_list_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # 13 added "OG0001642" "OG0002611" "OG0004344" "OG0007118" "OG0011865" "OG0011926" "OG0012919" "OG0013878" "OG0014276" "OG0016483" "OG0017158" "OG0017983" "OG0018491"

AIG1_XP_ALL_gff_GIMAP_species_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_species$protein_id))
length(AIG1_XP_ALL_gff_GIMAP_species_list)
#AIG1_XP_ALL_gff_GIMAP_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_XP_ALL_gff_GIMAP_species_list, collapse="|"), i))),]
length(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup)
#12

setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup) # 2 not found "OG0013109" "OG0016155"
setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 5 added

#### USE CG AND CV SPECIFIC IAP AND GIMAP HMMER/INTERPROSCAN LISTS TO PULL OUT ORTHOGROUPS ###
BIR_XP_gff_CG_CV <- rbind(BIR_XP_gff_CG, BIR_XP_gff_CV)
AIG1_XP_ALL_gff_GIMAP_CG_CV <- rbind(AIG1_XP_ALL_gff_GIMAP_CG, AIG1_XP_ALL_gff_GIMAP_CV)

# Use CV and CG only lists to pull out orthogroups
BIR_XP_gff_CG_CV_list <- as.list(unique(BIR_XP_gff_CG_CV$protein_id))
#BIR_XP_gff_CG_CV_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_CG_CV_list, collapse="|"), i))),]
length(BIR_XP_gff_CG_CV_lookup$Orthogroup)
# 26 orthogroups, 11 are added when you include all the mollusc proteins

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_CG_CV_lookup$Orthogroup) #  "OG0003807" "OG0006831" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
setdiff(BIR_XP_gff_CG_CV_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # "OG0004344" "OG0011865" "OG0012919" "OG0013878" "OG0017983"

AIG1_CDD_GIMAP_only_CV_CG_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_CG_CV$protein_id))
length(AIG1_CDD_GIMAP_only_CV_CG_list)
#AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_CDD_GIMAP_only_CV_CG_list, collapse="|"), i))),]
length(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup)
# 9 orthogroups,

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup) # "OG0013109" "OG0016155"
setdiff(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 0"OG0007435" "OG0014793"

#### GET ALL MOLLUSCS IAP ORTHOGROUP HITS #####

## IAP genes 
# Get full list of proteins for each species by transposing and uniting
# Transpose the rows and column 
BIR_XP_gff_species_list_lookup_transpose <- t(BIR_XP_gff_species_list_lookup)
class(BIR_XP_gff_species_list_lookup_transpose) # matrix
BIR_XP_gff_species_list_lookup_transpose <- as.data.frame(BIR_XP_gff_species_list_lookup_transpose)
# unite all columns into one column 
BIR_XP_gff_species_list_lookup_transpose_united <- unite(BIR_XP_gff_species_list_lookup_transpose, full_protein_list, sep=",")
# remove NAs
BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub("NA,", "",BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub(",NA", "",BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list <- gsub("NA", "", BIR_XP_gff_species_list_lookup_transpose_united$full_protein_list)
# Put all into single vector for annot and export to make tree
# Concatenate each into single vector
BIR_XP_gff_species_list_lookup_transpose_united_all <- BIR_XP_gff_species_list_lookup_transpose_united %>% summarise(combined =paste(full_protein_list, collapse=","))
BIR_XP_gff_species_list_lookup_transpose_united_all_col <- data.frame(protein_id = unlist(strsplit(as.character(BIR_XP_gff_species_list_lookup_transpose_united_all$combined), ",")))
# trimws and remove orthogroups
BIR_XP_gff_species_list_lookup_united_all_col <- BIR_XP_gff_species_list_lookup_transpose_united_all_col[-c(1:35),1]
BIR_XP_gff_species_list_lookup_united_all_col <- trimws(BIR_XP_gff_species_list_lookup_united_all_col, which="left")

# How many XPs identified?
length(BIR_XP_gff_species_list_lookup_united_all_col) # 601
length(BIR_XP_gff_species_list) # 791

# Are there any duplicated proteins found in orthogroups?
BIR_XP_gff_species_list_lookup_united_all_col[duplicated(BIR_XP_gff_species_list_lookup_united_all_col)] # 0 duplicated

# Were all the original proteins found? - NO
setdiff(BIR_XP_gff_species_list,  BIR_XP_gff_species_list_lookup_united_all_col) # 272 proteins are present in original HMMER list that are not in the Orthogroup list
length(setdiff(BIR_XP_gff_species_list_lookup_united_all_col, BIR_XP_gff_species_list)) # 82 are added in that are not in the original orthogroup list

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
BIR_XP_gff_species_list_lookup_united_all_col_df <- as.data.frame(BIR_XP_gff_species_list_lookup_united_all_col)
colnames(BIR_XP_gff_species_list_lookup_united_all_col_df)[1]<-"protein_id"
BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- left_join(BIR_XP_gff_species_list_lookup_united_all_col_df, All_molluscs_CDS_gff[,c("protein_id","gene","product")])
BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- left_join(BIR_XP_gff_species_list_lookup_united_all_col_df_genes, All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene","product")])
# are any still unfilled?
BIR_XP_gff_species_list_lookup_united_all_col_df_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
# remove duplicates
BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes[!duplicated(BIR_XP_gff_species_list_lookup_united_all_col_df_genes[,c("gene","locus_tag")]),]
# how many genes
nrow(BIR_XP_gff_species_list_lookup_united_all_col_df_genes ) # 308

## Write out to table the gene list to use for gathering sequences for MAFFT 
# Collapse the locus tag and gene columns into one
BIR_XP_gff_species_list_lookup_united_all_col_df_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes %>% 
        unite(gene_locus_tag, c("gene","locus_tag"))
BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "NA_")
BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "_NA")

write.table(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_mollusc_orthogroups_gene_locus_tag_list.txt",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

### GIMAP genes
# Get full list of proteins for each species by transposing and uniting
# Transpose the rows and column 
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose <- t(AIG1_XP_ALL_gff_GIMAP_species_list_lookup)
class(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose) # matrix
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose)
# unite all columns into one column 
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united <- unite(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose, full_protein_list, sep=",")
# remove NAs
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub("NA,", "",AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub(",NA", "",AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list <- gsub("NA", "", AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united$full_protein_list)
# Put all into single vector for annot and export to make tree
# Concatenate each into single vector
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united %>% summarise(combined =paste(full_protein_list, collapse=","))
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all_col <- data.frame(protein_id = unlist(strsplit(as.character(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all$combined), ",")))
# trimws and remove orthogroups
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_transpose_united_all_col[-c(1:12),1]
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col <- trimws(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col, which="left")

# How many XPs identified?
length(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col) # 438
length(AIG1_XP_ALL_gff_GIMAP_species_list ) # 403

# Are there any duplicated proteins found in orthogroups?
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col[duplicated(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col)] #0 duplicated

# Were all the original proteins found? - NO
setdiff(AIG1_XP_ALL_gff_GIMAP_species_list,  AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col) # 132 proteins are present in original HMMER list that are not in the Orthogroup list
length(setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col, AIG1_XP_ALL_gff_GIMAP_species_list)) # 167 are added in that are not in the original orthogroup list

## Find genes in original HMMER list for each protein NOTE THAT Elysia chlorotica, Lottia gigantea ONLY HAVE LOCUS TAGS AND NOT GENES
AIG1_XP_ALL_gff_GIMAP_species_join <- left_join(unique(AIG1_XP_ALL_gff_GIMAP_species[,c("protein_id","product","Species")]), All_molluscs_CDS_gff[,c("protein_id","gene")])
AIG1_XP_ALL_gff_GIMAP_species_join <- left_join(AIG1_XP_ALL_gff_GIMAP_species_join[,c("protein_id","product","gene","Species")], All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene")])
AIG1_XP_ALL_gff_GIMAP_species_join <- unique(AIG1_XP_ALL_gff_GIMAP_species_join)
# remove duplicates
AIG1_XP_ALL_gff_GIMAP_species_genes <- AIG1_XP_ALL_gff_GIMAP_species_join[!duplicated(AIG1_XP_ALL_gff_GIMAP_species_join[,c("gene","locus_tag")]),]

# How many total genes or locus tags?
nrow(AIG1_XP_ALL_gff_GIMAP_species_genes) # 252

## Find genes in full Orthogroup list and compare
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df <- as.data.frame(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col)
colnames(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df)[1]<-"protein_id"
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- left_join(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df, All_molluscs_CDS_gff[,c("protein_id","gene","product")])
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- left_join(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes, All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene","product")])
# are any still unfilled?
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
# remove duplicates
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[!duplicated(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[,c("gene","locus_tag")]),]
# how many genes
nrow(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes ) # 308

## Write out to table the gene list to use for gathering sequences for MAFFT 
# Collapse the locus tag and gene columns into one
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes %>% 
  unite(gene_locus_tag, c("gene","locus_tag"))
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "NA_")
AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag <- str_remove(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, "_NA")

write.table(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_mollusc_orthogroups_gene_locus_tag_list.txt",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

###### COMPARATIVE STATISTICS BETWEEN HMMER, ANNOTATION, AND ORTHOFINDER #####

### IAP ###
## Were all the original HMMER/Interproscan genes found in the Orthofinder results?
IAP_Orthogroup_missing_genes <- BIR_XP_gff_species_genes[!(BIR_XP_gff_species_genes$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]
IAP_Orthogroup_missing_genes <- IAP_Orthogroup_missing_genes %>% filter(!is.na(gene)) %>% distinct(gene)
length(setdiff(BIR_XP_gff_species_genes$locus_tag,  BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag)) # 2  locus tag added
nrow(IAP_Orthogroup_missing_genes) #120
# 120 genes + 2 locus tags is 122 genes
  # Mostly missing PREDICTED genes, partial genes, and uncharacterized protein genes. All missing are Lottia gigantea and Elysia chlorotica genes

## Were any genes added in Orthogroups that weren't in HMMER?
IAP_Orthogroup_added_genes <- BIR_XP_gff_species_list_lookup_united_all_col_df_genes[!(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag %in% BIR_XP_gff_species_genes$gene),]
# need to add in line to take into account locus tag vs gene here 
setdiff(BIR_XP_gff_species_list_lookup_united_all_col_df_genes$locus_tag,  BIR_XP_gff_species_genes$locus_tag) # 0 new locus tags added

## Are all the CV and CG IAP genes found from the genome in the Orthogroup search ?
View(Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]) # 7 CG genes were missed by Orthofinder that were annotated in genome
View(Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% BIR_XP_gff_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]) # 12 CV genes were missed by Orthofinder that were annotated in genome

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
GIMAP_Orthogroup_missing_genes <- AIG1_XP_ALL_gff_GIMAP_species_genes[!(AIG1_XP_ALL_gff_GIMAP_species_genes$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),]
GIMAP_Orthogroup_missing_genes <- GIMAP_Orthogroup_missing_genes %>% filter(!is.na(gene)) %>% distinct(gene)
length(setdiff(AIG1_XP_ALL_gff_GIMAP_species_genes$locus_tag,  AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag)) # 3  locus tag added
nrow(GIMAP_Orthogroup_missing_genes) #65
# 65 genes + 3 locus tags is 68 genes

## Were any genes added in Orthogroups that weren't in HMMER?
GIMAP_Orthogroup_added_genes <- AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes[!(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag %in% AIG1_XP_ALL_gff_GIMAP_species_genes$gene),]
# need to add in line to take into account locus tag vs gene here 
setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$locus_tag,  AIG1_XP_ALL_gff_GIMAP_species_genes$locus_tag) # 0 new locus tags added

## Are all the CV and CG IAP genes found from the genome in the Orthogroup search ?
GIMAP_missing_CG<- Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),] #  CG genes were missed by Orthofinder that were annotated in genome
length(unique(GIMAP_missing_CG$gene)) #35
GIMAP_missing_CV <- Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% AIG1_XP_ALL_gff_GIMAP_species_list_lookup_united_all_col_df_genes$gene_locus_tag),] #  CV genes were missed by Orthofinder that were annotated in genome
length(unique(GIMAP_missing_CV$gene)) #67

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

### RUN PROTEIN SEQUENCES IN MAFFT AND RAXML ####
 # RAxML and HMMER full results used

### IDENTICAL PROTEINS REMOVED BY CD-HIT ###
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


## Parse the CD-HIT cluster file
# Followed code from this site: https://rpubs.com/rmurdoch/cdhit_to_mapping_file
AIG_seq_rm_dup_clustering <- read.txt("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa.clstr")
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

#### PLOT IAP TREE ####
# Helpful online tutorial regarding tool: https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
# Tree data vignette https://yulab-smu.github.io/treedata-book/faq.html#different-x-labels-for-different-facet-panels

#### PLOT GIMAP TREE ####
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

# Plot circular tree
GIMAP_raxml_treedata_circular_product <- ggtree(GIMAP_raxml_treedata, layout="circular", aes(color=Species), branch.length = "none") + 
  geom_tiplab2(aes(label=product,angle=angle), size =2.2, offset=.5) + # geom_tiplab2 flips the labels correctly
  theme(legend.position = "right", legend.text = element_text(face = "italic")) + xlim(-130,130)  

GIMAP_raxml_treedata_circular_product + scale_color_discrete(name = "Species", labels = c("Aplysia californica", 
                             "Biomphalaria glabrata", "Crassostrea gigas", "Crassostrea virginica","Elysia chlorotica","Lottia gigantea","Mizuhopecten yessoensis",
                             "Pomacea canaliculata","NA"))
#Figure out how to change the colors later
                          
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

### Export GIMAP Protein Domains as BED file ###


#### INVESTIGATE POTENTIAL GENE ARTIFACTS ####

# LOAD BED FILES AND HAPLOTIG FILES FROM JON PURITZ
  ## Notes: I created BED files for the locations of all Cvirginica genes run with RAxML in my tree for both IAP and GIMAP. Jon used those coordinates to pull out 
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

# Is their overlap with the haplotigs file?
Cvir_GIMAP_haplomerger_haplotigs <- Cvir_GIMAP_meanCov_CD_Hit_95_length[Cvir_GIMAP_meanCov_CD_Hit_95_length$start %in% Cvir_haplotigs,]
  # 0 in the overlap
Cvir_IAP_haplomerger_haplotigs <- Cvir_IAP_meanCov_CD_Hit_95_length[Cvir_IAP_meanCov_CD_Hit_95_length$start %in% Cvir_haplotigs,]
  # 0 exact gene overlaps, need to check if my genes hit to any of these ranges 

# No exact overlaps with gene coordinates, check if it is within range
Cvir_IAP_meanCov_CD_Hit_95_length[inrange(start, Cvir_haplotigs$start, Cvir_haplotigs$end)]


#### PLOT DOMAIN STRUCTURE ####
# trying the package Sushi
library(Sushi)





#### PLOT ORTHOFINDER SPECIES TREE WITH GENE COUNTS ####
Mollusc_Species_Tree_text <-"((Octopus_bimaculoides:0.0710909,Octopus_sinensis:0.056727)N1:0.21781,((Mizuhopecten_yessoensis:0.315015,(Crassostrea_gigas:0.0955031,C_virginica:0.0982277)N5:0.236348)N3:0.0835452,(Lottia_gigantea:0.31253,(Pomacea_canaliculata:0.34807,(Elysia_chlorotica:0.303751,(Biomphalaria_glabrata:0.296022,Aplysia_californica:0.248891)N8:0.0608488)N7:0.129889)N6:0.0520687)N4:0.0492055)N2:0.21781)N0;"
Mollusc_Species_Tree <- read.newick(text=Mollusc_Species_Tree_text)

Mollusc_Species_tibble <- as.tibble(Mollusc_Species_Tree)

# Join with GIMAP gene counts (in case I want to add to plotting later)
colnames(Mollusc_Species_tibble)[4] <- "Species"
Mollusc_Species_tibble_GIMAP_genes <- left_join(Mollusc_Species_tibble, AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag_count)
colnames(Mollusc_Species_tibble_GIMAP_genes)[5] <- "GIMAP gene Count"
colnames(Mollusc_Species_tibble_GIMAP_genes)[4] <- "label"

Mollusc_Species_Treedata <- as.treedata(Mollusc_Species_tibble_GIMAP_genes)

# write out species and genus
genus <- c('Octopus', 'Octopus','Mizuhopecten','Crassostrea','Crassostrea','Lottia','Pomacea','Elysia','Biomphalaria','Aplysia')
species <- c('bimaculoides','sinensis','yessoensis','gigas','virginica','gigantea','canaliculata','chlorotica',
             'glabrata','californica')
d <- data.frame(label = Mollusc_Species_Tree$tip.label,
                genus = genus,
                species = species
                )

# Change tip labels to be properly formatted
lb = get.tree(Mollusc_Species_Tree)$tip.label
d = data.frame(label=lb, label2 = paste(lb))

# Plot species tree only
Mollusc_Species_Tree <- ggtree(Mollusc_Species_Treedata, branch.length = "none") %<+% d +
  geom_tiplab(align=TRUE, aes(label=paste0('italic(', genus,')~italic(', species, ')')), parse=T) + # italicize species labels 
  ggtitle("Species Tree Orthogroup Analysis") +  xlim(0,17)


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


