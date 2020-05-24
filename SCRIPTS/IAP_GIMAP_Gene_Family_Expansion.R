#### R script to identify and test gene family expansion
# Erin Roberts, 2020
# PhD Candidate University of Rhode Island 

#### Load packages ####
library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree) # install the dev version to get the get.tree function
library(phylotools)
library(treeio)
library(ggimage)
library(plyr)
library(tidyverse)
library(tidytext)
library(rtracklayer)
library(rentrez)
library(data.table)
library(chopper)

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
Cvir_gff_IAP_family_XP <- Cvir_gff_IAP_family_XP[!is.na(Cvir_gff_IAP_family_XP$protein_id),]
Cgig_gff_IAP_family <- C_gig_rtracklayer[grepl("inhibitor of apoptosis", C_gig_rtracklayer$product, ignore.case=TRUE) | grepl("XIAP", C_gig_rtracklayer$product, ignore.case=TRUE)|
                                           grepl("baculoviral", C_gig_rtracklayer$product, ignore.case=TRUE),]
Cgig_gff_IAP_family_XP <- Cgig_gff_IAP_family[!duplicated(Cgig_gff_IAP_family$protein_id),]
Cgig_gff_IAP_family_XP <- Cgig_gff_IAP_family_XP[!is.na(Cgig_gff_IAP_family_XP$protein_id),]

Cvir_gff_GIMAP_family <- C_vir_rtracklayer[grepl("IMAP", C_vir_rtracklayer$product, ignore.case=TRUE) | grepl("immune-associated nucleotide-binding protein", C_vir_rtracklayer$product, ignore.case=TRUE),]
Cvir_gff_GIMAP_family_XP <- Cvir_gff_GIMAP_family[!duplicated(Cvir_gff_GIMAP_family$protein_id),]
Cvir_gff_GIMAP_family_XP <- Cvir_gff_GIMAP_family_XP[!is.na(Cvir_gff_GIMAP_family_XP$protein_id),]

Cgig_gff_GIMAP_family <- C_gig_rtracklayer[grepl("IMAP", C_gig_rtracklayer$product, ignore.case=TRUE) | grepl("immune-associated nucleotide-binding protein", C_gig_rtracklayer$product, ignore.case=TRUE),]
Cgig_gff_GIMAP_family_XP <- Cgig_gff_GIMAP_family[!duplicated(Cgig_gff_GIMAP_family$protein_id),]
Cgig_gff_GIMAP_family_XP <- Cgig_gff_GIMAP_family_XP[!is.na(Cgig_gff_GIMAP_family_XP$protein_id),]

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

CV_CG_IAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_IAP_list, collapse="|"), i))),]
length(CV_CG_IAP_list_lookup$Orthogroup)
# 27 orthogroups

CV_CG_GIMAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(CV_CG_GIMAP_list, collapse="|"), i))),]
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

#Count IAP genes across species to compare with Lu et al. 2020 paper
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

## EXPORT GENE LISTS PER SPECIES TO RUN IN MAFFT AND RAXML
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

# Export gene lists by species (not needed since I'm going batch Entrez)
by(BIR_XP_gff_species_gene_locus_tag, BIR_XP_gff_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/BIR_Gene_list_", i$Species[1], ".txt"), 
quote = FALSE,col.names = FALSE, row.names=FALSE))

by(AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag, AIG1_XP_ALL_gff_GIMAP_species_gene_locus_tag$Species, FUN=function(i) write.table(i$gene, 
paste0("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/Gene_lists_by_species/AIG_GIMAP_Gene_list_", i$Species[1], ".txt"), 
quote = FALSE,col.names = FALSE, row.names=FALSE))


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

### FOR GIMAP GENES, PROCEEDING FORWARD WITH USING THE ORTHOGROUPS TO FIND MISSED PROTEINS 

##################### ORTHOGROUP ANALYSIS ###############################

#### USE FULL IAP AND GIMAP LISTS TO PULL OUT ALL MOLLUSC ORTHOGROUPS ####
BIR_XP_gff_species_list <- as.list(unique(BIR_XP_gff_species$protein_id))
BIR_XP_gff_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_species_list, collapse="|"), i))),]
length(BIR_XP_gff_species_list_lookup$Orthogroup)
# 35 orthogroups (got one extra orthogroup when setting hmmer to eval -3)

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_species_list_lookup$Orthogroup) # 4 missed "OG0003807" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
setdiff(BIR_XP_gff_species_list_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # 13 added "OG0001642" "OG0002611" "OG0004344" "OG0007118" "OG0011865" "OG0011926" "OG0012919" "OG0013878" "OG0014276" "OG0016483" "OG0017158" "OG0017983" "OG0018491"

AIG1_XP_ALL_gff_GIMAP_species_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_species$protein_id))
length(AIG1_XP_ALL_gff_GIMAP_species_list)
AIG1_XP_ALL_gff_GIMAP_species_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_XP_ALL_gff_GIMAP_species_list, collapse="|"), i))),]
length(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup)
#12

setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup) # 2 not found "OG0013109" "OG0016155"
setdiff(AIG1_XP_ALL_gff_GIMAP_species_list_lookup$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 5 added

#### USE CG AND CV SPECIFIC IAP AND GIMAP HMMER/INTERPROSCAN LISTS TO PULL OUT ORTHOGROUPS ###
BIR_XP_gff_CG_CV <- rbind(BIR_XP_gff_CG, BIR_XP_gff_CV)
AIG1_XP_ALL_gff_GIMAP_CG_CV <- rbind(AIG1_XP_ALL_gff_GIMAP_CG, AIG1_XP_ALL_gff_GIMAP_CV)

# Use CV and CG only lists to pull out orthogroups
BIR_XP_gff_CG_CV_list <- as.list(unique(BIR_XP_gff_CG_CV$protein_id))
BIR_XP_gff_CG_CV_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_gff_CG_CV_list, collapse="|"), i))),]
length(BIR_XP_gff_CG_CV_lookup$Orthogroup)
# 26 orthogroups, 11 are added when you include all the mollusc proteins

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_IAP_list_lookup$Orthogroup, BIR_XP_gff_CG_CV_lookup$Orthogroup) #  "OG0003807" "OG0006831" "OG0015932" "OG0016100" "OG0018222" "OG0020281"
setdiff(BIR_XP_gff_CG_CV_lookup$Orthogroup, CV_CG_IAP_list_lookup$Orthogroup) # "OG0004344" "OG0011865" "OG0012919" "OG0013878" "OG0017983"

AIG1_CDD_GIMAP_only_CV_CG_list <- as.list(unique(AIG1_XP_ALL_gff_GIMAP_CG_CV$protein_id))
length(AIG1_CDD_GIMAP_only_CV_CG_list)
AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_CDD_GIMAP_only_CV_CG_list, collapse="|"), i))),]
length(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup)
# 9 orthogroups,

#Compare to list from original orthogroup search using only proteins in annotation
setdiff(CV_CG_GIMAP_list_lookup$Orthogroup,AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup) # "OG0013109" "OG0016155"
setdiff(AIG1_CDD_GIMAP_only_CV_CG_list_lookup_CV_CG$Orthogroup, CV_CG_GIMAP_list_lookup$Orthogroup) # 0"OG0007435" "OG0014793"

#### GET ALL MOLLUSCS IAP ORTHOGROUP HITS #####

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
BIR_CDD_BIR_genes <- left_join(unique(BIR_CDD_BIR[,c("protein_id","product","species")]), All_molluscs_CDS_gff[,c("protein_id","gene")])
BIR_CDD_BIR_genes <- left_join(BIR_CDD_BIR_genes[,c("protein_id","product","gene","species")], All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene")])
  # are any still unfilled?
#BIR_CDD_BIR_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
# remove duplicates
BIR_CDD_BIR_genes <- BIR_CDD_BIR_genes[!duplicated(BIR_CDD_BIR_genes[,c("gene","locus_tag")]),]

# How many total genes or locus tags?
nrow(BIR_CDD_BIR_genes) # 377

## Find genes in full Orthogroup list and compare
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df <- as.data.frame(BIR_CDD_BIR_list_lookup_transpose_united_all_col)
colnames(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df)[1]<-"protein_id"
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes <- left_join(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df, All_molluscs_CDS_gff[,c("protein_id","gene","product")])
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes <- left_join(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes, All_molluscs_CDS_gff[,c("protein_id","locus_tag","gene","product")])
# are any still unfilled?
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes %>% filter(is.na(gene) & is.na(locus_tag)) %>% View() # none
# remove duplicates
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes <- BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes[!duplicated(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes[,c("gene","locus_tag")]),]
# how many genes
nrow(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes ) # 304

## Write out to table the gene list to use for gathering sequences for MAFFT 
# Collapse the locus tag and gene columns into one
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes <- BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes %>% 
        unite(gene_locus_tag, c("gene","locus_tag"))
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene_locus_tag, "NA_")
BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene_locus_tag <- str_remove(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene_locus_tag, "_NA")

write.table(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene_locus_tag, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_IAP_mollusc_orthogroups_gene_locus_tag_list.txt",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

## Were all the original HMMER/Interproscan genes found in the Orthofinder results?
IAP_Orthogroup_missing_genes <- BIR_CDD_BIR_genes[!(BIR_CDD_BIR_genes$gene %in% BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene),]
setdiff(BIR_CDD_BIR_genes$locus_tag,  BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$locus_tag) # 1  locus tag added
# 108 genes are missing what are the missing genes?
  # Mostly missing PREDICTED genes, partial genes, and uncharacterized protein genes. Perhaps Orthofinder was too stringent?

## Were any genes added that weren't in HMMER?
IAP_Orthogroup_added_genes <- BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes[!(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene %in% BIR_CDD_BIR_genes$gene),]
setdiff(BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$locus_tag,  BIR_CDD_BIR_genes$locus_tag) # two new locus tags added

## Are all the CV and CG IAP genes found from the genome in the Orthogroup search ?
Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene),] # 10 CG genes were missed by Orthofinder that were annotated in genome
Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% BIR_CDD_BIR_list_lookup_transpose_united_all_col_df_genes$gene),] # 14 CV genes were missed by Orthofinder that were annotated in genome

## Were any genes added by HMMER that were not in the genome? 
BIR_CDD_BIR_genes_CG <- BIR_CDD_BIR_genes %>% filter(species=="Crassostrea gigas")
BIR_CDD_BIR_genes_CG[!(BIR_CDD_BIR_genes_CG$gene %in% Cgig_gff_IAP_family_XP$gene),] # 8 uncharacterized Loci genes were added by HMMER that were not annotated in genome

BIR_CDD_BIR_genes_CV <- BIR_CDD_BIR_genes %>% filter(species=="Crassostrea virginica")
BIR_CDD_BIR_genes_CV[!(BIR_CDD_BIR_genes_CV$gene %in% Cvir_gff_IAP_family_XP$gene),] # 14 uncharacterized Loci genes were added by HMMER that were not annotated in genome

## Are all the CV and CG IAP genes found from the genome in the HMMER/Interproscan search ?
Cgig_gff_IAP_family_XP[!(Cgig_gff_IAP_family_XP$gene %in% BIR_CDD_BIR_genes$gene),] # 4 CG genes were missed by HMMER that were annotated in genome
    #LOC105333301 # only the zinc finger domain
    #LOC105336740 # only the RingUbox domain 
    #LOC105338773 # only has the UBCc no BIR
    #LOC105338774 # only has UBCc domains no BIR
Cvir_gff_IAP_family_XP[!(Cvir_gff_IAP_family_XP$gene %in% BIR_CDD_BIR_genes$gene),] # 6 CV genes were missed by HMMER that were annotated in genome
    #LOC111136287 # only the pfam zinc ring finger domain        
    #LOC111101682 # only the RingUbox domain no BIR repeat           
    #LOC111100802 # has the RING-HC_BIRC2_3_7; RING finger, HC subclass, found in apoptosis protein c-IAP1, c-IAP2, livin, and similar proteins but not the BIR repeat        
    #LOC111104430 # no BIR repeat only RingUbox in the NCBI domains            
    #LOC111109770 # no BIR repeat only RingUbox in the NCBI domains            
    #LOC111106726 # only has the RINGUbox           

### RUN GENE SEQUENCES IN MAFFT AND RAXML ####
#convert fasta file format to phylip format


#### PLOT DOMAIN STRUCTURE ####


#### PLOT IAP TREE ####


#### PLOT GIMAP TREE ####




#### PLOT SPECIES TREE  ####
# try ggtree instead! 
#https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html


# Helpful online tutorial regarding tool: https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
# Tree data vignette https://yulab-smu.github.io/treedata-book/faq.html#different-x-labels-for-different-facet-panels

Species_Tree_text <- "((Caenorhabditis_elegans:0.72042,((Amphibalanus_amphitrite:0.432549,Penaeus_vannamei:0.373059)0.265388:0.0455497,(Daphnia_magna:0.393008,Drosophila_melanogaster:0.468793)0.246603:0.0407988)0.398481:0.0643439)0.239408:0.0202884,((((Acropora_millepora:0.158152,(Orbicella_faveolata:0.153909,Pocillopora_damicornis:0.121481)0.736211:0.0511311)0.907274:0.119056,(Actinia_tenebrosa:0.158645,Exaiptasia_pallida:0.178692)0.944045:0.115837)0.929656:0.186001,(((Danio_rerio:0.184049,(Mus_musculus:0.0471454,Homo_sapiens:0.0387577)0.974021:0.142751)0.966827:0.167685,Ciona_intestinalis:0.500739)0.251399:0.0407473,(((Strongylocentrotus_purpuratus:0.298998,Acanthaster_planci:0.286896)0.847722:0.0814753,Saccoglossus_kowalevskii:0.336273)0.470024:0.0436738,(Branchiostoma_floridae:0.0977097,Branchiostoma_belcheri:0.0921079)0.963629:0.242379)0.151479:0.0255123)0.0743405:0.0289695)0.149081:0.041054,(Priapulus_caudatus:0.401133,(Lingula_anatina:0.309775,((Octopus_bimaculoides:0.0552811,Octopus_sinensis:0.0418086)0.97482:0.306353,((Mizuhopecten_yessoensis:0.244137,(Crassostrea virginica:0.0766964,Crassostrea_gigas:0.0711655)0.980815:0.185791)0.768585:0.0563199,(Lottia_gigantea:0.228622,(((Aplysia_californica:0.1813,Biomphalaria_glabrata:0.222248)0.377298:0.0526096,Elysia_chlorotica:0.226705)0.808553:0.104453,Pomacea_canaliculata:0.262457)0.546763:0.0388776)0.347322:0.0279315)0.26299:0.0240144)0.419265:0.0420681)0.421663:0.0401548)0.131495:0.02577)0.239408:0.0202884);"
Species_Tree <- read.newick(text=Species_Tree_text)

# Plot tree
# write out species and genus so I can 
genus <- c('Caenorhabditis','Amphibalanus','Penaeus','Daphnia','Drosophila','Acropora','Orbicella',
           'Pocillopora','Actinia','Exaiptasia','Danio','Mus','Homo','Ciona','Strongylocentrotus',
           'Acanthaster','Saccoglossus','Branchiostoma','Branchiostoma','Priapulus','Lingula',
           'Octopus','Octopus','Mizuhopecten','Crassostrea','Crassostrea','Lottia','Aplysia',
           'Biomphalaria','Elysia','Pomacea')
species <- c('elegans','amphitrite','vannamei','magna','melanogaster','millepora','faveolata',
             'damicornis','tenebrosa','pallida','rerio','musculus','sapiens','intestinalis',
             'purpuratus','planci','kowalevskii','floridae','belcheri','caudatus','anatina',
             'bimaculoides','sinensis','yessoensis','virginica','gigas','gigantea','californica',
             'glabrata','chlorotica','canaliculata')
d <- data.frame(label = Species_Tree$tip.label, genus = genus,
                species = species)
# add phlyopic image
#g <- ggimage::phylopic_uid(Species_Tree$tip.label)

ggtree(Species_Tree, branch.length = "none") %<+% d + #%<+% g + # make it a cladogram with branch.length = "none" rather than plotting by evolutionary distance
  geom_tiplab(align=TRUE, aes(label=paste0('italic(', genus, 
                                           ')~italic(', species, ')')), 
              parse=T) + # italicize species labels 
  #geom_tiplab(aes(image=uid), geom="phylopic", offset=2.5) +
  ggtitle("Species Tree Orthogroup Analysis") + 
  geom_hilight(node=52, fill="light green") + xlim(0,17)

#geom_hilight(node=33, fill="purple") +  
#geom_hilight(node=39, fill="blue") +    
#geom_hilight(node=44, fill="dark green") +    
#geom_text(aes(label=node), hjust=-.3) 
# Edit colors on tree after looking at the results for groups for gene family expansion 
# add blocks for protostome, deuterostomes, lophotrochozoa

# Change tip labels to be properly formatted
lb = get.tree(Species_Tree)$tip.label
d = data.frame(label=lb, label2 = paste(lb))
ggtree(Species_Tree) %<+% d + geom_tiplab(aes(label=label2))


#### SCRAP CODE ####

