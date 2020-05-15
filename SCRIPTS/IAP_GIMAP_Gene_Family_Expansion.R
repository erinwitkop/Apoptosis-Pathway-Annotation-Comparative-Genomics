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

#### Import Genomes and Annotations for each species in order to facilitate loookup #####
#load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
#load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")

#### Load March_25th Orthogroup Analysis of 31 species from Orthogroup.tsv file ####
# Load tsv

Orthogroups <- read_tsv("/Volumes/My Passport for Mac/OrthoFinder_3_25_2020_Bluewaves_Backup/Results_Mar25/Orthogroups/Orthogroups.tsv",
                        col_names = c("Orthogroup","Elysia_chlorotica",	"Amphibalanus_amphitrite","Drosophila_melanogaster",
                                      "Homo_sapiens",	"Mus_musculus",	"Danio_rerio", "Aplysia_californica",	"Strongylocentrotus_purpuratus","Caenorhabditis_elegans",
                                      "Saccoglossus_kowalevskii",	"Branchiostoma_floridae",	"Ciona_intestinalis",	"Crassostrea_gigas",	"Lottia_gigantea",
                                      "Biomphalaria_glabrata", "Priapulus_caudatus","Lingula_anatina","Octopus_bimaculoides",	"Exaiptasia_pallida",	"Branchiostoma_belcheri",
                                      "Acanthaster_planci",	"Crassostrea_virginica", "Orbicella_faveolata", "Mizuhopecten_yessoensis",
                                      "Pomacea_canaliculata",	"Pocillopora_damicornis",	"Penaeus_vannamei",	"Daphnia_magna","Acropora_millepora","Octopus_sinensis",	"Actinia_tenebrosa"))

#### IAP ANALYSIS ####

## Load XPs from all genomes that were found in both HMM and Interproscan on bluewaves 
BIR_XP_list <- read.table(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/BIR_hmmsearch_XP_seq.fa_1_BIR_repeat_sort_unique.tsv")
BIR_XP_list <- as.data.frame(BIR_XP_list)
AIG1_XP_gff <- readGFF(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/AIG1_hmmsearch_XP_seq.fa_1_AIG_Coil.gff3")
AIG1_XP_gff <- as.data.frame(AIG1_XP_gff)

## Filter out AIG Interproscan results for proteins that have a line with coil and a line with AIG domain and are GIMAP
AIG1_XP_gff_coil <- AIG1_XP_gff %>% filter(source=="Coils")
AIG1_XP_gff_coil <- AIG1_XP_gff_coil[!duplicated(AIG1_XP_gff_coil$seqid),]
AIG1_XP_gff_AIG1 <- AIG1_XP_gff %>% filter(source=="Pfam")
AIG1_XP_gff_AIG1 <- AIG1_XP_gff_AIG1[!duplicated(AIG1_XP_gff_AIG1$seqid),]

# merge the two to get matches in both
AIG1_XP_gff_GIMAP <- inner_join(AIG1_XP_gff_coil,AIG1_XP_gff_AIG1, by = "seqid")

### GET LIST OF JUST CV AND CGIG HITS FOR IAP AND GIMAP ###
#IAP
colnames(BIR_XP_list)[1] <- "protein_id"
BIR_XP_list_CV <- left_join(BIR_XP_list, C_vir_rtracklayer[,c("protein_id","product","gene")])
BIR_XP_list_CV <- drop_na(BIR_XP_list_CV)
BIR_XP_list_CV <- unique(BIR_XP_list_CV)
nrow(BIR_XP_list_CV) # 164 unique proteins after the HMM and Interproscan
BIR_XP_list_CV %>% group_by(gene) %>% dplyr::summarise(gene_count=n()) # 75 genes

BIR_XP_list_CG <- left_join(BIR_XP_list, C_gig_rtracklayer[,c("protein_id","product","gene")])
BIR_XP_list_CG <- drop_na(BIR_XP_list_CG)
BIR_XP_list_CG <- unique(BIR_XP_list_CG)
BIR_XP_list_CG %>% group_by(gene) %>% dplyr::summarise(gene_count=n()) # 39 genes
nrow(BIR_XP_list_CG) # 73 unique proteins after the HMM and Interproscan
BIR_XP_combined <- rbind(BIR_XP_list_CV, BIR_XP_list_CG)

#GIMAP
AIG1_XP_gff_GIMAP <- as.data.frame(AIG1_XP_gff_GIMAP$seqid)
colnames(AIG1_XP_gff_GIMAP)[1] <- "protein_id"
AIG1_XP_gff_GIMAP_CV <- left_join(AIG1_XP_gff_GIMAP, C_vir_rtracklayer[,c("protein_id","product","gene")])
AIG1_XP_gff_GIMAP_CV <- drop_na(AIG1_XP_gff_GIMAP_CV)
AIG1_XP_gff_GIMAP_CV <- unique(AIG1_XP_gff_GIMAP_CV)
nrow(AIG1_XP_gff_GIMAP_CV) # 109 unique proteins after the HMM and Interproscan
AIG1_XP_gff_GIMAP_CV %>% group_by(gene) %>% dplyr::summarise(gene_count=n()) # 59 genes

AIG1_XP_gff_GIMAP_CG <- left_join(AIG1_XP_gff_GIMAP, C_gig_rtracklayer[,c("protein_id","product","gene")])
AIG1_XP_gff_GIMAP_CG <- drop_na(AIG1_XP_gff_GIMAP_CG)
AIG1_XP_gff_GIMAP_CG <- unique(AIG1_XP_gff_GIMAP_CG)
AIG1_XP_gff_GIMAP_CG %>% group_by(gene) %>% dplyr::summarise(gene_count=n()) # 39 genes
nrow(AIG1_XP_gff_GIMAP_CG) # 30 unique proteins after the HMM and Interproscan

# combine lists
AIG1_XP_gff_GIMAP_combined <- rbind(AIG1_XP_gff_GIMAP_CV, AIG1_XP_gff_GIMAP_CG)

### USE FULL IAP AND GIMAP LISTS TO PULL OUT ORTHOGROUPS ###
BIR_XP_list <- as.list(BIR_XP_list$V1)
BIR_XP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_list, collapse="|"), i))),]
BIR_XP_list_lookup$Orthogroup
# 148 orthogroups

AIG1_XP_gff_GIMAP <- as.list(AIG1_XP_gff_GIMAP$seqid)
length(AIG1_XP_gff_GIMAP)
AIG1_XP_gff_GIMAP_list_lookup <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_XP_gff_GIMAP, collapse="|"), i))),]
AIG1_XP_gff_GIMAP_list_lookup$Orthogroup

### USE CG AND CV SPECIFIC IAP AND GIMAP LISTS TO PULL OUT ORTHOGROUPS ###
# Use CV and CG only lists to pull out orthogroups
BIR_XP_combined_list <- as.list(BIR_XP_combined$protein_id)
BIR_XP_list_lookup_CV_CG <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(BIR_XP_combined_list, collapse="|"), i))),]
BIR_XP_list_lookup$Orthogroup
# 15 orthogroups

AIG1_XP_gff_GIMAP_list <- as.list(AIG1_XP_gff_GIMAP_combined$protein_id)
length(AIG1_XP_gff_GIMAP_list)
AIG1_XP_gff_GIMAP_list_lookup_CV_CG <- Orthogroups[apply(Orthogroups, 1, function(i) any(grepl(paste(AIG1_XP_gff_GIMAP_list, collapse="|"), i))),]
AIG1_XP_gff_GIMAP_list_lookup_CV_CG$Orthogroup
#


## EXAMINE FULL ORTHOGROUP HITS THAT ARE NOT NA FOR C. VIR AND C. GIG ###
BIR_XP_list_lookup_CV <- BIR_XP_list_lookup %>% filter(!is.na(Crassostrea_virginica))
length(BIR_XP_list_lookup_CV$Orthogroup) # 85

BIR_XP_list_lookup_CG <- BIR_XP_list_lookup %>% filter(!is.na(Crassostrea_gigas))
length(BIR_XP_list_lookup_CG$Orthogroup) # 86

AIG1_XP_gff_GIMAP_list_lookup

## EXAMINE CV AND CG ORTHOGROUP HITS ####

## OUTPUT  LIST OF ALL XPS (REGARDLESS OF ORTHOGROUPS) SO I CAN PULL OUT ALL PROTEIN SEQUENCES IN BLUEWAVES
# Concatenate all into single vector
IAP_ALL_XP_lookup_all <- unite(IAP_ALL_XP_lookup, col = "all", sep = ",")
IAP_ALL_XP_lookup_all_sep <- str_split(IAP_ALL_XP_lookup_all, pattern = ",")
str(IAP_ALL_XP_lookup_all_sep) #$ list of 1
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep)
colnames(IAP_ALL_XP_lookup_all_sep)[1] <- "list"
# remove NAs 
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep[!grepl("NA", IAP_ALL_XP_lookup_all_sep$list),])
colnames(IAP_ALL_XP_lookup_all_sep)[1] <- "list"
# remove row 1
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep[-1,1])

# trimws
IAP_ALL_XP_lookup_all_sep <- trimws(IAP_ALL_XP_lookup_all_sep[,1], which="left")
length(IAP_ALL_XP_lookup_all_sep) # 18206
class(IAP_ALL_XP_lookup_all_sep) # character

# write out to file so I can get all of the sequences in bluewaves 
write.table(IAP_ALL_XP_lookup_all_sep, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IAP/IAP_ALL_XP_lookup_all_sep.txt",  
            row.names=FALSE, quote= FALSE, col.names = FALSE)

# transferring file to bluewaves to fetch sequences with `fetch_all_IAP_seq.sh` to get sequences

#### RUN MAFFT ALIGNMENT AND MAKE TREES IN BLUEWAVES ####
#Mafft chosen because of multi-threading capabilities


## OUTPUT  LIST OF ALL XPS (REGARDLESS OF ORTHOGROUPS) SO I CAN PULL OUT ALL PROTEIN SEQUENCES IN BLUEWAVES
# Concatenate all into single vector
IAP_ALL_XP_lookup_all <- unite(IAP_ALL_XP_lookup, col = "all", sep = ",")
IAP_ALL_XP_lookup_all_sep <- str_split(IAP_ALL_XP_lookup_all, pattern = ",")
str(IAP_ALL_XP_lookup_all_sep) #$ list of 1
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep)
colnames(IAP_ALL_XP_lookup_all_sep)[1] <- "list"
# remove NAs 
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep[!grepl("NA", IAP_ALL_XP_lookup_all_sep$list),])
colnames(IAP_ALL_XP_lookup_all_sep)[1] <- "list"
# remove row 1
IAP_ALL_XP_lookup_all_sep <- as.data.frame(IAP_ALL_XP_lookup_all_sep[-1,1])

# trimws
IAP_ALL_XP_lookup_all_sep <- trimws(IAP_ALL_XP_lookup_all_sep[,1], which="left")
length(IAP_ALL_XP_lookup_all_sep) # 18206
class(IAP_ALL_XP_lookup_all_sep) # character

# write out to file so I can get all of the sequences in bluewaves 
write.table(IAP_ALL_XP_lookup_all_sep, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IAP/IAP_ALL_XP_lookup_all_sep.txt",  
            row.names=FALSE, quote= FALSE, col.names = FALSE)

# transferring file to bluewaves to fetch sequences with `fetch_all_IAP_seq.sh` to get sequences

## CHECK FOR DUPLICATE PROTEIN HITS IN TWO SEPARATE ORTHOGROUPS TO ASSESS NEED FOR MERGE ##
# Parse and trimws
IAPC <- strsplit(IAP_ALL_XP_lookup$Crassostrea_virginica, split = ",")
IAPC_parse <-data.frame(Orthogroup = rep(IAP_ALL_XP_lookup$Orthogroup, sapply(IAPC, length)), protein_id = unlist(IAPC))
# trimws for checking list later 

IAPC_parse$protein_id <- trimws(IAPC_parse$protein_id, which=c("left"))

IAPG <- strsplit(IAP_ALL_XP_lookup$Crassostrea_gigas, split = ",")
IAPG_parse <-data.frame(Orthogroup = rep(IAP_ALL_XP_lookup$Orthogroup, sapply(IAPG, length)), protein_id  = unlist(IAPG))
# trimws for checking list later 
IAPG_parse$protein_id <- trimws(IAPG_parse$protein_id, which=c("left"))

# check for duplicated protein names
IAPC_parse[duplicated(IAPC_parse),] # 0 duplicated
IAPG_parse[duplicated(IAPG_parse$protein_id),] # some duplicated NA rows
# No orthogroups need to be combined

## CHECK FOR NON MAPPED PROTEINS FROM INITIAL CVIR CGIG ##
# Were any C_vir or C_gig TLR annotated proteins NOT mapped to Orthogroups?
IAP_CV_notmapped <- left_join(IAP_XP_CV, IAPC_parse, by ="protein_id") %>% filter(is.na(Orthogroup)) %>% View() # Two IAP 7's were not mapped, baculoviral IAP repeat-containing protein 7-A-like
IAP_CG_notmapped <- left_join(IAP_XP_CG, IAPG_parse, by ="protein_id") %>% filter(is.na(Orthogroup)) %>% View() # All were mapped

## ANNOTATE PROTEINS IN ALL C. VIR AND C. GIG ORTHOGROUPS

# Rbind the parsed files for each species by Orthogroup
IAPC_parse$species <- "C_virginica"
IAPG_parse$species <- "C_gigas"

IAP_CV_parsed_annot <- left_join(IAPC_parse, C_vir_rtracklayer[,c("gene","product","protein_id")])
IAP_CV_parsed_annot <- IAP_CV_parsed_annot[!duplicated(IAP_CV_parsed_annot$protein_id),]

IAP_CG_parsed_annot <- left_join(IAPG_parse, C_gig_rtracklayer[,c("gene","product","protein_id")])
IAP_CG_parsed_annot <- IAP_CG_parsed_annot[!duplicated(IAP_CG_parsed_annot$protein_id),]

# bind together the datatables 
IAP_CV_CG_parsed <- rbind(IAP_CV_parsed_annot, IAP_CG_parsed_annot)

# Remove orthogroups without any IAP proteins
IAP_CV_CG_parsed_IAP_only <- IAP_CV_CG_parsed %>% group_by(Orthogroup) %>% 
  
  
  ## WHAT ARE THE ANNOTATIONS FOR THE PROTEINS IN EACH ORTHOGROUP?
  
  # LOAD ALL XP ANNOTATIONS FROM BLUEWAVES FOR JOINING
  # Load annotations from bluewaves
  IAP_ALL_XP_annot <- readGFF(file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IAP/IAP_XP_proteins_all_genomes.txt")
IAP_ALL_XP_annot <- as.data.frame(IAP_ALL_XP_annot)
IAP_ALL_XP_annot <- IAP_ALL_XP_annot[!duplicated(IAP_ALL_XP_annot$Name),]

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

