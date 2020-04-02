#### R script to identify and test gene family expansion across 31 genes
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
library(tidyverse)
library(tidytext)
library(rtracklayer)
library(rentrez)

#### Upload Genome Annotations for each species in order to facilitate loookup #####
C_vir_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)
C_gig_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/GCF_000297895.1_oyster_v9_genomic.gff")
C_gig_rtracklayer <- as.data.frame(C_gig_rtracklayer)

#### Plot Species Tree for Figure 2 ####
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
d = data.frame(label=lb, label2 = paste("AA", substring(lb, 1, 5)))
ggtree(Species_Tree) %<+% d + geom_tiplab(aes(label=label2))

###### Identify C. vir and C. gig XPs for expanded gene families to search for ####

## List of expanded gene families to search for first 
#  TLR
#  IAP
#  GIMAP
#  programmed Cell death proteins
#  calpains
#  TNF receptor associated factor
#  IP3R
#  mitogen activated protein kinase kinase kinase
#  MyD88
#  cAMP responsive element binding protein
#  caspase
#  IFI44
#  tumor necrosis factor receptor
#  cathepsin 

# get list of XPs from those genes by using the genome
## TLR
TLR <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "toll-like receptor")
TLR_XP_CG <- unique(TLR$gene)
TLR_XP_CG <- as.data.frame(TLR_XP_CG)
colnames(TLR_XP_CG)[1] <- "gene"
TLR_XP_CG <- left_join(TLR_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
TLR_XP_CG <- TLR_XP_CG[!duplicated(TLR_XP_CG),]  
TLR_XP_CG <-  TLR_XP_CG %>% filter(!is.na(protein_id))

TLRCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "toll-like receptor")
TLR_XP_CV <- unique(TLRCV$gene)
TLR_XP_CV <- as.data.frame(TLR_XP_CV)
colnames(TLR_XP_CV)[1] <- "gene"
TLR_XP_CV <- left_join(TLR_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
TLR_XP_CV <- TLR_XP_CV[!duplicated(TLR_XP_CV),]
TLR_XP_CV <-  TLR_XP_CV %>% filter(!is.na(protein_id))

TLR <- rbind(TLR_XP_CG,TLR_XP_CV)
# lookup in gff3 for each species
#write.table(TLR$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/TLR_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## IAP 
IAP <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "baculoviral IAP")
IAP_XP_CG <- unique(IAP$gene)
IAP_XP_CG <- as.data.frame(IAP_XP_CG)
colnames(IAP_XP_CG)[1] <- "gene"
IAP_XP_CG <- left_join(IAP_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
IAP_XP_CG <- IAP_XP_CG %>% filter(!is.na(protein_id)) 

IAPCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "baculoviral IAP")
IAP_XP_CV <- unique(IAPCV$gene)
IAP_XP_CV <- as.data.frame(IAP_XP_CV)
colnames(IAP_XP_CV)[1] <- "gene"
IAP_XP_CV <- left_join(IAP_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
IAP_XP_CV <- IAP_XP_CV[!duplicated(IAP_XP_CV),]
IAP_XP_CV <- IAP_XP_CV %>% filter(!is.na(protein_id))

IAP <- rbind(IAP_XP_CG,IAP_XP_CV)
# lookup in gff3 for each species
# write.table(IAP$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IAP_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## GIMAP 
GIMAP <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "GTPase IMAP family member")
GIMAP_XP_CG <- unique(GIMAP$gene)
GIMAP_XP_CG <- as.data.frame(GIMAP_XP_CG)
colnames(GIMAP_XP_CG)[1] <- "gene"
GIMAP_XP_CG <- left_join(GIMAP_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
GIMAP_XP_CG <- GIMAP_XP_CG[!duplicated(GIMAP_XP_CG),]  
GIMAP_XP_CG  <- GIMAP_XP_CG %>% filter(!is.na(protein_id))

GIMAPCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "GTPase IMAP family member")
GIMAP_XP_CV <- unique(GIMAPCV$gene)
GIMAP_XP_CV <- as.data.frame(GIMAP_XP_CV)
colnames(GIMAP_XP_CV)[1] <- "gene"
GIMAP_XP_CV <- left_join(GIMAP_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
GIMAP_XP_CV <- GIMAP_XP_CV[!duplicated(GIMAP_XP_CV),]
GIMAP_XP_CV  <- GIMAP_XP_CV %>% filter(!is.na(protein_id))

GIMAP <- rbind(GIMAP_XP_CG,GIMAP_XP_CV)
# lookup in gff3 for each species
#write.table(GIMAP$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/GIMAP_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## PCDC
PCDC <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "programmed cell death protein")
PCDC_XP_CG <- unique(PCDC$gene)
PCDC_XP_CG <- as.data.frame(PCDC_XP_CG)
colnames(PCDC_XP_CG)[1] <- "gene"
PCDC_XP_CG <- left_join(PCDC_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
PCDC_XP_CG <- PCDC_XP_CG[!duplicated(PCDC_XP_CG),]  
PCDC_XP_CG <- PCDC_XP_CG %>% filter(!is.na(protein_id))

PCDCCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "programmed cell death protein")
PCDC_XP_CV <- unique(PCDCCV$gene)
PCDC_XP_CV <- as.data.frame(PCDC_XP_CV)
colnames(PCDC_XP_CV)[1] <- "gene"
PCDC_XP_CV <- left_join(PCDC_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
PCDC_XP_CV <- PCDC_XP_CV[!duplicated(PCDC_XP_CV),]
PCDC_XP_CV <- PCDC_XP_CV %>% filter(!is.na(protein_id))

PCDC <- rbind(PCDC_XP_CG,PCDC_XP_CV)
# lookup in gff3 for each species
#write.table(PCDC$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/PCDC_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## Calpain
calpain <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "calpain")
calpain_XP_CG <- unique(calpain$gene)
calpain_XP_CG <- as.data.frame(calpain_XP_CG)
colnames(calpain_XP_CG)[1] <- "gene"
calpain_XP_CG <- left_join(calpain_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
calpain_XP_CG <- calpain_XP_CG[!duplicated(calpain_XP_CG),]  
calpain_XP_CG <- calpain_XP_CG %>% filter(!is.na(protein_id))

calpainCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "calpain")
calpain_XP_CV <- unique(calpainCV$gene)
calpain_XP_CV <- as.data.frame(calpain_XP_CV)
colnames(calpain_XP_CV)[1] <- "gene"
calpain_XP_CV <- left_join(calpain_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
calpain_XP_CV <- calpain_XP_CV[!duplicated(calpain_XP_CV),]
calpain_XP_CV <- calpain_XP_CV %>% filter(!is.na(protein_id))

calpain <- rbind(calpain_XP_CG,calpain_XP_CV)
# lookup in gff3 for each species
#write.table(calpain$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/calpain_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## TRAF
TRAF <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "TNF receptor-associated factor")
TRAF_XP_CG <- unique(TRAF$gene)
TRAF_XP_CG <- as.data.frame(TRAF_XP_CG)
colnames(TRAF_XP_CG)[1] <- "gene"
TRAF_XP_CG <- left_join(TRAF_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
TRAF_XP_CG <- TRAF_XP_CG[!duplicated(TRAF_XP_CG),]  
TRAF_XP_CG <- TRAF_XP_CG %>% filter(!is.na(protein_id))

TRAFCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "TNF receptor-associated factor")
TRAF_XP_CV <- unique(TRAFCV$gene)
TRAF_XP_CV <- as.data.frame(TRAF_XP_CV)
colnames(TRAF_XP_CV)[1] <- "gene"
TRAF_XP_CV <- left_join(TRAF_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
TRAF_XP_CV <- TRAF_XP_CV[!duplicated(TRAF_XP_CV),]
TRAF_XP_CV <- TRAF_XP_CV %>% filter(!is.na(protein_id))

TRAF <- rbind(TRAF_XP_CG,TRAF_XP_CV)
# lookup in gff3 for each species
#write.table(TRAF$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/TRAF_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## IP3R
IP3R <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "inositol 1,4,5-trisphosphate receptor")
IP3R_XP_CG <- unique(IP3R$gene)
IP3R_XP_CG <- as.data.frame(IP3R_XP_CG)
colnames(IP3R_XP_CG)[1] <- "gene"
IP3R_XP_CG <- left_join(IP3R_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
IP3R_XP_CG <- IP3R_XP_CG[!duplicated(IP3R_XP_CG),]  
IP3R_XP_CG <- IP3R_XP_CG %>% filter(!is.na(protein_id))

IP3RCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "inositol 1,4,5-trisphosphate receptor")
IP3R_XP_CV <- unique(IP3RCV$gene)
IP3R_XP_CV <- as.data.frame(IP3R_XP_CV)
colnames(IP3R_XP_CV)[1] <- "gene"
IP3R_XP_CV <- left_join(IP3R_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
IP3R_XP_CV <- IP3R_XP_CV[!duplicated(IP3R_XP_CV),]
IP3R_XP_CV <- IP3R_XP_CV %>% filter(!is.na(protein_id))

IP3R <- rbind(IP3R_XP_CG,IP3R_XP_CV)
# lookup in gff3 for each species
#write.table(IP3R$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IP3R_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## MAP3K
MAP3K <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "mitogen-activated protein kinase kinase kinase")
MAP3K_XP_CG <- unique(MAP3K$gene)
MAP3K_XP_CG <- as.data.frame(MAP3K_XP_CG)
colnames(MAP3K_XP_CG)[1] <- "gene"
MAP3K_XP_CG <- left_join(MAP3K_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
MAP3K_XP_CG <- MAP3K_XP_CG[!duplicated(MAP3K_XP_CG),]  
MAP3K_XP_CG <- MAP3K_XP_CG %>% filter(!is.na(protein_id))

MAP3KCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "mitogen-activated protein kinase kinase kinase")
MAP3K_XP_CV <- unique(MAP3KCV$gene)
MAP3K_XP_CV <- as.data.frame(MAP3K_XP_CV)
colnames(MAP3K_XP_CV)[1] <- "gene"
MAP3K_XP_CV <- left_join(MAP3K_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
MAP3K_XP_CV <- MAP3K_XP_CV[!duplicated(MAP3K_XP_CV),]
MAP3K_XP_CV <- MAP3K_XP_CV %>% filter(!is.na(protein_id))

MAP3K <- rbind(MAP3K_XP_CG,MAP3K_XP_CV)
# lookup in gff3 for each species
#write.table(MAP3K$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/MAP3K_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## MyD88 
MyD88 <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "myeloid differentiation primary response protein MyD88")
MyD88_XP_CG <- unique(MyD88$gene)
MyD88_XP_CG <- as.data.frame(MyD88_XP_CG)
colnames(MyD88_XP_CG)[1] <- "gene"
MyD88_XP_CG <- left_join(MyD88_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
MyD88_XP_CG <- MyD88_XP_CG[!duplicated(MyD88_XP_CG),]  
MyD88_XP_CG <- MyD88_XP_CG %>% filter(!is.na(protein_id))

MyD88CV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "myeloid differentiation primary response protein MyD88")
MyD88_XP_CV <- unique(MyD88CV$gene)
MyD88_XP_CV <- as.data.frame(MyD88_XP_CV)
colnames(MyD88_XP_CV)[1] <- "gene"
MyD88_XP_CV <- left_join(MyD88_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
MyD88_XP_CV <- MyD88_XP_CV[!duplicated(MyD88_XP_CV),]
MyD88_XP_CV <- MyD88_XP_CV %>% filter(!is.na(protein_id))

MyD88 <- rbind(MyD88_XP_CG,MyD88_XP_CV)
# lookup in gff3 for each species
#write.table(MyD88$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/MyD88_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

#CREB
CREB <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "cAMP-responsive element")
CREB_XP_CG <- unique(CREB$gene)
CREB_XP_CG <- as.data.frame(CREB_XP_CG)
colnames(CREB_XP_CG)[1] <- "gene"
CREB_XP_CG <- left_join(CREB_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
CREB_XP_CG <- CREB_XP_CG[!duplicated(CREB_XP_CG),]  
CREB_XP_CG <- CREB_XP_CG %>% filter(!is.na(protein_id))

CREBCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "cAMP-responsive element")
CREB_XP_CV <- unique(CREBCV$gene)
CREB_XP_CV <- as.data.frame(CREB_XP_CV)
colnames(CREB_XP_CV)[1] <- "gene"
CREB_XP_CV <- left_join(CREB_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
CREB_XP_CV <- CREB_XP_CV[!duplicated(CREB_XP_CV),]
CREB_XP_CV <- CREB_XP_CV %>% filter(!is.na(protein_id))

CREB<- rbind(CREB_XP_CG,CREB_XP_CV)
# lookup in gff3 for each species
#write.table(CREB$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/CREB_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

##CASP
CASP <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "caspase-")
CASP_XP_CG <- unique(CASP$gene)
CASP_XP_CG <- as.data.frame(CASP_XP_CG)
colnames(CASP_XP_CG)[1] <- "gene"
CASP_XP_CG <- left_join(CASP_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
CASP_XP_CG <- CASP_XP_CG[!duplicated(CASP_XP_CG),]  
CASP_XP_CG<- CASP_XP_CG %>% filter(!is.na(protein_id))

CASPCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "caspase-")
CASP_XP_CV <- unique(CASPCV$gene)
CASP_XP_CV <- as.data.frame(CASP_XP_CV)
colnames(CASP_XP_CV)[1] <- "gene"
CASP_XP_CV <- left_join(CASP_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
CASP_XP_CV <- CASP_XP_CV[!duplicated(CASP_XP_CV),]
CASP_XP_CV<- CASP_XP_CV %>% filter(!is.na(protein_id))

CASP<- rbind(CASP_XP_CG,CASP_XP_CV)
# lookup in gff3 for each species
#write.table(CASP$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/CASP_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## IFI44
IFI44 <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "interferon-induced protein 44")
IFI44_XP_CG <- unique(IFI44$gene)
IFI44_XP_CG <- as.data.frame(IFI44_XP_CG)
colnames(IFI44_XP_CG)[1] <- "gene"
IFI44_XP_CG <- left_join(IFI44_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
IFI44_XP_CG <- IFI44_XP_CG[!duplicated(IFI44_XP_CG),]  
IFI44_XP_CG<- IFI44_XP_CG %>% filter(!is.na(protein_id))

IFI44CV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "interferon-induced protein 44")
IFI44_XP_CV <- unique(IFI44CV$gene)
IFI44_XP_CV <- as.data.frame(IFI44_XP_CV)
colnames(IFI44_XP_CV)[1] <- "gene"
IFI44_XP_CV <- left_join(IFI44_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
IFI44_XP_CV <- IFI44_XP_CV[!duplicated(IFI44_XP_CV),]
IFI44_XP_CV<- IFI44_XP_CV %>% filter(!is.na(protein_id))

IFI44<- rbind(IFI44_XP_CG,IFI44_XP_CV)
# lookup in gff3 for each species
#write.table(IFI44$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/IFI44_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## TNFR
TNFR <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "^tumor necrosis factor receptor")
TNFR_XP_CG <- unique(TNFR$gene)
TNFR_XP_CG <- as.data.frame(TNFR_XP_CG)
colnames(TNFR_XP_CG)[1] <- "gene"
TNFR_XP_CG <- left_join(TNFR_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
TNFR_XP_CG <- TNFR_XP_CG[!duplicated(TNFR_XP_CG),]  
TNFR_XP_CG<- TNFR_XP_CG %>% filter(!is.na(protein_id))

TNFRCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "^tumor necrosis factor receptor")
TNFR_XP_CV <- unique(TNFRCV$gene)
TNFR_XP_CV <- as.data.frame(TNFR_XP_CV)
colnames(TNFR_XP_CV)[1] <- "gene"
TNFR_XP_CV <- left_join(TNFR_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
TNFR_XP_CV <- TNFR_XP_CV[!duplicated(TNFR_XP_CV),]
TNFR_XP_CV<- TNFR_XP_CV %>% filter(!is.na(protein_id))

TNFR<- rbind(TNFR_XP_CG,TNFR_XP_CV)
# lookup in gff3 for each species
#write.table(TNFR$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/TNFR_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

## Cathepsin
Cathepsin <- C_gig_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "cathepsin")
Cathepsin_XP_CG <- unique(Cathepsin$gene)
Cathepsin_XP_CG <- as.data.frame(Cathepsin_XP_CG)
colnames(Cathepsin_XP_CG)[1] <- "gene"
Cathepsin_XP_CG <- left_join(Cathepsin_XP_CG, C_gig_rtracklayer[,c("gene","protein_id","product")])
Cathepsin_XP_CG <- Cathepsin_XP_CG[!duplicated(Cathepsin_XP_CG),]  
Cathepsin_XP_CG <- Cathepsin_XP_CG %>% filter(!is.na(protein_id))

CathepsinCV <- C_vir_rtracklayer_apop_product_final_product_joined %>% filter(apoptosis_names_query == "cathepsin")
Cathepsin_XP_CV <- unique(CathepsinCV$gene)
Cathepsin_XP_CV <- as.data.frame(Cathepsin_XP_CV)
colnames(Cathepsin_XP_CV)[1] <- "gene"
Cathepsin_XP_CV <- left_join(Cathepsin_XP_CV, C_vir_rtracklayer[,c("gene","protein_id","product")])
Cathepsin_XP_CV <- Cathepsin_XP_CV[!duplicated(Cathepsin_XP_CV),]
Cathepsin_XP_CV <- Cathepsin_XP_CV %>% filter(!is.na(protein_id))

Cathepsin<- rbind(Cathepsin_XP_CG,Cathepsin_XP_CV)
# lookup in gff3 for each species
#write.table(Cathepsin$protein_id, file ="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/OrthoFinder_DATA/OrthoFinder_Analysis/Results_Mar25/Cathepsin_CV_CG_XP.txt", row.names=FALSE, quote= FALSE, col.names = FALSE)

#### Search for proteins in the orthogroup.tsv file ####
# Load tsv

Orthogroups <- read_tsv("/Volumes/My Passport for Mac/OrthoFinder_3_25_2020_Bluewaves_Backup/Results_Mar25/Orthogroups/Orthogroups.tsv",
                        col_names = c("Orthogroup","Elysia_chlorotica",	"Amphibalanus_amphitrite","Drosophila_melanogaster",
                                      "Homo_sapiens",	"Mus_musculus",	"Danio_rerio", "Aplysia_californica",	"Strongylocentrotus_purpuratus","Caenorhabditis_elegans",
                                      "Saccoglossus_kowalevskii",	"Branchiostoma_floridae",	"Ciona_intestinalis",	"Crassostrea_gigas",	"Lottia_gigantea",
                                      "Biomphalaria_glabrata", "Priapulus_caudatus","Lingula_anatina","Octopus_bimaculoides",	"Exaiptasia_pallida",	"Branchiostoma_belcheri",
                                      "Acanthaster_planci",	"Crassostrea_virginica", "Orbicella_faveolata", "Mizuhopecten_yessoensis",
                                      "Pomacea_canaliculata",	"Pocillopora_damicornis",	"Penaeus_vannamei",	"Daphnia_magna","Acropora_millepora","Octopus_sinensis",	"Actinia_tenebrosa"))

# Grep lines based on protein lists above 
TLR_XP_CV_list <- as.list(TLR_XP_CV$protein_id)
TLR_XP_CG_list <- as.list(TLR_XP_CG$protein_id)
IAP_XP_CV_list <- as.list(IAP_XP_CV$protein_id )
IAP_XP_CG_list <- as.list(IAP_XP_CG$protein_id )
GIMAP_XP_CV_list <- as.list(GIMAP_XP_CV$protein_id )
GIMAP_XP_CG_list <- as.list(GIMAP_XP_CG$protein_id )
PCDC_XP_CV_list <- as.list(PCDC_XP_CV$protein_id )
PCDC_XP_CG_list <- as.list(PCDC_XP_CG$protein_id )
calpain_XP_CV_list <- as.list(calpain_XP_CV$protein_id )
calpain_XP_CG_list <- as.list(calpain_XP_CG$protein_id )
TRAF_XP_CV_list <- as.list(TRAF_XP_CV$protein_id )
TRAF_XP_CG_list <- as.list(TRAF_XP_CG$protein_id )
IP3R_XP_CV_list <- as.list(IP3R_XP_CV$protein_id )
IP3R_XP_CG_list <- as.list(IP3R_XP_CG$protein_id )
MAP3K_XP_CV_list <- as.list(MAP3K_XP_CV$protein_id )
MAP3K_XP_CG_list <- as.list(MAP3K_XP_CG$protein_id )
MyD88_XP_CV_list <- as.list(MyD88_XP_CV$protein_id )
MyD88_XP_CG_list <- as.list(MyD88_XP_CG$protein_id )
CREB_XP_CV_list <- as.list(CREB_XP_CV$protein_id)
CREB_XP_CG_list <- as.list(CREB_XP_CG$protein_id)
IFI44_XP_CV_list <- as.list(IFI44_XP_CV$protein_id)
IFI44_XP_CG_list <- as.list(IFI44_XP_CG$protein_id)
TNFR_XP_CV_list <- as.list(TNFR_XP_CV$protein_id)
TNFR_XP_CG_list <- as.list(TNFR_XP_CG$protein_id)
Cathepsin_XP_CV_list <- as.list(Cathepsin_XP_CV$protein_id)
Cathepsin_XP_CG_list <- as.list(Cathepsin_XP_CG$protein_id)

#### TLR analysis ####
TLR_CVlookup <- Orthogroups[grepl(paste(TLR_XP_CV_list,collapse="|"), Orthogroups$Crassostrea_virginica, ignore.case = TRUE),]
# Are any XPs shared
# Hit to four orthogroups: "OG0000019" "OG0012301" "OG0012607" "OG0015175"
TLR_CVlookup$Orthogroup
TLR_CGlookup <- Orthogroups[grepl(paste(TLR_XP_CG_list,collapse="|"), Orthogroups$Crassostrea_gigas, ignore.case = TRUE),]
# Hit to 6 orthgroups: "OG0000019" "OG0005331" "OG0013333" "OG0018123" "OG0019699" "OG0021942"

# Do any proteins hit to multiple orthogroups? Parse and check for duplicates
# Parse and trimws
TLRC <- strsplit(TLR_CVlookup$Crassostrea_virginica, split = ",")
TLRC_parse <-data.frame(Orthogroup = rep(TLR_CVlookup$Orthogroup, sapply(TLRC, length)), protein_id = unlist(TLRC))
# trimws for checking list later 
TLRC_parse$protein_id <- trimws(TLRC_parse$protein_id, which=c("left"))

TLRG <- strsplit(TLR_CGlookup$Crassostrea_gigas, split = ",")
TLRG_parse <-data.frame(Orthogroup = rep(TLR_CGlookup$Orthogroup, sapply(TLRG, length)), protein_id  = unlist(TLRG))
# trimws for checking list later 
TLRG_parse$protein_id <- trimws(TLRG_parse$protein_id, which=c("left"))

# check for duplicated protein names
TLRC_parse[duplicated(TLRC_parse),] # 0 duplicated
TLRG_parse[duplicated(TLRG_parse$protein_id),] # 0 duplicated
# No orthogroups need to be combined

# Were any C_vir or C_gig TLR annotated proteins NOT mapped to Orthogroups?
TLR_CV_notmapped <- left_join(TLR_XP_CV, TLRC_parse, by ="protein_id") %>% filter(is.na(Orthogroup)) %>% View() # Three were not mapped to Orthogroups (need to check these and add to counts list)
TLR_CG_notmapped <- left_join(TLR_XP_CG, TLRG_parse, by ="protein_id") %>% filter(is.na(Orthogroup)) %>% View() # All were mapped

# Full join to get full list of TLR orthogroups
TLR_Orthogroups <- full_join(TLR_CVlookup,TLR_CGlookup)
ncol(TLR_Orthogroups) # 32 
nrow(TLR_Orthogroups) # 9 total orthogroups
TLR_Orthogroups$Orthogroup # "OG0000019" "OG0012301" "OG0012607" "OG0015175" "OG0005331" "OG0013333" "OG0018123" "OG0019699" "OG0021942"

# Get full list of proteins for each species by transposing and uniting
# Transpose the rows and column 
TLR_Orthogroups_transpose <- t(TLR_Orthogroups)
class(TLR_Orthogroups_transpose) # matrix
TLR_Orthogroups_transpose <- as.data.frame(TLR_Orthogroups_transpose)
# unite all columns into one column 
TLR_Orthogroups_transpose_united <- unite(TLR_Orthogroups_transpose[,c(1:9)], full_protein_list, sep=",")

# Lookup gene LOC for each protein XP by using downloaded genome 

# Run each Orthogroup in Batch Entrez to see how they were annotated 
# Write each orthogroup to file
combined_OG0000019 <- unite(TLR_Orthogroups[1,], combined_OG0000019, 2:32, sep=',')

# "OG0012301" 
combined_OG0012301 <- unite(TLR_Orthogroups[2,], combined_OG0012301, 2:32, sep=',')

# also contains leucine-rich repeats and immunoglobulin-like domains protein 

# "OG0012607" 
combined_OG0012607 <- unite(TLR_Orthogroups[3,], combined_OG0012607, 2:32, sep=',')
# "OG0015175" 
combined_OG0015175<- unite(TLR_Orthogroups[4,], combined_OG0015175, 2:32, sep=',')
# "OG0005331" 
combined_OG0005331 <- unite(TLR_Orthogroups[5,], combined_OG0005331, 2:32, sep=',')
# "OG0013333" 
combined_OG0013333<- unite(TLR_Orthogroups[6,], combined_OG0013333, 2:32, sep=',')
# "OG0018123"
combined_OG0018123<- unite(TLR_Orthogroups[7,], combined_OG0018123, 2:32, sep=',')
# "OG0019699" 
combined_OG0019699 <- unite(TLR_Orthogroups[8,], combined_OG0019699, 2:32, sep=',')
# "OG0021942"
combined_OG0021942 <- unite(TLR_Orthogroups[9,], combined_OG0021942, 2:32, sep=',')
# looking at output from NCBI as to which TLRs are in this orthogroup 
# Orthogroup doesn't just contain TLR-annotated proteins, some mucins, uncharacterized proteins, etc. How do I weed these out? Do I prune for only those that are annotated

# Checking how many TLR XPs were assigned for each species
Elysia_chlorotica_TLR  <- unite(TLR_Orthogroups[,2], Elysia_chlorotica_TLR, sep=',')
Elysia_chlorotica_TLR <- Elysia_chlorotica_TLR %>% filter(Elysia_chlorotica_TLR != "NA")
ECTLR <- strsplit(Elysia_chlorotica_TLR$Elysia_chlorotica_TLR, split = ",")
Elysia_chlorotica_TLR_parse <-data.frame(protein_id = rep(Elysia_chlorotica_TLR$Elysia_chlorotica_TLR, sapply(ECTLR, length)), protein_id  = unlist(ECTLR))
# trimws for checking list later 
Elysia_chlorotica_TLR <- Elysia_chlorotica_TLR_parse$protein_id.1
summary(Elysia_chlorotica_TLR) # 19 in TLR orthogroups

apply(TLR_Orthogroups[,c(2:32)], 2, unite, na.rm =TRUE)


