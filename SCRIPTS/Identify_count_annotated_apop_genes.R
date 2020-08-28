#### R script to identify and count annotated apoptosis genes in the C. virginica and C. gigas genomes
# Erin Roberts, 2020
# PhD Candidate University of Rhode Island 

#### Load Packages ####

library(tidyverse)
library(rtracklayer)
library(sqldf)
library(bbplot)
library(viridis) #for colors

#### Import Genome Annotation ####

# Import gff file with rtracklayer for C vir
C_vir_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

# Load list of apoptosis names
Apoptosis_names_list <- c('bcl-2-related protein A1',
'apoptosis-inducing factor 1',
'akt1',
'RAC-alpha serine/threonine-protein kinase-',
'RAC-gamma serine/threonine-protein kinase-',
'methylthioribulose-1-phosphate dehydratase',
'tumor necrosis factor ligand superfamily member 10',
'tumor necrosis factor superfamily member 12',
'cell death regulator Aven',
'BCL2 associated agonist of cell death',
'BAG family molecular chaperone regulator',
'bcl-2 homologous antagonist/killer',
'apoptosis regulator BAX',
'bcl-2-like protein 2',
'bcl-2-like protein 1',
'bcl2 modifying factor',
'Bax inhibitor 1',
'BH3 interacting domain death agonist',
'Bcl-2 interacting killer',
'bcl-2 interacting protein BIM',
'Bik-like killer protein',
'Bcl-2 related ovarian killer',
'CASP8 and FADD like apoptosis regulator',
'transcription factor AP-1',
'caspase activity and apoptosis inhibitor 1',
'DNA fragmentation factor subunit beta',
'adenylate cyclase',
'caspase-1-like',
'caspase-2',
'caspase-3-like',
'caspase-7',
'caspase-8',
'caspase-9',
'caspase-10',
'caspase-11',
'caspase-6',
'caspase-4',
'caspase-5',
'cell division cycle and apoptosis regulator protein 1',
'CD151 antigen',
'protein BTG1',
'interferon alpha',
'caspase activity and apoptosis inhibitor 1',
'baculoviral IAP repeat-containing protein',
'cAMP-responsive element',
'cytochrome c-like',
'death-associated inhibitor of apoptosis',
'tumor necrosis factor receptor superfamily member',
'endonuclease G, mitochondrial',
'FAS-associated death domain protein',
'fas apoptotic inhibitory molecule 1',
'fas cell surface death receptor',
'GTPase IMAP family member',
'harakiri',
'baculoviral IAP repeat-containing protein',
'DNA fragmentation factor subunit alpha',
'interferon-induced protein 44',
'NF-kappa-B inhibitor alpha',
'NF-kappa-B inhibitor epsilon',
'inositol 1,4,5-trisphosphate receptor',
'stress-activated protein kinase JNK',
'lipopolysaccharide-induced tumor necrosis factor-alpha',
'induced myeloid leukemia cell differentiation protein Mcl-1-like',
'mitogen-activated protein kinase 1',
'dual specificity mitogen-activated protein kinase kinase 1',
'mitogen-activated protein kinase kinase kinase',
'transcriptional regulator Myc-A',
'myeloid differentiation primary response protein MyD88',
'phorbol-12-myristate-13-acetate-induced protein 1',
'nuclear factor NF-kappa-B p105 subunit',
'nuclear factor NF-kappa-B p100 subunit',
'transcription factor p65',
'RELB proto-oncogene',
'NF-kB subunit',
'reticuloendotheliosis oncogene',
'anti-apoptotic protein NR13',
'nuclear mitotic apparatus protein 1',
'dynamin-like 120 kDa protein',
'cyclin-dependent kinase 5 activator 1',
'cellular tumor antigen p53',
'programmed cell death protein',
'p53 and DNA damage-regulated protein 1',
'phosphatidylinositol 3-kinase',
'putative inhibitor of apoptosis',
'cAMP-dependent protein kinase',
'protein kinase C delta type',
'protein kinase C iota type',
'BCL2 binding component 3',
'cdc42 homolog',
'ras-like GTP-binding protein rho',
'rho-related GTP-binding protein RhoE-like',
'ras-related C3 botulinum toxin substrate 1',
'rho-related protein racA',
'mitochondrial Rho GTPase 1',
'receptor-interacting serine/threonine-protein kinase 1',
'receptor-interacting serine/threonine-protein kinase 4',
'diablo homolog, mitochondrial',
'toll-like receptor',
'tumor necrosis factor',
'lymphotoxin-alpha',
'CD40 ligand',
'TNFRSF1A associated via death domain',
'TNF receptor-associated factor',
'E3 ubiquitin-protein ligase XIAP',
'netrin receptor',
'neurotrophic receptor tyrosine kinase 1',
'sonic hedgehog receptor',
'receptor-interacting serine/threonine-protein kinase',
'mixed lineage kinase domain',
'heat shock protein',
'E3 ubiquitin-protein ligase CHIP',
'tumor necrosis factor alpha-induced protein 3',
'protein phosphatase 1B',
'aurora kinase A',
'glutathione peroxidase 4',
'gasdermin',
'poly \\[ADP-ribose]\\ polymerase 1-like',
'macrophage migration inhibitory factor',
'hexokinase-1',
'Raf-1 protooncogene serine/threonine kinase',
'elastase, neutrophil expressed',
'cathepsin',
'PRKC apoptosis WT1 regulator protein',
'apoptosis-stimulating of p53 protein 1',
'apoptosis-stimulating of p53 protein 2',
'apoptosis inhibitor 5',
'apoptotic chromatin condensation inducer in the nucleus',
'high mobility group box 1',
'ceramide synthase',
'cyclic AMP-responsive element-binding protein',
'cell death-inducing p53-target protein 1',
'TP53-binding protein 1',
'p53-induced death domain-containing protein 1',
'death domain-containing protein CRADD',
'p63',
'p73',
 'interferon regulatory factor',
 'stimulator of interferon genes',
'interleukin 17-like protein',
'interleukin-17',
'interleukin-1 receptor-associated kinase 4',
'interleukin-1 receptor-associated kinase 1-binding protein',
'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
'sterile alpha and TIR motif-containing protein',
'tyrosine-protein kinase JAK2',
'signal transducer and activator of transcription',
'serine-protein kinase ATM',
'MAP kinase-activating death domain protein',
'death domain-associated protein 6',
'leucine-rich repeat and death domain-containing protein',
'serine/threonine-protein kinase/endoribonuclease IRE1',
'eukaryotic translation initiation factor 2-alpha kinase 3',
'growth arrest and DNA damage-inducible protein GADD45',
'calpain',
'tyrosine-protein phosphatase non-receptor type 13',
'HTRA2',
'3-phosphoinositide-dependent protein kinase 1',
'dronc',
'pyrin',
'proto-oncogene c-Rel',
'leukocyte elastase inhibitor',
'protein patched homolog 1',
'cyclic AMP-dependent transcription factor ATF-4',
'Siva',
'immune-associated nucleotide-binding protein')

# Make this list more generic so it hits the previous ones subset and can be used to count genes in gene families
Apoptosis_names_df <- data.frame(product=c(# removing this because duplicated with bcl-2 'bcl-2-related protein A1',
'^apoptosis-inducing factor 1',
'nuclear apoptosis-inducing factor 1',
'akt1',
'RAC-alpha',
'RAC-gamma',
'methylthioribulose-1-phosphate dehydratase',
'^tumor necrosis factor receptor',
'^tumor necrosis factor ligand',
'tumor necrosis factor alpha-induced protein',
'lipopolysaccharide-induced tumor necrosis factor-alpha',
'cell death regulator Aven',
'bcl-2',
'BAX',
'BAG family molecular chaperone regulator',
'BH3 interacting domain death agonist',
'Bik-like killer protein',
'CASP8 and FADD like apoptosis regulator',
'transcription factor AP-1',
'adenylate cyclase',
'caspase-',
'caspase activity and apoptosis inhibitor 1',
'cell division cycle and apoptosis regulator protein 1',
'CD151 antigen',
'protein BTG1',
'baculoviral IAP',
'E3 ubiquitin-protein ligase XIAP',
'cAMP-responsive element',
'cytochrome c',
'endonuclease G, mitochondrial',
'FAS-associated death domain protein',
'fas apoptotic inhibitory molecule 1',
'fas cell surface death receptor',
'GTPase IMAP family member',
'harakiri',
'DNA fragmentation factor',
'interferon-induced protein 44',
'interferon alpha-inducible protein 27',
'NF-kappa-B',
'inositol 1,4,5-trisphosphate receptor',
'stress-activated protein kinase JNK',
'induced myeloid leukemia cell differentiation protein Mcl-1',
'mitogen-activated protein kinase kinase kinase',
'mitogen-activated protein kinase-',
'transcriptional regulator Myc-A',
'myeloid differentiation primary response protein MyD88',
'phorbol-12-myristate-13-acetate-induced protein 1',
'transcription factor p65',
'RELB proto-oncogene',
'reticuloendotheliosis oncogene',
'anti-apoptotic protein NR13',
'nuclear mitotic apparatus protein 1',
'dynamin-like 120 kDa protein',
'cyclin-dependent kinase 5 activator 1',
'cellular tumor antigen p53',
'programmed cell death protein',
'p53 and DNA damage-regulated protein 1',
'phosphatidylinositol 3-kinase',
'cAMP-dependent protein kinase',
'protein kinase C',
'BCL2 binding component 3',
'cdc42 homolog',
'ras-related',
'rho-related',
'ras-like',
'receptor-interacting serine/threonine-protein kinase',
'diablo homolog, mitochondrial',
'toll-like receptor',
'lymphotoxin-alpha',
'CD40 ligand',
'TNFRSF1A associated via death domain',
'TNF receptor-associated factor',
'netrin receptor',
'neurotrophic receptor tyrosine kinase 1',
'sonic hedgehog receptor',
'mixed lineage kinase domain',
'heat shock protein',
'E3 ubiquitin-protein ligase CHIP',
'protein phosphatase 1B',
'aurora kinase A',
'glutathione peroxidase 4',
'gasdermin',
'poly \\[ADP-ribose]\\ polymerase 1',
'macrophage migration inhibitory factor',
'hexokinase-1',
'Raf-1 protooncogene serine/threonine kinase',
'elastase, neutrophil expressed',
'cathepsin',
'PRKC apoptosis WT1 regulator protein',
'apoptosis-stimulating of p53 protein',
'apoptosis inhibitor 5',
'apoptotic chromatin condensation inducer in the nucleus',
"high mobility group box 1",
"ceramide synthase",
'inhibitor of apoptosis',
'cyclic AMP-responsive element-binding protein',
'cell death-inducing p53-target protein 1',
'TP53-binding protein 1',
'p53-induced death domain-containing protein 1',
'death domain-containing protein CRADD',
'p63',
'p73',
'interferon regulatory factor',
'stimulator of interferon genes',
'interleukin 17-like protein',
'interleukin-17',
'interleukin-1 receptor-associated kinase 4',
'interleukin-1 receptor-associated kinase 1-binding protein',
'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
'sterile alpha and TIR motif-containing protein',
'tyrosine-protein kinase JAK2',
'signal transducer and activator of transcription',
'serine-protein kinase ATM',
'MAP kinase-activating death domain protein',
'death domain-associated protein 6',
'leucine-rich repeat and death domain-containing protein',
'serine/threonine-protein kinase/endoribonuclease IRE1',
'eukaryotic translation initiation factor 2-alpha kinase 3',
'growth arrest and DNA damage-inducible protein GADD45',
'calpain',
'tyrosine-protein phosphatase non-receptor type 13',
'HTRA2',
'3-phosphoinositide-dependent protein kinase 1',
'dronc',
'pyrin',
'proto-oncogene c-Rel',
'leukocyte elastase inhibitor',
'protein patched homolog 1',
'cyclic AMP-dependent transcription factor ATF-4',
'dual specificity mitogen-activated protein kinase kinase 1',
'mitogen-activated protein kinase 1',
'mitochondrial Rho GTPase 1',
'Siva',
'immune-associated nucleotide-binding protein',
'uncharacterized')) # all the uncharacterized are IAPs 

## Note: August 5th, 2020: removed IAP search terms from list 
# Which terms to remove?
#View(BIR_XP_gff_species)
BIR_XP_gff_CG_CV %>% filter(Species == "Crassostrea_virginica") %>% distinct(product) %>% arrange(product) %>% View()

# IAP terms to remove: 
# 'baculoviral IAP repeat-containing protein',
# 'baculoviral IAP',
# 'death-associated inhibitor of apoptosis',
# 'E3 ubiquitin-protein ligase XIAP',
# 'putative inhibitor of apoptosis',
IAP_remove <- c('baculoviral IAP repeat-containing protein',
                'death-associated inhibitor of apoptosis',
                'E3 ubiquitin-protein ligase XIAP',
                'baculoviral IAP',
                'inhibitor of apoptosis')

Apoptosis_names_list <- setdiff(Apoptosis_names_list, IAP_remove)
# don't repeat for Apoptosis_names_df

#### Grep Apoptosis protein names in genome files ####

# C virginica 
C_vir_rtracklayer_filtered <- C_vir_rtracklayer %>% filter(type == "mRNA")
C_vir_rtracklayer_apop_product <- C_vir_rtracklayer_filtered [grepl(paste(Apoptosis_names_list,collapse="|"), 
                                                                    C_vir_rtracklayer_filtered$product, ignore.case = TRUE),]

## August 5th, 2020 join in with manually curated, haplotig collapsed list of C.virginica IAPs so I'm searching for the correct IAPs
# Use dataframe created in IAP_GIMAP_Gene_Family_Expansion.R (lines 229-233)
nrow(BIR_XP_gff_CV_uniq_XP_XM) # 164
class(BIR_XP_gff_CV_uniq_XP_XM$transcript_id)
BIR_XP_gff_CV_uniq_XP_XM$ID <- as.character(BIR_XP_gff_CV_uniq_XP_XM$ID)
# get row matches in the full mRNA dataframe so I can simply rbind the IAP matches to the rest of the identified IAP genes 
BIR_XP_gff_CV_uniq_XP_XM_apop_product <- C_vir_rtracklayer_filtered[C_vir_rtracklayer_filtered$ID %in% BIR_XP_gff_CV_uniq_XP_XM$ID,] # Join by ID

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and NOT other kDA protein names so I can keep all heat shock proteins
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product[ 
  !grepl("dual specificity protein phosphatase 1B-like", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
  !grepl("kunitz-type", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("complement C1q", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mannosyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("S-acyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("activator of 90 kDa", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("zinc finger protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
    !grepl("DDB_G0272098",C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("WD repeat-containing protein WRAP73", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("tumor protein p63-regulated gene 1-like protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mitogen-activated protein kinase kinase kinase kinase 3", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE),] # get rid of MAP4K

## Add curated IAPs to the apop final
C_vir_rtracklayer_apop_product_final <- rbind(C_vir_rtracklayer_apop_product_final, BIR_XP_gff_CV_uniq_XP_XM_apop_product)

nrow(C_vir_rtracklayer_apop_product_final) # 1312
#View(C_vir_rtracklayer_apop_product_final[26][!duplicated(C_vir_rtracklayer_apop_product_final$product),]) # NOT VALID TO REMOVE BASED ON DUPLICATED PRODUCT, ONLY DUPLICATED TRANSCRIPT

# remove rows with duplicated transcript_ids
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product_final[!duplicated(C_vir_rtracklayer_apop_product_final$transcript_id),]
nrow(C_vir_rtracklayer_apop_product_final) # 1290

# Export table as text file
write.table(C_vir_rtracklayer_apop_product_final, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/ PAPER DRAFTS AND OUTLINES/Additional:Supplementary_Data_files_for_paper/Supplementary_Table_3_C_virginica.txt",
            row.names = FALSE)

### Checked genes were added correctly and confirmed with previous results by comparing DF merged list to my previous research
# of which genes are there in my "Updated_APOPTOSIS_GENES, DOMAINS....".xlsx sheet., then checked by "Supplementary Table 2. C. virginica apoptosis.xlsx" table that I have
# Checking to make sure all YES proteins are actually present 

# How many unique gene LOC are in this list (this is going in summary paragraph for paper)
C_vir_rtracklayer_apop_product_final_LOC <- C_vir_rtracklayer_apop_product_final[!duplicated(C_vir_rtracklayer_apop_product_final$gene),]
nrow(C_vir_rtracklayer_apop_product_final_LOC )#676

### Add Gene family Key with gene family names in Apoptosis_names_df  ####
# Subset apoptosis grepl data frame and Apoptosis_names_df and check nrows
C_vir_rtracklayer_apop_product_final_product <- C_vir_rtracklayer_apop_product_final[,c("Name","product","gene","transcript_id")]
Apoptosis_names_df$rownames <- rownames(Apoptosis_names_df)
nrow(C_vir_rtracklayer_apop_product_final_product) #1290

# Make new index
idx2 <- sapply(Apoptosis_names_df$product, grep, C_vir_rtracklayer_apop_product_final_product$product)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))

# create df_merged which has duplicated product column for gene name hits, first column is the index name
df_merged <- cbind(Apoptosis_names_df[unlist(idx1),,drop=F], C_vir_rtracklayer_apop_product_final_product[unlist(idx2),,drop=F])
head(df_merged)
df_merged[is.na(df_merged),]
nrow(df_merged) # 1290
View(df_merged)
# change column names
colnames(df_merged)[1] <- "apoptosis_names_query" 
# subset rownames column
df_merged <- df_merged[,-2]
head(df_merged)

# Check df_merged apoptosis gene family name assignments here
setdiff(C_vir_rtracklayer_apop_product_final_product$transcript_id,df_merged$transcript_id) # none missing
setdiff(df_merged$transcript_id,C_vir_rtracklayer_apop_product_final_product$transcript_id)
df_merged[duplicated(df_merged$transcript_id),]

# check for NAs (meaning duplicate term hits)
df_merged %>% filter(is.na(product)) %>% View() # NO NA's yay! 
df_merged %>% filter(is.na(apoptosis_names_query)) %>% View() # NO NA's yay!!
nrow(df_merged) #1290 ( no duplicates)

# check that all transcript_ids originally found in search are still there
setdiff(C_vir_rtracklayer_apop_product_final_product$transcript_id, df_merged$transcript_id) # no differences, none misssing
setdiff(df_merged$transcript_id, C_vir_rtracklayer_apop_product_final_product$transcript_id) # no differences

# merge new index with old
C_vir_rtracklayer_apop_product_final_product_joined <- full_join(C_vir_rtracklayer_apop_product_final_product, df_merged) # Joining, by = c("Name", "product", "gene", "transcript_id")
#View(C_vir_rtracklayer_apop_product_final_product_joined)
nrow(C_vir_rtracklayer_apop_product_final_product_joined) #1290 rows PERFECT

# recode  mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog, transcript variant X1 (X2,X3)
class(C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query) # factor
C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query <- as.character(C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query)
C_vir_rtracklayer_apop_product_final_product_joined <- within(C_vir_rtracklayer_apop_product_final_product_joined, apoptosis_names_query[ apoptosis_names_query == 'mitogen-activated protein kinase kinase kinase' & product == 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog, transcript variant X1'] <- 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog')
C_vir_rtracklayer_apop_product_final_product_joined <- within(C_vir_rtracklayer_apop_product_final_product_joined, apoptosis_names_query[ apoptosis_names_query == 'mitogen-activated protein kinase kinase kinase' & product == 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog, transcript variant X2'] <- 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog')
C_vir_rtracklayer_apop_product_final_product_joined <- within(C_vir_rtracklayer_apop_product_final_product_joined, apoptosis_names_query[ apoptosis_names_query == 'mitogen-activated protein kinase kinase kinase' & product == 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog, transcript variant X3'] <- 'mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog')
nrow(C_vir_rtracklayer_apop_product_final_product_joined) # 1290

# Recode Gene Family Members to be all included together if they have several different names
C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query <- recode(C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query, 
      # the RHO superfamily to all be included together                                                                                 
      "cdc42 homolog" = "Rho superfamily", "ras-related"= "Rho superfamily", "ras-like"="Rho superfamily", "rho-related"= "Rho superfamily",
      # IAP repeat containing and putative inhibitor of apoptosis to be together in IAP family
      "inhibitor of apoptosis" = "Inhibitor of apoptosis", "IAP" = "Inhibitor of apoptosis", "uncharacterized" = "Inhibitor of apoptosis",
      # Protein BTG1
      "protein BTG1" = "B-cell translocation gene 1",
      # fas associated death domain
      "FAS-associated death domain protein" = "fas-associated death domain protein", 
      # Nf-kappa B family includes transcription factor p65 (aka RELA) and RELB proto-oncogene 
      "transcription factor p65"="NF-kappa-B",
      "RELB proto-oncogene "="NF-kappa-B",
      #recode Bcl-2 family
      "BAX"="bcl-2",
      #recode TNF family 
      "CD40 ligand"="tumor necrosis factor")

# NOTE : RAC-alpha and RAC-gamma are not part of the Rho superfamily
head(C_vir_rtracklayer_apop_product_final_product_joined)

### How many transcripts in each gene?
C_vir_rtracklayer_apop_product_final_product_joined_by_gene <- C_vir_rtracklayer_apop_product_final_product_joined %>%
  group_by(gene) %>% summarize(number_transcripts_per_gene = n() )
nrow(C_vir_rtracklayer_apop_product_final_product_joined_by_gene) #656 genes
# Match with names of gene family members
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_product <- left_join(C_vir_rtracklayer_apop_product_final_product_joined_by_gene, C_vir_rtracklayer[,c("gene","product")], by = "gene")
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_product <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_product[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_product$gene),]
View(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_product)

### How many transcripts in each overall gene family (combined across all genes in the family)?
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_family <- C_vir_rtracklayer_apop_product_final_product_joined %>%
  group_by(apoptosis_names_query) %>% summarize(number_transcripts_per_family = n() )

## How many genes in gene family 
#Match that gene ID with the gene name, (just name without the transcript variant)
# separate "transcript variant" by final comma
C_vir_rtracklayer_apop_product_final_product_joined_split_product <- separate(C_vir_rtracklayer_apop_product_final_product_joined, product, into = c("gene_name", "transcript_name"), sep =", transcript variant")

# subset for unique gene
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique <- unique(C_vir_rtracklayer_apop_product_final_product_joined_split_product[,c("gene","gene_name")])

# join just the gene (LOC) so the product name doesn't create duplicate columns
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name <- left_join(C_vir_rtracklayer_apop_product_final_product_joined_by_gene, C_vir_rtracklayer_apop_product_final_product_joined_split_product[,c("gene","gene_name","apoptosis_names_query")], by = c("gene"))
#remove duplicates
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name$gene),]

## Count number of genes in gene family
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>% summarize(gene_family_members = n())
View(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families)

## Number of genes with same gene name
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(gene_name) %>% summarize(gene_duplicates = n())
#View(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates)

# avg transcripts per gene family
C_vir_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>%
        summarize(avg_transcripts_per_family = mean(number_transcripts_per_gene))

## Join genes per gene family and transcripts for gene family into one table
C_vir_genes_transcripts_per_gene_family <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families,  C_vir_rtracklayer_apop_product_final_product_joined_by_gene_family)
C_vir_genes_transcripts_per_gene_family_transcripts_per_gene <- full_join(C_vir_genes_transcripts_per_gene_family, C_vir_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family)
View(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene)

# Save output to file
save(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates, C_vir_genes_transcripts_per_gene_family, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_vir_annotation_gene_family_info.RData")

#### Repeat Analysis using Crassostrea gigas ####

# Apoptosis names in C. gigas genome are lowercase like C.virginica
Apoptosis_names_list_CG <- c('bcl-2-related protein A1',
                          'apoptosis-inducing factor 1',
                          'akt1',
                          'RAC-alpha serine/threonine-protein kinase-',
                          'RAC-gamma serine/threonine-protein kinase-',
                          'methylthioribulose-1-phosphate dehydratase',
                          'tumor necrosis factor ligand superfamily member 10',
                          'tumor necrosis factor superfamily member 12',
                          'cell death regulator Aven',
                          'BCL2 associated agonist of cell death',
                          'BAG family molecular chaperone regulator',
                          'bcl-2 homologous antagonist/killer',
                          'apoptosis regulator BAX',
                          'bcl-2-like protein 2',
                          'bcl-2-like protein 1',
                          'bcl2 modifying factor',
                          'Bax inhibitor 1',
                          'BH3 interacting domain death agonist',
                          'Bcl-2 interacting killer',
                          'bcl-2 interacting protein BIM',
                          'Bik-like killer protein',
                          'Bcl-2 related ovarian killer',
                          'CASP8 and FADD like apoptosis regulator',
                          'transcription factor AP-1',
                          'caspase activity and apoptosis inhibitor 1',
                          'DNA fragmentation factor subunit beta',
                          'adenylate cyclase',
                          'caspase-1-like',
                          'caspase-2',
                          'caspase-3-like',
                          'caspase-7',
                          'caspase-8',
                          'caspase-9',
                          'caspase-10',
                          'caspase-11',
                          'caspase-6',
                          'caspase-4',
                          'caspase-5',
                          'cell division cycle and apoptosis regulator protein 1',
                          'CD151 antigen',
                          'B-cell translocation gene 1',
                          'interferon alpha',
                          'caspase activity and apoptosis inhibitor 1',
                          'baculoviral IAP repeat-containing protein',
                          'cAMP-responsive element',
                          'cytochrome c',
                          'death-associated inhibitor of apoptosis',
                          'tumor necrosis factor receptor superfamily member',
                          'endonuclease G, mitochondrial',
                          'FAS-associated death domain protein',
                          'fas apoptotic inhibitory molecule 1',
                          'fas cell surface death receptor',
                          'GTPase IMAP family member',
                          'harakiri',
                          'baculoviral IAP repeat-containing protein',
                          'DNA fragmentation factor subunit alpha',
                          'interferon-induced protein 44',
                          'NF-kappa-B inhibitor alpha',
                          'NF-kappa-B inhibitor epsilon',
                          'inositol 1,4,5-trisphosphate receptor',
                          'stress-activated protein kinase JNK',
                          'lipopolysaccharide-induced tumor necrosis factor-alpha',
                          'induced myeloid leukemia cell differentiation protein Mcl-1',
                          'mitogen-activated protein kinase kinase kinase',
                          'mitogen-activated protein kinase 1',
                          'dual specificity mitogen-activated protein kinase kinase 1',
                          'mitogen-activated protein kinase kinase kinase 7',
                          'transcriptional regulator Myc-A',
                          'myeloid differentiation primary response protein MyD88',
                          'phorbol-12-myristate-13-acetate-induced protein 1',
                          'nuclear factor NF-kappa-B p105 subunit',
                          'nuclear factor NF-kappa-B p100 subunit',
                          'transcription factor p65',
                          'RELB proto-oncogene',
                          'NF-kB subunit',
                          'reticuloendotheliosis oncogene',
                          'anti-apoptotic protein NR13',
                          'nuclear mitotic apparatus protein 1',
                          'dynamin-like 120 kDa protein',
                          'cyclin-dependent kinase 5 activator 1',
                          'cellular tumor antigen p53',
                          'programmed cell death protein',
                          'p53 and DNA damage-regulated protein 1',
                          'phosphatidylinositol 3-kinase',
                          'putative inhibitor of apoptosis',
                          'cAMP-dependent protein kinase',
                          'protein kinase C delta type',
                          'protein kinase C iota type',
                          'BCL2 binding component 3',
                          'cdc42 homolog',
                          'ras-like GTP-binding protein rho',
                          'rho-related GTP-binding protein RhoE',
                          'ras-related C3 botulinum toxin substrate 1',
                          'rho-related protein racA',
                          'mitochondrial Rho GTPase 1',
                          'receptor-interacting serine/threonine-protein kinase 1',
                          'receptor-interacting serine/threonine-protein kinase 4',
                          'diablo homolog, mitochondrial',
                          'toll-like receptor',
                          'tumor necrosis factor',
                          'lymphotoxin-alpha',
                          'CD40 ligand',
                          'TNFRSF1A associated via death domain',
                          'TNF receptor-associated factor',
                          'E3 ubiquitin-protein ligase XIAP',
                          'netrin receptor',
                          'neurotrophic receptor tyrosine kinase 1',
                          'sonic hedgehog receptor',
                          'receptor-interacting serine/threonine-protein kinase',
                          'mixed lineage kinase domain',
                          'heat shock protein',
                          'E3 ubiquitin-protein ligase CHIP',
                          'tumor necrosis factor alpha-induced protein 3',
                          'protein phosphatase 1B',
                          'aurora kinase A',
                          'glutathione peroxidase 4',
                          'gasdermin',
                          'poly \\[ADP-ribose]\\ polymerase 1-like',
                          'macrophage migration inhibitory factor',
                          'hexokinase-1',
                          'Raf-1 protooncogene serine/threonine kinase',
                          'elastase, neutrophil expressed',
                          'cathepsin',
                          'PRKC apoptosis WT1 regulator protein',
                          'apoptosis-stimulating of p53 protein 1',
                          'apoptosis-stimulating of p53 protein 2',
                          'apoptosis inhibitor 5',
                          'apoptotic chromatin condensation inducer in the nucleus',
                          "high mobility group box 1",
                          'ceramide synthase',
                          'cyclic AMP-responsive element-binding protein',
                          'cell death-inducing p53-target protein 1',
                          'TP53-binding protein 1',
                          'p53-induced death domain-containing protein 1',
                          'death domain-containing protein CRADD',
                          'p63',
                          'p73',
                          'interferon regulatory factor',
                          'stimulator of interferon genes',
                          'interleukin 17-like protein',
                          'interleukin-17',
                          'interleukin-1 receptor-associated kinase 4',
'interleukin-1 receptor-associated kinase 1-binding protein',
'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
'sterile alpha and TIR motif-containing protein',
'tyrosine-protein kinase JAK2',
'signal transducer and activator of transcription',
'serine-protein kinase ATM',
'MAP kinase-activating death domain protein',
'death domain-associated protein 6',
'leucine-rich repeat and death domain-containing protein',
'serine/threonine-protein kinase/endoribonuclease IRE1',
'eukaryotic translation initiation factor 2-alpha kinase 3',
'growth arrest and DNA damage-inducible protein GADD45',
'calpain',
'tyrosine-protein phosphatase non-receptor type 13',
'HTRA2',
'3-phosphoinositide-dependent protein kinase 1',
'dronc',
'pyrin',
'proto-oncogene c-Rel',
'leukocyte elastase inhibitor',
'protein patched homolog 1',
'cyclic AMP-dependent transcription factor ATF-4',
'Siva',
'immune-associated nucleotide-binding protein')

# Make this list more generic so it hits the previous ones subset and can be used to count genes in gene families
Apoptosis_names_df_CG <- data.frame(product=c(     # removing this because duplicated with bcl-2 'bcl-2-related protein A1',
                                           '^apoptosis-inducing factor 1',
                                           'nuclear apoptosis-inducing factor 1',
                                           'akt1',
                                           'RAC-alpha',
                                           'RAC-gamma',
                                           'methylthioribulose-1-phosphate dehydratase',
                                           '^tumor necrosis factor receptor',
                                           '^tumor necrosis factor ligand',
                                           '^tumor necrosis factor$',# dollar sign indicates the end of the string
                                           'tumor necrosis factor alpha-induced protein',
                                           'lipopolysaccharide-induced tumor necrosis factor-alpha',
                                           'cell death regulator Aven',
                                           'bcl-2',
                                           'BAX',
                                           'BAG family molecular chaperone regulator',
                                           'BH3 interacting domain death agonist',
                                           'Bik-like killer protein',
                                           'CASP8 and FADD like apoptosis regulator',
                                           'transcription factor AP-1',
                                           'adenylate cyclase',
                                           'caspase-',
                                           'caspase activity and apoptosis inhibitor 1',
                                           'cell division cycle and apoptosis regulator protein 1',
                                           'CD151 antigen',
                                           'B-cell translocation gene 1',
                                           'baculoviral IAP',
                                           'E3 ubiquitin-protein ligase XIAP',
                                           'cAMP-responsive element',
                                           'cytochrome c',
                                           'endonuclease G, mitochondrial',
                                           'FAS-associated death domain protein',
                                           'fas-associated death domain protein', #as with Heat shock protein beta will need to merge these
                                           'fas apoptotic inhibitory molecule 1',
                                           'fas cell surface death receptor',
                                           'GTPase IMAP family member',
                                           'harakiri',
                                           'DNA fragmentation factor',
                                           'interferon-induced protein 44',
                                           'interferon alpha-inducible protein 27',
                                           'NF-kappa-B',
                                           'inositol 1,4,5-trisphosphate receptor',
                                           'stress-activated protein kinase JNK',
                                           'induced myeloid leukemia cell differentiation protein Mcl-1',
                                           'mitogen-activated protein kinase kinase kinase',
                                           'mitogen-activated protein kinase-',
                                           'transcriptional regulator Myc-A',
                                           'myeloid differentiation primary response protein MyD88',
                                           'phorbol-12-myristate-13-acetate-induced protein 1',
                                           'transcription factor p65',
                                           'RELB proto-oncogene',
                                           'reticuloendotheliosis oncogene',
                                           'anti-apoptotic protein NR13',
                                           'nuclear mitotic apparatus protein 1',
                                           'dynamin-like 120 kDa protein',
                                           'cyclin-dependent kinase 5 activator 1',
                                           'cellular tumor antigen p53',
                                           'programmed cell death protein',
                                           'p53 and DNA damage-regulated protein 1',
                                           'phosphatidylinositol 3-kinase',
                                           'cAMP-dependent protein kinase',
                                           'protein kinase C',
                                           'BCL2 binding component 3',
                                           'cdc42 homolog',
                                           'ras-related',
                                           'rho-related',
                                           'ras-like',
                                           'receptor-interacting serine/threonine-protein kinase',
                                           'diablo homolog, mitochondrial',
                                           'toll-like receptor',
                                           'lymphotoxin-alpha',
                                           'CD40 ligand',
                                           'TNFRSF1A associated via death domain',
                                           'TNF receptor-associated factor',
                                           'netrin receptor',
                                           'neurotrophic receptor tyrosine kinase 1',
                                           'sonic hedgehog receptor',
                                           'mixed lineage kinase domain',
                                           'heat shock protein',
                                           'Heat shock protein', #will need to recode these later to be in the same group, case sensitivity caused problem
                                           'E3 ubiquitin-protein ligase CHIP',
                                           'protein phosphatase 1B',
                                           'aurora kinase A',
                                           'glutathione peroxidase 4',
                                           'gasdermin',
                                           'poly \\[ADP-ribose]\\ polymerase 1',
                                           'macrophage migration inhibitory factor',
                                           'hexokinase-1',
                                           'Raf-1 protooncogene serine/threonine kinase',
                                           'elastase, neutrophil expressed',
                                           'cathepsin',
                                           'PRKC apoptosis WT1 regulator protein',
                                           'apoptosis-stimulating of p53 protein',
                                           'apoptosis inhibitor 5',
                                           'apoptotic chromatin condensation inducer in the nucleus',
                                           "high mobility group box 1",
                                           "ceramide synthase",
                                           'inhibitor of apoptosis',
                                           'cyclic AMP-responsive element-binding protein',
                                           'cell death-inducing p53-target protein 1',
                                           'TP53-binding protein 1',
                                           'p53-induced death domain-containing protein 1',
                                           'death domain-containing protein CRADD',
                                           'p63',
                                           'p73',
                                           'interferon regulatory factor',
                                           'stimulator of interferon genes',
                                           'interleukin 17-like protein',
                                           'interleukin-17',
                                           'interleukin-1 receptor-associated kinase 4',
                                           'interleukin-1 receptor-associated kinase 1-binding protein',
                                           'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
                                           'sterile alpha and TIR motif-containing protein',
                                           'tyrosine-protein kinase JAK2',
                                           'signal transducer and activator of transcription',
                                           'serine-protein kinase ATM',
                                           'MAP kinase-activating death domain protein',
                                           'death domain-associated protein 6',
                                           'leucine-rich repeat and death domain-containing protein',
                                           'serine/threonine-protein kinase/endoribonuclease IRE1',
                                           'eukaryotic translation initiation factor 2-alpha kinase 3',
                                           'growth arrest and DNA damage-inducible protein GADD45',
                                           'calpain',
                                           'tyrosine-protein phosphatase non-receptor type 13',
                                           'HTRA2',
                                           '3-phosphoinositide-dependent protein kinase 1',
                                           'Dronc',
                                           'pyrin',
                                           'proto-oncogene c-Rel',
                                           'leukocyte elastase inhibitor',
                                           'protein patched homolog 1',
                                           'cyclic AMP-dependent transcription factor ATF-4',
                                           'dual specificity mitogen-activated protein kinase kinase 1',
                                           'mitogen-activated protein kinase 1',
                                           'mitochondrial Rho GTPase 1',
                                           'Siva',
                                           'immune-associated nucleotide-binding protein',
                                           'uncharacterized')) # all the uncharacterized will be IAP
## Note: August 5th, 2020: removed IAP search terms from list 
# Which terms to remove?
#View(BIR_XP_gff_species)
BIR_XP_gff_CG_CV %>% filter(Species == "Crassostrea_gigas") %>% distinct(product) %>% arrange(product) %>% View()

# IAP terms to remove: 
# 'baculoviral IAP repeat-containing protein',
# 'baculoviral IAP',
# 'death-associated inhibitor of apoptosis',
# 'E3 ubiquitin-protein ligase XIAP',
# 'putative inhibitor of apoptosis',
IAP_remove <- c('baculoviral IAP repeat-containing protein',
                'death-associated inhibitor of apoptosis',
                'E3 ubiquitin-protein ligase XIAP',
                'baculoviral IAP',
                'inhibitor of apoptosis')

Apoptosis_names_list_CG <- setdiff(Apoptosis_names_list_CG, IAP_remove)
# don't repeat for Apoptosis_names_df

# Import gff file, using new version of genome annotation
  C_gig_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/GCF_000297895.1_oyster_v9_genomic.gff")
  C_gig_rtracklayer <- as.data.frame(C_gig_rtracklayer)
  
#### Grep Apoptosis protein names in genome files ####
  
## Note unlike C. vir file, C. gigas file has gene names in the "product" column, and still includes comma with transcript info for some
  # than a comman and the transcript name
  # All transcripts from the same gene share the LOC ID in the "gene" column, though the Name column differs

# C gigas, filter for rows that have NA gene, this will also keep all the lines with transcript information. Filter out gbkey CDS 
  # because I only care about the number of transcripts and genes
C_gig_rtracklayer_filtered <- C_gig_rtracklayer %>% filter(type =="mRNA")

# only genes have the description, the mRNA id's will match with the ID
C_gig_rtracklayer_apop_product <- C_gig_rtracklayer_filtered[grepl(paste(Apoptosis_names_list_CG,collapse="|"), 
                                                          C_gig_rtracklayer_filtered$product, ignore.case = TRUE),]
nrow(C_gig_rtracklayer_apop_product) #884
  
## August 5th, 2020 join in with manually curated, haplotig collapsed list of C.virginica IAPs so I'm searching for the correct IAPs
# Use dataframe created in IAP_GIMAP_Gene_Family_Expansion.R (lines 229-233)
nrow(BIR_XP_gff_CG_uniq_XP_XM) # 74
# get row matches in the full mRNA dataframe so I can simply rbind the IAP matches to the rest of the identified IAP genes 
BIR_XP_gff_CG_uniq_XP_XM_apop_product <- C_gig_rtracklayer_filtered[C_gig_rtracklayer_filtered$transcript_id %in% BIR_XP_gff_CG_uniq_XP_XM$transcript_id,] # Join by transcript_id

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and other NOT kDA protein names so I can keep heat shock proteins in both
C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product[
# !grepl("dual specificity protein phosphatase 1B-like", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & not present in list
!grepl("kunitz-type", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("mannosyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("S-acyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
  !grepl("complement C1q", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("cytochrome c oxidase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("cytochrome c-type heme lyase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("cytochrome c1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("activator of 90 kDa", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
!grepl("phosphatidylinositol 3-kinase 1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
!grepl("phosphatidylinositol 3-kinase 2", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
!grepl("DDB_G0272098",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
!grepl("WD repeat-containing protein WRAP73", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
!grepl("tumor protein p63-regulated gene 1-like protein", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
  !grepl("mitogen-activated protein kinase kinase kinase kinase 3", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
  !grepl("mitogen-activated protein kinase kinase kinase kinase 4", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE),]
nrow(C_gig_rtracklayer_apop_product_final) #781
  
C_gig_rtracklayer_apop_product_final <- rbind(C_gig_rtracklayer_apop_product_final, BIR_XP_gff_CG_uniq_XP_XM_apop_product )

# remove duplicated rows
C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product_final[!duplicated(C_gig_rtracklayer_apop_product_final$transcript_id),]
nrow(C_gig_rtracklayer_apop_product_final) # 844 unique transcripts

# Export table as text file
write.table(C_gig_rtracklayer_apop_product_final, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/ PAPER DRAFTS AND OUTLINES/Additional:Supplementary_Data_files_for_paper/Supplementary_Table_4_C_gigas.txt",
            row.names = FALSE)

### Checked genes were added correctly and confirmed with previous results by comparing DF merged list to my previous research
# of which genes are there in my "Updated_APOPTOSIS_GENES, DOMAINS....".xlsx sheet.

# How many unique gene LOC are in this list (this is going in summary paragraph for paper)
C_gig_rtracklayer_apop_product_final_LOC <- C_gig_rtracklayer_apop_product_final[!duplicated(C_gig_rtracklayer_apop_product_final$gene),]
nrow(C_gig_rtracklayer_apop_product_final_LOC ) #511

### Add Gene family Key with gene family names in Apoptosis_names_df ####
  
# Subset apoptosis grepl data frame and Apoptosis_names_df and check nrows
C_gig_rtracklayer_apop_product_final_product <- C_gig_rtracklayer_apop_product_final[,c("Name","product","gene","transcript_id")]
Apoptosis_names_df_CG$rownames <- rownames(Apoptosis_names_df_CG)
nrow(C_gig_rtracklayer_apop_product_final_product) #564
unique(C_gig_rtracklayer_apop_product_final_product$product) #564

# Make new index
idx2 <- sapply(Apoptosis_names_df_CG$product, grep, C_gig_rtracklayer_apop_product_final_product$product)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))

# create df_merged which has duplicated product column for gene name hits, first column is the index name
df_merged_CG <- cbind(Apoptosis_names_df_CG[unlist(idx1),,drop=F], C_gig_rtracklayer_apop_product_final_product[unlist(idx2),,drop=F])
head(df_merged_CG)
# change column names
colnames(df_merged_CG)[1] <- "apoptosis_names_query" 
# subset rownames column
df_merged_CG <- df_merged_CG[,-2]
head(df_merged_CG)
df_merged_CG[is.na(df_merged_CG),]
nrow(df_merged_CG) # 844 (none duplicated, missing two)
View(df_merged_CG)

# Check df_merged_CG apoptosis gene family name assignments here
setdiff(C_gig_rtracklayer_apop_product_final_product$transcript_id,df_merged_CG$transcript_id) # "XM_011452283.2" "XM_011420990.2"
setdiff(df_merged_CG$transcript_id,C_gig_rtracklayer_apop_product_final_product$transcript_id)
df_merged_CG[duplicated(df_merged_CG$transcript_id),]

# check for NAs (meaning duplicate term hits)
df_merged_CG %>% filter(is.na(product)) %>% View() # NO NA's yay! 
df_merged_CG %>% filter(is.na(apoptosis_names_query)) %>% View() # NO NA's yay!!
nrow(df_merged_CG) #844

# merge new index with old
C_gig_rtracklayer_apop_product_final_product_joined <- full_join(C_gig_rtracklayer_apop_product_final_product, df_merged_CG) # Joining, by = c("Name", "product", "gene", "transcript_id")
#View(C_vir_rtracklayer_apop_product_final_product_joined)
nrow(C_gig_rtracklayer_apop_product_final_product_joined) #844 rows same as before

# Recode Gene Family Members to be all included together if they have several different names
C_gig_rtracklayer_apop_product_final_product_joined$apoptosis_names_query <- recode(C_gig_rtracklayer_apop_product_final_product_joined$apoptosis_names_query, 
     # the RHO superfamily to all be included together   
     "cdc42 homolog" = "Rho superfamily", "ras-related"= "Rho superfamily", "ras-like"="Rho superfamily", "rho-related"= "Rho superfamily",
     # IAP repeat containing and putative inhibitor of apoptosis to be together in IAP family
     "inhibitor of apoptosis" = "Inhibitor of apoptosis", "IAP" = "Inhibitor of apoptosis", "uncharacterized" = "Inhibitor of apoptosis",
    # recode heat shock proteins due to capitalization
      "Heat shock protein"="heat shock protein",
    # fas associated death domain
     "FAS-associated death domain protein" = "fas-associated death domain protein",
    # Nf-kappa B family includes transcription factor p65 (aka RELA) and RELB proto-oncogene 
    "transcription factor p65"="NF-kappa-B",
      "RELB proto-oncogene "="NF-kappa-B",
    "BAX"="bcl-2",
    "CD40 ligand"="tumor necrosis factor")

# NOTE : RAC-alpha and RAC-gamma are not part of the Rho superfamily, do not recode
  
head(C_gig_rtracklayer_apop_product_final_product_joined)

### How many transcripts in each gene?
C_gig_rtracklayer_apop_product_final_product_joined_by_gene <- C_gig_rtracklayer_apop_product_final_product_joined %>%
  group_by(gene) %>% summarize(number_transcripts_per_gene = n() )
nrow(C_gig_rtracklayer_apop_product_final_product_joined_by_gene) #511
# Match with names of gene family members
class(C_gig_rtracklayer_apop_product_final_product_joined_by_gene$gene)
class(C_gig_rtracklayer$gene)
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_product <- left_join(C_gig_rtracklayer_apop_product_final_product_joined_by_gene, C_gig_rtracklayer[,c("gene","product")], by = "gene") %>% filter(!is.na(product))
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_product <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_product[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_product$gene),]

View(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_product)

### How many transcripts in each overall gene family (combined across all genes in the family)?
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_family <- C_gig_rtracklayer_apop_product_final_product_joined %>%
  group_by(apoptosis_names_query) %>% summarize(number_transcripts_per_family = n() )

## How many genes in gene family 
#Match that gene ID with the gene name, (just name without the transcript variant)
# separate "transcript variant" by final comma
C_gig_rtracklayer_apop_product_final_product_joined_split_product <- separate(C_gig_rtracklayer_apop_product_final_product_joined, product, into = c("gene_name", "transcript_name"), sep =", transcript variant")

# subset for unique gene
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique <- unique(C_gig_rtracklayer_apop_product_final_product_joined_split_product[,c("gene","gene_name")])

# join just the gene (LOC) so the product name doesn't create duplicate columns
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name <- left_join(C_gig_rtracklayer_apop_product_final_product_joined_by_gene, C_gig_rtracklayer_apop_product_final_product_joined_split_product[,c("gene","gene_name","apoptosis_names_query")], by = c("gene"))
#remove duplicates
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name$gene),]

## Count number of genes in gene family
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>% summarize(gene_family_members = n())
View(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families)

## Number of genes with same gene name
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(gene_name) %>% summarize(gene_duplicates = n())
View(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates)

# avg transcripts per gene family
C_gig_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>%
  summarize(avg_transcripts_per_family = mean(number_transcripts_per_gene))

## Join genes per gene family and transcripts for gene family into one table
C_gig_genes_transcripts_per_gene_family <- full_join(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families,  C_gig_rtracklayer_apop_product_final_product_joined_by_gene_family)
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene <- full_join(C_gig_genes_transcripts_per_gene_family, C_gig_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family)
View(C_gig_genes_transcripts_per_gene_family_transcripts_per_gene)

# Save output to file
save(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates, C_gig_genes_transcripts_per_gene_family, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_annotation_gene_family_info.RData")

#### Combine C_vir and C_gig tables ######

### Make table of presence and absence of specific genes (not gene family) ###
# unique gene name lists for both species
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]

# rename gene column for joining
colnames(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name)[1] <- "C_vir_gene_LOC"
colnames(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name)[1] <- "C_gig_gene_LOC"

#remove "-like" from  gene name for correct name joining using regex
# this isn't exactly working
# https://datascience.stackexchange.com/questions/8922/removing-strings-after-a-certain-character-in-a-given-text 
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\like$")
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\like$")
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\-$")
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\-$")

C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name),]
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name),]

nrow(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name) # 226
nrow(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name) # 224
View(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name)

# https://stackoverflow.com/questions/50861626/removing-dot-from-the-end-of-string
# full join and places where names don't match will get an NA
combined_gene_name_yes_no_table <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name, C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name,
                                             by = "gene_name")
#remove duplicate gene names
combined_gene_name_yes_no_table_unique <- combined_gene_name_yes_no_table[!duplicated(combined_gene_name_yes_no_table$gene_name),]
write.table(combined_gene_name_yes_no_table_unique, file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/C_vir_C_gig_Apoptosis_Pathway_Annotation_Data/combined_gene_name_yes_no_table_unique.txt")
nrow(combined_gene_name_yes_no_table_unique) #275
#join on the pathway descriptions for the molecules 
# load in table with pathway descriptions for each in Excel file "Gene_name_pathway_key.csv"

# Still a few discrepancies in the gene names, namely when they are named in C. virginica "* homolog"
# "mitogen-activated protein kinase kinase kinase 7-interacting protein 3 homolog" is "TGF-beta-activated kinase 1 and MAP3K7-binding protein 3" - need to add one to the tally of shared

# Load gene pathway key with protein aliases curated in excel
# August 19th, 2020: fixed this table and added more specific pathway designations and the newly identified IAPs
gene_name_pathway_key_merged <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/C_vir_C_gig_Apoptosis_Pathway_Annotation_Data/Gene_name_pathway_key.csv", head=TRUE)

# Genes in C vir and not in C gig (all the C_gig_gene_LOC NA's)
c_vir_not_c_gig <- combined_gene_name_yes_no_table_unique %>% filter(is.na(C_gig_gene_LOC))
nrow(c_vir_not_c_gig) # 51
# Genes in C gig and not in C vir (all the C_vir_gene_LOC NA's)
c_gig_not_c_vir <- combined_gene_name_yes_no_table_unique %>% filter(is.na(C_vir_gene_LOC))
nrow(c_gig_not_c_vir) # 49
# shared in both 
shared_apoptosis_gene_names <- na.omit(combined_gene_name_yes_no_table_unique)
nrow(shared_apoptosis_gene_names) # 175
write.table(shared_apoptosis_gene_names , file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/C_vir_C_gig_Apoptosis_Pathway_Annotation_Data/shared_apoptosis_gene_names.txt")

# number of uniquely names proteins in C_vir for paper
C_vir_all <- combined_gene_name_yes_no_table_unique %>% filter(!is.na(C_vir_gene_LOC)) %>% select(C_vir_gene_LOC, gene_name) 
nrow(C_vir_all) #226
# number of uniquely names proteins in C_gig for paper
C_gig_all <- combined_gene_name_yes_no_table_unique %>% filter(!is.na(C_gig_gene_LOC)) %>% select(C_gig_gene_LOC, gene_name) 
nrow(C_gig_all) #224 

### Create Table of C. gig and C. vir gene names, keeping the -"like" to facilitate joining pathway designations in WGCNA ####
# Export table with the "-like" not removed (use this to present in publication)
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org <- C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org <- C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]

# rename gene column for joining
colnames(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org)[1] <- "C_vir_gene_LOC"
colnames(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org)[1] <- "C_gig_gene_LOC"

combined_gene_name_org_yes_no_table <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org, C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name_org,
                                             by = "gene_name")
#remove duplicate gene names
combined_gene_name_org_yes_no_table_unique <- combined_gene_name_org_yes_no_table [!duplicated(combined_gene_name_org_yes_no_table $gene_name),]
colnames(combined_gene_name_org_yes_no_table_unique)[2] <- "gene_name_original"

# Put the product with like removed in another column
combined_gene_name_org_yes_no_table_unique$gene_name <-  stringr::str_remove(combined_gene_name_org_yes_no_table_unique$gene_name_original, "\\like$")
combined_gene_name_org_yes_no_table_unique$gene_name <-  stringr::str_remove(combined_gene_name_org_yes_no_table_unique$gene_name, "\\-$")

# full_join with pathway table without the -like products so that pathway can be added
combined_gene_name_org_yes_no_table_unique_pathway_joined <- full_join(combined_gene_name_org_yes_no_table_unique, gene_name_pathway_key_merged[,c("gene_name","Pathway","Sub_pathway")], by = "gene_name") %>% arrange(desc(gene_name)) 
  # don't join by product name so that the manually edited uncharacterized IAP names would stay

write.table(combined_gene_name_org_yes_no_table_unique_pathway_joined[-4], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/C_vir_C_gig_Apoptosis_Pathway_Annotation_Data/combined_gene_name_org_yes_no_table_unique_pathway_joined.txt")
write.csv(combined_gene_name_org_yes_no_table_unique_pathway_joined[-4], file="/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/C_vir_C_gig_Apoptosis_Pathway_Annotation_Data/combined_gene_name_org_yes_no_table_unique_pathway_joined.csv")

nrow(combined_gene_name_org_yes_no_table_unique_pathway_joined) #405
#join on the pathway descriptions joined
## Use this version of the table, with the original product names to join in WGCNA 

### 2. Make combined table with gene family statistics
# investigate differences in apoptosis_names_df and apoptosis_names_df_CG 
setdiff(Apoptosis_names_df$product, Apoptosis_names_df_CG$product) #"protein BTG1" "dronc" 
setdiff(Apoptosis_names_df_CG$product, Apoptosis_names_df$product) # [1] "^tumor necrosis factor$"             "B-cell translocation gene 1"         "fas-associated death domain protein"
#[4] "Heat shock protein"                  "Dronc"    

C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_table <- C_vir_genes_transcripts_per_gene_family_transcripts_per_gene
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_table <- C_gig_genes_transcripts_per_gene_family_transcripts_per_gene
colnames(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_table)[2:4] <- c("C_vir_gene_family_members", "C_vir_number_transcripts_per_family", 
                                                                                 "C_vir_avg_transcripts_per_family")
colnames(C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_table)[2:4] <- c("C_gig_gene_family_members", "C_gig_number_transcripts_per_family", 
                                                                                 "C_gig_avg_transcripts_per_family")
# join the tables
gene_family_statistics_joined <- full_join(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_table, C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_table)

# Reorder columns so they are next to each other
gene_family_statistics_joined <- gene_family_statistics_joined[,c(1,2,5,3,6,4,7)]


### Make plot of gene family statistics ####

# color palette options
library(ggsci) # for publication colors
library(wesanderson)
names(wes_palettes)
Darjeeling1_palette <- wes_palette("Darjeeling1", 5, type = "discrete")
# helpful R colors website for publication graphics:https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
# scale_color_npg() scale_fill_npg() nature publishing group colors, pal_npg() to list pallette
# wes anderson color pallettes 
# Steve Zissou color pallette https://www.color-hex.com/color-palette/32550


# Add species column first
C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot  <- C_vir_genes_transcripts_per_gene_family_transcripts_per_gene
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot  <- C_gig_genes_transcripts_per_gene_family_transcripts_per_gene

C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot$Species <- "Crassostrea virginica"
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot$Species <- "Crassostrea gigas"

# find non shared families, format ab$a[!(ab$a %in% ab$b)]
C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query[!(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query %in% C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query)]
# "anti-apoptotic protein NR13" "lymphotoxin-alpha"     
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query[!(C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query %in% C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot$apoptosis_names_query)]
# cellular tumor antigen p53,    diablo homolog, mitochondrial, high mobility group box 1 

# created joined table
gene_family_statistics_joined_plot <- full_join(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot, C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot)

# In order to graph shift table into long format rather than wide
gene_family_statistics_joined_plot_long <- gather(gene_family_statistics_joined_plot, key = "gene_family_stat", value = "Count", c(-apoptosis_names_query,-Species))

# Subset plot for those that are actually real gene families of interest (some that I'm calling gene famlies right now are just gene names)
gene_family_statistics_joined_plot_long_gene_family_members <- gene_family_statistics_joined_plot_long %>% filter(gene_family_stat == "gene_family_members")

# remove non-shared gene families
gene_family_statistics_joined_plot_long_gene_family_members_shared <- gene_family_statistics_joined_plot_long_gene_family_members %>% filter(apoptosis_names_query !="anti-apoptotic protein NR13") %>%
                                                                                                                                      filter(apoptosis_names_query !="lymphotoxin-alpha") %>%
                                                                                                                                       filter(apoptosis_names_query != "protein BTG1") %>%
    filter(apoptosis_names_query != "cellular tumor antigen p53") %>% filter(apoptosis_names_query != "diablo homolog, mitochondrial") %>% filter(apoptosis_names_query != "high mobility group box 1")
  
# subset for gene families with less than one gene in them 
gene_family_statistics_joined_final_reduced <- gene_family_statistics_joined_plot_long_gene_family_members_shared %>% filter(Count >= 2)

## bar graph version
ggplot(gene_family_statistics_joined_plot_long_gene_family_members_shared, aes(x=apoptosis_names_query, y =Count, fill = Species)) + 
geom_col(position = "dodge") + theme(axis.text.x = element_text(hjust=0.5,angle = 45))  + bbc_style()
  
## heatmap version
# cut the continues Count variable in discrete bins
gene_family_statistics_joined_plot_long_gene_family_members_shared_binned <- gene_family_statistics_joined_plot_long_gene_family_members_shared %>%
mutate(countfactor = cut(Count, breaks=c(0,5,10,20,40,60,80,100,120,140, max(Count,na.rm=TRUE)),   # create a new variable from count
        labels=c("0-5","5-10","10-20","20-40","40-60","60-80","80-100",
        "100-120","120-140",">140"))) %>%
mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor)))) # change level order

# helpful link for pretty ggplot heatmaps: https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/
# get viridis color pallette of hexcodes
viridis(10, alpha=1)

# plot heatmap of all the gene families
gene_family_cvir_cgig_heatmap <- ggplot(gene_family_statistics_joined_plot_long_gene_family_members_shared_binned, aes(x=Species, y =apoptosis_names_query, fill = countfactor)) +
                    geom_tile(colour="white",size=0.25) +    #add border white colour of line thickness 0.25
                    scale_y_discrete(expand=c(0,0)) + #remove extra space
                    scale_fill_manual(values=c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                           "#B4DE2CFF", "#FDE725FF"),na.value = "grey90")+
                    guides(fill=guide_legend(title="Number of Genes in Gene Family"))+
                    labs(x="",y="",title="Gene Family Members in C. gigas and C. virginica")+
                    theme_grey(base_size=9) 
                    
# subset heatmap for the more interesting gene families I care about

# extra plot formatting code for reference 
#+ xlab("Cell Type") +
#  ylab("Percent Hemocytes") +
#  ggtitle("Percent Caspase 3/7 Active Hemocytes") +
#  theme(panel.background=element_blank(),
#        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
#        text=element_text(family="serif",size=16), 
#        axis.title.y=element_text(family="serif",size=16),
#        axis.title.x=element_text(family="serif",size=16),
#        legend.key=element_rect(fill=NA)) + 
#  theme(text=element_text(size=16))  +
#  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
#  theme(axis.text.x = element_text(size=16)) +
#  theme(legend.text = element_text(size=16)) +
#  scale_x_discrete(labels=c("casp_apop_combined_granular"="Granular", "casp_apop_combined_agranular"="Agranular")) +
#  scale_fill_manual(name="Cell Type", labels=c("Notched Control", "Dermo Injected"), values = c("#e08c67","#6388ca")) 

#### EXPORT DATA FRAMES FOR USE IN DESEQ AND WGCNA ####

save(C_gig_rtracklayer, C_vir_rtracklayer, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
save(C_gig_rtracklayer_apop_product_final, C_vir_rtracklayer_apop_product_final, combined_gene_name_org_yes_no_table_unique_pathway_joined, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")


