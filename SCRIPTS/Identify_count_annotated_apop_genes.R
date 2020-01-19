# R script to identify and count annotated apoptosis genes in the genome

#### Load Packages ####

library(tidyverse)
library(rtracklayer)
library(sqldf)

#### Import Genome Annotation and add GO terms ####

# Import gff file with rtracklayer for C vir
C_vir_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

# Load list of apoptosis names
Apoptosis_names_list <- c('bcl-2-related protein A1',
'apoptosis-inducing factor 1',
'akt1',
'RAC-alpha',
'RAC-gamma serine/threonine-protein kinase-',
'APAF-1 interacting protein',
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
'caspase-14-like',
'cell division cycle and apoptosis regulator protein 1',
'CD151 antigen',
'protein BTG1',
'interferon alpha',
'caspase activity and apoptosis inhibitor 1',
'baculoviral IAP repeat-containing protein',
'cAMP-responsive element modulator',
'cytochrome c-like',
'death-associated inhibitor of apoptosis',
'tumor necrosis factor receptor superfamily member 25',
'tumor necrosis factor receptor superfamily member 10A',
'tumor necrosis factor receptor superfamily member 10B',
'endonuclease G, mitochondrial',
'FAS-associated death domain protein',
'fas apoptotic inhibitory molecule 1',
'tumor necrosis factor receptor superfamily member 6',
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
'mitogen-activated protein kinase kinase kinase 1-like',
'mitogen-activated protein kinase 1',
'dual specificity mitogen-activated protein kinase kinase 1',
'mitogen-activated protein kinase kinase kinase 7-like',
'MYC proto-oncogene',
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
'Early 35 kDa protein',
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
'tumor necrosis factor receptor superfamily member',
'TNFRSF1A associated via death domain',
'TNF receptor-associated factor',
'E3 ubiquitin-protein ligase XIAP',
'netrin receptor',
'neurotrophic receptor tyrosine kinase 1',
'sonic hedgehog receptor',
'peptidyl-prolyl cis-trans isomerase A-like',
'peptidyl-prolyl cis-trans isomerase B-like',
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
'interleukin 17-like protein',
'calreticulin',
'serine protease inhibitor',
"high mobility group box 1",
"pyrin",
'ceramide synthase',
'fatty acid synthase')

# Make this list more generic so it hits the previous ones subset and can be used to count genes in gene families
Apoptosis_names_df <- data.frame(product=c(# removing this because duplicated with bcl-2 'bcl-2-related protein A1',
'apoptosis-inducing factor 1',
'akt1',
'RAC',
'APAF-1 interacting protein',
'tumor necrosis factor',
'cell death regulator Aven',
'bcl-2',
'BAX',
'BAG family molecular chaperone regulator',
'BH3 interacting domain death agonist',
'Bik-like killer protein',
'CASP8 and FADD like apoptosis regulator',
'transcription factor AP-1',
'DNA fragmentation factor subunit beta',
'adenylate cyclase',
'caspase',
'cell division cycle and apoptosis regulator protein 1',
'CD151 antigen',
'protein BTG1',
'interferon alpha',
'IAP',
'cAMP-responsive element modulator',
'cytochrome c',
'endonuclease G, mitochondrial',
'FAS-associated death domain protein',
'fas apoptotic inhibitory molecule 1',
'fas cell surface death receptor',
'GTPase IMAP family member',
'harakiri',
'DNA fragmentation factor subunit alpha',
'interferon-induced protein',
'NF-kappa-B',
'inositol 1,4,5-trisphosphate receptor',
'stress-activated protein kinase JNK',
'induced myeloid leukemia cell differentiation protein Mcl-1',
'mitogen-activated protein kinase',
'MYC proto-oncogene',
'myeloid differentiation primary response protein MyD88',
'phorbol-12-myristate-13-acetate-induced protein 1',
'transcription factor p65',
'RELB proto-oncogene',
'reticuloendotheliosis oncogene',
'anti-apoptotic protein NR13',
'nuclear mitotic apparatus protein 1',
'dynamin-like 120 kDa protein',
'Early 35 kDa protein',
'cellular tumor antigen p53',
'programmed cell death protein',
'p53 and DNA damage-regulated protein 1',
'phosphatidylinositol 3-kinase',
'cAMP-dependent protein kinase',
'protein kinase C delta',
'protein kinase C iota',
'BCL2 binding component 3',
'cdc42 homolog',
'ras-',
'Rho',
'racA',
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
'peptidyl-prolyl cis-trans isomerase',
'receptor-interacting serine/threonine-protein kinase',
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
'interleukin',
'calreticulin',
'serine protease inhibitor',
"high mobility group box 1",
"pyrin",
"ceramide synthase",
'inhibitor of apoptosis',
'fatty acid synthase'))

#### Grep Apoptosis protein names in genome files ####

# C virginica 
C_vir_rtracklayer <- C_vir_rtracklayer %>% filter(type == "mRNA")
C_vir_rtracklayer_apop_product <- C_vir_rtracklayer[grepl(paste(Apoptosis_names_list,collapse="|"), 
                        C_vir_rtracklayer$product, ignore.case = TRUE),]

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and other kDA protein names
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product[!grepl("complement C1q", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
                                 !grepl("activator of 90 kDa", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
                          !grepl("60 kDA", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
                            !grepl("10 kDA", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
  !grepl("dual specificity protein phosphatase 1B-like", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
  !grepl("kunitz-type", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mannosyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("S-acyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE),]

### Checked genes were added correctly and confirmed with previous results by comparing DF merged list to my previous research
# of which genes are there in my "Updated_APOPTOSIS_GENES, DOMAINS....".xlsx sheet., then checked by "Supplementary Table 2. C. virginica apoptosis.xlsx" table that I have
# saved in my Qualifying Exam table to check my numbers for C. virginica

### Add Gene family Key with gene family names in Apoptosis_names_df  ####

# Subset apoptosis grepl data frame and Apoptosis_names_df and check nrows
C_vir_rtracklayer_apop_product_final_product <- C_vir_rtracklayer_apop_product_final[,c("Name","product","gene","transcript_id")]
Apoptosis_names_df$rownames <- rownames(Apoptosis_names_df)
nrow(C_vir_rtracklayer_apop_product_final_product) #1091
unique(C_vir_rtracklayer_apop_product_final_product$product) #467

# Make new index
idx2 <- sapply(Apoptosis_names_df$product, grep, C_vir_rtracklayer_apop_product_final_product$product)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))

# create df_merged which has duplicated product column for gene name hits, first column is the index name
df_merged <- cbind(Apoptosis_names_df[unlist(idx1),,drop=F], C_vir_rtracklayer_apop_product_final_product[unlist(idx2),,drop=F])
head(df_merged)

# change column names
colnames(df_merged)[1] <- "apoptosis_names_query" 
# subset rownames column
df_merged <- df_merged[,-2]
head(df_merged)

# check for NAs (meaning duplicate term hits)
df_merged %>% filter(is.na(product)) %>% View() # NO NA's yay! 
df_merged %>% filter(is.na(apoptosis_names_query)) %>% View() # NO NA's yay!!
nrow(df_merged) #1116 (duplicate rows were added! BAD)

# remove duplicated transcript_id rows, meaning that two search terms hit to it
df_merged <- df_merged[!duplicated(df_merged$transcript_id),]
nrow(df_merged) #1091 duplicates removed

# check that all transcript_ids originally found in search are still there
setdiff(C_vir_rtracklayer_apop_product_final_product$transcript_id, df_merged$transcript_id) # no differences, none misssing

# merge new index with old
C_vir_rtracklayer_apop_product_final_product_joined <- full_join(C_vir_rtracklayer_apop_product_final_product, df_merged) # Joining, by = c("Name", "product", "gene", "transcript_id")
#View(C_vir_rtracklayer_apop_product_final_product_joined)
nrow(C_vir_rtracklayer_apop_product_final_product_joined) #1091 rows PERFECT

# Recode members of the RHO superfamily to all be included together and IAP repeat containing and putative inhibitor of apoptosis to be together in IAP family
C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query <- recode(C_vir_rtracklayer_apop_product_final_product_joined$apoptosis_names_query, 
  "cdc42 homolog" = "Rho superfamily", "racA" = "Rho superfamily", "Rho "= "Rho superfamily", "ras-"="Rho superfamily", "RAC" = "Rho superfamily",
    "inhibitor of apoptosis" = "Inhibitor of apoptosis", "IAP" = "Inhibitor of apoptosis")
head(C_vir_rtracklayer_apop_product_final_product_joined)

### How many transcripts in each gene?
C_vir_rtracklayer_apop_product_final_product_joined_by_gene <- C_vir_rtracklayer_apop_product_final_product_joined %>%
  group_by(gene) %>% summarize(number_transcripts_per_gene = n() )

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
View(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates)

# avg transcripts per gene family
C_vir_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>%
        summarize(avg_transcripts_per_family = mean(number_transcripts_per_gene))

## Join genes per gene family and transcripts for gene family into one table
C_vir_genes_transcripts_per_gene_family <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families,  C_vir_rtracklayer_apop_product_final_product_joined_by_gene_family)
C_vir_genes_transcripts_per_gene_family_transcripts_per_gene <- full_join(C_vir_genes_transcripts_per_gene_family, C_vir_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family)
View(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene)

#### Repeat Analysis using Crassostrea gigas ####

# Apoptosis names in C. gigas genome are lowercase like C.virginica
Apoptosis_names_list_CG <- c('bcl-2-related protein A1',
                          'apoptosis-inducing factor 1',
                          'akt1',
                          'RAC-alpha',
                          'RAC-gamma serine/threonine-protein kinase-',
                          'APAF-1 interacting protein',
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
                          'caspase-14-like',
                          'cell division cycle and apoptosis regulator protein 1',
                          'CD151 antigen',
                          'B-cell translocation gene 1',
                          'interferon alpha',
                          'caspase activity and apoptosis inhibitor 1',
                          'baculoviral IAP repeat-containing protein',
                          'cAMP-responsive element modulator',
                          'cytochrome c',
                          'death-associated inhibitor of apoptosis',
                          'tumor necrosis factor receptor superfamily member 25',
                          'tumor necrosis factor receptor superfamily member 10A',
                          'tumor necrosis factor receptor superfamily member 10B',
                          'endonuclease G, mitochondrial',
                          'FAS-associated death domain protein',
                          'fas apoptotic inhibitory molecule 1',
                          'tumor necrosis factor receptor superfamily member 6',
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
                          'mitogen-activated protein kinase kinase kinase 1,',
                          'mitogen-activated protein kinase 1',
                          'dual specificity mitogen-activated protein kinase kinase 1',
                          'mitogen-activated protein kinase kinase kinase 7',
                          'MYC proto-oncogene',
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
                          'Early 35 kDa protein',
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
                          'tumor necrosis factor receptor superfamily member',
                          'TNFRSF1A associated via death domain',
                          'TNF receptor-associated factor',
                          'E3 ubiquitin-protein ligase XIAP',
                          'netrin receptor',
                          'neurotrophic receptor tyrosine kinase 1',
                          'sonic hedgehog receptor',
                          'peptidyl-prolyl cis-trans isomerase A',
                          'peptidyl-prolyl cis-trans isomerase B',
                          'receptor-interacting serine/threonine-protein kinase',
                          'mixed lineage kinase domain',
                          'heat shock protein',
                          'E3 ubiquitin-protein ligase CHIP',
                          'tumor necrosis factor alpha-induced protein 3',
                          'protein phosphatase 1B',
                          'aurora kinase A',
                          'glutathione peroxidase 4',
                          'gasdermin',
                          'poly \\[ADP-ribose]\\ polymerase 1 ',
                          'poly \\[ADP-ribose]\\ polymerase 1-',
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
                          'interleukin 17-like protein',
                          'calreticulin',
                          'serine protease inhibitor',
                          "high mobility group box 1",
                          "pyrin",
                          'ceramide synthase',
                          'fatty acid synthase')

# Make this list more generic so it hits the previous ones subset and can be used to count genes in gene families
Apoptosis_names_df_CG <- data.frame(product=c(     # removing this because duplicated with bcl-2 'bcl-2-related protein A1',
                                           'apoptosis-inducing factor 1',
                                           'akt1',
                                           'RAC',
                                           'APAF-1 interacting protein',
                                           'tumor necrosis factor',
                                           'cell death regulator Aven',
                                           'bcl-2',
                                           'BAX',
                                           'BAG family molecular chaperone regulator',
                                           'BH3 interacting domain death agonist',
                                           'Bik-like killer protein',
                                           'CASP8 and FADD like apoptosis regulator',
                                           'transcription factor AP-1',
                                           'DNA fragmentation factor subunit beta',
                                           'adenylate cyclase',
                                           'caspase',
                                           'cell division cycle and apoptosis regulator protein 1',
                                           'CD151 antigen',
                                           'B-cell translocation gene 1',
                                           'interferon alpha',
                                           'IAP',
                                           'cAMP-responsive element modulator',
                                           'cytochrome c',
                                           'endonuclease G, mitochondrial',
                                           'FAS-associated death domain protein',
                                           'fas-associated death domain protein', #as with Heat shock protein beta will need to merge these
                                           'fas apoptotic inhibitory molecule 1',
                                           'fas cell surface death receptor',
                                           'GTPase IMAP family member',
                                           'harakiri',
                                           'DNA fragmentation factor subunit alpha',
                                           'interferon-induced protein',
                                           'NF-kappa-B',
                                           'inositol 1,4,5-trisphosphate receptor',
                                           'stress-activated protein kinase JNK',
                                           'induced myeloid leukemia cell differentiation protein Mcl-1',
                                           'mitogen-activated protein kinase',
                                           'MYC proto-oncogene',
                                           'myeloid differentiation primary response protein MyD88',
                                           'phorbol-12-myristate-13-acetate-induced protein 1',
                                           'transcription factor p65',
                                           'RELB proto-oncogene',
                                           'reticuloendotheliosis oncogene',
                                           'anti-apoptotic protein NR13',
                                           'nuclear mitotic apparatus protein 1',
                                           'dynamin-like 120 kDa protein',
                                           'Early 35 kDa protein',
                                           'cellular tumor antigen p53',
                                           'programmed cell death protein',
                                           'p53 and DNA damage-regulated protein 1',
                                           'phosphatidylinositol 3-kinase',
                                           'cAMP-dependent protein kinase',
                                           'protein kinase C delta',
                                           'protein kinase C iota',
                                           'BCL2 binding component 3',
                                           'cdc42 homolog',
                                           'ras-',
                                           'Rho',
                                           'racA',
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
                                           'peptidyl-prolyl cis-trans isomerase',
                                           'receptor-interacting serine/threonine-protein kinase',
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
                                           'interleukin',
                                           'calreticulin',
                                           'serine protease inhibitor',
                                           "high mobility group box 1",
                                           "pyrin",
                                           "ceramide synthase",
                                           'inhibitor of apoptosis',
                                           'fatty acid synthase'))


# Import gff file, using new version of genome annotation
  C_gig_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_GIGAS_PIPELINE/Sequences_Genomes_IDLists/GCF_000297895.1_oyster_v9_genomic.gff")
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
  nrow(C_gig_rtracklayer_apop_product) #797
  
  # Terms to remove
  # remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and other kDA protein names
  C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product[!grepl("complement C1q", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                           !grepl("activator of 90 kDa", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
   !grepl("60 kDA", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
   !grepl("10 kDA", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
     !grepl("75 kDA", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & # for the heat shock proteins
   # !grepl("dual specificity protein phosphatase 1B-like", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & not present in list
   !grepl("kunitz-type", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
   !grepl("mannosyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
   !grepl("S-acyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
     !grepl("cytochrome c oxidase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
   !grepl("cytochrome c-type heme lyase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
     !grepl("cytochrome c1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
     !grepl("interferon alpha-inducible protein 27-like protein 2",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE),]
  nrow(C_gig_rtracklayer_apop_product_final) #686

### Checked genes were added correctly and confirmed with previous results by comparing DF merged list to my previous research
# of which genes are there in my "Updated_APOPTOSIS_GENES, DOMAINS....".xlsx sheet.

### Add Gene family Key with gene family names in Apoptosis_names_df ####
  
# Subset apoptosis grepl data frame and Apoptosis_names_df and check nrows
C_gig_rtracklayer_apop_product_final_product <- C_gig_rtracklayer_apop_product_final[,c("Name","product","gene","transcript_id")]
Apoptosis_names_df_CG$rownames <- rownames(Apoptosis_names_df_CG)
nrow(C_gig_rtracklayer_apop_product_final_product) #686
unique(C_gig_rtracklayer_apop_product_final_product$product) #418

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

# check for NAs (meaning duplicate term hits)
df_merged_CG %>% filter(is.na(product)) %>% View() # NO NA's yay! 
df_merged_CG %>% filter(is.na(apoptosis_names_query)) %>% View() # NO NA's yay!!
nrow(df_merged_CG) #693 (duplicate rows were added! BAD)

# remove duplicated transcript_id rows, meaning that two search terms hit to it
#investigate duplicates
CG_duplicate <- df_merged_CG[duplicated(df_merged_CG$transcript_id),]
df_merged_CG_no_dup <- df_merged_CG[!duplicated(df_merged_CG$transcript_id),]
nrow(df_merged_CG_no_dup) #685, duplicates removed (one line disappeared?)

# check that all transcript_ids originally found in search are still there
setdiff(C_gig_rtracklayer_apop_product_final_product$transcript_id, df_merged_CG_no_dup$transcript_id) # none missing

# merge new index with old
C_gig_rtracklayer_apop_product_final_product_joined <- full_join(C_gig_rtracklayer_apop_product_final_product, df_merged_CG_no_dup) # Joining, by = c("Name", "product", "gene", "transcript_id")
#View(C_vir_rtracklayer_apop_product_final_product_joined)
nrow(C_gig_rtracklayer_apop_product_final_product_joined) #686 rows same as before

# Recode members of the RHO superfamily to all be included together and IAP repeat containing and putative inhibitor of apoptosis to be together in IAP family
C_gig_rtracklayer_apop_product_final_product_joined$apoptosis_names_query <- recode(C_gig_rtracklayer_apop_product_final_product_joined$apoptosis_names_query, 
                                                                                    "cdc42 homolog" = "Rho superfamily", "racA" = "Rho superfamily", "Rho "= "Rho superfamily", "ras-"="Rho superfamily", "RAC" = "Rho superfamily",
                                                                                    "inhibitor of apoptosis" = "Inhibitor of apoptosis", "IAP" = "Inhibitor of apoptosis","Heat shock protein"="heat shock protein")


head(C_gig_rtracklayer_apop_product_final_product_joined)

### How many transcripts in each gene?
C_gig_rtracklayer_apop_product_final_product_joined_by_gene <- C_gig_rtracklayer_apop_product_final_product_joined %>%
  group_by(gene) %>% summarize(number_transcripts_per_gene = n() )
nrow(C_gig_rtracklayer_apop_product_final_product_joined_by_gene) #446

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
#View(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families)

## Number of genes with same gene name
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(gene_name) %>% summarize(gene_duplicates = n())
#View(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates)

# avg transcripts per gene family
C_gig_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family <- C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final %>% group_by(apoptosis_names_query) %>%
  summarize(avg_transcripts_per_family = mean(number_transcripts_per_gene))

## Join genes per gene family and transcripts for gene family into one table
C_gig_genes_transcripts_per_gene_family <- full_join(C_gig_rtracklayer_apop_product_final_product_joined_by_gene_name_final_families,  C_gig_rtracklayer_apop_product_final_product_joined_by_gene_family)
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene <- full_join(C_gig_genes_transcripts_per_gene_family, C_gig_rtracklayer_apop_product_final_product_joined_avg_transcripts_per_gene_family)
View(C_gig_genes_transcripts_per_gene_family_transcripts_per_gene)

#### Combine C_vir and C_gig tables ######

### 1.Make table of presence and absence of specific genes (not gene family)
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

# https://stackoverflow.com/questions/50861626/removing-dot-from-the-end-of-string
# full join and places where names don't match will get an NA
combined_gene_name_yes_no_table <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name, C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name,
                                             by = "gene_name")
#remove duplicate gene names
combined_gene_name_yes_no_table_unique <- combined_gene_name_yes_no_table[!duplicated(combined_gene_name_yes_no_table$gene_name),]

# Still a few discrepancies in the gene names, namely when they are named in C. virginica "* homolog"
  # Need to mention in methods that I rmemoved "like" from the end of names
  # C_ virginica: lipopolysaccharide-induced tumor necrosis factor-alpha factor homolog, C_gig: lipopolysaccharide-induced tumor necrosis factor-alpha
  # C_vir: putative transcription factor p65 homolog, C_gig: transcription factor p65 homolog
  # C_vir: macrophage migration inhibitory factor homolog, C_gig: macrophage migration inhibitory factor

# Genes in C vir and not in C gig (all the C_gig_gene_LOC NA's)
c_vir_not_c_gig <- combined_gene_name_yes_no_table_unique %>% filter(is.na(C_gig_gene_LOC))

### 2. Make combined table with gene family statistics
# investigate differences in apoptosis_names_df and apoptosis_names_df_CG 
setdiff(Apoptosis_names_df$product, Apoptosis_names_df_CG$product)

# Only differences are 
  #[1]  "protein BTG1" (because in C gigas its called B cell translocation gene 1, yes these are the same, see Liu et al., 2017) 
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

### 3.  Make plot of gene family statistics
# Add species column first
C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot  <- C_vir_genes_transcripts_per_gene_family_transcripts_per_gene
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot  <- C_gig_genes_transcripts_per_gene_family_transcripts_per_gene

C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot$Species <- "Crassostrea virginica"
C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot$Species <- "Crassostrea gigas"

gene_family_statistics_joined_plot <- full_join(C_vir_genes_transcripts_per_gene_family_transcripts_per_gene_plot, C_gig_genes_transcripts_per_gene_family_transcripts_per_gene_plot)

# In order to graph shift table into long format rather than wide
gene_family_statistics_joined_plot_long <- gather(gene_family_statistics_joined_plot, key = "gene_family_stat", value = "Count", c(-apoptosis_names_query,-Species))

# Subset plot for those that are actually real gene families of interest (some that I'm calling gene famlies right now are just gene names)
gene_family_statistics_joined_plot_long_gene_family_members <- gene_family_statistics_joined_plot_long %>% filter(gene_family_stat == "gene_family_members")

# subset into two plots for those that are 0 -10 and then 10-100 
gene_family_statistics_joined_plot_long_gene_family_members_0_10 <- 
  
  
ggplot(gene_family_statistics_joined_plot_long_gene_family_members, aes(x=apoptosis_names_query, y =Count, fill = Species)) + 
  geom_col(position = "dodge") + theme(axis.text.x = element_text(hjust=0.5,angle = 45)) + coord_flip()



+ xlab("Cell Type") +
  ylab("Percent Hemocytes") +
  ggtitle("Percent Caspase 3/7 Active Hemocytes") +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),panel.border=element_rect(fill=NA), 
        text=element_text(family="serif",size=16), 
        axis.title.y=element_text(family="serif",size=16),
        axis.title.x=element_text(family="serif",size=16),
        legend.key=element_rect(fill=NA)) + 
  theme(text=element_text(size=16))  +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100))+
  theme(axis.text.x = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  scale_x_discrete(labels=c("casp_apop_combined_granular"="Granular", "casp_apop_combined_agranular"="Agranular")) +
  scale_fill_manual(name="Cell Type", labels=c("Notched Control", "Dermo Injected"), values = c("#e08c67","#6388ca")) 

