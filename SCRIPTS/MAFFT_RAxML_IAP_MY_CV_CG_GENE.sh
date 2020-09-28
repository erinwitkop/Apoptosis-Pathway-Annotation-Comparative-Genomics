#!/bin/bash
#SBATCH -t 400:00:00
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_IAP_tree_MY_CV_CG_9_28_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_IAP_tree_error_MY_CV_CG_9_28_2020

# note, to get multithreading to work with mpiexec, remove #SBATCH --export=NONE
echo "START $(date)"

# Set paths needed
M=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT
R=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/RAxML_2020
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/GFFs_All_genomes_comb
F=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/Full_IAP_GIMAP_gene_prot_lists_for_MAFFT


# Load MAFFT first - run the following code locally
module load MAFFT/7.453-GCC-8.3.0-with-extensions

### Generate alignments of all protein sequences using MAFFT ###
echo "Start MAFFT mollusc IAP gene"
cd /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT/
mafft --quiet --thread 20 --auto /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT/BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene.fa > /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT/BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene_MSA.fa
echo "done MAFFT all mollusc IAP"

# --auto selects an appropriate algorithm strategy option based on the size of the data
# don-t use phyllip out because I want the full fasta output. This can be used in RAxML
#  --thread each node on the cluster has 20 threads
# adding --quiet resolves issue where it was throwing errors 

echo "MAFFT done $(date)"

# Unload MAFFT and load RAxML
module purge
module load RAxML/8.2.10-goolf-2016b-hybrid-avx

### Perform ML search and rapid bootstrapping with one command ###
cd $R/
echo "Start RAxML IAP"
mpiexec -npernode 5 raxmlHPC-HYBRID-AVX -T 4 -s $M/BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene_MSA.fa -n BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene_MSA_RAxML -m PROTGAMMAAUTO -x 12345 -p 12345 -f a -N autoMRE
echo "DONE RAxML IAP gene "

# Notes: Option -T does not have any effect with the sequential or parallel MPI version.
# -s is the sequence file name
# -n is the outputFileName
# -m is the substitution model
# -x rapidBootstrapRandomNumberSeed
# -p parsimonyRandomSeed
# -f a is rapid Bootstrap analysis and search for best-scoring ML tree in one program run
# -N specifies the number of alternative runs on distinct starting trees to run.
    # -N: In combination with ­ combined with ­-f a -­x a rapid BS search and thereafter a thorough ML search on the original alignment.
    # -N: if you want to use the bootstrapping criteria rather than a set number, specifiy autoMRE with the -x or -b option
    # autoMRE applies the bookstrap convergence criterion
