#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_GIMAP_5_22_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_GIMAP_error_5_22_2020

echo "START $(date)"

# Set paths needed
M=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT
R=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/RAxML_2020
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/GFFs_All_genomes_comb
F=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/Full_IAP_GIMAP_gene_prot_lists_for_MAFFT

# Fetch protein sequences using list output from HMMER and Interproscan
# script grab sequences of HMMER output from $O/All_mollusc_prot.faa to run through Inteproscan

array1=($(cat $F/AIG_GIMAP_HMMER_Interpro_XP_list_all.txt))
for i in ${array1[@]}; do
	sed -n "/${i}/,/^>/p" $O/All_mollusc_prot.faa | sed '$d' >> $F/AIG_GIMAP_HMMER_Interpro_XP_list_all.fa
	echo "done"
done

echo "fetch sequences done $(date)"

## Remove duplicates sequences that are identical and keep longest
module load CD-HIT/4.8.1-foss-2018b
cd-hit -G 1 -c 1.0 -t 1 -i $F/AIG_GIMAP_HMMER_Interpro_XP_list_all.fa -o $F/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa
echo "done rm dup $(date)"

# Load MAFFT first
module purge
module load MAFFT/7.453-GCC-8.3.0-with-extensions

### Generate alignments of all protein sequences using MAFFT ###
echo "Start MAFFT all mollusc GIMAP"
cd $M/
mafft --auto --thread 20 $F/AIG_GIMAP_HMMER_Interpro_XP_list_all_rm_dup.fa > $M/AIG_GIMAP_HMMER_Interpro_XP_list_all_MSA.fa
echo "done MAFFT all mollusc GIMAP"

# --auto selects an appropriate algorithm strategy option based on the size of the data
# don-t use phyllip out because I want the full fasta output. This can be used in RAxML
#  --thread each node on the cluster has 20 threads

echo "MAFFT done $(date)"

# Unload MAFFT and load RAxML
module purge
module load RAxML/8.2.10-goolf-2016b-mpi-avx

### Perform ML search and rapid bootstrapping with one command ###
cd $R/
echo "Start RAxML GIMAP"
mpiexec raxmlHPC-MPI-AVX -s $M/AIG_GIMAP_HMMER_Interpro_XP_list_all_MSA.fa  -n AIG_GIMAP_HMMER_Interpro_XP_list_all_MSA_RaxML -m PROTGAMMAAUTO -x 12345 -p 12345 -f a -N autoMRE
echo "DONE RAxML GIMAP"

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
