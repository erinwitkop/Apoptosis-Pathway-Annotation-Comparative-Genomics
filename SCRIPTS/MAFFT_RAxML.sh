#!/bin/bash
#SBATCH -t 400:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_5_15_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MAFFT_RAxML_error_5_15_2020

echo "START $(date)"

# Set paths needed
IAP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/OrthoFinder_Data_Analysis/Results_Mar25/IAP
GIMAP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/OrthoFinder_Data_Analysis/Results_Mar25/GIMAP
M=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT
R=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/RAxML_2020

# Load MAFFT first
module load MAFFT/7.310-foss-2016b-with-extensions

# Fetch protein sequences using list from Orthogroups in R
# script grab sequences of HMMER output from $O/All_genomes_prot.faa to run through Inteproscan
array1=($(cat $IAP/BIR_IAP_mollusc_orthogroups_XP_list.txt))
for i in ${array1[@]}; do
	sed -n "/${i}/,/^>/p" $O/All_genomes_prot.faa | sed '$d' >> $IAP/BIR_IAP_mollusc_orthogroups_XP_seq.fa
	echo "done"
done

array2=($(cat $GIMAP/AIG1_GIMAP_mollusc_orthogroups_XP_list.txt))
for i in ${array2[@]}; do
	sed -n "/${i}/,/^>/p" $O/All_genomes_prot.faa | sed '$d' >> $GIMAP/AIG1_GIMAP_mollusc_orthogroups_XP_seq.fa
	echo "done"
done

echo "fetch sequences done $(date)"

# Generate alignments of all protein sequences using MAFFT
echo "Start MAFFT IAP"
mafft --auto --phylipout --thread 20 $IAP/ > $M/BIR_IAP_mollusc_orthogroups_MSA.phy
echo "done MAFFT IAP"
echo "Start MAFFT GIMAP"
mafft --auto --phylipout --thread 20 $GIMAP/ > $M/AIG_GIMAP_mollusc_orthogroups_MSA.phy
echo "done MAFFT GIMAP"

# --auto selects an appropriate algorithm strategy option based on the size of the data
# --phylipout will output in phylip format that I can use to make my RAxML tree
#  --thread each node on the cluster has 20 threads

echo "MAFFT done $(date)"

# Unload RAxML and load RAxML
module purge
module load RAxML/8.2.10-goolf-2016b-mpi-avx

# Perform ML search and rapid bootstrapping with one command
mpiexec raxmlHPC-MPI-AVX -T 20 -s $M/BIR_IAP_mollusc_orthogroups_MSA.phy -n $R/BIR_IAP_mollusc_orthogroups_RAxML -m PROTGAMMAAUTO -x 12345 -p 12345 -f a -N autoMRE
mpiexec raxmlHPC-MPI-AVX -T 20 -s $M/AIG_GIMAP_mollusc_orthogroups_MSA.phy  -n $R/AIG_GIMAP_mollusc_orthogroup_RAxML -m PROTGAMMAAUTO -x 12345 -p 12345 -f a -N autoMRE

# -s is the sequence file name
# -n is the outputFileName
# -m is the substitution model
# -x rapidBootstrapRandomNumberSeed
# -p parsimonyRandomSeed
# -T 20 uses 20 threads
# -f a is rapid Bootstrap analysis and search for best-scoring ML tree in one program run
# -N specifies the number of alternative runs on distinct starting trees to run.
    # -N: In combination with ­ combined with ­-f a -­x a rapid BS search and thereafter a thorough ML search on the original alignment.
    # -N: if you want to use the bootstrapping criteria rather than a set number, specifiy autoMRE with the -x or -b option
    # autoMRE applies the bookstrap convergence criterion

echo "RAxML done $(date)"
