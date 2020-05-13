#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_seq_5_13_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_seq_error_5_13_2020

echo "START $(date)"

# Set paths needed
H=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/HMMER
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020

# script grab sequences of HMMER output from $O/All_genomes_prot.faa to run through Inteproscan
grep -w -A 1 -Ff $H/BIR_hmmsearch_XP_list.txt $O/All_genomes_prot.faa --no-group-separator > $H/BIR_hmmsearch_XP_seq.fa
grep -w -A 1 -Ff $H/AIG1_hmmsearch_XP_list.txt $O/All_genomes_prot.faa --no-group-separator > $H/AIG1_hmmsearch_XP_seq.fa


echo "STOP $(date)"
