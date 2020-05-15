#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_seq_5_15_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_seq_error_5_15_2020

echo "START $(date)"

# Set paths needed
H=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/HMMER
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/GFFs_All_genomes_comb

# script grab sequences of HMMER output from $O/All_mollusc_prot.faa to run through Inteproscan
array1=($(cat $H/BIR_hmmsearch_XP_list.txt))
for i in ${array1[@]}; do
	sed -n "/${i}/,/^>/p" $O/All_mollusc_prot.faa | sed '$d' >> $H/BIR_hmmsearch_XP_seq.fa
	echo "done"
done

array2=($(cat $H/AIG1_hmmsearch_XP_list.txt))
for i in ${array2[@]}; do
	sed -n "/${i}/,/^>/p" $O/All_mollusc_prot.faa | sed '$d' >> $H/AIG1_hmmsearch_XP_seq.fa
	echo "done"
done

echo "STOP $(date)"
