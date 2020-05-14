#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Interproscan_out_5_13_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Interproscan_error_5_13_2020

echo "START $(date)"

# Set paths needed
H=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/HMMER
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020

# Load module
module load InterProScan/5.33-72.0-foss-2018b

# Run InterProScan
interproscan.sh -i $H/BIR_hmmsearch_XP_seq.fa -d $H/ -f GFF3,TSV
interproscan.sh -i $H/AIG1_hmmsearch_XP_seq.fa -d $H/ -f GFF3,TSV

# -i is the input data
# -b is the output file base
# -d is the output directory, the output filenames are the same as the input filename
# -f is formats

echo "DONE $(date)"
