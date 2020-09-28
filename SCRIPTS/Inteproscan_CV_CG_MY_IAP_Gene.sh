#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Gene_Interproscan_out_9_28_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Gene_Interproscan_error_9_28_2020

echo "START $(date)"

# Set paths needed
M=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MAFFT
I=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/Interproscan

# Load module
module load InterProScan/5.44-79.0-foss-2018b

# Run InterProScan
interproscan.sh -t n -i $M/BIR_XP_gff_species_join_haplotig_collapsed_CV_CG_MY_Gene.fa -d $I/ -f GFF3,TSV

# this should translate the input first using EMBOSS getorf tool

# -i is the input data
# -t is the seqtype, labeling here as nucleotide
# -b is the output file base
# -d is the output directory, the output filenames are the same as the input filename
# -f is formats

echo "DONE $(date)"
