#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HMMER_out_5_13_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HMMER_error_5_13_2020

echo "START $(date)"

# Set paths needed
H=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/HMMER
O=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020

# Load HMMER
module load HMMER/3.2.1-foss-2018b

#Step 1: build a profile HMM with hmmbuild
#input file as Stockholm or FASTA alignments
#It expects Stockholm by default. To read aligned FASTA files, which HMMER calls “afa” format,
#specify --informat afa on the command line of any program that reads an input alignment

#Use first line of code if in mfasta format
hmmbuild --informat afa $H/BIR.hmm $H/PF00653_full_BIR_alignment.fa
hmmbuild --informat afa $H/AIG1.hmm $H/PF04548_full_AIG1_alignment.fa

#Search sequence database of all the protein sequences with hmmsearch
#hmmsearch accepts any FASTA file as input. It also accepts EMBL/Uniprot text format.
#It will automatically determine what format your file is in; you don’t have to say.

hmmsearch -E 0.00001 --tblout $H/BIR_hmmsearch_tbl.out --domtblout $H/BIR_hmmsearch_domtbl.out $H/BIR.hmm $O/All_genomes_prot.faa
hmmsearch -E 0.00001 --tblout $H/AIG1_hmmsearch_tbl.out --domtblout $H/AIG1_hmmsearch_domtbl.out $H/AIG1.hmm $O/All_genomes_prot.faa

# -E sets the evalue cutoff , selecting E value cutoff of 10^-5 as my cutoff

# Multithreading? HMMER already uses all available processors

echo "DONE $(date)"
