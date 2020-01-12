#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_output_BAD
#SBATCH -e /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_error_BAD
#SBATCH -D /data3/marine_diseases_lab/erin/HMMER_analysis

#-D submits the start path
echo "START $(date)"

module load HMMER/3.1b2-foss-2016b
F=/data3/marine_diseases_lab/erin/HMMER_analysis

#Step 1: build a profile HMM with hmmbuild
#input file as Stockholm or FASTA alignments
#It expects Stockholm by default. To read aligned FASTA files, which HMMER calls “afa” format, 
#specify --informat afa on the command line of any program that reads an input alignment

#Use first line of code if in mfasta format
hmmbuild --informat afa $F/BAD.hmm $F/BAD_conserved_domains.fasta
#hmmbuild $F/siglec.hmm $F/Ig_conserved_Cgsiglec.fasta

#Search sequence database with hmmsearch
#hmmsearch accepts any FASTA file as input. It also accepts EMBL/Uniprot text format. 
#It will automatically determine what format your file is in; you don’t have to say. 

hmmsearch $F/BAD.hmm /data3/marine_diseases_lab/shared/GCF_002022765.2_C_virginica-3.0_protein.faa > BAD_CV_search.out
