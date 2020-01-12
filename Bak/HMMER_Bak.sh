#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_output_Bak
#SBATCH -e /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_error_Bak
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
#hmmbuild --informat afa $F/Akt.hmm $F/Akt_conserved_domains.fasta
hmmbuild $F/Bak.hmm $F/Bak_SMARTBLAST_nocvir.out

#Search sequence database with hmmsearch
#hmmsearch accepts any FASTA file as input. It also accepts EMBL/Uniprot text format. 
#It will automatically determine what format your file is in; you don’t have to say. 

hmmsearch $F/Bak.hmm /data3/marine_diseases_lab/shared/GCF_002022765.2_C_virginica-3.0_protein.faa > Bak_CV_search.out
