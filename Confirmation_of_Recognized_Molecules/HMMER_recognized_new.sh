#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_output_recognized
#SBATCH -e /data3/marine_diseases_lab/erin/HMMER_analysis/HMM_error_recognized
#SBATCH -D /data3/marine_diseases_lab/erin/HMMER_analysis

#-D submits the start path
echo "START $(date)"

module load HMMER/3.1b2-foss-2016b
F=/data3/marine_diseases_lab/erin/HMMER_analysis/recognized_molecules

#Step 1: build a profile HMM with hmmbuild
#all are in fasta format so include afa
# the full text file with each conserved domain and the fasta file is called "recognized_CDD_list.txt"
#  cut -f1 recognized_CDD_list.txt | sed 's/^/"/;s/$/"/'

array=("NFkB" "p38\ MAPK" "p38\ MAPK" "PDCD" "PDCD" "PDRP" "PI3" "PKA" "PKC" "RhoGTP" "RIP" "sAC" "TLR" "TNFR" "TRAF" "TRAF" "TNF")
array2=("IPT_NKKappaB_CDD.fasta" "SIK_p38_CDD.fasta" "p38_MAPK.fasta" "IgV_PD1.fasta" "zf-MYND_PDCD_CDD.fasta" "Prefoldin_alpha.fasta" "PIK3_CDD.fasta" "STKc_PKA_CDD.fasta" "STKc_nPKC_delta_CDD.fasta" "Cdc42_CDD.fasta" "Death_RIP1 _CDD.fasta" "AcyC_CDD.fasta" "TIR_TLR_CDD.fasta" "TNFRSF1A_CDD.fasta" "zf_TRAF_CDD.fasta" "MATH_TRAF6_CDD.fasta" "TNF_CDD.fasta")
for ((i=0;i<${#array[@]};++i)); do
  hmmbuild --informat afa $F/${array[i]}/${array[i]}.hmm $F/${array[i]}/${array2[i]}
  echo "done ${array[i]} $(date)"
done


#Search sequence database with hmmsearch
#hmmsearch accepts any FASTA file as input. It also accepts EMBL/Uniprot text format.
#It will automatically determine what format your file is in; you don’t have to say.
for ((i=0;i<${#array[@]};++i)); do
  hmmsearch $F/${array[i]}/${array[i]}.hmm /data3/marine_diseases_lab/shared/GCF_002022765.2_C_virginica-3.0_protein.faa > $F/${array[i]}/${array[i]}.out
  echo "done ${array[i]} $(date)"
done

echo "done $(date)"
