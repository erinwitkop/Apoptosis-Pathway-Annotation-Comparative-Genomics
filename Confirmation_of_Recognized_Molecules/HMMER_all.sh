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

array=("AIF" "Bag" "Bax" "Bcl-w" "CAAP" "CAD" "Casp3_6_7_9_2" "Casp3_6_7_9_2" "Casp3_6_7_9_2" "Casp3_6_7_9_2" "CCAR" "CD151" "ceramide" "ceramide" "CgBTG1" "c-jun" "CREM" "cyto-c" "DIAP" "DIAP" "endo_G" "FADD" "FADD" "FAIM" "GPCR" "GPCR" "IAP" "ICAD" "IFI44" "IkB" "IP3R" "JNK" "LITAF" "MEKK:MAPK" "MyD88" "NFkB" "p38\ MAPK" "p38\ MAPK" "PDCD" "PDCD" "PDRP" "PI3" "PKA" "PKC" "RhoGTP" "RIP" "sAC" "TLR" "TNFR" "TRAF" "TRAF" "TNF")
array2=("AIF-C_CDD.fasta" "BAG1_CDD.fasta" "BCL2_CDD.fasta" "Bclw_CDD.fasta" "CAAP1_CDD.fasta" "CIDE_CAD_CDD.fasta" "CARD_CASP2_CDD.fasta" "DED_Caspase8_r1.fasta" "d00032_CAS8_CAD.fasta" "CARD_CASP9_CDD.fasta" "DBC1_CCAR.fasta" "CD151_like_LEL.fasta" "TLC_ceramide.fasta" "homeodomain_ceramide.fasta" "btg_CDD.fasta" "JUN_CDD.fasta" "bZIP_CREB1.fasta" "cytochrome_c_CDD.fasta" "BIR_DIAP.fasta" "UBA_IAPs_CDD_DIAP.fasta" "NUC_CDD.fasta" "DEATH_FADD.fasta" "DED_FADD_CDD.fasta" "FAIM1_CDD.fasta" "7tmB2_GPR133.fasta" "GPS_CDD.fasta" "BIR_repeat.fasta" "CIDE_N_ICAD_CDD.fasta" "IFI44.fasta" "ANK_CDD(IkB).fasta" "Ins145_P3_CDD.fasta" "STKc_JNK_CDD.fasta" "LITAF_CDD.fasta" "MEKK1_CDD.fasta" "MyD88_CDD.fasta" "IPT_NKKappaB_CDD.fast" "SIK_p38_CDD.fasta" "p38_MAPK.fasta" "IgV_PD1.fasta" "zf-MYND_PDCD_CDD.fasta" "Prefoldin_alpha.fasta" "PIK3_CDD.fasta" "STKc_PKA_CDD.fasta" "STKc_nPKC_delta_CDD.fasta" "Cdc42_CDD.fasta" "Death_RIP1 _CDD.fasta" "AcyC_CDD.fasta" "TIR_TLR_CDD.fasta" "TNFRSF1A_CDD.fasta" "zf_TRAF_CDD.fasta" "MATH_TRAF6_CDD.fasta" "TNF_CDD.fasta")

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
