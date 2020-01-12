#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_apop_clustalw_tcoffee/tcoffee_default_siglec_heminectin_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_apop_clustalw_tcoffee/tcoffee_default_siglec_heminectin_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_apop_clustalw_tcoffee

#-D submits the start path
echo "START $(date)"

module load tcoffee/11.00.ddc7141-foss-2016b
F=/data3/marine_diseases_lab/erin/CV_apop_clustalw_tcoffee

t_coffee -seq $F/Cvir_heminectin_conserved.fasta -outfile $F/Cvir_heminectin_conserved.out

echo "STOP $(date)"