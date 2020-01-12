#!/bin/bash

#-D submits the start path
#have to use the command tee for redirection of sed output
array1=($(cat GIMAP_XP_list.txt))

for i in ${array1[@]}; do
	sed -n "/${i}/,/^>/p" ~/Downloads/GCF_002022765.2_C_virginica-3.0_rna.fna | head -n-1 >> GIMAP_XP.fasta 
	echo "done"
done

echo "STOP $(date)"
