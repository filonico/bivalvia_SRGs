#!/bin/bash

for i in $(tail -n +2 13_distribution_divergence/models_perOrthogroup_Rformatted.tsv | awk -F "\t" '{print $1}'); do
	
	if [[ $i =~ ^OG.* ]]; then

		LINK="$(realpath 09_orthogroup_alignments_withoutSgloAmar/"$i")"
		
		ln -s "$LINK" 13_distribution_divergence/01_input_alignments/
	
	else
		LINK="$(realpath 11_SRG_alignments/"$i")"
		
		ln -s "$LINK" 13_distribution_divergence/01_input_alignments/
		
	fi
done
