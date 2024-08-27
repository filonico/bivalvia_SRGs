#!/bin/bash

mkdir 01d_FINAL_dataset/01_PROTEOMES/01_orthofinder_input/

# create symlink to species to be included in the orthofinder analysis
while read j; do
	
	OF="$(echo $j | awk -F " " '{print $3}')"
	
	if [ "$OF" == "yes" ]; then
		
		SpID="$(echo $j | awk -F " " '{print $2}')"
		LINK="$(realpath 01d_FINAL_dataset/01_PROTEOMES/"$SpID"*)"
		
		ln -t 01d_FINAL_dataset/01_PROTEOMES/01_orthofinder_input/ -s "$LINK"
		
	fi

done <00_input/complete_dataset_withTaxonomy.tsv
