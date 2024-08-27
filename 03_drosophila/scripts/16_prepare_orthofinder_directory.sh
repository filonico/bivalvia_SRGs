#!/bin/bash

mkdir 01_datasets/01_PROTEOMES/01_orthofinder_input/

# create symlink to species to be included in the orthofinder analysis
while read j; do
	
	SpID="$(echo $j | awk -F " " '{print $2}')"
	
	if [ "$SpID" != "Agam" ]; then
		
		LINK="$(realpath 01_datasets/01_PROTEOMES/"$SpID"*)"
		
		ln -t 01_datasets/01_PROTEOMES/01_orthofinder_input/ -s "$LINK"
		
	fi

done <00_input/drosophila_genome_toDownload.tsv
