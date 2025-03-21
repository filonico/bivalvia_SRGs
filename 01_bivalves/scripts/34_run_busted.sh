#!/bin/bash

for i in 14_selection_analysis/01_input/*busted.txt; do
	
	GENE="$(basename "${i%.*}")"
	PREFIX="$(echo 14_selection_analysis/"$GENE")"
		
	hyphy busted --alignment 14_selection_analysis/01_input/dmrt_disco_tree3_busted.txt --save-fit "$PREFIX".fit --output "$PREFIX".json | tee -a "$PREFIX".md

done
