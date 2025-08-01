#!/bin/bash

for i in 15_selection_analysis/01_input/*busted.txt; do
	
	GENE="$(basename "${i%.*}")"
	PREFIX="$(echo 15_selection_analysis/"$GENE"_leaves)"
		
	hyphy busted --alignment $i --branches Leaves --save-fit "$PREFIX".fit --output "$PREFIX".json | tee -a "$PREFIX".md

done
