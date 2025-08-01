#!/bin/bash

for i in 15_selection_analysis/01_input/*relax.txt; do
	
	GENE="$(basename "${i%.*}")"
	PREFIX="$(echo 15_selection_analysis/"$GENE"_leavesTest)"
		
	hyphy relax --alignment $i --test Foreground --save-fit "$PREFIX".fit --output "$PREFIX".json | tee -a "$PREFIX".md

done
