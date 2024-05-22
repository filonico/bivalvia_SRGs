#!/bin/bash

for i in 11_SRG_alignments/*trim.fna; do

	LINK="$(realpath "$i")"

	ln -s "$LINK" 14_selection_analysis/01_input/"$(basename $i)"

done

for i in 14_selection_analysis/01_input/*fna; do
	
	iqtree2 -s $i -m MFP -T AUTO -bb 1000

done
