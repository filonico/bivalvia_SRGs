#!/bin/bash

for i in 14_selection_analysis/01_input/*fna; do
	
	GENE="$(basename $i | sed -E 's/_aligned.+$//')"
	OUTNAME="$(echo 14_selection_analysis/01_input/"$GENE"_busted.txt)"
	
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i |\
		tail -n +2 | sed -E 's/TAA\-\-/-----/; s/TGA\-\-/-----/; s/TAA$//; s/TGA$//' > "$OUTNAME"
		
	cat "$i".treefile >> "$OUTNAME"

done
