#!/bin/bash

for i in 15_selection_analysis/01_input/*fna; do
	
	GENE="$(basename $i | sed -E 's/_aligned.+$//')"
	
	# generate busted files
	BUSTED_OUTNAME="$(echo 15_selection_analysis/01_input/"$GENE"_busted.txt)"
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i |\
		tail -n +2 |\
		sed -E 's/TAA\-\-/-----/; s/TAG\-\-/-----/; s/TGA\-\-/-----/; s/TAG$//; s/TAA$//; s/TGA$//' > "$BUSTED_OUTNAME"
		
	cat "$i".treefile >> "$BUSTED_OUTNAME"

	# generate relax files
	RELAX_OUTNAME="$(echo 15_selection_analysis/01_input/"$GENE"_relax.txt)"
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i |\
		tail -n +2 |\
		sed -E 's/TAA\-\-/-----/; s/TAG\-\-/-----/; s/TGA\-\-/-----/; s/TAG$//; s/TAA$//; s/TGA$//' > "$RELAX_OUTNAME"
	
	gotree labels -i "$i".treefile |\
	 awk '{print $1"\t"$1"{Foreground}"}' > "$i".labels

	gotree rename -i "$i".treefile -m "$i".labels >> "$RELAX_OUTNAME"

done
