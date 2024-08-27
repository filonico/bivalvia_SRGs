#!/bin/bash

echo Running disco on Dmrt, Sox and Fox gene families...

# create a new directory where to store trees without outgroups
mkdir trees_TMP
cp 06_possvm_orthology/dmrt_ALL_aligned_trim04.faa.treefile.ortholog_groups.newick trees_TMP/
gotree prune -i 06_possvm_orthology/sox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick "Neur_XP.062685912.1.OUT | NA | NA |" "Neur_CCD57795.1.OUT | NA | NA |" > trees_TMP/sox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick
gotree prune -i 06_possvm_orthology/fox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick "Acas_XP.004368148.1.OUT | NA | NA |" "Acas_XP.004333268.1.OUT | NA | NA |" > trees_TMP/fox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick

for i in trees_TMP/*newick; do

	GENE="$(basename $i | sed -E 's/_.+$//')"
	OUTPUT_FILE="$(echo 10_SRG_decomposition/"$GENE"_disco_decomposed.txt)"
	
	# keep trees with at least the 50% of species (9 out of 17)
	python3 compiled_softwares/disco.py -i $i -o "$OUTPUT_FILE" -m 9 -d "_" -v --keep-labels > 10_SRG_decomposition/"$GENE"_disco_decomposed.log 2> /dev/null

	# if the produced file is empty, skip to the next iteration
	if [ ! -s "$OUTPUT_FILE" ]; then
		continue
	fi

	counter_tree=1

	if [ -f 10_SRG_decomposition/"$GENE"_disco_conversion.tsv ]; then
		rm -rf 10_SRG_decomposition/"$GENE"_disco_conversion.tsv
	fi
	
	# write a separate file for each disco sub-tree
	while read tree; do

		DISCOTREE="$GENE"_disco_tree"$counter_tree".nwk
		
		echo $tree | sed -E 's/ [^:]*//g' > 10_SRG_decomposition/"$DISCOTREE"

		# write a file with annotation for every decomposed tree
		ANN=$(echo $tree | gotree labels -i - | sed -E 's/ //g' | awk -F "|" '{print $2"_"$3}' | sort -u | tr $'\n' ',' | sed -E 's/,$//')
		echo -e "$DISCOTREE"$'\t'"$ANN" >> 10_SRG_decomposition/"$GENE"_disco_conversion.tsv

		counter_tree=$((counter_tree+1))
	
	done <"$OUTPUT_FILE"

	echo -e $'\t'"$GENE": done!
	
done

rm -rf trees_TMP

echo DONE
