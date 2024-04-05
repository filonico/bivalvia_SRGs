#!/bin/bash

TOT_TREES=$(find 07_orthofinder/Results_splitted_hogs/Gene_Trees/ -iname "*txt" | xargs ls -lrt | wc -l)
counter_file=0

echo Running disco on "$TOT_TREES" OrthoFinder gene trees...

for i in 07_orthofinder/Results_splitted_hogs/Resolved_Gene_Trees/*txt; do

	OG="$(basename $i | sed -E 's/_.+$//')"
	OUTPUT_FILE="$(echo 08_orthogroup_decomposition_alignment/"$OG"_disco_decomposed.txt)"

	# keep trees with at least the 50% of species (16 out of 32)
	python3 compiled_softwares/disco.py -i $i -o "$OUTPUT_FILE" -d "_" -m 18 -v --keep-labels > 08_orthogroup_decomposition_alignment/"$OG"_disco_decomposed.log 2> /dev/null

	# if the produced file is empty, skip to the next iteration
	if [ ! -s "$OUTPUT_FILE" ]; then
		continue
	fi

	counter_tree=1

	# write a separate file for each disco sub-tree
	while read tree; do

		echo $tree > 08_orthogroup_decomposition_alignment/"$OG"_disco_tree"$counter_tree".nwk

		counter_tree=$((counter_tree+1))
	
	done <"$OUTPUT_FILE"
	
	counter_file=$((counter_file+1))

	if ! (( "$counter_file" % 1000 )); then
		echo -e $'\t' "$counter_file" analysed
	fi
done

echo DONE
