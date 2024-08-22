#!/bin/bash

# model selection
mkdir -p 12_model_selection/01_random_trees

MODELS="$(tail -n +2 00_input/aa_subst_models.tsv | awk -F "\t" '{print $1}' | tr $'\n' ',' | sed -E 's/,$//')"

echo "Running IQTREE ModelFinder on trimmed amino acids..."
for i in 09_orthogroup_alignments_withoutSgloAmar/*trim.fna; do

	FILENAME="$(basename "${i%%.*}")".faa

	iqtree2 -s "${i%%.*}".faa -m TESTONLY --mset "$MODELS" -T 15 --prefix 12_model_selection/"$FILENAME"

done

for i in 11_SRG_alignments/*trim.fna; do

	FILENAME="$(basename "${i%%.*}")".faa

	iqtree2 -s "${i%%.*}".faa -m TESTONLY --mset "$MODELS" -T 15 --prefix 12_model_selection/"$FILENAME"

done

# random trees
for i in $(find 11_SRG_alignments/*trim.faa 09_orthogroup_alignments_withoutSgloAmar/*trim.faa -maxdepth 1 -type f | shuf | head -200); do
	LINK="$(realpath $i)"
	
	ln -s "$LINK" 12_model_selection/01_random_trees/

done

iqtree2 -S 12_model_selection/01_random_trees -m TESTNEW --mset $MODELS -T AUTO
