#!/bin/bash

mkdir 12_model_selection

#echo "Preparing directory for IQTREE..."
#mkdir tmp_aa

# copy trimmed SRG amino acids into a separate folder
#for i in 09_orthogroup_alignments/*trim.fna; do cp "${i%%.*}".faa tmp_aa/; done
#for i in 11_SRG_alignments/*trim.fna; do cp "${i%%.*}".faa tmp_aa/; done

MODELS="$(tail -n +2 00_input/aa_subst_models.tsv | awk -F "\t" '{print $1}' | tr $'\n' ',' | sed -E 's/,$//')"

echo "Running IQTREE ModelFinder on trimmed amino acids..."
for i in 09_orthogroup_alignments/*trim.fna; do

	FILENAME="$(basename "${i%%.*}")".faa

	iqtree2 -s "${i%%.*}".faa -m TESTONLY --mset "$MODELS" -T 8 --prefix 12_model_selection/"$FILENAME"

done

for i in 11_SRG_alignments/*trim.fna; do 

	FILENAME="$(basename "${i%%.*}")".faa

        iqtree2 -s "${i%%.*}".faa -m TESTONLY --mset "$MODELS" -T 8 --prefix 12_model_selection/"$FILENAME"

done

#rm -rf tmp
