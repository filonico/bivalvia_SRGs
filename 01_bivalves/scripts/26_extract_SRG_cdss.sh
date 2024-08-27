#!/bin/bash

# generate a unique fasta file with all the nucleotide sequences of input species
cat 01d_FINAL_dataset/02_CDSs/*fna > tmp.fna

for i in 10_SRG_decomposition/*nwk; do
	
	LIST_FILE="${i/nwk/ls}"
	FILE_BASENAME="$(basename ${i%%.*})"
	FNA_OUTNAME="$(echo 11_SRG_alignments//"$FILE_BASENAME".fna)"

	# extract species IDs of bivalves included in orthofinder analysis
	grep yes 00_input/complete_dataset_withTaxonomy.tsv | awk -F "\t" '{print $2}' > tmp.ls

	# extract tips for each decomposed OG
	gotree labels -i $i | grep -f tmp.ls - | grep -Ev "Sglo|Amar" > $LIST_FILE 2> /dev/null

	FINAL_SPECIES="$(wc -l $LIST_FILE)"

#	if ( $FINAL_SPECIES

	# extract fasta files
	python3 scripts/08_extract_sequences_from_fasta.py -l $LIST_FILE -f tmp.fna -o $FNA_OUTNAME

done

rm -rf tmp*
