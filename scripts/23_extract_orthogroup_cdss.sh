#!/bin/bash

# generate a unique fasta file with all the nucleotide sequences of input species
cat 01d_FINAL_dataset/02_CDSs/*fna > tmp.fna

for i in 08_orthogroup_decomposition/*nwk; do
	
	LIST_FILE="${i/nwk/ls}"
	FILE_BASENAME="$(basename ${i%%.*})"
	FNA_OUTNAME="$(echo 09_orthogroup_alignments/"$FILE_BASENAME".fna)"

	# extract tips for each decomposed OG
	gotree labels -i $i > $LIST_FILE 2> /dev/null

	# extract fasta files
	python3 scripts/08_extract_sequences_from_fasta.py -l $LIST_FILE -f tmp.fna -o $FNA_OUTNAME

done

rm -rf tmp*
