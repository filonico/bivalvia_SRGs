#!/bin/bash


for i in 01_datasets/GC*/*noIso*gff; do
	
	DIR=$(dirname $i) &&
	spID=$(basename $i | awk -F "_" '{print $2}') &&
	OUT_LIST=$(echo "$DIR"/"$spID"_header_noIso.ls) &&
	
	# generate a list of the survived CDSs after agat removal of isoforms
	grep -P '\t'CDS'\t' $i | sed -E "s/^.+cds-/"$spID"_/; s/;Parent.+$//; s/_/./2; s/-.+$//" | sort -u > $OUT_LIST &&
	
	# extract survived sequences from CDS fastas	
	python3 scripts/08_extract_sequences_from_fasta.py -l $OUT_LIST -f "$DIR"/*rn.fna -o "$DIR"/"$spID"_cds_noIso.fna

	# extract survived sequences from protein fastas
        python3 scripts/08_extract_sequences_from_fasta.py -l $OUT_LIST -f "$DIR"/*rn.faa -o "$DIR"/"$spID"_protein_noIso.faa


done
