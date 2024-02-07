#!/bin/bash

for i in 00_input/*stk; do
	
	GENE=$(echo $i | awk -F "." '{print $2}') &&
		
	mkdir -p 03_family_identification/"$GENE" &&
	
	# build hmm profile from pfam stockholm alignment
	hmmbuild "${i::-3}"hmm $i &&
		
	for j in 01d_FINAL_dataset/01_PROTEOMES/*faa; do
		
		SP=$(basename $j | awk -F "_" '{print $1}') &&
		OUTPREFIX=$(echo 03_family_identification/"$GENE"/"$SP"_"$GENE") &&
		
		# for each gene, extract sequences from each species proteome
		hmmsearch -o "$OUTPREFIX".out --tblout "$OUTPREFIX".tblout --domtblout "$OUTPREFIX".domtblout -A "$OUTPREFIX"_dom_aligned.stk \
			--cpu 10 --incE 0.00001 --incdomE 0.00001 "${i::-3}"hmm $j &&
		
		# extract sequences IDs for identified sequences
		grep "#=GS" "$OUTPREFIX"*stk | awk -F " " '{print $NF}' | sort -u > "$OUTPREFIX".ls

		# extract complete aminoacid and nucleotide sequences of identified genes
		touch "$OUTPREFIX".faa
		python3 scripts/08_extract_sequences_from_fasta.py -l "$OUTPREFIX".ls -f $j -o "$OUTPREFIX".faa

		touch "$OUTPREFIX".fna
		python3 scripts/08_extract_sequences_from_fasta.py -l "$OUTPREFIX".ls -f 01d_FINAL_dataset/02_CDSs/"$SP"*fna -o "$OUTPREFIX".fna
	
	done
done
