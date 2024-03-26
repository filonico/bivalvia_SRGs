#!/bin/bash

for gene in {dmrt,sox,fox}; do
	
	FASTA_out=$(echo 05_family_phylogeny/"$gene"_ALL.faa) &&
	HMM=$(ls 00_input/*"$gene"*hmm) &&

	# retrieve the fasta file with all the genes that have survived both panther and CDD annotation
	cat 03_family_identification/"$gene"/*faa > "$gene".faa.tmp &&
		
	python3 scripts/08_extract_sequences_from_fasta.py \
		-l 04_family_annotation_filtering/"$gene"/"$gene"_PANTHER_vs_CDD_annotation_survivedBoth.ls \
		-f "$gene".faa.tmp \
		-o 05_family_phylogeny/"$gene"_ALL.faa &&
		
	rm -f *tmp &&

	# add outgroups to sox and fox fasta files
	cat 00_input/"$gene"*faa >> 05_family_phylogeny/"$gene"_ALL.faa

	# align each fasta file using clustal omega
	# REQUIRES: conda_envs/alignments_env.yml

	echo Aligning "$gene" genes with clustal omega...

	clustalo -i $FASTA_out --hmm-in="$HMM" --outfmt=fa --threads=10 -o "${FASTA_out::-4}"_aligned.faa -v -l tmp.log --output-order=tree-order &&

	rm -r tmp.log

	# trim each alignment uing trimal
	# REQUIRES: conda_envs/alignments_env.yml

	echo Trimming the "$gene" gene alignment with trimAl...

	trimal -in "${FASTA_out::-4}"_aligned.faa -out "${FASTA_out::-4}"_aligned_trim04.faa -gt 0.4
	
#	if [[ $gene == dmrt ]]; then
#		trimal -in "${FASTA_out::-4}"_aligned.faa -out "${FASTA_out::-4}"_aligned_trim04.faa -gt 0.4
#	elif [[ $gene == sox ]]; then
#		trimal -in "${FASTA_out::-4}"_aligned.faa -out "${FASTA_out::-4}"_aligned_trim035.faa -gt 0.35
#	else
#		trimal -in "${FASTA_out::-4}"_aligned.faa -out "${FASTA_out::-4}"_aligned_trim03.faa -gt 0.3
#	fi

done
