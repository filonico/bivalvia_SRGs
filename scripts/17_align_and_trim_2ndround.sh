#!/bin/bash

for i in 05_family_phylogeny/*reduced.faa; do

	GENE=$(basename $i | sed -E 's/_.+$//')	
	FASTA_out="${i/reduced/reduced_aligned}" &&
	HMM=$(ls 00_input/*"$GENE"*hmm) &&

	# align each fasta file using clustal omega
	# REQUIRES: conda_envs/alignments_env.yml
	echo Aligning "$GENE" genes with clustal omega...

	clustalo -i $i --hmm-in="$HMM" --outfmt=fa --threads=10 -o "$FASTA_out" -v -l tmp.log --output-order=tree-order &&

	rm -r tmp.log

	# trim each alignment uing trimal
	# REQUIRES: conda_envs/alignments_env.yml
	echo Trimming the "$GENE" gene alignment with trimAl...

	trimal -in "$FASTA_out" -out "${FASTA_out::-4}"_trim04.faa -gt 0.4	

done
