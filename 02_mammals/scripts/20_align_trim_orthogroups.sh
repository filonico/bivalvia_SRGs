#!/bin/bash

if [ -f 09_orthogroup_alignments/discarded_orthogroups_after_trimming.ls ]; then
	rm -f 09_orthogroup_alignments/discarded_orthogroups_after_trimming.ls
fi

if [ -f 09_orthogroup_alignments/align_trim_orthogroups.log ]; then
	rm -f 09_orthogroup_alignments/align_trim_orthogroups.log
fi
	

while read j; do
	
	OUT_DIR="09_orthogroup_alignments/"
	ALN_OUT=""$OUT_DIR""$j"_aligned"
	TRIM_OUT=""$ALN_OUT"_trim"
	LOG_OUT=""$TRIM_OUT".log"
	
	echo -e Dealing with "$j"...$"\n"$"\t"Aligning with CLustalW...
	
	# align nucleotide sequences on the basis of amino acid alignemnt with CLustalW
	translatorx -i "$OUT_DIR""$j".fna -o "$ALN_OUT" -p C -t T 2> "$LOG_OUT"

	echo -e $"\t"Trimming with TrimAl...

	# trim amino acid alignemnt
	trimal -in "$ALN_OUT".aa_ali.fasta -out "$TRIM_OUT".faa -gt 0.4 -resoverlap 0.5 -seqoverlap 50 2>> "$LOG_OUT"

	# trim nucleotide alignment
	trimal -in "$ALN_OUT".aa_ali.fasta -backtrans "$OUT_DIR""$j".fna -out "$TRIM_OUT".fna -gt 0.4 -resoverlap 0.5 -seqoverlap 50 2>> "$LOG_OUT"

	# check if after trimming the number of species is still greater than of equal to 16 (half); else, remove the trimmed file and keep track
	NUM_SP="$(grep -c ">" "$TRIM_OUT".fna)"

	echo -e $"\t"$NUM_SP survived after trimming...

	if [ "$NUM_SP" -lt 16 ]; then

		echo "$j" >> "$OUT_DIR"discarded_orthogroups_after_trimming.ls
		
		rm -f "$TRIM_OUT"*

		echo -e $"\t"Orthogroup removed...

	fi

	echo -e DONE$"\n"

done <09_orthogroup_alignments/decomposed_orthogroups_tokeep.ls | tee -a 09_orthogroup_alignments/align_trim_orthogroups.log
