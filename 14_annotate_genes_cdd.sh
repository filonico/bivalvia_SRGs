#!/bin/bash

pantherScore="./compiled_softwares/pantherScore2.2/pantherScore2.2_update.pl"

# PLEASE edit with the full path of your panther- and cdd-databases
PANTHERdb="/DATABIG/filipponicolini/DATABASES/PANTHER18.0"
CDDdb="/DATABIG/filipponicolini/DATABASES/CDD/Cdd"

for i in 03_family_identification/*; do

        GENE=$(basename $i) &&
        mkdir -p 04_family_annotation_filtering/"$GENE" &&

        for j in "$i"/*faa; do
                FILENAME=$(basename $j) &&

                # annotate sequences with the panther scoring system
#                $pantherScore -l $PANTHERdb -E 0.00001 -D B -V -i $j -o 04_family_annotation_filtering/"$GENE"/"${FILENAME::-4}"_panther.tsv -c 12 -n

                # annotate sequences with cdd 
                # REQUIRES: conda_env/alignments_env.yml
                rpsblast -query $j -db $CDDdb -max_hsps 1 -evalue 1e-7 -num_threads 8 \
                        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
                        -out 04_family_annotation_filtering/"$GENE"/"${FILENAME::-4}"_cdd.tsv

        done
done
