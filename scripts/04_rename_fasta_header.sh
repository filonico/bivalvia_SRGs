#!/bin/bash

for i in 01_datasets/GC*/*{cds,protein}*; do

        OUTNAME=${i%.*} &&
        EXT=${i##*.} &&
        spID=$(echo $i | awk -F "_" '{print $4}') &&

        # convert multi-line fasta into one-line fasta
        bash scripts/05_onelinerizer.sh $i |\

        # rename header
        sed -E "/^>/ s/^>.+_cds_/>/; s/P_/P./; s/_.+$//; s/ .+$//; s/^>/>"$spID"_/" > "$OUTNAME"_rn."$EXT"

        gzip -9 $i

done
