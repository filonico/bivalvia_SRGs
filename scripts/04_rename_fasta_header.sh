#!/bin/bash

for i in 01a_assemblies_NCBI/GC*/*{cds,protein}*; do

        OUTNAME=${i%.*} &&
        EXT=${i##*.} &&
        spID=$(basename $i | awk -F "_" '{print $2}') &&

        # convert multi-line fasta into one-line fasta
        bash scripts/05_onelinerizer.sh $i |\

        # rename header
        sed -E "s/^>.+_cds_/>/; s/P_/P./; s/_.+$//; s/ .+$//; s/^>/>"$spID"_/" $i > "$OUTNAME"_rn."$EXT"

        gzip -9 $i

done
