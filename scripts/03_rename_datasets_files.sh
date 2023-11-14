#!/bin/bash

for i in 01a_assemblies_NCBI/GC*/{cds,genomic,protein,rna}*; do

        FILENAME=$(echo $i | awk -F "/" '{print $2}') &&
        OUTPUTNAME=$(echo "$(dirname $i)"/"$FILENAME"_"$(basename $i)" | sed -E 's/_from_genomic//') &&

        mv $i $OUTPUTNAME

done
