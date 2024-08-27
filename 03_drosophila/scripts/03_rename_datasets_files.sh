#!/bin/bash

for i in 01_datasets/GC*/{cds,genomic,protein}*; do

        FILENAME=$(echo $i | awk -F "/" '{print $2}') &&
        OUTPUTNAME=$(echo "$(dirname $i)"/"$FILENAME"_"$(basename $i)" | sed -E 's/_from_genomic//') &&

        mv $i $OUTPUTNAME

done
