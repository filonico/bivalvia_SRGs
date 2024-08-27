#!/bin/bash

for i in 01_datasets/*zip; do

        ACC=$(basename $i | sed -E 's/_.+.zip//; s/\./_/') &&

        unzip $i -d "${i::-4}" &&
        mv "${i::-4}"/ncbi_dataset/data/"$ACC"/* "${i::-4}" &&
        rm -r "${i::-4}"/ncbi_dataset/ $i
done
