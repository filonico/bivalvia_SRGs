#!/bin/bash

yes | for i in 01a_assemblies_NCBI/*zip; do

        ACC=$(basename $i | sed -E 's/_.+.zip//; s/\./_/') &&

        unzip $i -d "${i::-4}" &&
        mv "${i::-4}"/ncbi_dataset/data/"$ACC"/* "${i::-4}" &&
        rm -r "${i::-4}"/ncbi_dataset/ $i
done
