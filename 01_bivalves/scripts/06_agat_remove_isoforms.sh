#!/bin/bash

for i in 01a_assemblies_NCBI/GC*/*gff; do

        agat_sp_keep_longest_isoform.pl -gff $i -o "${i::-4}"_noIso.gff &&
        gzip -9 $i

done
