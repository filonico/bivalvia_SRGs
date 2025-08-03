#!/bin/bash

# get gene annotation and busted p-values
join -j 1 -a 1 -a 2 -t $'\t' -e NA -o auto \
    <(cat 10_SRG_decomposition/*_disco_conversion.tsv | sed -E 's/\.nwk//' | LC_COLLATE=C sort -k1 -t $'\t') \
    <(grep p-value 15_selection_analysis/*busted*json | sed -E 's/^.+\///; s/_busted.+:/\t/' | LC_COLLATE=C sort -k1 -t $'\t') |\
    sed -E '1i disco_tree\tannotation\tbusted_pval' \
    > 15_selection_analysis/busted_results.tsv

# get relax p-values and K parameter
join -j 1 -a 1 -a 2 -t $'\t' -e NA -o auto \
    <(grep -E "intensification parameter" 15_selection_analysis/*relax*json | sed -E 's/^.+\///; s/_relax.+:/\t/; s/\,$//' | LC_COLLATE=C sort -k1 -t $'\t') \
    <(grep -E "p-value" 15_selection_analysis/*relax*json | sed -E 's/^.+\///; s/_relax.+:/\t/; s/\,$//' | LC_COLLATE=C sort -k1 -t $'\t') |\
    sed -E '1i disco_tree\trelax_Kparam\trelax_pval' \
    > 15_selection_analysis/relax_results.tsv

# get 1 table
join -j 1 -a 1 -a 2 -t $'\t' -e NA -o auto \
    <(LC_COLLATE=C sort -k1 -t $'\t' 15_selection_analysis/busted_results.tsv) \
    <(LC_COLLATE=C sort -k1 -t $'\t' 15_selection_analysis/relax_results.tsv) | sed -E 's/\t/ /g; s/\./,/g'