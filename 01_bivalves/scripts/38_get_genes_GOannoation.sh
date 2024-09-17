#!/bin/bash

echo -e "\nGenerating the tsv file with selected genes per orthogroup..."

for i in 09_orthogroup_alignments_withoutSgloAmar/*trim.faa; do
    OG="$(basename "$i" | sed -E 's/_aligned.+$//')"
    for species in Pmax Cgig Hbia Tsqu Sgra; do
        if grep -q "$species" "$i"; then
            echo -e "$OG\t$(grep "$species" "$i" | sed -E 's/>//')"
            break
        fi
    done
done > 14_GO_enrichment/genes_for_GOannotation.tsv

echo -e "\tDone!"

cat 01d_FINAL_dataset/01_PROTEOMES/{Pmax,Cgig,Hbia,Tsqu,Sgra}*faa > proteome.tmp
awk -F "\t" '{print $2}' 14_GO_enrichment/genes_for_GOannotation.tsv > list.tmp

echo -e "\nGenerating the fasta file..."

python3 scripts/08_extract_sequences_from_fasta.py -l list.tmp -f proteome.tmp -o 14_GO_enrichment/genes_for_GOannotation.faa > /dev/null

echo -e "\tDone!"

rm -rf *tmp
