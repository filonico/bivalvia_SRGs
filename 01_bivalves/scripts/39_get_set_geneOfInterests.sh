#!/bin/bash

# NOTE THAT THE GO ANNOTATION HAS TO BE PERFORMED MANUALLY ON OMA

# create table for the genes of the 2 upper quantiles
grep OG 13_distribution_divergence/distance_median_values_quant.tsv |\
	grep -E "1$|2$" |\
	sed -E 's/^.+\///; s/\t.+$//; s/_align.+$//' |\
	grep -f - 14_GO_enrichment/genes_for_GOannotation.tsv |\
	awk -F "\t" '{print $2}' |\
	grep -f - 14_GO_enrichment/geneUniverse_GOannotation.tsv > 14_GO_enrichment/genesOfInterest2quants_GOannotation.tsv


# create table for the genes of the 1 upper quantiles
grep OG 13_distribution_divergence/distance_median_values_quant.tsv |\
        grep -E "1$" |\
        sed -E 's/^.+\///; s/\t.+$//; s/_align.+$//' |\
        grep -f - 14_GO_enrichment/genes_for_GOannotation.tsv |\
        awk -F "\t" '{print $2}' |\
        grep -f - 14_GO_enrichment/geneUniverse_GOannotation.tsv > 14_GO_enrichment/genesOfInterest1quants_GOannotation.tsv
