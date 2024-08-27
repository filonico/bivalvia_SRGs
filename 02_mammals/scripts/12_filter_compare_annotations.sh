#!/bin/bash

# get PANTHER annotation of dmrt, sox and fox genes
# save results with IDs and descriptions to a tsv file
echo "Retrieving Panther annotation of Dmrt, Sox and Fox genes..."

cat 04_family_annotation_filtering/dmrt/*panther* |\
	awk -F "\t" '{print $2"\t"$3}' |\
	sort -u > 04_family_annotation_filtering/dmrt/dmrt_PANTHER_annotation.tsv

cat 04_family_annotation_filtering/sox/*panther* |\
	awk -F "\t" '{print $2"\t"$3}' |\
	sort -u |\
	grep -Ei "sox|sex-determining region y" > 04_family_annotation_filtering/sox/sox_PANTHER_annotation.tsv

cat 04_family_annotation_filtering/fox/*panther* |\
	awk -F "\t" '{print $2"\t"$3}' |\
	sort -u |\
	grep -Ei "fox|forkhead" > 04_family_annotation_filtering/fox/fox_PANTHER_annotation.tsv

# get CDD annotation of dmrt, sox and fox genes
# save results with IDs and descriptions to a tsv file
# REQUIRES: conda_env/ncbi_env.yml
echo "Retrieving CDD annotation of Dmrt, Sox and Fox genes..."

cat 04_family_annotation_filtering/dmrt/*cdd* |\
	awk '{print $2}' |\
	awk -F "|" '{print $NF}' |\
	sort -u |\
	esummary -db cdd | xtract -pattern DocumentSummary -element Id,Subtitle |\
	grep -Ei "double|DM" | grep -Evi "CUE|DMRTA|C1|C2" | sed -E 's/^/gnl|CDD|/' > 04_family_annotation_filtering/dmrt/dmrt_CDD_annotation.tsv

cat 04_family_annotation_filtering/sox/*cdd* |\
	awk '{print $2}' |\
	awk -F "|" '{print $NF}' |\
	sort -u |\
	esummary -db cdd | xtract -pattern DocumentSummary -element Id,Subtitle |\
	grep -Ei "sox|sry" | sed -E 's/^/gnl|CDD|/' > 04_family_annotation_filtering/sox/sox_CDD_annotation.tsv

cat 04_family_annotation_filtering/fox/*cdd* |\
	awk '{print $2}' |\
	awk -F "|" '{print $NF}' |\
	sort -u |\
	esummary -db cdd | xtract -pattern DocumentSummary -element Id,Subtitle |\
	grep -i "fox" | grep -Eiv "fha|coiled|transactivation" | sed -E 's/^/gnl|CDD|/' > 04_family_annotation_filtering/fox/fox_CDD_annotation.tsv

# filter hmm results on the basis of panther annotations
echo "Filtering results on the basis of Panther annotation..."

for gene in {dmrt,sox,fox}; do
	for i in 04_family_annotation_filtering/"$gene"/*panther.tsv; do
		awk -F "\t" '{print $1}' 04_family_annotation_filtering/"$gene"/"$gene"_PANTHER_annotation.tsv |\
			grep -f - $i > "${i::-4}"_filtered.tsv
	done
done

# filter rpsblast results on the basis of cdd annotations
echo "Filtering results on the basis of CDD annotation..."

for gene in {dmrt,sox,fox}; do
	for i in 04_family_annotation_filtering/"$gene"/*cdd.tsv; do
		awk -F "\t" '{print $1}' 04_family_annotation_filtering/"$gene"/"$gene"_CDD_annotation.tsv |\
			grep -f - $i > "${i::-4}"_filtered.tsv
	done
done

# check genes that survived both panther and annotation filtering
echo "Checking genes that have survived both Panther and CDD filtering..."

for gene in {dmrt,fox,sox}; do
	comm -1 -2 \
		<(cat 04_family_annotation_filtering/"$gene"/*panther_filtered* | awk -F "\t" '{print $1}' | sort -u) \
		<(cat 04_family_annotation_filtering/"$gene"/*cdd_filtered* | awk -F "\t" '{print $1}' | sort -u) \
		> 04_family_annotation_filtering/"$gene"/"$gene"_PANTHER_vs_CDD_annotation_survivedBoth.ls
done

# check genes that have not survived either panther or cdd annotation filtering
echo "Checking genes that have not survived either Panther or CDD filtering..."

for gene in {dmrt,fox,sox}; do
	diff --suppress-common-lines -y \
		<(cat 04_family_annotation_filtering/"$gene"/*panther_filtered* | awk -F "\t" '{print $1}' | sort -u) \
		<(cat 04_family_annotation_filtering/"$gene"/*cdd_filtered* | awk -F "\t" '{print $1}' | sort -u) \
		> 04_family_annotation_filtering/"$gene"/"$gene"_PANTHER_vs_CDD_annotation_notSurvived.txt
done
