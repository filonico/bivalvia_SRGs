#!/bin/bash

mkdir 06_possvm_orthology

####################
#    ROOT TREES    #
####################

# generate the list of outgroup sequences for sox and fox gene families
echo -e $'\n'Generating lists of outgroup sequences for Sox and Fox gene families and rooting trees...
for i in 05_family_phylogeny/{sox,fox}*ALL.faa; do GENE="$(basename $i | awk -F "_" '{print $1}')"; grep -i OUT $i | sed -E 's/^>//' > 06_possvm_orthology/"$GENE"_outgroups.ls; done

gotree reroot outgroup -i 05_family_phylogeny/sox_ALL_aligned_trim04.faa.treefile -l 06_possvm_orthology/sox_outgroups.ls -o 05_family_phylogeny/sox_ALL_aligned_trim04.faa_rooted.treefile
gotree reroot outgroup -i 05_family_phylogeny/fox_ALL_aligned_trim04.faa.treefile -l 06_possvm_orthology/fox_outgroups.ls -o 05_family_phylogeny/fox_ALL_aligned_trim04.faa_rooted.treefile


####################################################
#    GENERATE INPUT ANNOTATION FILES FOR POSSVM    #
####################################################

## generate the table with gene IDs and names that posvm will use to annotate orthogroups
# dmrt
echo -e $'\n'Generating annotation tables...
while read j; do

	SPECIES="$(echo $j | sed -E 's/_.+$//')"
	GENE="$(echo $j | sed -E 's/^.+_//')"

	echo $j$'\t'$(zgrep --no-filename $GENE 01_datasets/GCF*"$SPECIES"/*"$SPECIES"*faa.gz)

done <04_family_annotation_filtering/dmrt/dmrt_PANTHER_vs_CDD_annotation_survivedBoth.ls |\

	sed -E 's/>.+sex.+factor /dmrt/; s/ PREDICTED://; s/ .+$//; s/-li.+$//; s/,$//' > 06_possvm_orthology/dmrt_key_all.tsv

echo -e $'\t'Dmrt genes: done!


# sox
while read j; do

	SPECIES="$(echo $j | sed -E 's/_.+$//')"
	GENE="$(echo $j | sed -E 's/^.+_//')"

	echo $j$'\t'$(zgrep --no-filename $GENE 01_datasets/GCF*"$SPECIES"/*"$SPECIES"*faa.gz)

done <04_family_annotation_filtering/sox/sox_PANTHER_vs_CDD_annotation_survivedBoth.ls |\

	sed -E 's/>.+SOX/sox/; s/>.+ion Y/sry/; s/>.+CH/CH/; s/-like.+$//; s/ .+$//; s/,$//' > 06_possvm_orthology/sox_key_all.tsv

echo -e $'\t'Sox genes: done!


# fox
while read j; do

	SPECIES="$(echo $j | sed -E 's/_.+$//')"
	GENE="$(echo $j | sed -E 's/^.+_//')"

	echo $j$'\t'$(zgrep --no-filename $GENE 01_datasets/GCF*"$SPECIES"/*"$SPECIES"*faa.gz)

done <04_family_annotation_filtering/fox/fox_PANTHER_vs_CDD_annotation_survivedBoth.ls |\

	sed -E 's/>.+forkhead box protein /fox/; s/>.+hepatocyte nuclear factor /hnf/; s/>.+forkhead box /fox/; s/-alpha/a/; s/-beta/b/; s/-gamma/g/; s/,.+$//; s/-.+$//; s/ .+$//' > 06_possvm_orthology/fox_key_all.tsv

echo -e $'\t'Fox genes: done!
echo -e Done$'\n'


####################
#    RUN POSSVM    #
####################

echo Running possvm...

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/dmrt_ALL_aligned_trim04.faa.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/dmrt_key_all.tsv -refsps Hsap,Mmus,Emax,Oana -method mclw -itermidroot 10 -skipprint -clean_gene_names > /dev/null 2> 06_possvm_orthology/dmrt.possvm.log

echo -e $'\t'Dmrt genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/sox_ALL_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/sox_key_all.tsv -refsps Hsap,Mmus,Emax,Oana -method mclw -outgroup 06_possvm_orthology/sox_outgroups.ls -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/sox.possvm.log

echo -e $'\t'Sox genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/fox_ALL_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/fox_key_all.tsv -refsps Hsap,Mmus,Emax,Oana -method mclw -outgroup Acas -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/fox.possvm.log

echo -e $'\t'Fox genes: done!

echo -e Done$'\n'
