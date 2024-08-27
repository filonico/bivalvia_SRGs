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

	sed -E 's/>.+sex.+dmd/dmd/; s/>.+unch.+protein /NA /; s/>.+doublesex- and mab-3-related transcription factor /dmrt-/; s/>.+doublesex-Mab related /dmrt-/; s/>.+doublesex/dsx/; s/>.+abnormal /mab-/; s/ .+$//; s/,$//' > 06_possvm_orthology/dmrt_key_all.tsv

echo -e $'\t'Dmrt genes: done!


# sox
while read j; do

	SPECIES="$(echo $j | sed -E 's/_.+$//')"
	GENE="$(echo $j | sed -E 's/^.+_//')"

	echo $j$'\t'$(zgrep --no-filename $GENE 01_datasets/GCF*"$SPECIES"/*"$SPECIES"*faa.gz)

done <04_family_annotation_filtering/sox/sox_PANTHER_vs_CDD_annotation_survivedBoth.ls |\

	sed -E 's/>.+myb/myb-like-Q /; s/>.+dichaete/dichaete /; s/>.+uncharacterized /NA /; s/>.+sox/sox/i; s/>.+Hmx/hmx /; s/Sox box protein 14/sox-14/i; s/Sox box protein 15/sox-15/i; s/>.+box /sox-/; s/>.+hormone receptor /hormone-receptor-/; s/>.+ Eaf/eaf/; s/>.+mef/mef-/; s/>.+GATA/GATA-zinc-finger-10 /; s/>.+sem/sem/; s/>.+TPR/tpr/; s/>.+alpha-/alpha-protein-kinease-1 /; s/>.+sem/sem/; s/soxNeuro/sox-N/; s/>.+fusilli/fusilli/; s/>.+sterile/female-sterile /; s/>.+kinase/serine\/threonine-kinase /; s/ .+$//; s/,$//; s/-like$//; s/sox21/sox-21/; s/-B/b/' > 06_possvm_orthology/sox_key_all.tsv

echo -e $'\t'Sox genes: done!


# fox
while read j; do

	SPECIES="$(echo $j | sed -E 's/_.+$//')"
	GENE="$(echo $j | sed -E 's/^.+_//')"

	echo $j$'\t'$(zgrep --no-filename $GENE 01_datasets/GCF*"$SPECIES"/*"$SPECIES"*faa.gz)

done <04_family_annotation_filtering/fox/fox_PANTHER_vs_CDD_annotation_survivedBoth.ls |\

	sed -E 's/>.+biniou/biniou /; s/>.+forkhead box protein /fox-/; s/>.+uncharacterized /NA /; s/>.+FD/fd-/; s/>.+crocodile/crocodile /; s/>.+adhesin/adhesin /; s/>.+Mark-A/mark-A/i; s/>.+slp/slp-/; s/>.+homeobox protein /homeobox-/; s/>.+mei4/mei-4/; s/>.+mucin-19/mucin-19/; s/>.+SPT20/spt-20/; s/>.+pakB/pak-B/; s/>.+fork head transcription factor 1/fkh/; s/>.+alpha-protein kinase 1/alpha-protein-kinase-1/; s/>.+forkhead box C1/fox-C1/; s/>.+receptor 4/hormone-receptor-4/; s/>.+HCM1/hcm-1/; s/>.+pfl2/pfl-2/; s/fork-head transcriptional regulator /fox-/; s/>.+hepatocyte nuclear factor 3-beta/hnf-3b/; s/>.+hyphally regulated cell wall protein 3/hyphally-cell-wall-protein-3/; s/>.+myb-like protein P/myb-P/; s/>.+myosin-G heavy chain/myosin-G-heavy-chain/; s/>.+fox-/fox-/; s/>.+forkhead box /fox-/; s/>.+96C/fd-96C/; s/>.+3F/fox-3F/; s/>.+sloppy paired /slp-/; s/>.+jumeau/jumeau/; s/>.+sub-group O/fox-O/; s/>.+102C/fd-102C/; s/>.+59A/fd-59A/; s/>.+19B/fd-19B/; s/>.+serine-rich repeat protein/serine-rich-repeat-protein/; s/>.+female sterile/female-sterile/; s/>.+GATA zinc finger domain-containing protein 7/GATA-zinc-finger-protein7/; s/>.+fork head/fkh /; s/-like//; s/ .+$//; s/,//' > 06_possvm_orthology/fox_key_all.tsv

echo -e $'\t'Fox genes: done!
echo -e Done$'\n'


####################
#    RUN POSSVM    #
####################

echo Running possvm...

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/dmrt_ALL_aligned_trim04.faa.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/dmrt_key_all.tsv -refsps Dhyd,Dpse,Dsuz,Dmel -method mclw -itermidroot 10 -skipprint -clean_gene_names > /dev/null 2> 06_possvm_orthology/dmrt.possvm.log

echo -e $'\t'Dmrt genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/sox_ALL_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/sox_key_all.tsv -refsps Dhyd,Dpse,Dsuz,Dmel -method mclw -outgroup 06_possvm_orthology/sox_outgroups.ls -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/sox.possvm.log

echo -e $'\t'Sox genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/fox_ALL_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 06_possvm_orthology/fox_key_all.tsv -refsps Dhyd,Dpse,Dsuz,Dmel -method mclw -outgroup Acas -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/fox.possvm.log

echo -e $'\t'Fox genes: done!

echo -e Done$'\n'
