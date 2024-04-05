#!/bin/bash

mkdir 06_possvm_orthology


####################
#    ROOT TREES    #
####################

# generate the list of outgroup sequences for sox and fox gene families
echo -e $'\n'Generating lists of outgroup sequences for Sox and Fox gene families and rooting trees...
for i in 05_family_phylogeny/{sox,fox}*ALL.faa; do GENE="$(basename $i | awk -F "_" '{print $1}')"; grep -i OUT $i | sed -E 's/^>//' > 06_possvm_orthology/"$GENE"_outgroups.ls; done

gotree reroot outgroup -i 05_family_phylogeny/sox_ALL_reduced_aligned_trim04.faa.treefile -l 06_possvm_orthology/sox_outgroups.ls -o 05_family_phylogeny/sox_ALL_reduced_aligned_trim04.faa_rooted.treefile
gotree reroot outgroup -i 05_family_phylogeny/fox_ALL_reduced_aligned_trim04.faa.treefile -l 06_possvm_orthology/fox_outgroups.ls -o 05_family_phylogeny/fox_ALL_reduced_aligned_trim04.faa_rooted.treefile


####################
#    RUN POSSVM    #
####################

echo Running possvm...

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/dmrt_ALL_reduced_aligned_trim04.faa.treefile -o 06_possvm_orthology/ -r 00_input/dmrt_reference.tsv -method mclw -itermidroot 10 -skipprint -clean_gene_names > /dev/null 2> 06_possvm_orthology/dmrt.possvm.log

echo -e $'\t'Dmrt genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/sox_ALL_reduced_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 00_input/sox_reference_outgroup.tsv -method mclw -outgroup Neur -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/sox.possvm.log

echo -e $'\t'Sox genes: done!

python3 compiled_softwares/possvm.py -i 05_family_phylogeny/fox_ALL_reduced_aligned_trim04.faa_rooted.treefile -o 06_possvm_orthology/ -r 00_input/fox_reference_outgroup.tsv -method mclw -outgroup Acas -skipprint -skiproot -clean_gene_names > /dev/null 2> 06_possvm_orthology/fox.possvm.log

echo -e $'\t'Fox genes: done!

echo -e Done$'\n'

