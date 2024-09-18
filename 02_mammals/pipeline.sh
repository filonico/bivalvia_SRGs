#!/bin/bash


#########################################
#     Downoad files and rename them     #
#########################################

# retrieve assembly features from NCBI
# REQUIRES: conda_envs/ncbi_env.yml
python3 scripts/01_download_genomes_with_datasets.py -i 00_input/mammalia_genome_toDownload.tsv -spID -d compiled_softwares/datasets -f cds,protein,gff3 -o 01_datasets

# unzip output file and change output directory structure to have just protein, gff, cds and README
bash scripts/02_change_datasets_structure.sh

# rename downloaded genome assembly features
bash scripts/03_rename_datasets_files.sh


################################
#     Rename fasta headers     #
################################

# rename header of each fasta to match the following structure: ">spID_geneID"
bash scripts/04_rename_fasta_header.sh


###########################
#     Remove isoforms     #
###########################

# remove isoforms from gffs using agat
# REQUIRES: conda_envs/gff_manipulation_env.yml
bash scripts/06_agat_remove_isoforms.sh

# remove isoforms from CDS and protein fasta files
# REQUIRES: conda_envs/gff_manipulation_env.yml
bash scripts/07_remove_isoforms_from_fasta.sh

# move reduced fasta files into a dedicated folder
mkdir -p 01_datasets/{01_PROTEOMES,02_CDSs}
mv 01_datasets/GC*/*noIso*fna 01_datasets/02_CDSs/
mv 01_datasets/GC*/*noIso*faa 01_datasets/01_PROTEOMES/

# remove renamed fasta files containig isoforms
rm -rf 01_datasets/GC*/*rn.{fna,faa}

# gzip gff files without isoforms
gzip -9 01_datasets/GC*/*noIso*gff


############################
#     BUSCO assessment     #
############################

# REQUIRES: conda_env/busco_env.yml

mkdir 02_busco

bash scripts/09_run_busco.sh

tar -cvf - 02_busco | gzip -9 - > 02_busco.tar.gz
rm -rf 02_busco/


##########################################################
#     Extract ad analyze dmrt, sox and fox candidates    #
##########################################################

# extract genes from each species starting from the pfam stockholm alignment
# REQUIRES: conda_env/alignments_env.yml
bash scripts/10_get_genes_from_proteomes.sh

# annotate extracted genes with both panther and cdd
# REQUIRES: conda_envs/ncbi_env.yml
# REQUIRES: CDD "little endian" database (originally downloaded from: https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz)
# REQUIRES: panther HMM database (originally downloaded from: http://data.pantherdb.org/ftp/hmm_scoring/18.0/PANTHER18.0_hmmscoring.tgz)
# REQUIRES: edit the script with the full path of your panther- and cdd-databases
bash scripts/11_annotate_genes.sh

# filter genes on the basis of their annotation and compare Panther and CDD results
# NB: only genes surviving both Panther and CDD annotation filtering are kept for subsequent analysis
# REQUIRES: conda_envs/ncbi_env.yml
bash scripts/12_filter_compare_annotations.sh


########################################################
#     Phylogeny-based orthology inference (Possvm)     #
########################################################

mkdir 05_family_phylogeny

# align and trim each gene family
# REQUIRES: conda_env/alignments_env.yml
bash scripts/13_align_and_trim.sh

# infere phylogenetic tree of each gene family
# REQUIRES: conda_envs/phylogeny_env.yml
bash scripts/14_phylo_inference.sh

# infere orthology relationships of dmrt, sox and fox genes from ML trees
# REQUIRES: conda_envs/possvm.yml 
bash scripts/15_possvm_orthology.sh


#######################################################
#     MCL-based orthology inference (OrthoFinder)     #
#######################################################

# prepare the directory of input files to be given to OrthoFinder (i.e., exclude Ggal)
bash scripts/16_prepare_orthofinder_directory.sh

# run OrthoFinder
# REQUIRES: conda_envs/orthofinder_env.yml
bash scripts/17_run_of.sh

##################################################
#     Orthogroup decomposition and alignment     #
##################################################

mkdir 08_orthogroup_decomposition

# run disco on orthogroup gene trees
# REQUIRES: conda_envs/disco_env.yml
bash scripts/18_run_disco.sh

mkdir 09_orthogroup_alignments

# retrieve CDSs fasta file for each decomposed orthogroup
# REQUIRES: conda_envs/disco_env.yml
scripts/19_extract_orthogroup_cdss.sh

# generate a table with genes per each decomposed orthogroup
for i in 09_orthogroup_alignments/*fna; do OG="$(basename ${i%.*})"; GENES="$(grep ">" $i | sed -E 's/^>//' | tr '\n' ',' | sed -E 's/,$//')"; echo -e $OG$'\t'$GENES; done > 09_orthogroup_alignments/decomposed_orthogroups.tsv

# select decomposed orthogroups to keep for alignments (i.e., orthogroups that do not contain Dmrt, Sox and Fox genes)
grep -h ">" 05_family_phylogeny/*ALL.faa | grep -v OUT | sed -E 's/^>//' | grep -v -f - 09_orthogroup_alignments/decomposed_orthogroups.tsv | awk -F "\t" '{print $1}' > 09_orthogroup_alignments/decomposed_orthogroups_tokeep.ls

sed -Ei 's/TAG$//; s/TGA$//; s/TAA$//' 09_orthogroup_alignments/*fna

# align and trim orthogorups
# REQUIRES: conda_envs/alignments_env.yml
bash scripts/20_align_trim_orthogroups.sh


###########################################
#     SRG decomposition and alignment     #
###########################################

mkdir 10_SRG_decomposition

# run disco on SRG trees
# REQUIRES: conda_envs/disco_env.yml
bash scripts/21_run_disco_SRGs.sh

mkdir 11_SRG_alignments

# retrieve CDSs fasta file for each decomposed orthogroup
# REQUIRES: conda_envs/disco_env.yml
bash scripts/22_extract_SRG_cdss.sh

# generate a table with genes per each decomposed orthogroup
for i in 11_SRG_alignments/*fna; do OG="$(basename ${i%.*})"; GENES="$(grep ">" $i | sed -E 's/^>//' | tr '\n' ',' | sed -E 's/,$//')"; echo -e $OG$'\t'$GENES; done > 11_SRG_alignments/decomposed_orthogroups.tsv

sed -Ei 's/TAG$//; s/TGA$//; s/TAA$//' 11_SRG_alignments/*fna

# align and trim orthogorups
# REQUIRES: conda_envs/alignments_env.yml
bash scripts/23_align_trim_SRGs.sh


###################################################
#     Model selection on orthogroups and SRGs     #
###################################################

# run ModelFinder (from IQTREE) on decomposed orthogroups and SRGs
# REQUIRES: conda_envs/phylogeny_env.yml
bash scripts/24_modelSelection.sh


#############################################################
#     Compute the distribution of amino acid divergence     #
#############################################################

mkdir -p 13_distribution_divergence/01_input_alignments

# create a file ith the selected substitution model per each orthogroup
grep Best 12_model_selection/*log | sed -E 's/^.+\///; s/\.log.+: /\t/; s/\+.+$//; s/ .+$//; 1i alignment\tmodel' > 13_distribution_divergence/models_perOrthogroup.tsv

# substitute model names to match names accepted by dist.ml in the subsequent R script
python3 scripts/25_ReDictio.py 13_distribution_divergence/models_perOrthogroup.tsv 00_input/aa_subst_models.tsv new
mv 13_distribution_divergence/models_perOrthogroup.tsv_replaced 13_distribution_divergence/models_perOrthogroup_Rformatted.tsv

# prepare directory with symlinks to alignments
bash scripts/26_prepare_R_directory.sh
sed -Ei '/disco/ s/^/13_distribution_divergence\/01_input_alignments\//' 13_distribution_divergence/models_perOrthogroup_Rformatted.tsv

# plot the distribution of amino acid divergence
# REQUIRES: conda_envs/R_env.yml
Rscript scripts/27_compute_median_distances.R


##############################
#     Selection analysis     #
##############################

mkdir -p 14_selection_analysis/01_input

# infer ML trees on nucleotide alignments of SRGs
# REQUIRES: conda_envs/phylogeny_env.yml
bash scripts/28_compute_MLtrees_SRGs.sh

# prepare busted input
bash scripts/29_prepare_hyphy_files.sh

# run BUSTED
# REQUIRES: conda_envs/selection_env.yml
bash scripts/30_run_busted.sh


########################
#     Plot results     #
########################

# supp fig S4A: compilation
mkdir 06_possvm_orthology/01_plot_occurrence
if [[ -f 06_possvm_orthology/01_plot_occurrence/possvm_orthology_all_withOUT.tsv ]]; then rm -f 06_possvm_orthology/01_plot_occurrence/possvm_orthology_all.tsv; fi; for i in 06_possvm_orthology/*.ortholog_groups.csv; do GENE="$(basename $i | sed -E 's/_.+$//')"; tail -n +2 $i | sed -E "s%OG%"$GENE"_OG%" | awk -F "\t" '{print $1"\t"$2"\t"$4}' | sed -E 's/_[^\t]+//' >> 06_possvm_orthology/01_plot_occurrence/possvm_orthology_all_withOUT.tsv; done
Rscript scripts/31_plot_occurrence.R

# figure 4A: diversity
mkdir 13_distribution_divergence/02_plot_diversity/
Rscript scripts/32_plot_diversity.R

# supp fig S1-3: SRG trees
mkdir 06_possvm_orthology/02_plot_trees
Rscript scripts/33_plot_trees.R


#########################
#     GO enrichment     #
#########################

# create a tsv table with one gene per orthogroup
bash scripts/34_get_genes_GOannoation.sh

# NOTE THAT HERE IT IS NECESSARY TO MANUALLY UPLOAD THE "genes_for_GOannotation.faa" TO THE OMA WEB SERVER, IN ORDER TO OBTIAN THE GO ANNOTATION
# RESULTS OF THE GO ANNOTATION CAN BE FOUND IN intermediate_results/09_GO_enrichment/geneUniverse_GOannotation.tsv.gz

# get the files with the list of genes of interest (genes from the upper quantiles) 
# REQUIRES: conda_envs/possvm_env.yml
bash scripts/35_get_set_geneOfInterests.sh

# perform GO enrichment
Rscript scripts/36_topgGO_enrichment.R 14_GO_enrichment/geneUniverse_GOannotation.tsv 14_GO_enrichment/genesOfInterest1quants_GOannotation.tsv 14_GO_enrichment/genesOfInterest1quants_
Rscript scripts/36_topgGO_enrichment.R 14_GO_enrichment/geneUniverse_GOannotation.tsv 14_GO_enrichment/genesOfInterest2quants_GOannotation.tsv 14_GO_enrichment/genesOfInterest2quants_
