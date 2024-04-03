#!/bin/bash


###########################################################
#     Downoad files and rename them (NCBI assemblies)     #
###########################################################

# retrieve assembly features from NCBI
# REQUIRES: conda_envs/ncbi_env.yml
python3 scripts/01_download_genomes_with_datasets.py -i 00_input/bivalve_genomes_toDownload.tsv -spID -d /home/PERSONALE/filippo.nicolini6/SOFTWARES/datasets -o 01a_assemblies_NCBI/

# unzip output file and change output directory structure to have just features from the genome assembly
bash scripts/02_change_datasets_structure.sh

# rename downloaded genome assembly features
bash scripts/03_rename_datasets_files.sh


######################################################################
#     Rename fasta headers and remove isoforms (NCBI assemblies)     #
######################################################################

# rename header of each fasta to match the following structure: ">spID_geneID"
bash scripts/04_rename_fasta_header.sh

# remove isoforms from gffs using agat
# REQUIRES: conda_envs/gff_manipulation_env.yml
bash scripts/06_agat_remove_isoforms.sh

# remove isoforms form CDS and protein fastas
bash scripts/07_remove_isoforms_from_fasta.sh

# move reduced fasta files into a dedicated folder
mkdir -p 01a_assemblies_NCBI/{01_PROTEOMES,02_CDSs}
mv 01a_assemblies_NCBI/GC*/*noIso*fna 01a_assemblies_NCBI/02_CDSs/
mv 01a_assemblies_NCBI/GC*/*noIso*faa 01a_assemblies_NCBI/01_PROTEOMES/

# remove renamed fasta files containig isoforms
rm -rf 01a_assemblies_NCBI/GC*/*rn.{fna,faa}


#################################################################################
#     Rename fasta headers and remove isoforms (other-resources assemblies)     #
#################################################################################

### ATTENTION: ADD THE DOWNLOADED GENOMES IN SOME WAY ###

# extract assembly files
# note that they have been download from different resources
tar -xf 00_input/assemblies_otherResources.tar.gz
mv assemblies_otherResources 01b_assemblies_otherResources

# rename header of each fasta to match the following structure: ">spID_geneID"
# NB: here fasta headers were edited one-by-one because of format inconsistence
bash scripts/09_rename_fasta_header_otherAssemblies.sh

# remove isoforms from gffs using agat
# REQUIRES: conda_envs/gff_manipulation_env.yml
bash scripts/10_agat_remove_isoforms_otherAssemblies.sh

# remove isoforms from CDS and protein fastas
bash scripts/11_remove_isoforms_from_fasta_otherAssemblies.sh

# move reduced fasta files into a dedicated folder
mkdir -p 01b_assemblies_otherResources/{01_PROTEOMES,02_CDSs}
mv 01b_assemblies_otherResources/*/*noIso*fna 01b_assemblies_otherResources/02_CDSs/
mv 01b_assemblies_otherResources/*/*noIso*faa 01b_assemblies_otherResources/01_PROTEOMES/

# remove renamed fasta files containig isoforms
rm -rf 01b_assemblies_otherResources/*/*rn.{fna,faa}


######################################
#     Prpare transcriptome files     #
######################################

# see 01c_transcriptomes/pipeline.sh


################################
#     Run BUSCO on genomes     #
################################

mkdir 02_busco

bash scripts/12_run_busco.sh


##################################################################################
#     Aggregate results, then extract ad analyze dmrt, sox and fox candidates    #
##################################################################################

mkdir -p 01d_FINAL_dataset/{01_PROTEOMES,02_CDSs}

cp 01*/01_PROTEOMES/*faa 01d_FINAL_dataset/01_PROTEOMES/
cp 01c_transcriptomes/01_transdecoder/01_PROTEOMES_final/*faa 01d_FINAL_dataset/01_PROTEOMES/

cp 01*/02_CDSs/*fna 01d_FINAL_dataset/02_CDSs/
cp 01c_transcriptomes/01_transdecoder/02_CDSs_final/*fna 01d_FINAL_dataset/02_CDSs/

# extract genes from each species starting from the pfam stockholm alignment
# REQUIRES: conda_env/alignments_env.yml
bash scripts/13_get_genes_from_proteomes.sh

# annotate extracted genes with both panther and cdd
# REQUIRES: conda_envs/ncbi_env.yml
# REQUIRES: CDD "little endian" database (originally downloaded from: https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz)
# REQUIRES: panther HMM database (originally downloaded from: http://data.pantherdb.org/ftp/hmm_scoring/18.0/PANTHER18.0_hmmscoring.tgz)
# REQUIRES: edit the script with the full path of your panther- and cdd-databases
bash scripts/14_annotate_genes.sh

# filter genes on the basis of their annotation and compare Panther and CDD results
# NB: only genes surviving both Panther and CDD annotation filtering are kept for subsequent analysis
# REQUIRES: conda_envs/ncbi_env.yml
bash scripts/15_filter_compare_annotations.sh


########################################################
#     Phylogeny-based orthology inference (Possvm)     #
########################################################

mkdir 05_family_phylogeny

# align and trim each gene family
# REQUIRES: conda_env/alignments_env.yml
bash scripts/16_align_and_trim.sh

# alignments were then manually inspected to remove sequences with poorly aligning domains
# sequeces were then de-aligned and re-aligned
bash scripts/17_align_and_trim_2ndround.sh

# infere phylogenetic tree of each gene family
# REQUIRES: conda_envs/phylogeny_env.yml
bash scripts/18_phylo_inference.sh

# infere orthology relationships of dmrt, sox and fox genes from ML trees
# REQUIRES: conda_envs/possvm.yml 
bash scripts/19_possvm_orthology.sh


#######################################################
#     MCL-based orthology inference (OrthoFinder)     #
#######################################################

# prepare the directory of input files to be given to OrthoFinder
bash scripts/20_prepare_orthofinder_directory.sh

# run OrthoFinder
# REQUIRES: conda_envs/orthofinder_env.yml
python3 compiled_softwares/OrthoFinder_source/orthofinder_update.py -f 01d_FINAL_dataset/01_PROTEOMES/01_orthofinder_input/ -t 15 -a 4 -S diamond_ultra_sens --fewer-files -y -o 07_orthofinder/ -n splitted_hogs -X
