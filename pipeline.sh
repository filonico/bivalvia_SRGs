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
