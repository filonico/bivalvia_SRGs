#!/bin/bash


#########################################
#     Downoad files and rename them     #
#########################################

# retrieve assembly features from NCBI
python3 scripts/01_download_genomes_with_datasets.py -i 00_input/bivalve_genomes_toDownload.tsv -spID -d /home/PERSONALE/filippo.nicolini6/softwares/datasets -o 01_datasets

# unzip output file and change output directory structure to have just features from the genome assembly
bash scripts/02_change_datasets_structure.sh

# rename downloaded genome assembly features
bash scripts/03_rename_datasets_files.sh

################################
#     Rename fasta headers     #
################################

# rename header of each fasta to match the following structure: ">spID_geneID"
bash scripts/04_rename_fasta_header.sh

#####################################
#     Remove isoforms from gffs     #
#####################################

# remove isoforms using agat and the information in gff files
bash scripts/06_agat_remove_isoforms.sh
