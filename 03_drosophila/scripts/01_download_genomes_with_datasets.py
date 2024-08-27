#!/bin/env python3

# Given a file with NCBI genome accession numbers, this script is intented to download the corresponding genome features from NCBI.
#
# If species ID are to be included in the name of output directory, please provide a tsv file with:
#   * accession numbers on the first column;
#   * species identifiers on the second column;
#   * additional columns are accepted but will not be used.
#
# If species ID are not to be included, please provide a file with genome assembly accesion number on the first column.
#
# REQUIRED SOFTWARES: datasets
#
# The final output would be a directory structured as follow:
#
# ./
# └── your_output_dir/
#     ├── GCX.XXXXXXXX1.X[_spID1].zip
#     ├── GCX.XXXXXXXX2.X[_spID2].zip
#     ...
#     └── GCX.XXXXXXXXN.X[_spIDN].zip
#
# Written by:   Filippo Nicolini
# Last updated: 12/10/2023


import subprocess, argparse, sys, os


############################################
#     Defining arguments of the script     #
############################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Download genome assembly features from NCBI through the datasets utility")

# Define some options/arguments/parameters
parser.add_argument("-i", "--input_file",
                    required = True,
                    help = "Path to input file. If a tsv, note that accession numbers should be on the first column and speciesID on the second.")

parser.add_argument("-d", "--datasets_path",
                    required = True,
                    help = "Full path to the \"datasets\" executable.")

parser.add_argument("-spID", "--species_ID",
                    action = "store_true",
                    help = "Include the species ID of the species in the name of output directory. (default False)",
                    default = False)

parser.add_argument("-f", "--features",
                    help = "Features of the assembly to be downloaded. Please provided a comma separated string. (defaul: genome,rna,protein,cds,gff3,gtf,gbff,seq-report)",
                    default = "genome,rna,protein,cds,gff3,gtf,gbff,seq-report")

parser.add_argument("-o", "--output_dir",
                    help = "Name of the output directory (default \"00_datasets\")",
                    default = "00_datasets")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


##########################################################################
#     Defining functions to download reads and perform quality check     #
##########################################################################

# Function to retrieve genome assembly feature

def download_assembly_feature(accession,exe_path,features,output):
    try:
        datasets_process = subprocess.run(f"{exe_path} download genome accession {accession} --filename {output} --include {features}",
        shell = True,
        capture_output = True,
        text = True)

    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)

    
##############################
#     Reading input file     #
##############################

# Create output direcotory
if not os.path.isdir(args.output_dir):
    print()
    print(f"Creating output directory in {args.output_dir}/")
    subprocess.run(f"mkdir {args.output_dir}", shell = True)

# Read in genome assembly accession list and store into a list object
print()
print(f"Reading {args.input_file}...")

acc_list = []
with open(args.input_file) as input_acc:
    for line in input_acc.readlines():
        acc_list.append(line.strip().split("\t"))

print(f"    {len(acc_list)} genome assembly accession numbers found.")
print("    " + args.features.replace(",", ", ") + " will be downloaded for each of them.")
print()

# Check if correct arguments have been passed to the script
if args.species_ID and len(acc_list[0]) == 1:
    print("Species ID is to be included in the name of output directory, but only a list has been provided.\n"
          "Please provide a tsv file or do not include species ID. See help message.")
    exit()

#############################################
#     Download genome assembly features     #
#############################################

for assembly in acc_list:
    if len(assembly) == 1:
        output_name = args.output_dir + "/" + assembly[0].replace("_", ".") + ".zip"
    else:
        output_name = args.output_dir + "/" + assembly[0].replace("_", ".") + "_" + assembly[1] + ".zip"

    print(f"-- {assembly[0]} --")
    print("  Retrieving selected features...")
    download_assembly_feature(assembly[0], args.datasets_path, args.features, output_name)
    print("Done")
    print()
