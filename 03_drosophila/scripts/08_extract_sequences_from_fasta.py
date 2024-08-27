#!/bin/env python3

# Given a list of fasta headers, this script is intented to extract the corresponding sequences from another fasta.
#
# It returns:
#   * a fasta file with extracted sequences according to the provided list;
#   * a list file with sequences that could not be retrieved in the fasta file. 
#
#
# Written by:   Filippo Nicolini
# Last updated: 18/10/2023

import subprocess, argparse, sys, os
from Bio import SeqIO


##########################################
#     Define arguments of the script     #
##########################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Extract sequences from a fasta file given a list of headers.")

# Define some options/arguments/parameters
parser.add_argument("-l", "--list",
                    required = True,
                    help = "The list of fasta headers to be extracted.")

parser.add_argument("-f", "--fasta",
                    required = True,
                    help = "The fasta file from which to extract sequences.")

parser.add_argument("-o", "--output",
                    help = "Name of the output fasta file.",
                    default = "./subsampled_fasta.fasta")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


############################
#     Define functions     #
############################

# Function to extract sequences from fasta according to a list

def extract_sequences(list, fasta):
    sequences_extracted = {}
    sequences_not_found = []

    for sequence in list:
        if sequence in fasta and sequence:
            sequences_extracted[sequence] = fasta[sequence]
        else:
            sequences_not_found.append(sequence)

    return sequences_extracted, sequences_not_found
    

#------------------------------------------------------------------------------------------


############################
#     Read input files     #
############################

# Read in header list
print()
print(f"Reading {args.list}...")

header_list = []
with open(args.list) as input_list:
    [header_list.append(line.strip()) for line in input_list.readlines()]

print(f"    Attempting to extract {len(header_list)} fasta sequences from {args.fasta}.")

# Read in fasta file and remove sequences with duplicated header
seen = {}
input_fasta = []
for sequence in SeqIO.parse(args.fasta, "fasta"):
    if str(sequence.id) not in seen:
        input_fasta.append(sequence)
        seen[(str(sequence.id))] = ""

input_fasta_dict = SeqIO.to_dict(input_fasta)


#################################################
#     Extract selected sequences from fasta     #
#################################################

sequences_extracted_dict,sequences_not_found_list = extract_sequences(header_list, input_fasta_dict)

print(f"    {len(sequences_extracted_dict)} "
      f"({round(len(sequences_extracted_dict)/len(header_list)*100, 2)}%) "
      "sequences were actually extracted.")

try:
    notFound_filename = os.path.splitext(args.output)[0] + "_notFound.ls"
except:
    notFound_filename = args.output + "_notFound.ls"

if not len(sequences_not_found_list) == 0:
    print(f"    {len(sequences_not_found_list)} "
          f"({round(len(sequences_not_found_list)/len(header_list)*100, 2)}%) "
          f"were not present in the fasta file. If you want to check them, see {notFound_filename}.")


##########################################
#     Write output subsampled fasta      #
##########################################

print("Creating output files...")

with open(args.output, 'w') as output_file:
    SeqIO.write(sequences_extracted_dict.values(), output_file, "fasta-2line")

with open(notFound_filename, 'w') as notFound_file:
    [notFound_file.write(f"{line}\n") for line in sequences_not_found_list]
