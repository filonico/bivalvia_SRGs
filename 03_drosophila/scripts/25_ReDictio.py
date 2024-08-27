#!/bin/env python3

"""
### ReDictio is a python script that can _Re_place strings in a file according to a _Dictio_nary. ###

Imagine you have the following fasta file:
>gene1
TTCCACTTGGTTGGGCCGGCTAGGCCTCTCTGCCCGGAGTTTCGGCGCAC
>gene2
TGCTGCCGACAGCCGGGCATTGTTTTAGGGGCGTTATTCGAGGGCACTCG
>gene3
TGCTGCCGACAGCCGGGCATTGTTTTAGGGGCGTTATTCGAGGGCACTCG
>gene4
GAGCTAACTTGTCGGGACCAGCCGGGGTAGTCATCGGGCTTATACAGCGA

and you want to automatically replace headers according to this other file (which we're calling a dictionary):
gene1   sox-7
gene2   dmrt2
gene3   tssk1
gene4   antp

ReDictio can handle this for you!

Written by: Filippo Nicolini
Last update: 08/05/2024

"""

import argparse, sys, re


##########################################
#     Define arguments of the script     #
##########################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "test")

# Define options/arguments/parameters
parser.add_argument("target_file",
                    help = "The file where you want to replace strings according to a dictionary. *DISCLAIMER*: please note that this script currently replaces ONLY EXTACT MACTHES, thus it does not replace substrings inside words.")

parser.add_argument("dictionary",
                    help = "The file that is going to be used as a dictionary. The file has to be a tsv file where the first column is the list of strings you want to replace, and the second the list of the corresponding replacement patterns. REMOVE any header from your dictionary file!")

parser.add_argument("inplace_newFile",
                    choices=["inplace", "new"],
                    help = "If you want to replace strings in the original file, type \"inplace\". If you want to keep the original file, type \"new\".")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


#######################
#     Actual code     #
#######################

# read in the dictionary file and create a dictionary
dictionary = {}
with open(args.dictionary) as dictionary_file:
    for line in dictionary_file.readlines():
        dictionary[line.split("\t")[0].strip()] = line.split("\t")[1].strip()

# read in the target file
with open(args.target_file, 'r') as target_file:
  filedata = target_file.read()

# replace according to dictionary
for key in dictionary:
    replacement = r"\b" + key + r"\b"
    filedata = re.sub(replacement, dictionary[key], filedata)

# replace in the original file if "inplace" is chosen
if args.inplace_newFile == "inplace":
    with open(args.target_file, 'w') as target_file:
        target_file.write(filedata)

# create a new file if "new" is chosen
else:
    outfilename = f"{args.target_file}_replaced"
    with open(outfilename, 'w') as outfile:
        [outfile.write(filedata)]