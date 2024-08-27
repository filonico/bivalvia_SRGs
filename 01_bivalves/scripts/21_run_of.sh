#!/bin/bash

orthofinder -f 01d_FINAL_dataset/01_PROTEOMES/01_orthofinder_input/ -t 15 -a 4 -S diamond_ultra_sens --fewer-files -y -o 07_orthofinder/ -n splitted_hogs -s 00_input/species_tree.nwk -X
