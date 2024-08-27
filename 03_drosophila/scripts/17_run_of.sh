#!/bin/bash

python3 compiled_softwares/OrthoFinder_source/orthofinder_update.py -f 01_datasets/01_PROTEOMES/01_orthofinder_input/ -t 15 -a 4 -S diamond_ultra_sens --fewer-files -y -o 07_orthofinder/ -n splitted_hogs -X
