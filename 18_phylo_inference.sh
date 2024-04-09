#!/bin/bash

# REQUIRES: conda_envs/phylogeny_env.yml
for i in 05_family_phylogeny/sox*reduced*trim04.faa; do
	
	iqtree2 -s "$i" -m TESTNEW -bb 1000 -nstop 200 -T AUTO --runs 5

done
