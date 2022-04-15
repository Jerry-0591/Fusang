#!/bin/bash

#module load R/3.6.2
#module load gcc

cd ./code
time python simulate_topology.py --num_of_topology 200000 --len_of_msa 1000
cd ../simulate_data
./indelible
cd ../code
time python extract_fasta_data.py
cp ../simulate_data/trees.txt ../label_file/trees.txt
time python gen_numpy.py
rm -r ../simulate_data