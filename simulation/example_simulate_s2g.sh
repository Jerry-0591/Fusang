#!/bin/bash

#module load R/3.6.2
#module load gcc

chmod 777 ./indelible
cd ./code
time python simulate_topology.py --num_of_topology 20 --taxa_num 7 --range_of_taxa_num '[7, 40]' \
--len_of_msa_upper_bound 10000 --len_of_msa_lower_bound 1200 --num_of_process 24 --distribution_of_internal_branch_length '[1, 0.5, 0.3]' \
--distribution_of_external_branch_length '[1, 0.5, 0.3]' --range_of_mean_pairwise_divergence '[0.03,0.3]' \
--range_of_indel_substitution_rate '[0.01,0.25]' --max_indel_length 50
cd ../simulate_data
./indelible
cd ../code
time python extract_fasta_data.py
cp ../simulate_data/trees.txt ../label_file/trees.txt
time python gen_numpy.py
rm -r ../simulate_data