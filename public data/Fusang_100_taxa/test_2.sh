#!/bin/bash
#SBATCH -J test_2
#SBATCH -p bioinfo_fat
#SBATCH -n 1 
#SBATCH -c 12
#SBATCH -o test_2.out
#SBATCH -e test_2.err
#SBATCH --mem-per-cpu=2G
module load raxml
raxmlHPC -m GTRGAMMA -p 12345 -s ./simulation/fasta_file/t_sim1_TRUE.fas -n a -T 12

