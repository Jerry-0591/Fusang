#!/bin/bash
#SBATCH -J test_1
#SBATCH -p bioinfo_fat
#SBATCH -n 1 
#SBATCH -c 12
#SBATCH -o test_1.out
#SBATCH -e test_1.err
#SBATCH --mem-per-cpu=2G
module load iqtree   
iqtree -s ./simulation/fasta_file/t_sim1_TRUE.fas -m GTR+G4+F -b 1 -nt 12
