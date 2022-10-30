#!/bin/bash
#SBATCH -J test_1
#SBATCH -p q_ai
#SBATCH -n 1 
#SBATCH -c 12
#SBATCH -o test_1.out
#SBATCH -e test_1.err
#SBATCH --mem-per-cpu=2G
#SBATCH --gres=gpu:1
source /GPFS/zhangli_lab_permanent/wangzhicheng/phylogenetic_reconstruction/bin/activate
module load cuda/10.1   
time sh z.sh