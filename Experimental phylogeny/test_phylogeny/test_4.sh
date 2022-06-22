#!/bin/bash
#SBATCH -J test_4
#SBATCH -p bioinfo_ai
#SBATCH -n 1 
#SBATCH -c 12
#SBATCH -o test_4.out
#SBATCH -e test_4.err
#SBATCH --mem-per-cpu=2G
#SBATCH --gres=gpu:1
source /GPFS/zhangli_lab_permanent/wangzhicheng/phylogenetic_reconstruction/bin/activate
module load cuda/10.1
python fusang.py --msa_dir a.fas --save_prefix test1 --beam_size 1 --sequence_type standard --window_coverage 1
