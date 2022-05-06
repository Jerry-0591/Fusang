import os
import shutil
import argparse
from multiprocessing import Process, Pool
import multiprocessing
import subprocess
from Bio import AlignIO

folder = '../simulate_data/'
folder_fasta = '../fasta_file/'

if not os.path.exists(folder_fasta):
    os.mkdir(folder_fasta)

parser = argparse.ArgumentParser('filter output msa length')
p_input = parser.add_argument_group("INPUT")
p_input.add_argument("--length", action="store", type=int, default=1e10, required=False)

args = parser.parse_args()
length = args.length

def get_msa_length(msa_dir):
	alignment = AlignIO.read(open(msa_dir), 'fasta')
	len_of_msa = len(alignment[0].seq)
	return len_of_msa

def extract(ele):
	if('.fas' in ele and 'TRUE' in ele):
		file = folder + ele
		file_fasta = folder_fasta + ele
		shutil.copy(file, file_fasta)

para_list = os.listdir('../simulate_data/')
pool = Pool(8)
pool.map(extract, para_list)
pool.close()
pool.join()

if not os.path.exists('../fasta_file/fail'):
    os.mkdir('../fasta_file/fail')

for i in range(1, len(os.listdir('../fasta_file/'))+1):
	msa_length = get_msa_length('../fasta_file/t_sim{}_TRUE.fas'.format(i))
	if msa_length > length:
		subprocess.call('mv ../fasta_file/t_sim{}_TRUE.fas ../fasta_file/fail/'.format(i), shell=True)

