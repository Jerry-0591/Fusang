import os
import shutil
from multiprocessing import Process, Pool
import multiprocessing

folder = '../simulate_data/'
folder_fasta = '../fasta_file/'

if not os.path.exists(folder_fasta):
    os.mkdir(folder_fasta)

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