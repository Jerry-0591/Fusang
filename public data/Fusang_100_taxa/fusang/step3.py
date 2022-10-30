import os
import subprocess

file_num = len(os.listdir('./tmp_fasta'))

for i in range(0, file_num):
    tmp_cmd = 'python fusang.py --msa_dir ./tmp_fasta/{}.fas --save_prefix {} --beam_size 1'.format(i,i)
    tmp_cmd_2 = 'cp ./dl_output/{}.txt ./subset-{}-outof-5.tre'.format(i, i+1)
    subprocess.call(tmp_cmd,shell=True)
    subprocess.call(tmp_cmd_2,shell=True)

tmp_cmd_3 = ''
tmp_cmd_3 += 'python njmerge.py -t'

for i in range(0, file_num):
    tmp_cmd_3 += ' subset-{}-outof-{}.tre'.format(i+1,file_num)
    
tmp_cmd_3 += ' -m distance.mat -x distance.mat_taxlist -o njmerge.tre'
subprocess.call(tmp_cmd_3,shell=True)    
