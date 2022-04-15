import os
import re
import numpy as np
import random
import pandas as pd
import ete3
from ete3 import Tree

folder_fasta = '../fasta_file/'
folder_numpy = '../numpy_file/'

if not os.path.exists(folder_numpy):
    os.mkdir(folder_numpy)

folder_numpy_seq = folder_numpy + 'seq/'
folder_numpy_label = folder_numpy + 'label/'
if not os.path.exists(folder_numpy_seq):
    os.mkdir(folder_numpy_seq)
if not os.path.exists(folder_numpy_label):
    os.mkdir(folder_numpy_label)

csv_data = pd.read_table("../label_file/trees.txt", skiprows=5,sep='\t',header=None)
file = list(csv_data[0])
topo = list(csv_data[8])
dic = {}

for i in range(0, len(file)):
    dic[file[i]] = topo[i]

def assign_label(str):
    t1 = Tree(str,format=5)
    label_0 = Tree('((0,1),2,3);')
    label_1 = Tree('((0,2),1,3);')
    label_2 = Tree('((0,3),1,2);')
    if t1.robinson_foulds(label_0, unrooted_trees=True)[0] == 0:
        return 0
    elif t1.robinson_foulds(label_1, unrooted_trees=True)[0] == 0:   
        return 1
    elif t1.robinson_foulds(label_2, unrooted_trees=True)[0] == 0:
        return 2

def get_numpy(folder_dir, fasta_dir):
    aln_file = folder_dir + fasta_dir
    aln = open(aln_file)
    dic = {'A':'0','T':'1','C':'2','G':'3','-':'4'}

    matrix_out=[]
    fasta_dic={}
    for line in aln:
        if line[0]==">":
            header=line[1:].rstrip('\n').strip()
            fasta_dic[header]=[]
        elif line[0].isalpha() or line[0]=='-':
            for base, num in dic.items():
                line=line[:].rstrip('\n').strip().replace(base,num)
            line=list(line)
            line=[int(n) for n in line]

            tmp_line = line+[4]*(2000-len(line))
            fasta_dic[header] += tmp_line[0:2000]

    taxa_block=[]
    for taxa in sorted(list(fasta_dic.keys())):
        taxa_block.append(fasta_dic[taxa.strip()])
    fasta_dic={}
    matrix_out.append(taxa_block)

    return np.array(matrix_out)

file_list = os.listdir(folder_fasta)
random.shuffle(file_list)
cnt = 0

for ele in file_list:
    tmp_file = ele.split('.')[0][:-5]
    current_label = dic[tmp_file] 
    current_label = np.array(assign_label(current_label))
    current_seq = get_numpy(folder_fasta, ele)
    seq_dir = folder_numpy_seq + str(cnt) + '.npy'
    label_dir = folder_numpy_label + str(cnt) + '.npy'
    np.save(seq_dir, current_seq)
    np.save(label_dir, current_label)
    print(current_label.shape, current_seq.shape)
    cnt += 1
