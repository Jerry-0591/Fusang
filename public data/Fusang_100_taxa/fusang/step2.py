import os
import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from Bio import AlignIO
aln = AlignIO.read('./t_sim1_TRUE.fas', 'fasta')

if not os.path.exists('./tmp_fasta'):
    os.mkdir('./tmp_fasta')

node_cnt = len(aln)

file_num = node_cnt // 20
mode = node_cnt % 20

for i in range(0,file_num):
    with open('./tmp_fasta/{}.fas'.format(i),'w') as f:
        for j in range(i*20,i*20+20):
            record = aln[j]
            f.write('>'+str(record.id)+'\n')
            f.write(str(record.seq)+'\n')

if mode != 0:
    with open('./tmp_fasta/{}.fas'.format(i),'w') as f:
        for j in range(file_num*20,node_cnt):
            record = aln[j]
            f.write('>'+str(record.id)+'\n')
            f.write(str(record.seq)+'\n')
