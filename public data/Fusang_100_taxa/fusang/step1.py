import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from Bio import AlignIO
aln = AlignIO.read('./t_sim1_TRUE.fas', 'fasta')

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

with open('./distance.mat','w') as f:
    f.write(str(len(aln))+'\n')

    for i in range(0, len(aln)):
        tmp_str = aln[i].id
        for j in range(0, len(aln)):
            tmp_str += ' ' 
            tmp_str += str(dm[i,j])

        f.write(tmp_str+'\n')

with open('./distance.mat_taxlist','w') as f:
    for i in range(0, len(aln)):
        f.write(aln[i].id+'\n')