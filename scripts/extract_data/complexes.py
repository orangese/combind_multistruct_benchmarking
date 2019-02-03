"""
Get names of complexes used in docking for all protein
"""
import os
import sys
sys.path.append('../../dock')
sys.path.append('../../2_ifp')
sys.path.append('../../score')
from containers import Dataset
from shared_paths import shared_paths

fname = '/scratch/PI/rondror/combind/bpp_outputs/structs.tsv'
with open(fname, 'w') as fp:
    for protein in os.listdir(shared_paths['data']):
        print(protein)
        if protein[0] == '.': continue
    
        data = Dataset(shared_paths, [protein])
        lm = data.proteins[protein].lm
        for ligand in lm.pdb[:20]:
            fp.write('\t'.join([protein, lm.st, ligand.split('_')[0]])+'\n')
