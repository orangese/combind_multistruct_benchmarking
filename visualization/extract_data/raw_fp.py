import numpy as np
import sys
import os

sys.path.append('../../1_dock/')
sys.path.append('../../2_ifp/')
sys.path.append('../../3_analyze/')
from pairs import LigPair
from shared_paths import shared_paths
from containers import Dataset


num_ligs = 20
features = ['sb2', 'mcss', 'hbond', 'pipi', 'contact']

with open('/scratch/PI/rondror/combind/bpp_outputs/fps.csv', 'w') as out:
    out.write('\t'.join(['Protein', 'Struct', 'Ligand1', 'Ligand2', 'Rank1', 'Rank2',
                         'RMSD1', 'RMSD2', 'GSCORE1',  'GSCORE2'] + features) + '\n')

    for protein in os.listdir(shared_paths['data']):
        print(protein)
        if protein[0] == '.': continue

        data = Dataset(shared_paths, [protein])
        lm = data.proteins[protein].lm
        docking = data.proteins[protein].docking[lm.st]
        ligs = lm.docked(lm.pdb)[:num_ligs]
        data.load({protein:ligs}, load_fp=True, load_mcss=False)

        for lig_name, ligand in docking.ligands.items():
            for rank, pose in enumerate(ligand.poses):
                fp = ';'.join(
                    ["{},{},{}".format(interaction, residue, score)
                     for (interaction, residue), score in pose.fp.items()])
                out.write('\t'.join(map(str,
                                        [protein, lm.st, lig_name, rank, pose.gscore, pose.emodel,
                                         pose.rmsd, fp]))+'\n')
