import numpy as np
import sys
import os

sys.path.append('../1_dock/')
sys.path.append('../2_ifp/')
sys.path.append('../3_analyze/')
from pairs import LigPair
from shared_paths import shared_paths
from containers import Dataset


num_ligs = 20
features = ['sb2', 'mcss', 'hbond', 'pipi', 'contact']

with open('/scratch/PI/rondror/combind/bpp_outputs/stats.csv', 'w') as fp:
    fp.write('\t'.join(['Protein', 'Struct', 'Ligand1', 'Ligand2', 'Rank1', 'Rank2',
                        'RMSD1', 'RMSD2', 'GSCORE1',  'GSCORE2'] + features) + '\n')

    for protein in os.listdir(shared_paths['data']):
        print(protein)
        if protein[0] == '.': continue

        data = Dataset(shared_paths, [protein])
        lm = data.proteins[protein].lm
        docking = data.proteins[protein].docking[lm.st]
        ligs = lm.docked(lm.pdb)[:num_ligs]
        data.load({protein:ligs}, load_fp=True, load_mcss=True)

        for i1, l1 in enumerate(ligs):
            for l2 in ligs[i1+1:]:
                ligpair = LigPair(docking.ligands[l1], docking.ligands[l2],
                                  features, lm.mcss, max_poses=100, normalize_fp=True)
                ligpair.init_pose_pairs()

                for rank1, rank2 in ligpair.pose_pairs:
                    g1, g2 = ligpair.get_gscores(rank1, rank2)
                    rmsd1, rmsd2 = ligpair.get_rmsds(rank1, rank2)
                    x = [ligpair.get_feature(feature, rank1, rank2) for feature in features]
                    
                    fp.write('\t'.join(map(str, [protein, lm.st, l1, l2, rank1, rank2, rmsd1, rmsd2,
                                                 g1, g2] + x)) + '\n')
