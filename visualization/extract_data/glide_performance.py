import os
import sys
from containers import Dataset
sys.path.append('../1_dock')
from shared_paths import shared_paths

fname = '/scratch/PI/rondror/combind/bpp_outputs/glide_performance{}.tsv'.format(shared_paths['docking'])
with open(fname, 'w') as fp:
    for protein in os.listdir(shared_paths['data']):
        print(protein)
        if protein[0] == '.': continue
        data = Dataset(shared_paths, [protein])
        lm = data.proteins[protein].lm
        print(lm.pdb)

        data.load({protein: lm.pdb}, {protein: [lm.st]}, load_mcss = False, load_fp = False)
        
        docking =  data.proteins[protein].docking[lm.st]
        for lig_name, ligand in docking.ligands.items():
            fp.write('\t'.join([
                protein,
                lm.st,
                lig_name,
                ','.join(map(lambda x: str(x.rmsd), ligand.poses)),
                ','.join(map(lambda x: str(x.gscore), ligand.poses))
            ]) + '\n')
