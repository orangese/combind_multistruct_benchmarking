import os
import pandas as pd
from utils import grouper
from glob import glob

queue = 'rondror'
group_size = 50

def _prep_ligands(ligands, root):   
    os.system('mkdir -p {}'.format(root))
    os.system('rm {}/*.sbatch'.format(root))
    os.system('rm {}/slurm*'.format(root))

    for i, ligs in enumerate(grouper(group_size, ligands)):
        batch_file = '{}/batch-{}.sbatch'.format(root, i)
        with open(batch_file,'w') as batch:
            batch.write('#!/bin/bash\n')
            batch.write('#SBATCH --chdir={}\n'.format(root))
            batch.write('#SBATCH -t 1:00:00\n')
            batch.write('#SBATCH -p {}\n'.format(queue))
            for name, smiles in ligs:
                batch.write('cd {0}; ligprep -WAIT -epik -ismi {0}.smi -omae {0}.mae; cd ..\n'.format(name))
                os.system('rm -rf {}/{}'.format(root, name))
                os.system('mkdir {}/{}'.format(root, name))

                with open('{0}/{1}/{1}.smi'.format(root, name), 'w') as fp:
                    fp.write('{} {}\n'.format(smiles, name))
        os.system('sbatch {}'.format(batch_file))

def prep_ligands(lm):
    ligands = lm.pdb.copy()
    ligands.update(lm.chembl)
    unfinished = []
    for name, (smiles, affinity) in ligands.items():
        if not name in lm.prepped:
            unfinished += [(name, smiles)]

    if unfinished:
        print('Processing {} ligands'.format(len(unfinished)))
        _prep_ligands(unfinished, lm.path('LIGANDS_ROOT'))
