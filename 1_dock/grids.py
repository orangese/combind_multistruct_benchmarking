import os
import sys

from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid

def make_grids():

    os.system('mkdir -p docking')
    os.system('mkdir -p docking/grids')

    for lig in os.listdir('structures/ligands') + ['3J5Q_lig']:
        pdb = lig.split('_')[0]
        
        if os.path.exists('docking/grids/{}/{}.zip'.format(pdb, pdb)): 
            continue
        if not os.path.exists('structures/proteins/{}_prot.mae'.format(pdb)):
            continue
        os.system('rm -rf docking/grids/{}'.format(pdb))        
        os.system('mkdir -p docking/grids/{}'.format(pdb))

        if lig == '3J5Q_lig': lig = '5IRX_lig.mae'
        st_2 = StructureReader('structures/ligands/{}'.format(lig)).next()
        c2 = get_centroid(st_2)
        x,y,z = c2[:3]
        #else:
        #    x,y,z = -21.5,5.5,-16.5

        with open('docking/grids/{}/{}.in'.format(pdb, pdb), 'w') as f:
            f.write('GRID_CENTER {},{},{}\n'.format(x,y,z))
            f.write('GRIDFILE {}.zip\n'.format(pdb))
            f.write('INNERBOX 15,15,15\n')
            f.write('OUTERBOX 32,32,32\n')
            f.write('RECEP_FILE ../../../structures/proteins/{}_prot.mae\n'.format(pdb))

        with open('docking/grids/{}/grid_in.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('$SCHRODINGER/glide -WAIT {}.in'.format(pdb))
        print 'making grid', pdb
        os.chdir('docking/grids/{}'.format(pdb))
        os.system('sbatch -p rondror -t 00:30:00 -o grid.out grid_in.sh')
        os.chdir('../../..')

