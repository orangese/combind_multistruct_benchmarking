import os
import sys

from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid

def make_grids():

    os.system('mkdir -p docking')
    os.system('mkdir -p docking/grids')

    for prot in os.listdir('structures/proteins'):
        pdb = prot.split('_')[0]
        
        if os.path.exists('docking/grids/{}/{}.zip'.format(pdb, pdb)): 
            continue
        if not (os.path.exists('structures/proteins/{}_prot.mae'.format(pdb))
                and os.path.exists('structures/ligands/{}_lig.mae'.format(pdb))):
            continue
        os.system('rm -rf docking/grids/{}'.format(pdb))        

        if pdb == '3J5Q':
            x,y,z = 13,-30,-33
        else:
            st_2 = next(StructureReader('structures/ligands/{}_lig.mae'.format(pdb)))
            c2 = get_centroid(st_2)
            x,y,z = c2[:3]
                      
        out_f = pdb  
        os.system('mkdir -p docking/grids/{}'.format(out_f))
        with open('docking/grids/{}/{}.in'.format(out_f,out_f), 'w') as f:
            f.write('GRID_CENTER {},{},{}\n'.format(x,y,z))
            f.write('GRIDFILE {}.zip\n'.format(out_f))
            f.write('INNERBOX 15,15,15\n')
            f.write('OUTERBOX 30,30,30\n')
            f.write('RECEP_FILE ../../../structures/proteins/{}_prot.mae\n'.format(pdb))

        with open('docking/grids/{}/grid_in.sh'.format(out_f), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('$SCHRODINGER/glide -WAIT {}.in'.format(out_f))
                        
        print('making grid', out_f)
        os.chdir('docking/grids/{}'.format(out_f))
        os.system('sbatch -p owners -t 00:30:00 -o grid.out grid_in.sh')
        os.chdir('../../..')
