import os
from glob import glob
from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid

GRID_IN = """
GRID_CENTER {x},{y},{z}
GRIDFILE {pdb}.zip
INNERBOX 15,15,15
OUTERBOX 30,30,30
RECEP_FILE ../../../structures/proteins/{pdb}_prot.mae
"""

CMD = ('sbatch -p rondror -t 00:30:00 -o grid.out'
       ' --wrap="glide -WAIT {infile}"'
       ' --chdir={wd}')

WD = 'docking/grids/{pdb}/'
INFILE = '{pdb}.in'
ZIPFILE = '{pdb}.zip'

PROTFILE = 'structures/proteins/{pdb}_prot.mae'
LIGFILE = 'structures/ligands/{pdb}_lig.mae'

def make_grids(pdb=None):
    if not glob(PROTFILE.format(pdb='*')):
        return

    if pdb is None:
        pdb = sorted(glob(PROTFILE.format(pdb='*')))[0]
        pdb = pdb.split('/')[-1].split('_')[0]

    if os.path.exists((WD+ZIPFILE).format(pdb=pdb)):
        return # Done.
    if not (os.path.exists(LIGFILE.format(pdb=pdb))
            and os.path.exists(PROTFILE.format(pdb=pdb))):
        return # Not ready.

    print('making grid', pdb)

    for path in glob(WD.format(pdb=pdb) + '*'):
        os.remove(path)
    os.makedirs(WD.format(pdb=pdb), exist_ok=True)

    st = next(StructureReader(LIGFILE.format(pdb=pdb)))
    c = get_centroid(st)
    x,y,z = c[:3]
        
    with open((WD+INFILE).format(pdb=pdb), 'w') as fp:
        fp.write(GRID_IN.format(x=x, y=y, z=z, pdb=pdb))

    os.system(CMD.format(infile=INFILE.format(pdb=pdb),
                         wd=WD.format(pdb=pdb)))
