import os
import sys

from schrodinger.structure import StructureReader, StructureWriter

os.chdir('../../data/CDK2')

ref_ligs = 'jain_files/TestRef.mol2'
dock_ligs = 'jain_files/TestMols.mol2'

for st in StructureReader(ref_ligs):
    pdb_id = st._getTitle().split('-')[0]

    st_path = 'structures/ligands/{}_lig.mae'.format(pdb_id)
    if os.path.exists(st_path): continue

    st._setTitle('{}_lig'.format(pdb_id))
    st_wr = StructureWriter(st_path)
    st_wr.append(st)
    st_wr.close()

for st in StructureReader(dock_ligs):
    pdb_id = st._getTitle().split('-')[0]

    st_path = 'ligands/raw_files/{}_lig.mae'.format(pdb_id)
    if os.path.exists(st_path): continue

    st._setTitle('{}_lig'.format(pdb_id))
    st_wr = StructureWriter(st_path)
    st_wr.append(st)
    st_wr.close()

