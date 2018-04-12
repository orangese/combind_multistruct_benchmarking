import os
import sys

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty
from schrodinger.structutils.analyze import evaluate_smarts

data = '/scratch/PI/rondror/jbelk/method/data'
prot = sys.argv[1]
os.chdir('{}/{}'.format(data, prot))

in_dir = 'ligands/mcss/mcss7'
out_dir = 'mcss/mcss_debug2'

os.system('mkdir -p {}'.format(out_dir))

#all_mcss = sorted([m.split('.')[0] for m in os.listdir(in_dir) if m.split('.')[-1] == 'mae'])
#all_mcss = [tuple(m.split('-')) for m in all_mcss if m[-3:] != '_in']
all_mcss = [('CHEMBL2426396_lig', 'CHEMBL2426384_lig'),('CHEMBL2426396_lig', 'CHEMBL2426392_lig')]
for l1,l2 in all_mcss:
    #if l1[:6] == 'CHEMBL' or l2[:6] == 'CHEMBL':
    #    continue
    if os.path.exists('{}/{}-{}.mae'.format(out_dir,l1,l2)):
        continue
    print l1, l2

    out = StructureWriter('{}/{}-{}.mae'.format(out_dir,l1,l2))
    for st_out in StructureReader('{}/{}-{}.mae'.format(in_dir,l1,l2)):
        smarts = _StructureProperty(st_out)['s_canvas_MCS_SMARTS']
        mcss = evaluate_smarts(st_out, smarts, unique_sets=True)
        for match in mcss:
            out.append(st_out.extract(match))
    out.close()


