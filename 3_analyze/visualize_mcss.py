import sys
import os
from matplotlib import cm
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.analyze import evaluate_smarts

from containers import LigandManager

sys.path.append('../1_dock')
from  MCSSController import MCSSController
from shared_paths import shared_paths

def color(i, total):
    rgba = cm.jet(i/total)
    return map(lambda x: int(256*x), rgba[:3])
    

MCSS_VIZ = '/scratch/PI/rondror/combind/bpp_outputs/mcss_viz'

l1, l2, prot = sys.argv[1:]
if l1 > l2: l1, l2 = l2, l1

lm = LigandManager(shared_paths, prot)
lm.mcss.load_mcss()

mcss = lm.mcss.MCSSs["{}-{}".format(l1, l2)]

D = "{}/{}-{}-{}".format(MCSS_VIZ, prot, l1, l2)

os.system('mkdir -p {}'.format(D))

# Write MCSS    
with open("{}/mcss.txt".format(D), 'w') as fp:
    fp.write(str(mcss)+'\n')
    fp.write(str(mcss.is_valid())+'\n')

# Write ligands to MAE file
l1_struct = next(StructureReader(lm.mcss.lig_template.format(l1)))
l2_struct = next(StructureReader(lm.mcss.lig_template.format(l2)))

with StructureWriter("{}/ligs.mae".format(D)) as writer:
    writer.append(l1_struct)
    writer.append(l2_struct)

# Write projected ASLs
l1_atom_idxss = [evaluate_smarts(l1_struct, smarts, unique_sets=True)
                 for smarts in mcss.smarts_l1]
l2_atom_idxss = [evaluate_smarts(l2_struct, smarts, unique_sets=True)
                 for smarts in mcss.smarts_l2]

for k, (l1_atom_idxs, l2_atom_idxs) in enumerate(zip(l1_atom_idxss, l2_atom_idxss)):
    j = 0
    for l1_atom_idx in l1_atom_idxs:
        for l2_atom_idx in l2_atom_idxs:
            with StructureWriter('{}/mcss_{}_{}.mae'.format(D, k, j)) as writer:
                for i, (l1_idx, l2_idx) in enumerate(zip(l1_atom_idx, l2_atom_idx)):
                    l1_struct.atom[l1_idx].setColorRGB(*color(i, len(l1_atom_idx)))
                    l2_struct.atom[l2_idx].setColorRGB(*color(i, len(l1_atom_idx)))
                l1_sub = l1_struct.extract(l1_atom_idx)
                l2_sub = l2_struct.extract(l2_atom_idx)
                writer.append(l1_sub)
                writer.append(l2_sub)
            j += 1
