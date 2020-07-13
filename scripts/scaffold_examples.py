"""
Pick out examples of where docking fails to produce any good poses for one
ligand, but succeeds for a ligand that shares a scaffold.

View the xtal + best docked poses for the ligands using the commands in the
"scaffold_examples_pymol.txt" file. First navigate to the root of the data
directory, then copy/paste the commands for the pair of interest.
"""

import os
from containers import Protein
import utils
import config

paths = {'CODE': os.path.realpath('../'),
         'DATA': '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems',
         'PDB': '{ROOT}/structures/pdb.csv'}
paths.update(config.PATHS)
paths = utils.resolve(paths)
params = config.STATS['rd1']

proteins = utils.get_proteins(paths, [])

data = []
for prot in proteins:
    protein = Protein(prot, params, paths)
    protein.lm.mcss.load_mcss()
    ligands = protein.lm.get_xdocked_ligands(20)
    protein.load_docking(ligands)
    for i, ligand1 in enumerate(ligands):
        for ligand2 in ligands[i+1:]:
            if protein.lm.mcss.get_mcss_size(ligand1, ligand2) > 0.5:
                rmsd1 = float('inf')
                pose1 = None
                for pose in protein.docking[ligand1].poses:
                    if pose.rmsd < rmsd1:
                        rmsd1 = pose.rmsd
                        pose1 = pose.rank
                        
                rmsd2 = float('inf')
                pose2 = None
                for pose in protein.docking[ligand2].poses:
                    if pose.rmsd < rmsd2:
                        rmsd2 = pose.rmsd
                        pose2 = pose.rank
                data += [(prot, ligand1, ligand2, protein.lm.st, rmsd1, rmsd2, pose1, pose2)]

with open('scaffold_examples.txt', 'w') as fp:
	fp.write('protein ligand1 ligand2 docking_struct rmsd1 rmsd2 best_pose1 best_pose2\n')
	for pair in data:
	    if min(pair[4:6]) < 2.0 and max(pair[4:6]) > 3.0:
	        fp.write(' '.join(map(str, pair))+'\n')

with open('scaffold_examples_pymol.txt', 'w') as fp:
	for pair in data:
	    if min(pair[4:6]) < 2.0 and max(pair[4:6]) > 3.0:
	        fp.write('load {0}/structures/ligands/{1}.mae\n'.format(*pair))
	        fp.write('load {0}/structures/ligands/{2}.mae\n'.format(*pair))
	        fp.write('load {0}/structures/ligands/{3}_lig.mae\n'.format(*pair))
	        fp.write('load {0}/structures/proteins/{3}_prot.mae\n'.format(*pair))
	        
	        if pair[-2] == 0:
	            rank = ''
	        elif pair[-2] < 10:
	            rank = '_0{}'.format(pair[-2])
	        else:
	            rank = '_{}'.format(pair[-2])
	        fp.write('load {0}/docking/confgen_es4/{1}-to-{3}/{1}-to-{3}_pv.maegz\n'.format(*pair))
	        fp.write('set_name {1}-to-{3}_pv.{1}{8}, {1}_dock\n'.format(*(pair+(rank,))))
	        fp.write('ungroup {1}_dock\n'.format(*pair))

	        if pair[-1] == 0:
	            rank = ''
	        elif pair[-1] < 10:
	            rank = '_0{}'.format(pair[-1])
	        else:
	            rank = '_{}'.format(pair[-1])
	        fp.write('load {0}/docking/confgen_es4/{2}-to-{3}/{2}-to-{3}_pv.maegz\n'.format(*pair))
	        fp.write('set_name {2}-to-{3}_pv.{2}{8}, {2}_dock\n'.format(*(pair+(rank,))))
	        fp.write('ungroup {2}_dock\n'.format(*pair))
	        
	        fp.write('delete *-to-*\n')
	        
	        fp.write('show cartoon, *prot*\n')
	        fp.write('as sticks, *lig*\n')
	        fp.write('hide everything, element H and not (element N+O+S extend 1)\n')

	        if pair[4] < 2.0:
	            fp.write('util.cbag {1}*\n'.format(*pair))
	            fp.write('util.cbay {2}*\n'.format(*pair))
	        else:
	            fp.write('util.cbag {2}*\n'.format(*pair))
	            fp.write('util.cbay {1}*\n'.format(*pair))
	        fp.write('util.cbaw {3}*\n'.format(*pair))
	        fp.write('zoom {3}_lig\n'.format(*pair))
	        fp.write('\n')
