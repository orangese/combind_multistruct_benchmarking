import os
import subprocess
from glob import glob
import pandas as pd

# # for i in *; do cd $i; sbatch -p rondror --wrap="combind featurize . docking/*/*_pv.maegz --processes 6" -J $i -t 06:00:00 -n 6; cd -; done;

# data_root = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020_v2'
# ligands_fname = 'pdb.smi'

# for protein in glob(data_root + '/*'):
# 	print(protein)
# 	cwd = os.getcwd()
# 	os.chdir(protein)

# 	ligands = pd.read_csv(ligands_fname, sep=' ')
# 	struct = sorted(ligands['ID'])[0]

# 	# Protein preparation

# 	# Ligand preparation

# 	# Docking

# 	self_docking_pv = ('docking/{struct}-to-{struct}/'
# 	                   '{struct}-to-{struct}_pv.maegz').format(struct=struct)
# 	self_docking_native_pv = self_docking_pv.replace('_pv', '_native_pv')
# 	struct_lig = ('structures/ligands/{struct}_lig.mae').format(struct=struct)

# 	# Filter-native poses for template structure.
# 	if (os.path.exists(self_docking_pv)
# 	    and os.path.exists(struct_lig)
# 	    and not os.path.exists(self_docking_native_pv)):
# 		filter_native = 'combind filter-native {} {}'
# 		filter_native = filter_native.format(self_docking_pv, struct_lig)
# 		print(filter_native)
# 		subprocess.call(filter_native, shell=True)

# 	# Featurization (not including native)
# 	docking = glob('docking/*/*_pv.maegz')
# 	docking.remove(self_docking_native_pv)
# 	featurize = 'combind featurize . ' + ' '.join(docking)
# 	print(featurize)
# 	subprocess.call(featurize, shell=True)
	
# 	# Featurization (to native)
# 	docking = glob('docking/*/*_pv.maegz')
# 	docking.remove(self_docking_pv)
# 	featurize = 'combind featurize . ' + ' '.join(docking)
# 	print(featurize)
# 	subprocess.call(featurize, shell=True)

# 	os.chdir(cwd)

import sys
sys.path.append('/home/users/jpaggi/combind')
from features.features import Features
from score.pose_prediction import PosePrediction
from score.statistics import read_stats
max_poses = 100


helpers = {}
with open('chembl/helpers/stats41-best_affinity_diverse.txt') as fp:
	for line in fp:
		line = line.strip()
		k, v = line.split(':')
		v = v.split(',')
		helpers[k] = v

def get_pv_name(ligand):
	pv = glob('docking/{}*/*_pv.maegz'.format(ligand))
	assert len(pv) <= 1, pv
	if pv:
		return pv[0]
	return None

pvs = []
for query, helper in helpers.items():
	_pvs = [get_pv_name(query)]
	for _helper in helper:
		pv = get_pv_name(_helper)
		if pv:
			_pvs += [pv]
	pvs += [_pvs]

# features = Features('.', ifp_version='rd1', shape_version='pharm_max',
#                     mcss_version='mcss16', max_poses=100)
# features.compute_single_features(pvs, processes=16)
# features.compute_pair_features(pvs, processes=16)

###############################################################################

features = 'mcss,hbond,saltbridge,contact'
features = features.split(',')

stats = read_stats('stats/', features)

for _pvs in pvs:
	protein = Features('.', ifp_version='rd1', mcss_version='mcss16', max_poses=max_poses)
	protein.load_features(pvs=_pvs, features=features)

	ligands = [_pv.split('/')[-1].split('.')[0] for _pv in _pvs]
	ps = PosePrediction(ligands, protein.raw, stats, [], features,
	                    max_poses, 1.0, -7.0)
	best_poses = ps.max_posterior(1000, 50)

	print(best_poses)

	out = 'scores/paper_dev/{}.csv'.format(ligands[0])

	with open(out, 'w') as fp:
	    fp.write('ID,POSE,PROB,COMBIND_RMSD,GLIDE_RMSD,BEST_RMSD\n')
	    for ligand in best_poses:
	        if 'rmsd' in protein.raw and ligand in protein.raw['rmsd']:
	            rmsds = protein.raw['rmsd'][ligand]
	            grmsd = rmsds[0]
	            crmsd = rmsds[best_poses[ligand]]
	            brmsd = min(rmsds)
	        else:
	            grmsd, crmsd, brmsd = None, None, None
	        fp.write(','.join(map(str, [ligand,
	                                    best_poses[ligand],
	                                    0.0,
	                                    crmsd, grmsd, brmsd]))+ '\n')
