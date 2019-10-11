import os

metric = 'tanimoto_exclude'
poses_equal = False
ligands_equal = False
max_poses = 100
native_poses = 100

for reference_poses in [1, 5, 100]:
      for stats_sd in [0.02, 0.03, 0.04, 0.05, 0.06, 0.07]:
            root = '/oak/stanford/groups/rondror/users/jpaggi/stats'
            name = '-'.join(map(str, [metric, poses_equal, ligands_equal, max_poses, native_poses, reference_poses, stats_sd]))
            fname = 'run.py'

            script = """import os
import sys
sys.path.append(os.environ['COMBINDHOME'])

from score.statistics import Statistics
from shared_paths import proteins, feature_defs, StringFunction

settings = {}'metric'         : "{}",
            'ifp_version'     : 'ifp5',
            'mcss_version'    : 'mcss16',
            'mcss_func'       : min,
            'mcss_rel_min'    : 0.5,
            'mcss_abs_min'    : 10,
            'docking_version' : 'confgen_es4',
            'pdb_order'       : 'First',
            'poses_equal'     : {},
            'ligands_equal'   : {},
            'native_thresh'   : 2.0,
            'n_ligs'          : 20,
            'max_poses'       : {},
            'native_poses'    : {},
            'reference_poses' : {},
            'stats_sd'        : {},
            'stats_points'    : 100,{}

Statistics(proteins, feature_defs, settings, path="{}/{}/{}-{}-{}.de")
""".format('{', metric, poses_equal, ligands_equal, max_poses, native_poses, reference_poses, stats_sd, '}', root, name, '{}', '{}', '{}')

            if os.path.exists('{}/{}'.format(root, name)): continue
            print(script)
            os.mkdir('{}/{}'.format(root, name))
            with open('{}/{}/{}'.format(root, name, fname), 'w') as fp:
                 fp.write(script)

            os.system('sbatch -p owners -t 04:00:00 -D {}/{} --wrap="$SCHRODINGER/run {}"'.format(root, name, fname))
