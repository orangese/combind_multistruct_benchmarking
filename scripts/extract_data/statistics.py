import sys
from containers import LigandManager
from statistics import Statistics

sys.path.append('../dock')
from shared_paths import shared_paths

name, proteins, features, num_ligs = sys.argv[1:]
proteins = proteins.split(',')
features = features.split(',')
num_ligs = int(num_ligs)

stats_ligs = {}
stats_st = {}
for p in proteins:
    lm = LigandManager(shared_paths, p)
    stats_ligs[p] = lm.docked(lm.pdb)[:num_ligs]
    stats_st[p] = lm.st
stats = Statistics(stats_ligs, stats_st, features)
stats.read(shared_paths['data'], shared_paths['stats'])

for k in features:
    out_f = '/scratch/PI/rondror/combind/bpp_outputs/stats/{}_{}.txt'.format(name, k)
    stats.write(out_f, k)
