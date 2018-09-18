import os
import sys
from containers import LigandManager
from statistics import gscore_statistics
sys.path.append('../1_dock/')
from shared_paths import shared_paths

n_ligs = 20
max_poses = 100
native_thresh = 2.0
points = 1000

if sys.argv[1] == 'scaled':
    scale_by_top = True
    sd = 0.4
    domain = (0, 10)
    out = '/scratch/PI/rondror/combind/bpp_outputs/glide_scaled.'
elif sys.argv[1] == 'unscaled':
    scale_by_top = False
    sd = 0.4
    domain = (-16, 1)
    out = '/scratch/PI/rondror/combind/bpp_outputs/glide.'
    

datasets = [d for d in sorted(os.listdir(shared_paths['data']))
            if d[0] != '.' and d[-3:] != 'old']

os.chdir(shared_paths['data'])
data = {}
for d in datasets:
    os.chdir(d)
    lm = LigandManager(shared_paths, d)
    data[d] = lm.docked(lm.pdb)[:n_ligs]
    os.chdir('..')

native, reference = gscore_statistics(data, max_poses, native_thresh, points,
                                      scale_by_top, sd, domain)


native.write(out+'native.de')
reference.write(out+'reference.de')

ratio = native.ratio(reference, prob = False)

ratio.write(out+'pnative.de')

