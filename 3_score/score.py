import sys
import os
import math

import numpy as np
import itertools as it

sys.path.append('../notebooks')
from load_data import *

from optimize import Scores

receptor = sys.argv[1]
structure = sys.argv[2]

output_dir = 'scoring_output2'

def start(receptor, structure):
    glide_ifp='xglide8'
    crystal_ifp='xcrystal8'
    w = [0,0,10,10,10,5,10,10,10,0.2,0,0,0]
    g = 27

    os.system('mkdir -p /scratch/PI/rondror/docking_data/{}/{}'.format(receptor, output_dir))
        
    (xcrystals, xglides) = load_data(receptor,
                                     structure,
                                     w=w,
                                     glide_ifp=glide_ifp,
                                     crystal_ifp=crystal_ifp)
    
    ligs = xglides.keys()
    best_rmsds = [np.min([xglides[l].poses[p].rmsd for p in range(min(25, len(xglides[l].poses.keys())))]) for l in ligs]
    filt_lig = sorted([l for i, l in enumerate(ligs) if l != structure and best_rmsds[i] <= 2])
    #filt_lig = sorted([l for l in xcrystals if l != structure])

    done = [False for l in filt_lig]
    done[0] = True
    done[1] = True

    already_done = load_combinations(receptor, structure, output_dir)
    n = {k:set(already_done.get(k, {}).keys()) for k in range(2,len(filt_lig))}

    total = {k:math.factorial(len(filt_lig))/( math.factorial(k)*math.factorial(len(filt_lig) - k) ) for k in range(len(filt_lig))}

    print 'Receptor {}, Structure {}'.format(receptor, structure)
    #print 'Writing output to {}'.format(output_dir)
    for k in n:
        print '{}: {} done of {} total'.format(k, len(n[k]), total[k])
        if len(n[k]) == total[k]:
            done[k] = True
    #return
    if False not in done: return

    all_scores = Scores(xglides, xcrystals, filt_lig, structure, 25, gscore_weight=g)
    all_pair_scores = all_scores.all_scores
        
    for k in range(2, len(filt_lig)+1):
        with open('/scratch/PI/rondror/docking_data/{}/{}/{}_struct_{}.txt'.format(receptor, output_dir, structure, k), 'a') as output_f:
            output_f.write('# {} g={} {} {}\n'.format(w, g, glide_ifp, crystal_ifp))
    print_cluster(all_scores, receptor, structure, len(filt_lig))
    
    while False in done:
        for k in range(len(filt_lig)):
            if done[k]: continue
            if len(n[k]) < total[k]/4:
                rand_combo = tuple(sorted(np.random.choice(len(filt_lig), size=k, replace=False)))
                while rand_combo in n[k]:
                    rand_combo = tuple(sorted(np.random.choice(len(filt_lig), size=k, replace=False)))
                n[k].add(rand_combo)
                combo_scores = Scores(xglides, xcrystals, [filt_lig[i] for i in rand_combo], 
                                      structure, 25, gscore_weight=g, all_scores = all_pair_scores)
                print_cluster(combo_scores, receptor, structure, k)
            else:
                for rand_combo in it.combinations(range(len(filt_lig)), k):
                    if tuple(rand_combo) in n[k]: continue
                    combo_scores = Scores(xglides, xcrystals, [filt_lig[i] for i in rand_combo], 
                                          structure, 25, gscore_weight=g, all_scores = all_pair_scores)
                    print_cluster(combo_scores, receptor, structure, k)
                done[k] = True

def print_cluster(scores, receptor, structure, k):
    with open('/scratch/PI/rondror/docking_data/{}/{}/{}_struct_{}.txt'.format(receptor, output_dir, structure, k), 'a') as output_f:
        output_f.write(','.join(scores.ligands) + '\n' )    
        output_f.write(','.join([str(int(scores.get_top_num(l, 1, glide=False))) for l in scores.ligands]) + '\n' )

start(receptor, structure)
