import os
import sys

import numpy as np

import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec

sys.path.append('../4_analyze')
from containers import Dataset
from score_query import ScoreQuery

def load_score_file(f_path, lig_objs):
    go = 0
    with open(f_path, 'r') as f:
        for line in f:
            line = line.strip().split(',')
            if line[0] == 'best_affinity':
                continue
            elif line[0][:4] == 'prot':
                pr = line[0].split('=')[1]
                st = line[1].split('=')[1]
                gdir = line[2].split('=')[1]
                idir = line[3].split('=')[1]
                lval = float(line[4].split('=')[1])
                nposes = int(line[5].split('=')[1])
                w_dict = eval(','.join([line[6].split('=')[1]]+line[7:]))
                go = 1
            elif go == 1:
                l_q = line[0]
                l_i = line[1:]
                sq = ScoreQuery(l_q, l_i, lig_objs, lval, nposes, w_dict)
                nei = {}
                go = 2
                #print len(l_i), lval, nposes, w_dict
            elif go == 2:
                p_list = [int(i) for i in line]
                #print f_path
                #print p_list
                #p_q = int(line[0])
                #p_i = [int(p) for p in line[1:]]
                nei[p_list[0]] = {l:lig_objs[l].poses[p_list[i]] for i,l in enumerate([l_q]+l_i)}
    sq.load(nei)
    return sq

export_script = '/scratch/PI/rondror/jbelk/prospective/code/3_analyze/export_cluster.py'
def export(data, cluster, cluster_name, receptor, struct=None, ligs=None, verbose=False, glide_dir='docking/glide12'):
    glide_dir = '{}/{}/{}'.format(data, receptor, glide_dir)
    out_dir = '{}/outputs'.format('/'.join(data.split('/')[:-1]))
    #os.system('mkdir -p {}'.format(out_dir))
    command = [cluster_name, glide_dir, out_dir]
    #if ligs is None: ligs = self.ligs
    for l,p in cluster.items():
        if struct:
            if l[-3:] == 'lig':
                command.append('{}-to-{}'.format(l, struct))
            else:
                command.append(l)
        else:
            command.append('{}-to-{}'.format(l, l))
        #if self.poses[l].true:
        #    command.append('T')
        #elif self.poses[l].refined:
        #    command.append('R')
        #else:
        command.append(str(p.rank))
    if verbose: print '$SCHRODINGER/run {} {}'.format(export_script, ' '.join(command))
    os.system('$SCHRODINGER/run {} {}'.format(export_script, ' '.join(command)))

def get_fp_vectors(c):
    all_interactions = {}
    fp_map = {}
    for l, p in c.items():
        fp_map[l] = {}
        for r, sc in p.fp.items():
            fp_map[l][r] = sc 
            all_interactions[r] = all_interactions.get( r, 0 ) + sc

    all_interactions = sorted(all_interactions.keys(), key=lambda x: -all_interactions[x])
    return all_interactions, {l: [fp_map[l].get(x, 0) for x in all_interactions] for l in c}

def get_fp_matrix(c, sorted_ligs, all_i=None, num_i=15):
    if all_i is None:
        all_i, fp_vectors = get_fp_vectors(c)
        
    fp_matrix = np.zeros(( len(sorted_ligs), num_i ))
    for i,l in enumerate(sorted_ligs):
        for j,r in enumerate(all_i[:num_i]):#range(num_i):
            #r = all_i[j]
            fp_matrix[i,j] = c[l].fp.get(r, 0)
    return all_i[:num_i], fp_matrix

def get_rmsd_matrix(c, sorted_ligs):
    fp_matrix = np.zeros( (len(sorted_ligs), 1) )
    for i, l in enumerate(sorted_ligs):
        if c[l].rmsd <= 2:
            fp_matrix[i,0] = 2
        elif c[l].rmsd <= 4:
            fp_matrix[i,0] = 5
        else:
            fp_matrix[i,0] = 0
    return fp_matrix

def show_side_by_side(c1, c2, sorted_ligs, t1='Cluster 1', t2='Cluster 2', num_i=15):
    i_key, fp_mat1 = get_fp_matrix(c1, sorted_ligs, num_i=num_i)
    i_key, fp_mat2 = get_fp_matrix(c2, sorted_ligs, all_i=i_key, num_i=num_i)
    fp_mat3 = fp_mat2 - fp_mat1

    heat_max = max(np.max(fp_mat1), np.max(fp_mat2))

    r_mat1 = get_rmsd_matrix(c1, sorted_ligs)
    r_mat2 = get_rmsd_matrix(c2, sorted_ligs)

    f, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(1, 5, gridspec_kw = {'width_ratios':[1,1]+[num_i]*3}, dpi=400)

    ax0.matshow(r_mat1, cmap='Set1', vmin=0, vmax=9)
    ax0.set_yticks(np.arange(len(sorted_ligs)), minor=False)
    ax0.set_yticklabels(sorted_ligs, minor=False, size=14)

    ax1.matshow(r_mat2, cmap='Set1', vmin=0, vmax=9)
    ax1.set_yticks(np.arange(len(sorted_ligs)), minor=False)

    interaction_heatmap(fp_mat1, sorted_ligs, i_key, m=heat_max, fig=f, ax=ax2)
    interaction_heatmap(fp_mat2, sorted_ligs, i_key, m=heat_max, fig=f, ax=ax3)
    interaction_heatmap(fp_mat3, sorted_ligs, i_key, fig=f, ax=ax4, difference=True)

    ax0.set_xticks(np.arange(1), minor=False)
    ax0.xaxis.tick_top()
    ax0.set_xticklabels(['RMSD 1'], minor=False, rotation = 'vertical', size=14)
    ax1.set_xticks(np.arange(1), minor=False)
    ax1.xaxis.tick_top()
    ax1.set_xticklabels(['RMSD 2'], minor=False, rotation = 'vertical', size=14)

    ax2.set_xlabel(t1, size=20)
    ax3.set_xlabel(t2, size=20)
    ax4.set_xlabel('Difference', size=20)

    ax1.set_yticklabels(['' for i in sorted_ligs])
    ax2.set_yticklabels(['' for i in sorted_ligs])
    ax3.set_yticklabels(['' for i in sorted_ligs])
    ax4.set_yticklabels(['' for i in sorted_ligs])

    plt.tight_layout()
    plt.show()

def interaction_heatmap(A, structs, res, fname='', m=None, fig=None, ax=None, difference=False):
    sq_size = 0.9
    if fig is None:
        fig, ax = plt.subplots()
        sq_size = 0.3

    fig.set_size_inches(sq_size*len(res), sq_size*len(structs), forward=True)

    def i_matrix(A, res, i):
        aa = np.zeros(A.shape)
        for j, r in enumerate(res):
            if r[0] == i: aa[:,j] = A[:,j]
            else: aa[:,j] = np.nan*A[:,j]
        return aa

    colors = {2:cm.Oranges, 3:cm.Oranges,11:cm.Greys, 13:cm.Oranges, 14:cm.Oranges, 7:cm.Purples,
              8:cm.Greens, 9:cm.Greens, 6:cm.Purples, 4:cm.Greens, 5:cm.Purples, 10:cm.Blues}

    if difference:
        max_abs = np.max(np.abs(A))
        ax.matshow(A, cmap='bwr', vmin=-max_abs, vmax=max_abs)
    else:
        if m is None: m = np.max(A)
        for i in range(20):
            ax.matshow(i_matrix(A, res, i), cmap=colors.get(i, cm.Blues), vmin=0, vmax=m)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]), minor=False)
    ax.set_yticks(np.arange(A.shape[0]), minor=False)
    ax.xaxis.tick_top()
    #print res
    xlabels = ['{}, {}-{}'.format(i, r.split('(')[1].split(')')[0], r.split(':')[1].split('(')[0]) for (i,r) in res]
    ax.set_xticklabels(xlabels, minor=False, rotation = 'vertical', size=14)
    ax.set_yticklabels(structs, minor=False, size=14)
