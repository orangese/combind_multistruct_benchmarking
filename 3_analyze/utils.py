import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def plot_docking(rmsds_list, title_list, plt_title=''):
    plt.figure(figsize=(6,6))
    plt.plot([2,2],[0,1],'--k')
    
    prop_ligands = np.cumsum([-1.0/len(rmsds_list[0])] + 
                             [1.0/len(rmsds_list[0]) for i in rmsds_list[0]] + 
                             [1.0/len(rmsds_list[0])])

    count_none = len([i for i in rmsds_list[0] if i is None])
    if count_none != 0:  
        print(plt_title, count_none, 'did not dock')
        frac_docked = 1 - float(count_none)/float(len(rmsds_list[0]))
        plt.plot([0,8],[frac_docked, frac_docked], '--k')
    for i, r in enumerate(rmsds_list):
        r = [j if j is not None else 100 for j in r]
        x = sorted([0] + r + [max(r)])
        plt.step(x, prop_ligands, label=title_list[i], linewidth=2)

    plt.gca().set_xlim([0,6])
    plt.gca().set_ylim([0,1])
    plt.xlabel('Error [RMSD, $\AA$]', size=20)
    plt.ylabel('Cumulative Proportion of Ligands', size=20)
    plt.title(plt_title, size=24)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=20)#, ncol=3)
    #plt.legend()

    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.gca().title.set_position([.5, 1.00])
    x0, x1 = plt.xlim()
    y0, y1 = plt.ylim()
    plt.gca().set_aspect(abs(x1-x0)/abs(y1-y0))
    
    plt.show()

def stats_hist(green_dist,blue_dist,normalize=True,logratio=True):

    if green_dist is None or blue_dist is None: return

    plt.rcParams['figure.figsize'] = (10,7.5)
    plt.rcParams['figure.dpi'] = 400 
    plt.rcParams['font.size'] = 20

    if normalize:
        xmin, xmax = -0.3, 1.3
        binwidth = 0.05
    else:
        allx = set(green_dist.prob.keys()+blue_dist.prob.keys())
        xmin, xmax = min(allx),max(allx)
        binwidth = (xmax-xmin)/100
    bins = np.arange(xmin, xmax + binwidth, binwidth)
    
    g_x = sorted(green_dist.prob.keys())#np.arange(xmin - std,xmax + std,std*0.01)#0.02)
    g_px = [green_dist.prob[x] for x in g_x]

    b_x = sorted(blue_dist.prob.keys())
    b_px = [blue_dist.prob[x] for x in b_x]

    plt.plot(b_x, b_px, linewidth=2, color='b')
    plt.plot(g_x, g_px, linewidth=2, color='g')

    plt.fill_between(b_x, 0, b_px, facecolor='b',alpha=0.2)
    plt.fill_between(g_x, 0, g_px, facecolor='g',alpha=0.2)

    plt.xlim(xmin, xmax)

    ax = plt.subplot(111)  
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().set_visible(False)#.tick_left()

    if normalize: plt.xticks([0,1])
    plt.tick_params(axis='x', direction='inout', length=13, width=3, pad=15)
    ax.axhline(linewidth=4, color='k')#, color="g")

    if logratio:
        ax2 = ax.twinx()
        x2 = [x for x in g_x if x >= 0 and x <= 1]
        y2 = [np.log10(green_dist.prob[x]/blue_dist.prob[x]) for x in x2]
        xmax2 = x2[np.argmax(y2)]
        ax2.plot(x2, y2, linewidth=2, color='k')
        ax2.plot([xmax2,xmax2],[1.1*min(y2),1.1*max(y2)],'--k')
        plt.gca().set_ylim([max(-2,min(y2)),1.1*max(y2)])
        ax2.spines["top"].set_visible(False)  
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.get_yaxis().set_visible(False)#.tick_left()
        #plt.gca().y

    plt.gca().set_xlim([xmin,xmax])
    plt.show()

def heatmap(A, row_vals, col_vals, red=10.25):
    fig, ax = plt.subplots()
    
    sq_size = 0.4
    fig.set_size_inches(sq_size*A.shape[0] + 2, sq_size*A.shape[1], forward=True)
    
    cmap = cm.jet
    cmap.set_under('black')

    heatmap = ax.pcolor(A, cmap=cmap, vmin=0, vmax=red)
    ax.set_xlim((0,A.shape[0]))
    ax.set_ylim((0,A.shape[1]))
    plt.colorbar(heatmap)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(A.shape[0]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #labels
    
    ax.set_xticklabels(col_vals, minor=False, rotation = 'vertical')
    ax.set_yticklabels(row_vals, minor=False)
    plt.show()

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
def export(data, cluster, cluster_name, receptor, struct=None, ligs=None, verbose=False, glide_dir='docking/glide12epik'):
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
    if verbose: print('$SCHRODINGER/run {} {}'.format(export_script, ' '.join(command)))
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

def show_side_by_side(c1, c2, sorted_ligs, t1='Cluster 1', t2='Cluster 2', num_i=15, size=1):
    i_key, fp_mat1 = get_fp_matrix(c1, sorted_ligs, num_i=num_i)
    i_key, fp_mat2 = get_fp_matrix(c2, sorted_ligs, all_i=i_key, num_i=num_i)
    fp_mat3 = fp_mat1 - fp_mat2

    heat_max = max(np.max(fp_mat1), np.max(fp_mat2))

    r_mat1 = get_rmsd_matrix(c1, sorted_ligs)
    r_mat2 = get_rmsd_matrix(c2, sorted_ligs)

    f, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(1, 5, gridspec_kw = {'width_ratios':[1,1]+[num_i]*3}, dpi=400)

    ax0.matshow(r_mat1, cmap='Set1', vmin=0, vmax=9)
    ax0.set_yticks(np.arange(len(sorted_ligs)), minor=False)
    ax0.set_yticklabels(sorted_ligs, minor=False, size=14)

    ax1.matshow(r_mat2, cmap='Set1', vmin=0, vmax=9)
    ax1.set_yticks(np.arange(len(sorted_ligs)), minor=False)

    interaction_heatmap(fp_mat1, sorted_ligs, i_key, m=heat_max, fig=f, ax=ax2, size=size)
    interaction_heatmap(fp_mat2, sorted_ligs, i_key, m=heat_max, fig=f, ax=ax3, size=size)
    interaction_heatmap(fp_mat3, sorted_ligs, i_key, fig=f, ax=ax4, difference=True, size=size)

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

def interaction_heatmap(A, structs, res, fname='', m=None, fig=None, ax=None, difference=False, size=1):
    sq_size = 0.9*size
    if fig is None:
        fig, ax = plt.subplots()
        sq_size = 0.3*size

    fig.set_size_inches(sq_size*A.shape[0], sq_size*A.shape[1], forward=True)

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

def feature_heatmap(A, structs, fname='', ma=None, mi=None, fig=None, ax=None, difference=False, size=1):
    sq_size = 0.9*size
    if fig is None:
        fig, ax = plt.subplots()
        sq_size = 0.3*size

    fig.set_size_inches(sq_size*A.shape[0], sq_size*A.shape[1], forward=True)

    colors = {2:cm.Oranges, 3:cm.Oranges,11:cm.Greys, 13:cm.Oranges, 14:cm.Oranges, 7:cm.Purples,
              8:cm.Greens, 9:cm.Greens, 6:cm.Purples, 4:cm.Greens, 5:cm.Purples, 10:cm.Blues}

    if difference:
        max_abs = np.max(np.abs(A))
        img = ax.matshow(A, cmap='bwr', vmin=-max_abs, vmax=max_abs)
    else:
        if ma is None: ma = np.max(A)
        if mi is None: mi = np.min(A)
        img = ax.matshow(A, cmap=cm.Greys, vmin=mi, vmax=ma)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]), minor=False)
    ax.set_yticks(np.arange(A.shape[0]), minor=False)
    ax.xaxis.tick_top()
    ax.set_xticklabels(structs, minor=False, rotation = 'vertical', size=14)
    ax.set_yticklabels(structs, minor=False, size=14)
    
    return img

def show_features(c1, mat1, c2, mat2, sorted_ligs, t1='Cluster 1', t2='Cluster 2', size=1,mi=None,ma=None):
    num_l = len(sorted_ligs)
    mat3 = mat1 - mat2

    r_mat1 = get_rmsd_matrix(c1, sorted_ligs)
    r_mat2 = get_rmsd_matrix(c2, sorted_ligs)

    f, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(1, 5, gridspec_kw = 
                                                {'width_ratios':[1,1]+[num_l]*2 + [num_l*1.1]}, dpi=400)

    ax0.matshow(r_mat1, cmap='Set1', vmin=0, vmax=9)
    ax0.set_yticks(np.arange(len(sorted_ligs)), minor=False)
    ax0.set_yticklabels(sorted_ligs, minor=False, size=14)

    ax1.matshow(r_mat2, cmap='Set1', vmin=0, vmax=9)
    ax1.set_yticks(np.arange(len(sorted_ligs)), minor=False)

    feature_heatmap(mat1, sorted_ligs, ma=ma, mi=mi, fig=f, ax=ax2, size=size)
    feature_heatmap(mat2, sorted_ligs, ma=ma, mi=mi, fig=f, ax=ax3, size=size)
    img = feature_heatmap(mat3, sorted_ligs, fig=f, ax=ax4, difference=True, size=size)

    ax0.set_xticks(np.arange(1), minor=False)
    ax0.xaxis.tick_top()
    ax0.set_xticklabels(['RMSD 1'], minor=False, rotation = 'vertical', size=14)
    ax1.set_xticks(np.arange(1), minor=False)
    ax1.xaxis.tick_top()
    ax1.set_xticklabels(['RMSD 2'], minor=False, rotation = 'vertical', size=14)

    ax2.set_xlabel(t1, size=20)
    ax3.set_xlabel(t2, size=20)
    ax4.set_xlabel('Difference', size=20)

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="8%", pad=0.05)
    plt.colorbar(img, cax=cax)

    ax1.set_yticklabels(['' for i in sorted_ligs])
    ax2.set_yticklabels(['' for i in sorted_ligs])
    ax3.set_yticklabels(['' for i in sorted_ligs])
    ax4.set_yticklabels(['' for i in sorted_ligs])

    plt.tight_layout()
    plt.show()
