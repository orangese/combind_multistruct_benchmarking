import sys
import os
import numpy as np
sys.path.append('../../')
from shared_paths import shared_paths
shared_paths['data'] = '/Users/jpaggi/Downloads/combind_data/bpp_data'

from score.density_estimate import DensityEstimate
from score.prob_opt import PredictStructs
from score.statistics import statistics
from containers import Protein
from score.scores import ScoreContainer

import matplotlib.pyplot as plt


def load(prot, stats, scores, fname, struct):
    root = '{}/{}/scores/{}'.format(shared_paths['data'], prot, scores)
    stats_root = '{}/{}/scores/{}'.format(shared_paths['data'], prot, stats)
    sc = ScoreContainer(root, stats_root, prot, struct)
    cluster = sc.read_results(fname)
    glide_cluster = {ligand: 0 for ligand in cluster}
    best_cluster = {}
    for name, ligand in sc.ps.ligands.items():
        best_rmsd, best_pose = float('inf'), -1
        for i, pose in enumerate(ligand.poses[:100]):
            if pose.rmsd < best_rmsd:
                best_rmsd = pose.rmsd
                best_pose = i
        if best_rmsd < 2:
            best_cluster[name] = best_pose
    
    return cluster, glide_cluster, best_cluster, sc

def stats_plot(sc):
    features = sc.settings['k_list']
    f, ax = plt.subplots(len(features), 2, figsize = (10, 4 * len(features)))

    for i, feature in enumerate(sc.settings['k_list']):
        nat = sc.ps.stats['native'][feature]
        ref = sc.ps.stats['reference'][feature]
        ax[i, 0].plot(nat.x, nat.fx, c = 'g')
        ax[i, 0].plot(ref.x, ref.fx, c = 'b')

        ratio = nat.ratio(ref)
        ax[i, 1].plot(ratio.x, np.log(nat.fx) - np.log(ref.fx))
        ax[i, 0].set_ylabel(feature)
    plt.show()
    
def performance_plot(cluster, glide_cluster, sc):
    f, ax = plt.subplots(figsize = (4, 4))
    ligands = list(cluster)
    glide = [sc.ps.ligands[ligand].poses[glide_cluster[ligand]].rmsd for ligand in ligands]
    combind = [sc.ps.ligands[ligand].poses[cluster[ligand]].rmsd for ligand in ligands]
    best = [min(pose.rmsd for pose in sc.ps.ligands[ligand].poses) for ligand in ligands]
    plt.scatter(glide, combind)
    for b, g, c, ligand in zip(best, glide, combind, ligands):
        if b >= 2.0:
            print(ligand, b)
            plt.scatter([g], [c], c = 'orange')
            
        if c > 8:
            print(ligand, g, c)

    for i, txt in enumerate(ligands):
        txt = txt.replace('_lig', '')
        ax.annotate(txt, (glide[i]+.1, combind[i]+.1))

    m = max(glide+combind)
    plt.ylim(0, m+.5)
    plt.xlim(0, m+.5)
    plt.plot([0, m+.5], [0, m+.5], c = 'grey')
    ax.set_aspect('equal')
    plt.xlabel('Glide')
    plt.ylabel('ComBind')
    plt.show()

def energy_plot(cluster1, cluster2, sc):
    ligands = [ligand for ligand in cluster1 if ligand in cluster2]
    for feature in sc.settings['k_list']:
        x1, like1 =  sc.ps.likelihood_and_feature_matrix(cluster1, feature, ligands)
        x2, like2 =  sc.ps.likelihood_and_feature_matrix(cluster2, feature, ligands)

        print('Cluster1 - Cluster2 = {} for {}'.format(np.sum(like1-like2), feature))
        f, ax = plt.subplots(1, 2, figsize = (10, 4))
        if feature != 'mcss':
            ax[0].imshow(x1-x2, vmin = -2, vmax = 2)
        else:
            ax[0].imshow(x1-x2)
        ax[0].set_yticks(range(len(ligands)))
        ax[0].set_xticks(range(len(ligands)))
        ax[0].set_yticklabels(ligands)
        ax[0].set_xticklabels(ligands, rotation = 'vertical')


        cb = ax[1].imshow(like1-like2, cmap = 'coolwarm', vmin = -1, vmax = 1)
        ax[1].set_yticks(range(len(ligands)))
        ax[1].set_xticks(range(len(ligands)))
        ax[1].set_yticklabels(ligands)
        ax[1].set_xticklabels(ligands, rotation = 'vertical')
        plt.show()
