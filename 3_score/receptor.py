import os
import sys

from matplotlib import pyplot as plt
import numpy as np

import load_data as load
import objective as obj

data = '/scratch/PI/rondror/docking_data'
export_script = '/share/PI/rondror/jbelk/combind/3_score/export_cluster.py'

class Cluster:
    def __init__(self, ligs, poses, all_data):
        self.ligs = ligs
        self.all_data = all_data

        self.our_data = {ligs[i]:all_data[ligs[i]][poses[i]] for i in range(len(self.ligs))}
        self.glide_data = {ligs[i]:all_data[ligs[i]][0] for i in range(len(self.ligs))}
        self.min_data = {ligs[i]:all_data[ligs[i]][p] for i, p in enumerate(self.min_rmsd_poses())}

    def min_rmsd_poses(self):
        return [np.argmin([self.all_data[l][i].rmsd for i in range(min(25, len(self.all_data[l])))]) for l in self.ligs]

    def score_breakdown(self):
        phys_score = obj.objective(self.our_data, -2.7, 0)
        pair_score = obj.objective(self.our_data, 0, 1)
 
        n = float(len(self.ligs))
        all_interactions = {}
        for i in range(len(self.ligs)):
            for j in range(i+1, len(self.ligs)):
                p1 = self.our_data[self.ligs[i]]
                p2 = self.our_data[self.ligs[j]]
                for r in p1.fp:
                    if r in p2.fp:
                        prod = np.multiply(np.sqrt(p1.fp[r]), np.sqrt(p2.fp[r]))#/n**2 # (n*(n-1))
                        if r not in all_interactions:
                            all_interactions[r] = prod
                        else:
                            all_interactions[r] = np.add(all_interactions[r], prod)
        return phys_score, pair_score, all_interactions

    def plot_final_rmsds(self, title=''):

        for i, d in [('min',self.min_data), ('glide',self.glide_data), ('us',self.our_data)]:
            rmsd_list = [d[l].rmsd for l in self.ligs]
            plt.plot(rmsd_list, marker='.', markersize=10, label='{}: {}'.format(i, str(np.mean(rmsd_list))[:4]))

        plt.legend()
        plt.gca().set_xticklabels(self.ligs, minor=False, rotation='vertical')
        plt.gca().set_xticks(np.arange(0,len(self.ligs),1))
        plt.gca().set_ylim([0,12])

        plt.title(title)
        plt.gca().set_xlabel('Ligands')
        plt.gca().set_ylabel('RMSD [A]')

        plt.show()

    def export(self, cluster_name, receptor, struct):
        glide_dir = '{}/{}/xglide'.format(data, receptor)
        out_dir = '{}/{}/cluster_outputs'.format(data, receptor)
        command = [cluster_name, glide_dir, out_dir]
        for l in self.ligs:
            command.extend(['{}_ligand-to-{}'.format(l, struct), str(self.our_data[l].rank)])
        os.system('$SCHRODINGER/run {} {}'.format(export_script, ' '.join(command)))

class Structure:
    def __init__(self, dataset, struct):
        self.dataset = dataset
        self.struct = struct

        self.all_poses = load.load_data(dataset, struct)
        self.all_ligs = sorted(self.all_poses.keys())
        self.all_clusters = self.load_clusters()

    def load_clusters(self):
        all_combos = load.load_combinations(self.dataset, self.struct)
        all_clusters = {}
        for k in all_combos:
            all_clusters[k] = []
            for l, p in all_combos[k].items():
                all_clusters[k].append(Cluster(l, p, self.all_poses))
        return all_clusters

    def biggest_cluster(self):
        max_ligs = max(self.all_clusters.keys())
        if max_ligs + 1 != len(self.all_ligs):
            print 'Warning: {} ligs in biggest cluster, {} total ligs'.format(max_ligs, len(self.all_ligs))
        if len(self.all_clusters[max_ligs]) != 1:
            print 'Warning: {} clusters with {} ligs found. {} total ligs'.format(len(self.all_clusters[max_ligs]), max_ligs, len(self.all_ligs))
        return self.all_clusters[max_ligs][0]

    def top_poses(self, l, n, glide=True):
        pass

    def view_interaction_heatmap(self):
        pass

    def view_interaction_difference(self, other_cluster):
        pass

    def view_objective_breakdown(self, cluster):
        pass

    def view_performance(self, n, ligs=None):
        pass

class Receptor:
    def __init__(self, dataset, structs):
        self.name = dataset
        self.true_poses = load.load_data(dataset, 'crystals_only')
        if structs == 'all':
            self.structs = {s: Structure(dataset, s) for s in os.listdir('{}/{}/grids'.format(data, dataset))}
        else:
            self.structs = {s: Structure(dataset, s) for s in structs}
