import os
import sys
import numpy as np
import pandas as pd

from containers import Protein
from score.prob_opt import PredictStructs
from score.density_estimate import DensityEstimate

from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

class ScoreContainer:
    """
    Manages execution of ComBind for a given protein and settings.
    """
    def __init__(self, paths, params, features, stats_root, prot,
                 struct=None, num_poses=100, alpha=1.0,
                 max_iterations=1000000, restart=500):

        # Initialize PredictStructs object.
        self.features = features
        self.num_poses = num_poses
        self.alpha = alpha
        self.native_thresh = params['native_thresh']
        self.max_iterations = max_iterations
        self.restart = restart
        self.stats = self.read_stats(stats_root)
        self.ps = PredictStructs(self.stats, self.features, self.num_poses, self.alpha)

        # Initialize docking data.
        self.predict_data = Protein(prot, params, paths)
        if struct is not None:
            self.predict_data.lm.st = struct

    def read_stats(self, stats_root):
        stats = {'native':{}, 'reference':{}}
        for dist in stats:
            for interaction in self.features:
                fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
                assert os.path.exists(fname), fname
                stats[dist][interaction] = DensityEstimate.read(fname)
        return stats

    def load_docking(self, ligands, xtal=[]):
        self.predict_data.load_docking(ligands, load_fp=True,
                                       load_mcss='mcss' in self.features,
                                       load_shape='shape' in self.features)

        for ligand in xtal:
            self.predict_data.docking[ligand].load_native_poses(True,
                                                                self.native_thresh)

            if not self.predict_data.docking[ligand].poses:
                print(('WARNING: Ligand {} '
                      'was specified as an XTAL ligand '
                      'but it has no correct poses. It is being removed from '
                      'the query list').format(ligand))
                ligands.remove(ligand)

        ligands = {lig: self.predict_data.docking[lig] for lig in ligands}
        self.ps.set_ligands(ligands, self.predict_data.lm.mcss, self.predict_data.lm.shape)

    def compute_results(self, queries, xtal=[]):
        self.load_docking(queries, xtal)
        best_cluster = self.ps.max_posterior(max_iterations=self.max_iterations,
                                             restart=self.restart)
        return best_cluster

    def write_results(self, cluster, fname):
        """
        Part 1:
        ligand id, combind rank, combind rmsd, glide rank, glide rmsd, best rank, best rmsd
        Part 2:
        combind score, glide score, best score (if all pdb, else 0)
        """
        with open(fname, 'w') as f:
            f.write('lig,combind_rank,combind_rmsd,glide_rank,glide_rmsd,best_rank,best_rmsd\n')
            glide_cluster, best_cluster = {}, {}
            for lig, combind_pose in sorted(cluster.items()):
                poses = self.predict_data.docking[lig].poses
                best_rmsd = float('inf')
                for i, pose in enumerate(poses[:self.num_poses]):
                    if pose.rmsd is not None and pose.rmsd < best_rmsd:
                        best_cluster[lig] = i
                        best_rmsd = pose.rmsd
                
                best_pose = best_cluster[lig] if lig in best_cluster else 0
                f.write(','.join(map(str, [lig,
                                           poses[combind_pose].rank, poses[combind_pose].rmsd,
                                           poses[0].rank, poses[0].rmsd,
                                           poses[best_pose].rank, poses[best_pose].rmsd
                                           ]))+'\n')
            f.write('combind={},glide={},best={}\n'.format(
                    self.ps.log_posterior(cluster),
                    self.ps.log_posterior(glide_cluster),
                    self.ps.log_posterior(best_cluster) if len(best_cluster) == len(cluster) else 0))

    def read_results(self, fname):
        cluster = {}
        with open(fname) as fp:
            fp.readline()
            for line in fp:
                if line[:7] == 'combind': continue

                lig, combind_pose, *rest = line.split(',')
                cluster[lig] = int(combind_pose)
        return cluster

################################################################################

    def plot_interactions(self, combind, feature):
        glide = {k:0 for k in combind}
        interactions = self.interactions([glide, combind], feature)
        ligands = sorted(combind)

        if not interactions:
            return

        f, ax = plt.subplots(1, 2, figsize=(10, 4))
        self.gel_plot(combind, feature, ligands, interactions=interactions, ax=ax[0])
        self.gel_plot(glide, feature, ligands, interactions=interactions, ax=ax[1])
        
        ax[1].set_yticklabels([])
        plt.tight_layout()
        plt.savefig(feature+'.pdf')
        plt.show()

    def likelihood_and_feature_matrix(self, pose_cluster, k, lig_order):
        """
        Returns the feature values and likelihood ratios, P(X|l)/P(X)
        for feature 'fname' for poses in 'pose_cluster'
        as len(pose_cluster) x len(pose_cluster) numpy arrays.
        """
        x                    = np.zeros((len(lig_order), len(lig_order)))
        log_likelihood_ratio = np.zeros((len(lig_order), len(lig_order)))

        for i, ligname1 in enumerate(lig_order):
            for j, ligname2 in enumerate(lig_order):
                if j <= i: continue
                x_k, p_x_native, p_x = self.ps._likelihoods_for_pair_and_single_feature(k, pose_cluster,
                                                                                        ligname1, ligname2)
                x[i, j] = x_k
                log_likelihood_ratio[i, j] = np.log(p_x_native) - np.log(p_x)
        return x, log_likelihood_ratio

    def interactions(self, pose_clusters, k):
        interactions = set()
        for pose_cluster in pose_clusters:
            _interactions = set([interaction
                                 for ligand, pose in pose_cluster.items()
                                 for interaction in self.ps.ligands[ligand].poses[pose].fp
                                 if (interaction[0] in self.feature_defs[k])])
            interactions = interactions.union(_interactions)
        interactions = sorted(interactions)
        X = 0
        for pose_cluster in pose_clusters:
            X += self.interaction_matrix(pose_cluster, interactions, list(pose_cluster.keys()))

        # Remove interactions that are never seen
        interactions = [interaction
                        for interaction, present in zip(interactions, X.sum(axis = 1) > 0)
                        if present]
        X = X[X.sum(axis = 1) > 0]

        # Sort by interacton frequency
        idx = np.argsort(X.sum(axis = 1))[::-1]
        return [interactions[i] for i in idx]

    def interaction_matrix(self, pose_cluster, interactions, lig_order):
        X = []
        for ligand in lig_order:
            X += [[]]
            pose = self.ps.ligands[ligand].poses[pose_cluster[ligand]]
            for interaction in interactions:
                X[-1] += [pose.fp[interaction] if interaction in pose.fp else 0]
        return np.array(X).T

    def gel_plot(self, cluster, k, lig_order, interactions=None,
                 divide=2, resname_size=14, ax=None, pretty=True):
        if interactions is None:
            interactions = self.interactions([cluster], k)
        X = self.interaction_matrix(cluster, interactions, lig_order)
        labels = [self._format_int(interaction) for interaction in interactions]

        # Y is X + empty cells seperating boxes.
        Y = np.zeros((X.shape[0]*divide-1, X.shape[1]*divide-1))
        Y[::divide, ::divide] = X
        if ax is None:
            f, ax = plt.subplots(figsize = (9, 5))

        vmax = np.max(Y)
        ax.imshow(Y, cmap='binary', aspect = 'auto', vmin=0, vmax=vmax)

        if pretty:
            self._pretty(lig_order, labels, divide, resname_size, ax)

    def _pretty(self, lig_order, resnames, divide, resname_size, ax):
        for spine in ax.spines.values():
            spine.set_visible(False)
        for axis in ['x', 'y']:
            ax.tick_params(axis=axis, which='both',bottom=False,
                               top=False, left=False, right=False)
        ax.set_yticks(range(0, divide*len(resnames), divide))
        ax.set_yticklabels(resnames)
        ax.set_xticks(range(0, divide*len(lig_order), divide))
        ax.set_xticklabels(lig_order, rotation='vertical')

    def _format_int(self, interaction):
        three_to_one = {
            'TYR': 'Y',
            'VAL': 'V',
            'ARG': 'R',
            'THR': 'T',
            'GLU': 'E',
            'SER': 'S',
            'PHE': 'F',
            'ALA': 'A',
            'MET': 'M',
            'ILE': 'I',
            'ASP': 'D',
            'GLN': 'Q',
            'ASN': 'N',
            'GLY': 'G',
            'PRO': 'P',
            'TRP': 'W',
            'LEU': 'L',
            'LYS': 'K',
            'CYS': 'C',
            'HIS': 'H'
        }
        name = interaction[1].split('(')[1][:-1]
        num = interaction[1].split('(')[0].split(':')[1]
        if name in three_to_one:
            name = three_to_one[name]
        
        feature = [feature for feature, codes in self.feature.items()
                   if interaction[0] in codes]

        if 'hbond' in feature:
            feature.remove('hbond')
        assert len(feature) == 1
        feature = feature[0]

        if 'hbond_' in feature:
            feature = feature.replace('hbond_', '')
        return '{}{}:{}'.format(name, num, feature)

def score(paths, params, features, stats_root, protein, queries,
          xtal=[], pose_fname='poses.sc', struct=None, num_poses=100, alpha=1.0,
          max_iterations=1000000, restart=500):
    sc = ScoreContainer(paths, params, features, stats_root, protein,
                        num_poses=num_poses, alpha=alpha, struct=struct,
                        max_iterations=max_iterations, restart=restart)

    if 'all' in queries:
        assert len(queries) == 1
        queries = sc.predict_data.lm.docked(list(sc.predict_data.lm.pdb.keys()))

    if 'cross' in queries:
        assert len(queries) == 1
        queries = sc.predict_data.lm.docked(list(sc.predict_data.lm.pdb.keys()))
        self_dock = sc.predict_data.lm.st + '_lig'
        if self_dock in queries:
            queries.remove(self_dock)

    cluster = sc.compute_results(queries, xtal)
    sc.write_results(cluster, pose_fname)
