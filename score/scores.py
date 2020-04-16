import os
import sys
import numpy as np

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
    def __init__(self, root, paths, feature_defs, stats_root, prot, struct):
        self.prot = prot
        self.struct = struct
        self.root = root

        self.feature_defs = feature_defs
        self.settings = self.read_settings()
        self.stats = self.read_stats(stats_root)
 
        self.predict_data = Protein(prot, self.settings['stats'], paths)
        features = {feature: feature_defs[feature] for feature in self.settings['k_list']}
        self.ps = PredictStructs({}, self.predict_data.lm.mcss,
                                 self.stats, features,
                                 self.settings['num_poses'],
                                 self.settings['alpha'])

    def read_settings(self):
        tr = {'chembl': False}
        with open('{}/settings.py'.format(self.root)) as f:
            for line in f:
                var, val = line.split('=')
                tr[var] = eval(val)
        return tr

    def read_stats(self, stats_root):
        stats = {'native':{}, 'reference':{}}
        for dist in stats:
            for interaction in self.settings['k_list']:
                fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
                assert os.path.exists(fname), fname
                stats[dist][interaction] = DensityEstimate.read(fname)
        return stats

    def compute_results_chembl(self, query):
        assert self.settings['chembl']
        chembl_ligs = self.predict_data.lm.get_helpers(query, self.settings['chembl_file'],
                                          num=self.settings['num_pred_chembl'],
                                          struct=self.struct)
        return self.compute_results([query]+chembl_ligs)

    def compute_results(self, queries):
        self.predict_data.load_docking(queries, load_fp = True,
                                       load_mcss = 'mcss' in self.settings['k_list'],
                                       st=self.struct)

        # Set ligands and optimize!
        self.ps.ligands = {lig: self.predict_data.docking[self.struct][lig]
                           for lig in queries}
        best_cluster = self.ps.max_posterior()
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
                poses = self.predict_data.docking[self.struct][lig].poses
                best_rmsd = float('inf')
                for i, pose in enumerate(poses[:self.settings['num_poses']]):
                    if pose.rmsd is not None and pose.rmsd < best_rmsd:
                        best_cluster[lig] = i
                        best_rmsd = pose.rmsd
                
                best_pose = best_cluster[lig] if lig in best_cluster else None
                f.write(','.join(map(str, [lig,
                                           combind_pose, poses[combind_pose].rmsd,
                                           0, poses[0].rmsd,
                                           best_pose, best_rmsd
                                           ]))+'\n')
            f.write('combind={},glide={},best={}\n'.format(
                    self.ps.log_posterior(cluster),
                    self.ps.log_posterior(glide_cluster),
                    self.ps.log_posterior(best_cluster) if len(best_cluster) == len(cluster) else 0))

    def read_results(self, fname):
        cluster = {}
        with open('{}/{}'.format(self.root, fname)) as fp:
            fp.readline()
            for line in fp:
                if line[:7] == 'combind': continue

                lig, combind_pose, *rest = line.split(',')
                cluster[lig] = int(combind_pose)

        self.predict_data.load_docking(list(cluster.keys()), load_fp = True,
                                       load_mcss = 'mcss' in self.settings['k_list'],
                                       st = self.struct)
        self.ps.ligands = {lig: self.predict_data.docking[self.struct][lig]
                           for lig in list(cluster.keys())}
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
        
        feature = [feature for feature, codes in self.feature_defs.items()
                   if interaction[0] in codes]

        if 'hbond' in feature:
            feature.remove('hbond')
        assert len(feature) == 1
        feature = feature[0]

        if 'hbond_' in feature:
            feature = feature.replace('hbond_', '')
        return '{}{}:{}'.format(name, num, feature)

def main(paths, feature_defs, stats_root, struct, protein, queries, plot=None):
    sc = ScoreContainer(os.getcwd(), paths, feature_defs,
                        stats_root, protein, struct)

    if sc.settings['chembl']:
        for query in queries:
            combind_cluster = sc.compute_results_chembl(query)
            fname = '{}/{}-to-{}.sc'.format(sc.root, query, struct)
            sc.write_results(combind_cluster, fname)
    else:
        # This needs to be before compute results as therein queries is mutated
        fname = '{}/{}.sc'.format(sc.root, 'pdb' if len(queries) > 1 else queries[0])
        combind_cluster = sc.compute_results(queries)
        sc.write_results(combind_cluster, fname)

        if plot:
            for feature in feature_defs:
                if feature != 'mcss':
                    sc.plot_interactions(combind_cluster, feature)
