import os
import sys
import numpy as np

import config
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
        self.params = config.STATS[self.settings['stats_version']]
 
        self.predict_data = Protein(prot, params, paths)
        
        if self.struct is None:
            self.struct = self.predict_data.lm.st

        features = {feature: feature_defs[feature]
                    for feature in self.settings['features']}

        self.ps = PredictStructs({}, self.predict_data.lm.mcss,
                                 self.stats, features,
                                 self.settings['num_poses'],
                                 self.settings['alpha'])

    def read_settings(self):
        tr = {'num_poses': 100,
              'alpha': 1.0,
              'stats_version': 'default',
              'features': ['hbond', 'sb', 'contact', 'mcss']}
        if os.path.exists(settings_file):
            with open('{}/settings.py'.format(self.root)) as f:
                for line in f:
                    var, val = line.split('=')
                    tr[var] = eval(val)
        return tr

    def read_stats(self, stats_root):
        stats = {'native':{}, 'reference':{}}
        for dist in stats:
            for interaction in self.settings['features']:
                fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
                assert os.path.exists(fname), fname
                stats[dist][interaction] = DensityEstimate.read(fname)
        return stats

    def compute_results(self, queries):
        self.predict_data.load_docking(queries, load_fp = True,
                                       load_mcss = 'mcss' in self.settings['features'],
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
                                       load_mcss = 'mcss' in self.settings['features'],
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

def screen(paths, feature_defs, stats_root, struct, protein, queries):
    sc = ScoreContainer(os.getcwd(), paths, feature_defs,
                        stats_root, protein, struct)

    if 'all' in queries:
        queries = sc.predict_data.lm.docked(list(sc.predict_data.lm.pdb.keys()))

    #pose_cluster = sc.read_results(cluster)
    pose_cluster = {'CHEMBL135076_lig': 0,
                    'CHEMBL85194_lig': 0}

    all_ligands = queries + list(pose_cluster.keys())


    sc.predict_data.load_docking(all_ligands, load_fp = True,
                                 #load_mcss = 'mcss' in sc.settings['features'],
                                 st=sc.struct)

    sc.ps.ligands = sc.predict_data.docking[sc.struct]

    affinities, gscores, cscores = [], [], []
    for query in queries:
        ligand = sc.predict_data.docking[sc.struct][query]
        affinities += [sc.predict_data.lm.pdb[query]['AFFINITY']]
        gscores += [ligand.poses[0].gscore]
        cscores += [-sc.ps.score_new_ligand(pose_cluster, ligand)]

    affinities = np.log(affinities)-9

    from scipy import stats

    def enrich(affinities, scores, p=0.1, cut=-7):
        N = int(p*len(affinities))
        print(N)
        sorted_affinities = affinities[np.argsort(scores)]
        top = np.mean(sorted_affinities[:N] < cut)
        overall = np.mean(sorted_affinities < cut)
        return  top / overall 

    rho, _ = stats.spearmanr(affinities, gscores)
    r, _ = stats.pearsonr(affinities, gscores)
    print(rho, r, enrich(affinities, gscores))
    plt.scatter(affinities, gscores, s=2, c='k')
    plt.savefig('gscore.pdf')
    plt.close()

    rho, _ = stats.spearmanr(affinities, cscores)
    r, _ = stats.pearsonr(affinities, cscores)
    print(rho, r, enrich(affinities, cscores))
    plt.scatter(affinities, cscores, s=2, c='k')
    plt.savefig('cscore.pdf')

def main(paths, feature_defs, stats_root, protein, queries,
         fname='pdb.sc', plot=None, struct=None):
    sc = ScoreContainer(os.getcwd(), paths, feature_defs,
                        stats_root, protein, struct)

    combind_cluster = sc.compute_results(queries)
    sc.write_results(combind_cluster, fname)

    if plot:
        for feature in feature_defs:
            if feature != 'mcss':
                sc.plot_interactions(combind_cluster, feature)
